# `src/ed_solver/frequency_class.f90`

## Overview

*   **Purpose:** This module defines a Fortran derived type, `freq_type`, and associated procedures to create, manage, and describe various types of frequency meshes. These meshes are essential for calculations involving frequency-dependent quantities such as Green's functions, self-energies, and other correlation functions.
*   **Role in Project:** It provides a standardized data structure for frequency grids within an Exact Diagonalization (ED) or DMFT solver. The module supports both Matsubara (imaginary) frequencies, typically used in finite-temperature quantum field theory calculations, and real frequencies, often used for spectral functions and including phenomenological or physically motivated broadening. It also incorporates functionality for distributing frequency points across MPI processes for parallel computations.

## Key Components

*   **`TYPE freq_type`**:
    A Fortran derived type to encapsulate the properties of a frequency mesh.
    *   `INTEGER :: Nw = 0`: The total number of discrete frequency points in the mesh.
    *   `INTEGER :: iwprint = 0`: The 1-based index of a specific frequency point that might be singled out for printing or detailed analysis (e.g., \(\omega_n=0\) or \(\omega=0\)).
    *   `CHARACTER(LEN=9) :: title = '\0'`: A string identifying the type of frequency mesh. Expected values include `FERMIONIC`, `BOSONIC` (for Matsubara frequencies), or `RETARDED`, `ADVANCED` (for real frequencies).
    *   `REAL(DP) :: beta = 0`: The inverse temperature \(1/(k_B T)\). This is relevant only for Matsubara frequency meshes.
    *   `REAL(DP) :: wmax = 0`: The maximum real frequency value for a real frequency grid.
    *   `REAL(DP) :: wmin = 0`: The minimum real frequency value for a real frequency grid.
    *   `REAL(DP) :: width = 0`: A base constant broadening factor \(\eta\) added to the imaginary part for real frequency grids (e.g., \(\omega + i\eta\)).
    *   `COMPLEX(DP), POINTER :: vec(:) => NULL()`: A 1D allocatable array storing the actual complex frequency values (\(i\omega_n\) for Matsubara, or \(\omega + i\eta_{eff}\) for real frequencies).
    *   `INTEGER, POINTER :: iwmin(:) => NULL(), iwmax(:) => NULL(), chunk(:) => NULL()`: Arrays used for distributing frequency points across MPI processes. For each process, `iwmin` and `iwmax` define the starting and ending global indices of its local frequency chunk, and `chunk` stores the size of this local segment.

*   **Module-Private Broadening Parameters:**
    *   `REAL(8) :: rdelta_frequ_eta1, rdelta_frequ_w0, rdelta_frequ_T, rdelta_frequ_eta2`: These module-level private variables store parameters for a sophisticated, frequency-dependent broadening scheme that can be applied to real-frequency grids. This scheme involves terms based on Fermi functions, allowing the effective broadening \(\eta_{eff}\) to vary with \(\omega\). They are set by `init_frequency_module`.

*   **`INTERFACE new_freq`**: A generic interface for constructing `freq_type` objects.
    *   `MODULE PROCEDURE new_rfreq_from_scratch(FREQ, Nw, wmin, wmax, width, sense)`: Creates a real frequency grid. `sense` should be `RETARDED` or `ADVANCED`. The imaginary part includes the base `width` plus contributions from the Fermi-function-based broadening scheme if its parameters (`rdelta_frequ_eta1`, etc.) are non-zero.
    *   `MODULE PROCEDURE new_ifreq_from_scratch(FREQ, Nw, beta, stat)`: Creates a Matsubara frequency grid. `stat` should be `FERMIONIC` (\(\omega_n = (2n+1)\pi/\beta\)) or `BOSONIC` (\(\omega_n = 2n\pi/\beta\)).
    *   `MODULE PROCEDURE new_freq_from_old(FREQOUT, FREQIN)`: Creates `FREQOUT` as a deep copy of an existing `freq_type` object `FREQIN`.

*   **Public Management Subroutines:**
    *   **`SUBROUTINE init_frequency_module(rdelta_frequ_eta1_, rdelta_frequ_w0_, rdelta_frequ_T_, rdelta_frequ_eta2_)`**: Initializes the module-level private parameters (`rdelta_frequ_eta1`, `rdelta_frequ_w0`, `rdelta_frequ_T`, optional `rdelta_frequ_eta2`) used for the advanced frequency-dependent broadening on the real axis.
    *   **`SUBROUTINE copy_frequency(FREQOUT, FREQIN)`**: Performs a deep copy from `FREQIN` to `FREQOUT`.
    *   **`SUBROUTINE delete_frequency(FREQ)`**: Destructor for `freq_type` objects, deallocating all pointer components (`vec`, `iwmin`, `iwmax`, `chunk`).
    *   **`SUBROUTINE read_raw_freq(FREQ, UNIT)` / `SUBROUTINE write_raw_freq(FREQ, UNIT)`**: For unformatted (binary) input/output of `freq_type` objects, allowing for efficient saving and loading of frequency mesh definitions.

*   **Internal Subroutines:**
    *   **`SUBROUTINE split_frequency(FREQ)`**: Called after a frequency mesh is generated. It distributes the `Nw` frequency points among `nproc` MPI processes by populating the `FREQ%iwmin`, `FREQ%iwmax`, and `FREQ%chunk` arrays using the `split` routine (from `mpi_mod.f90`).

## Important Variables/Constants

*   **Within `freq_type`:**
    *   `vec(:) :: COMPLEX(DP), POINTER`: The array containing the actual complex frequency values. This is the primary data produced and used by this module.
    *   `Nw :: INTEGER`: The number of points in the frequency mesh.
    *   `title :: CHARACTER(LEN=9)`: Specifies the type of frequency axis (`FERMIONIC`, `BOSONIC`, `RETARDED`, `ADVANCED`).
    *   `beta :: REAL(DP)`: Inverse temperature, essential for Matsubara frequencies.
    *   `wmin, wmax, width :: REAL(DP)`: Define the range and base broadening for real-axis frequencies.
*   **Module-private broadening parameters (`rdelta_frequ_eta1`, etc.)**: These allow for a flexible, frequency-dependent broadening to be added to real-axis frequencies, potentially mimicking physical effects or improving numerical stability.

## Usage Examples

```fortran
MODULE example_frequency_usage
  USE frequency_class
  USE genvar, ONLY: FERMIONIC, RETARDED, DP, imi ! For imi if manually checking values

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: setup_frequency_meshes

CONTAINS
  SUBROUTINE setup_frequency_meshes
    TYPE(freq_type) :: mats_fermion_freq, real_retarded_freq
    INTEGER :: num_mats_points, num_real_points
    REAL(DP) :: temperature_beta, omega_min, omega_max, eta_broadening

    num_mats_points = 1024
    temperature_beta = 100.0_DP

    num_real_points = 500
    omega_min = -5.0_DP
    omega_max =  5.0_DP
    eta_broadening = 0.05_DP

    ! Initialize module-level parameters for advanced broadening (optional, for real frequencies)
    CALL init_frequency_module(rdelta_frequ_eta1_=0.0_DP, rdelta_frequ_w0_=0.0_DP, &
     &                         rdelta_frequ_T_=1.0_DP) ! Using simple constant broadening here

    ! 1. Create a Matsubara frequency mesh for fermionic particles
    CALL new_freq(mats_fermion_freq, Nw=num_mats_points, beta=temperature_beta, stat=FERMIONIC)

    PRINT *, "Matsubara Frequency Mesh ("//TRIM(mats_fermion_freq%title)//"):"
    PRINT *, "  Number of points (Nw): ", mats_fermion_freq%Nw
    PRINT *, "  Beta: ", mats_fermion_freq%beta
    IF (mats_fermion_freq%Nw > 0) THEN
      PRINT *, "  First Matsubara freq (omega_0): ", mats_fermion_freq%vec(1)
      PRINT *, "  omega_print (iwprint=", mats_fermion_freq%iwprint, "): ", mats_fermion_freq%vec(mats_fermion_freq%iwprint)
    END IF

    ! 2. Create a real frequency mesh for retarded functions
    CALL new_freq(real_retarded_freq, Nw=num_real_points, wmin=omega_min, wmax=omega_max, &
     &            width=eta_broadening, sense=RETARDED)

    PRINT *, "Real Frequency Mesh ("//TRIM(real_retarded_freq%title)//"):"
    PRINT *, "  Number of points (Nw): ", real_retarded_freq%Nw
    PRINT *, "  w_min, w_max, base_width: ", real_retarded_freq%wmin, real_retarded_freq%wmax, real_retarded_freq%width
    IF (real_retarded_freq%Nw > 0) THEN
      PRINT *, "  First real freq point (omega_0 + i*eta_eff): ", real_retarded_freq%vec(1)
    END IF

    ! Clean up
    CALL delete_frequency(mats_fermion_freq)
    CALL delete_frequency(real_retarded_freq)

  END SUBROUTINE setup_frequency_meshes
END MODULE example_frequency_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `genvar`: Provides `DP` (double precision kind), character constants for frequency types (`ADVANCED`, `BOSONIC`, `FERMIONIC`, `RETARDED`), `imi` (imaginary unit), `iwprint` (default print index), `pi`, `pi2` (mathematical constants), `nproc` (number of MPI processes), and `messages3` (for debug message control).
    *   `linalg`: The `dexpc` (complex exponential) and `ramp` functions are used in `new_rfreq_from_scratch` and `new_ifreq_from_scratch` respectively.
    *   `mpi_mod`: The `split` subroutine is used by the internal `split_frequency` routine to determine how frequency points are distributed among MPI processes.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules.
    *   MPI (Message Passing Interface) if the code is compiled and run in a parallel environment, as `split_frequency` uses MPI information via `nproc` (from `genvar`, which is MPI-aware) to set up distributed frequency chunks.

*   **Interactions with other components:**
    *   **`correl_class.f90`**: This is the primary consumer of `frequency_class`. The `freq_type` is a component of `correl_type`. When a `correl_type` object (representing a Green's function, self-energy, etc.) is created, an appropriate `freq_type` object is also created and associated with it to define its frequency domain.
    *   **Solver Routines:** Any part of the ED or DMFT solver that computes or manipulates frequency-dependent quantities will rely on the `freq%vec` array from a `freq_type` object to know the specific frequency points.
    *   **Parallel Execution:** The `split_frequency` routine and the `iwmin`, `iwmax`, `chunk` members of `freq_type` are essential for parallel algorithms that distribute calculations over frequency points. Each MPI process can then operate on its local subset of frequencies.
    *   The `init_frequency_module` allows for setting global parameters that define a potentially sophisticated, energy-dependent broadening scheme for real-axis frequencies. This can be useful for more realistic simulations of spectral functions.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
