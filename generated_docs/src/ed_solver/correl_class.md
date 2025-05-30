# `src/ed_solver/correl_class.f90`

## Overview

*   **Purpose:** This module defines a Fortran derived type, `correl_type`, and a comprehensive suite of associated procedures for creating, managing, manipulating, and performing I/O operations on matrix-valued physical correlation functions. These functions are typically dependent on frequency, such as Green's functions, self-energies, or hybridization functions.
*   **Role in Project:** Provides a standardized and versatile data structure for handling various types of correlation functions that are central to Exact Diagonalization (ED) solvers, Dynamical Mean-Field Theory (DMFT), and other many-body electronic structure methods. The module supports different frequency domains (Matsubara or real axis), particle statistics (fermionic or bosonic), and allows for efficient storage and manipulation through masked matrix structures.

## Key Components

*   **`TYPE correl_type`**:
    A Fortran derived type to encapsulate a matrix-valued correlation function.
    *   `CHARACTER(LEN=100) :: title`: A descriptive title or name for the correlation function (e.g., "SelfEnergy_iw", "GreensFunction_rt").
    *   `INTEGER :: N`: The dimension of the square matrix (e.g., number of orbitals or sites).
    *   `INTEGER :: Nw`: The number of frequency points at which the function is defined.
    *   `CHARACTER(LEN=9) :: stat`: Particle statistics, typically `FERMIONIC` or `BOSONIC`.
    *   `TYPE(masked_cplx_matrix_type) :: MM`: An object from `masked_matrix_class_mod` that defines the mask and structure of non-zero elements for a single N x N complex matrix. This mask is applied at each frequency point.
    *   `COMPLEX(DP), POINTER :: fctn(:,:,:) => NULL()`: A 3D allocatable array of shape `(N, N, Nw)` that stores the actual complex values of the correlation function \(C_{ij}(\omega_k)\).
    *   `COMPLEX(DP), POINTER :: vec(:,:) => NULL()`: A 2D allocatable array of shape `(MM%MASK%nind, Nw)`. For each frequency point `iw`, `vec(:,iw)` stores the independent (non-zero, unique by symmetry if applicable) elements of the matrix `fctn(:,:,iw)` as a flat vector, according to the mask in `MM`.
    *   `TYPE(freq_type) :: freq`: An object from `frequency_class` that defines the frequency mesh (e.g., Matsubara frequencies, real frequency grid details).

*   **`INTERFACE new_correl`**: A generic interface for constructing `correl_type` objects.
    *   `MODULE PROCEDURE new_correl_from_scratch_rfreq(CORREL, title, N, Nw, wmin, wmax, width, STAT, IMASK, AB)`: Creates a `correl_type` object for real frequencies, initializing its structure and frequency mesh.
    *   `MODULE PROCEDURE new_correl_from_scratch_ifreq(CORREL, title, N, Nw, beta, STAT, IMASK, AB)`: Creates for Matsubara (imaginary) frequencies.
    *   `MODULE PROCEDURE new_correl_from_scratch(CORREL, title, N, Nw, STAT, IMASK, AB)`: A core routine that initializes the structure (N, Nw, title, stat, mask) without setting up the specific frequency array.
    *   `MODULE PROCEDURE new_correl_from_old(CORRELOUT, CORRELIN)`: Creates `CORRELOUT` as a deep copy of an existing `CORRELIN` object.

*   **Data Conversion Subroutines:**
    *   **`SUBROUTINE correl2vec(CORREL, extmat, vecext)`**: Flattens the matrix representation `CORREL%fctn` (or an optional external matrix `extmat`) into its vectorized form `CORREL%vec` (or optional `vecext`), using the sparsity mask defined in `CORREL%MM`. This is done for each frequency point.
    *   **`SUBROUTINE vec2correl(CORREL, CORRELEXT)`**: Reconstructs the full matrix `CORREL%fctn` from its vectorized representation `CORREL%vec` (or an optional external vector `CORRELEXT%vec`), using the mask in `CORREL%MM`.

*   **Lifecycle and I/O Subroutines:**
    *   **`SUBROUTINE copy_correl(CORRELOUT, CORRELIN)`**: Performs a deep copy from `CORRELIN` to `CORRELOUT`.
    *   **`SUBROUTINE delete_correl(CORREL)`**: Destructor for `correl_type` objects, deallocating all pointer components.
    *   **`SUBROUTINE write_raw_correl(CORREL, UNIT)` / `SUBROUTINE read_raw_correl(CORREL, UNIT, FILEIN)`**: For unformatted (binary) input/output of `correl_type` objects.
    *   **`SUBROUTINE write_correl(CORREL, FILEOUT)` / `SUBROUTINE read_correl(CORREL, FILEIN)`**: For formatted text input/output. `write_correl` also includes commented-out calls to plotting routines.
    *   **`SUBROUTINE glimpse_correl(CORREL, UNIT)`**: Prints a "glimpse" (the matrix at `CORREL%freq%iwprint`) of the correlation function to a specified unit.

*   **Manipulation and Analysis Subroutines:**
    *   **`SUBROUTINE diff_correl(diff, CORREL1, CORREL2)`**: Calculates a weighted difference metric `diff` between two `correl_type` objects (`CORREL1` and `CORREL2`). This is often used as an objective function in fitting procedures (e.g., fitting a model hybridization to a target). It includes options for frequency windowing (`window_hybrid`, `window_hybrid2`), frequency point limiting (`fit_nw`), and power-law weighting (`fit_weight_power`, `fit_shift`).
    *   **`SUBROUTINE average_correlations(G, Gm1, offdiag_also, mask, boundup)`**: Computes an averaged version of the correlation function `G`, storing the result in `Gm1`. The averaging is performed based on a provided `mask` matrix and can optionally include off-diagonal elements. It handles symmetries for fermionic/bosonic statistics and different frequency types (Matsubara/real) by applying appropriate transformations (e.g., `CONJG`, `mirror_array`).
    *   **`SUBROUTINE pad_correl(CORRELOUT, CORRELIN)`**: "Pads" `CORRELOUT` using data from `CORRELIN`. If `CORRELIN` is provided, it first copies the vectorized elements from `CORRELIN%vec` to corresponding elements in `CORRELOUT%vec` based on matching mask indices. Then, it expands `CORRELOUT%fctn` and seems to fill elements based on mask equality, ensuring that if `imat(i1,i2) == imat(ii1,ii2)`, then `fctn(i1,i2) = fctn(ii1,ii2)`.
    *   **`SUBROUTINE slice_correl(CORRELOUT, CORRELIN, rbounds, cbounds)`**: Extracts a sub-matrix (a "slice") defined by row bounds `rbounds` and column bounds `cbounds` from `CORRELIN` and stores it in `CORRELOUT`. The frequency data and other properties are copied.
    *   **`SUBROUTINE transform_correl(CORRELOUT, CORRELIN)`**: Performs a transformation \(C_{out}(z) = (-1)^F C_{in}(-z^*)\), where \(F=0\) for bosons and \(F=1\) for fermions. For real frequencies, \(z^* = z\), so it involves \(C_{out}(\omega) = (-1)^F C_{in}(-\omega)\), which is achieved using `mirror_array`. For Matsubara frequencies, \((-i\omega_n)^* = i\omega_n\), so it's \(C_{out}(i\omega_n) = (-1)^F C_{in}(-i\omega_n)\), effectively \(C_{out}(i\omega_n) = (-1)^F \text{conj}(C_{in}(i\omega_n))\) due to \(C(i\omega_n)^* = C(-i\omega_n)\).
    *   **`SUBROUTINE dump_chi0(G, beta)`**: Calculates a bare susceptibility \(\chi_0(\Omega_m) \sim \sum_{\omega_n} G(i\omega_n) G(i\omega_n + i\Omega_m)\) from a given Green's function `G`. Includes commented-out plotting calls.

## Important Variables/Constants

Within the `correl_type` derived type:

*   **`fctn(:,:,:) :: COMPLEX(DP), POINTER`**: This 3D array is the core data store for the correlation function values \(C_{ij}(\omega_k)\), where \(i,j\) are matrix indices and \(k\) is the frequency index.
*   **`vec(:,:) :: COMPLEX(DP), POINTER`**: The vectorized form of the independent matrix elements. For each frequency `iw`, `vec(:,iw)` holds the elements selected by `MM%MASK%ivec`. This form is particularly useful for optimization algorithms.
*   **`MM :: TYPE(masked_cplx_matrix_type)`**: An instance of `masked_cplx_matrix_type` that defines the sparsity structure (mask) of the `N x N` matrix. This mask determines which elements are stored and considered non-zero. `MM%mat` points to `fctn(:,:,iwprint)`.
*   **`freq :: TYPE(freq_type)`**: An object (from `frequency_class.f90`) that encapsulates all information about the frequency mesh on which the correlation function is defined (e.g., Matsubara frequencies, real frequency grid, temperature/beta).
*   **`stat :: CHARACTER(LEN=9)`**: A string indicating the particle statistics, typically `FERMIONIC` or `BOSONIC`. This affects transformations and symmetries.

## Usage Examples

```fortran
MODULE example_correl_usage
  USE correl_class
  USE frequency_class, ONLY: new_freq ! For detailed frequency setup if not using new_correl_*freq
  USE genvar, ONLY: FERMIONIC, BOSONIC, DP

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: test_correl_module

CONTAINS
  SUBROUTINE test_correl_module
    TYPE(correl_type) :: g_omega_fermionic
    TYPE(correl_type) :: chi_omega_bosonic
    INTEGER :: num_orbitals, num_frequencies
    REAL(DP) :: temperature_beta

    num_orbitals = 2
    num_frequencies = 256
    temperature_beta = 50.0

    ! Create a Fermionic Matsubara frequency correlation function (e.g., Green's function)
    CALL new_correl_from_scratch_ifreq(g_omega_fermionic, title="Green_Func_iw", &
     &                                 N=num_orbitals, Nw=num_frequencies, beta=temperature_beta, &
     &                                 STAT=FERMIONIC)
    ! g_omega_fermionic%fctn can now be filled with data, e.g., from an impurity solver

    ! Create a Bosonic Matsubara frequency correlation function (e.g., susceptibility)
    CALL new_correl_from_scratch_ifreq(chi_omega_bosonic, title="Susceptibility_iOmega", &
     &                                 N=num_orbitals, Nw=num_frequencies, beta=temperature_beta, &
     &                                 STAT=BOSONIC)
    ! chi_omega_bosonic%fctn can be filled

    ! Example: Convert matrix to vector form for g_omega_fermionic
    CALL correl2vec(g_omega_fermionic)
    ! Now g_omega_fermionic%vec contains the vectorized data

    ! Example: Write to file
    CALL write_correl(g_omega_fermionic, FILEOUT="green_function_output.dat")

    ! Clean up
    CALL delete_correl(g_omega_fermionic)
    CALL delete_correl(chi_omega_bosonic)

  END SUBROUTINE test_correl_module
END MODULE example_correl_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `frequency_class`: Essential for the `freq_type` component, which defines the frequency mesh and its properties.
    *   `genvar`: Provides global constants (`DP`, `FERMIONIC`, `BOSONIC`, `RETARDED`, `ADVANCED`, `pi2`), error/status variables (`istati`), and logging/MPI related variables (`log_unit`, `messages3`, `strongstop`).
    *   `globalvar_ed_solver`: Used for `istati` and various global parameters that control the behavior of `diff_correl` (e.g., `fit_nw`, `fit_shift`, `lambda_sym_fit`, `MASK_AVERAGE_SIGMA`, `weight_expo`, `window_hybrid`, `window_weight`).
    *   `masked_matrix_class_mod`: Crucial for the `masked_cplx_matrix_type` (`MM` component) and its associated procedures which handle the matrix structure, masking, and vectorization of individual frequency slices.
    *   `common_def`: For utility functions like `dump_message`, string conversion (`c2s`, `i2c`), and file operations (`open_safe`, `close_safe`, `skip_line`).
    *   `matrix`: For matrix utilities like `write_array` and `average_matrix`.
    *   `mesh`: The `mirror_array` routine is used in `average_correlations` and `transform_correl` for handling real-frequency data symmetries.
    *   `string5`: For the `get_unit` utility in file operations.
    *   `stringmanip` (used in `dump_chi0`): For `toString` function.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules. No explicit third-party external libraries are directly used by this module itself.

*   **Interactions with other components:**
    *   This module is central to any part of the codebase that deals with frequency-dependent quantities.
    *   **Impurity Solvers (e.g., ED solver):** Would use `correl_type` to receive hybridization functions as input and to output Green's functions and self-energies.
    *   **`bath_class_hybrid`:** This module likely uses `correl_type` extensively. For instance, `bath2hybrid` would populate a `correl_type` object representing the hybridization function. `hybrid2bath` would take a target `correl_type` (hybridization) and fit bath parameters to it, using `diff_correl` to measure the fit quality.
    *   **Post-processing tools:** Any tool analyzing results like spectral functions, susceptibilities, etc., would read data stored in files generated by `write_correl` or `write_raw_correl`.
    *   The vectorization routines (`correl2vec`, `vec2correl`) are key for interfacing with generic optimization algorithms that operate on flat 1D parameter vectors.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
