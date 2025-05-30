# `src/ed_solver/fermion_hilbert_class.f90`

## Overview

*   **Purpose:** This module defines Fortran derived types (`fermion_sector_type` and the internal `fermion_Hilbert_type`) and associated procedures to construct, manage, and provide access to fermionic Hilbert spaces. A key feature is the decomposition of the full Hilbert space into distinct sectors, typically characterized by a fixed number of particles (or a fixed total Sz quantum number).
*   **Role in Project:** It provides the foundational data structures for representing many-body basis states in a fermionic system, specifically for Exact Diagonalization (ED) calculations. The module pre-calculates and stores all possible Fock states (represented as integers) for a given number of orbitals, organizing them into these quantum number sectors. This pre-computation allows for efficient lookup (ranking) and manipulation of states, which is critical for constructing Hamiltonian matrices and applying operators within ED solvers.

## Key Components

*   **`TYPE fermion_sector_type`**:
    Represents a specific sector of the fermionic Hilbert space, characterized by a fixed number of particles or Sz.
    *   `LOGICAL :: SZ = .false.`: If true, this sector is defined by a total Sz value rather than particle number `npart`. `npart` might then be re-interpreted (e.g. related to \(2S_z + N_{sites}\)).
    *   `LOGICAL :: not_physical = .false.`: Set to true if `npart < 0` or `npart > norbs`, indicating an invalid sector.
    *   `INTEGER :: norbs = 0`: The total number of single-particle orbitals in the system.
    *   `INTEGER :: npart = 0`: The number of occupied orbitals (i.e., number of fermions) in this particular sector.
    *   `INTEGER :: nstates = 0`: The total number of states in the full Hilbert space for `norbs` orbitals (i.e., \(2^{norbs}\)). This is stored for context/consistency.
    *   `CHARACTER(LEN=100) :: title = '\0'`: A descriptive title for the sector, e.g., "Nocc = 3" or "Sz*2 = 1".
    *   `INTEGER :: dimen = 0`: The dimension of this specific sector, i.e., the number of unique many-body states that have `npart` fermions distributed among `norbs` orbitals (calculated as \(C(\text{norbs}, \text{npart})\)).
    *   `INTEGER, POINTER :: state(:) => NULL()`: An allocatable array storing the integer representations of all unique many-body basis states belonging to this sector. Each integer represents a Fock state (e.g., bit `i` is 1 if orbital `i` is occupied, 0 otherwise). This points to a slice of the global `Hilbert%sector(npart)%state` array.
    *   `INTEGER, POINTER :: rank(:) => NULL()`: A pointer to the global `Hilbert%rank` array. This array maps an integer representation of any state in the full Hilbert space to its 1-based index within the `state` array of its specific particle-number sector.
    *   `INTEGER, POINTER :: chunk(:) => NULL()`: Array storing the number of states in this sector managed by each MPI process when the Hilbert space is distributed.
    *   `INTEGER, POINTER :: istatemin(:) => NULL()`: Array storing the starting global index of states for each MPI process.
    *   `INTEGER, POINTER :: istatemax(:) => NULL()`: Array storing the ending global index of states for each MPI process.

*   **`TYPE fermion_Hilbert_type` (Module-private, global instance `Hilbert` with `SAVE` attribute)**:
    Manages the entire Hilbert space for a fixed `norbs`. It is constructed once per `norbs` value.
    *   `INTEGER :: norbs`: Total number of orbitals.
    *   `INTEGER(8) :: dimen`: Total dimension of the full Hilbert space (\(2^{norbs}\)).
    *   `INTEGER :: nsectors`: Number of sectors, equal to `norbs + 1` (for particle numbers 0 to `norbs`).
    *   `INTEGER :: max_dimen = 2**30 + 1`: A hardcoded limit for `dimen`. If `2**norbs` exceeds this, storing the full `rank` array is deemed too memory-intensive, and the program stops.
    *   `TYPE(fermion_sector_type), POINTER :: sector(:) => NULL()`: An array (indexed 0 to `norbs`) where each element is a `fermion_sector_type` object containing all states for that specific particle number.
    *   `INTEGER, POINTER :: rank(:) => NULL()`: A global array of size \(2^{norbs}\). For any state (integer representation `s`), `rank(s)` gives its 1-based index within the `state` array of the sector to which `s` belongs (i.e., `Hilbert%sector(npart_of_s)%state(rank(s)) == s`).

*   **`INTERFACE new_fermion_sector`**: A generic interface for creating `fermion_sector_type` objects.
    *   `MODULE PROCEDURE new_fermion_sector_from_scratch(SEC, npart, norbs, SZ)`: This is the primary constructor for users.
        *   If this is the first time it's called for a given `norbs` (or if `norbs` changes), it triggers the construction of the global `Hilbert` object by calling `new_fermion_Hilbert_space(Hilbert, norbs, SZ)`.
        *   It then sets up the output `SEC` by making its `state` and `rank` pointers point to the appropriate data within the global `Hilbert` object (e.g., `SEC%state => Hilbert%sector(npart)%state`).
        *   Calls `split_fermion_sector(SEC)` to define MPI distribution for this sector.
    *   `MODULE PROCEDURE new_fermion_sector_from_old(SECOUT, SECIN)`: Creates `SECOUT` by first calling `new_fermion_sector_from_scratch` and then `copy_fermion_sector`.

*   **Public Module Variables (Flags):**
    *   `LOGICAL :: FLAG_FORBID_SPLITTING = .false.`: If set to `.TRUE.`, prevents the distribution of Hilbert space sectors across MPI nodes, even if running in parallel.
    *   `LOGICAL :: HILBERT_SPACE_SPLITED_AMONG_NODES`: A status flag, set to `.TRUE.` if the Hilbert space sectors are currently being treated as distributed among MPI processes (depends on `FLAG_FORBID_SPLITTING`, `size2`, `no_mpi`).

*   **Public Management Subroutines:**
    *   **`SUBROUTINE copy_fermion_sector(SECOUT, SECIN)`**: Copies metadata and makes `SECOUT` pointers point to the same global `Hilbert` data as `SECIN`.
    *   **`SUBROUTINE delete_fermion_sector(SEC)`**: Nullifies pointers in `SEC` (`state`, `rank`) and deallocates MPI-specific arrays (`chunk`, `istatemin`, `istatemax`). It does *not* deallocate the global `Hilbert` data.
    *   **`SUBROUTINE read_raw_fermion_sector(SEC, UNIT, SZ)` / `SUBROUTINE write_raw_fermion_sector(sec, UNIT)`**: For unformatted (binary) I/O of sector definitions (primarily `norbs` and `npart` to allow recreation via `new_fermion_sector`).

*   **Internal Core Subroutines:**
    *   **`SUBROUTINE new_fermion_Hilbert_space(Hilbert, norbs, SZ)`**: This crucial routine builds the global `Hilbert` data structure for a given `norbs`.
        1.  Calculates total dimension \(2^{norbs}\). If too large ( > `max_dimen`), it stops.
        2.  Allocates `Hilbert%sector(0:norbs)`.
        3.  For each particle number `p` from 0 to `norbs`, calls `allocate_fermion_sector` to set dimensions and allocate `Hilbert%sector(p)%state`.
        4.  Allocates the global `Hilbert%rank(0:2**norbs - 1)` array.
        5.  Iterates through all possible integer states from 0 to \(2^{norbs}-1\). For each `state`:
            *   Determines its particle number `npart` (using `noccupied` from `fermion_ket_class`).
            *   Increments a counter for that `npart` sector (`dimen_check(npart)`).
            *   Stores the `state` in `Hilbert%sector(npart)%state(dimen_check(npart))`.
            *   Stores its rank: `Hilbert%rank(state) = dimen_check(npart)`.
        6.  Performs consistency checks on sector dimensions.
        7.  Makes each `Hilbert%sector(npart)%rank` point to the global `Hilbert%rank`.
    *   **`SUBROUTINE delete_fermion_Hilbert_space(Hilbert)`**: Destructor for the global `Hilbert` object, deallocating its `sector` array and `rank` array.
    *   **`SUBROUTINE allocate_fermion_sector(SEC, npart, norbs, SZ)`**: Sets metadata (`norbs`, `npart`, `title`), calculates sector dimension using `combination(norbs, npart)`, allocates `SEC%state`, and calls `split_fermion_sector`.
    *   **`SUBROUTINE split_fermion_sector(SEC)`**: Based on `HILBERT_SPACE_SPLITED_AMONG_NODES` flag (derived from MPI status and `FLAG_FORBID_SPLITTING`), it either sets the local chunk to the full sector size or calls `split` (from `mpi_mod`) to determine `SEC%chunk`, `SEC%istatemin`, `SEC%istatemax` for each MPI process.
    *   **`FUNCTION is_in_fermion_sector(state, SEC) AS LOGICAL`**: Checks if a given integer `state` has the correct particle number `SEC%npart` for the specified `SEC%norbs`.

## Important Variables/Constants

*   **Module-private `Hilbert :: TYPE(fermion_Hilbert_type), SAVE`**: This is the central, persistent data structure. It holds the pre-computed information (all states and their ranks within sectors) for a given number of orbitals (`norbs`). Once created for a specific `norbs`, it's reused by all `fermion_sector_type` instances that refer to that `norbs`.
*   **Within `fermion_sector_type`:**
    *   `state(:) :: INTEGER, POINTER`: Points to an array (part of the global `Hilbert` object) containing the integer representations of all Fock basis states belonging to this sector.
    *   `rank(:) :: INTEGER, POINTER`: Points to the global `Hilbert%rank` array, which provides a mapping from any integer state representation to its 1-based index within its specific sector's `state` array. This is crucial for efficient state lookup.
    *   `dimen :: INTEGER`: The dimension (number of states) of this specific sector.
    *   `chunk(:), istatemin(:), istatemax(:)`: Define how the states of this sector are distributed among MPI processes if parallel processing is enabled.

## Usage Examples

```fortran
MODULE example_fermion_hilbert_usage
  USE fermion_hilbert_class
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: explore_hilbert_sector

CONTAINS
  SUBROUTINE explore_hilbert_sector
    TYPE(fermion_sector_type) :: two_particle_sector
    INTEGER, PARAMETER :: num_orbitals_example = 4
    INTEGER, PARAMETER :: num_particles_example = 2
    INTEGER :: i, state_int, rank_in_sector
    LOGICAL :: is_sz_sector_example

    is_sz_sector_example = .FALSE. ! Using particle number conservation

    ! Create/initialize a sector for 2 particles in 4 orbitals.
    ! This will trigger the creation of the global 'Hilbert' object for 4 orbitals
    ! if it hasn't been created yet by a previous call with norbs=4.
    CALL new_fermion_sector(two_particle_sector, npart=num_particles_example, &
     &                      norbs=num_orbitals_example, SZ=is_sz_sector_example)

    PRINT *, "Sector Details:"
    PRINT *, "  Title: ", TRIM(two_particle_sector%title)
    PRINT *, "  Number of orbitals (norbs): ", two_particle_sector%norbs
    PRINT *, "  Number of particles (npart): ", two_particle_sector%npart
    PRINT *, "  Dimension of this sector: ", two_particle_sector%dimen

    IF (two_particle_sector%dimen > 0) THEN
      PRINT *, "  States in this sector (integer representation and rank within sector):"
      DO i = 1, two_particle_sector%dimen
        state_int = two_particle_sector%state(i)
        ! The rank pointer in fermion_sector_type points to the global Hilbert rank array.
        ! So, Hilbert%rank(state_int) gives the rank of state_int within its sector.
        rank_in_sector = two_particle_sector%rank(state_int)
        PRINT *, "    State (int): ", state_int, " (Binary: ", TO_BINARY(state_int, num_orbitals_example), ")", &
         &         " Rank in Sector: ", rank_in_sector
      END DO
    END IF

    ! Clean up the local sector structure (does not delete the global Hilbert data)
    CALL delete_fermion_sector(two_particle_sector)

    ! The global 'Hilbert' object for norbs=4 still exists due to SAVE attribute.
    ! To truly release it, one would need a separate call to a routine that calls
    ! delete_fermion_Hilbert_space(Hilbert) when norbs is expected to change or at program end.

  END SUBROUTINE explore_hilbert_sector

  ! Helper function to convert integer to binary string for printing
  FUNCTION TO_BINARY(val, num_bits) RESULT(bin_str)
    INTEGER, INTENT(IN) :: val, num_bits
    CHARACTER(LEN=num_bits) :: bin_str
    INTEGER :: i
    bin_str = ''
    DO i = num_bits - 1, 0, -1
      IF (BTEST(val, i)) THEN
        bin_str(num_bits-i:num_bits-i) = '1'
      ELSE
        bin_str(num_bits-i:num_bits-i) = '0'
      END IF
    END DO
  END FUNCTION TO_BINARY

END MODULE example_fermion_hilbert_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `common_def`: For string conversion utilities (`c2s`, `i2c`) and messaging (`dump_message`).
    *   `quantum_algebra`: The comment mentions `fermion`, but the direct usage for state properties (like `noccupied`) comes from `fermion_ket_class`.
    *   `genvar`: Provides MPI-related variables (`iproc`, `log_unit`, `messages3`, `no_mpi`, `nproc`, `size2`) and status variables (`istati`).
    *   `mpi_mod`: The `split` subroutine is used by `split_fermion_sector` for distributing states across MPI processes.
    *   `specialfunction` (used in `allocate_fermion_sector`): The `combination` function \(C(N,p)\) is used to calculate the dimension of a sector.
    *   `fermion_ket_class` (used in `new_fermion_Hilbert_space`): For `fermion_ket_type`, `new_ket`, and `noccupied` to determine the particle number of an integer-represented state.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules.
    *   MPI (Message Passing Interface) if the parallel features for distributing sectors (`split_fermion_sector`) are actively used.

*   **Interactions with other components:**
    *   This module is foundational for any ED solver component that needs to work with explicit many-body basis states.
    *   The `new_fermion_sector` procedure is the primary interface for other modules. When called, it ensures that the global `Hilbert` object (for the specified number of orbitals `norbs`) is created if it doesn't already exist. This `Hilbert` object pre-computes and stores all possible Fock states and their ranks.
    *   `fermion_sector_type` instances then merely point to the relevant segments of this global `Hilbert` data. This is memory efficient as the full state listing and rank array are stored only once per `norbs`.
    *   **Hamiltonian Construction (`h_class`):** Routines that build the Hamiltonian matrix would use `fermion_sector_type` to iterate over basis states of a given sector.
    *   **Operator Application (`apply_c`, `apply_ns`, `apply_p`):** These modules would use the `state` integer representation and the `rank` mapping from `fermion_sector_type` to apply operators and find the resulting states within the appropriate sectors.
    *   **ED Solvers (`ed_arpack`, `block_lanczos`):** The dimension of the sector (`SEC%dimen`) is passed to these solvers, and the resulting eigenvectors are represented in the basis defined by `SEC%state`.
    *   The `split_fermion_sector` routine, along with `chunk`, `istatemin`, `istatemax` members, facilitates the distribution of work (e.g., matrix-vector products) across MPI processes in parallel ED implementations.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
