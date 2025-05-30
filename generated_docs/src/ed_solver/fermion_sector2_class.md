# `src/ed_solver/fermion_sector2_class.f90`

## Overview

*   **Purpose:** This module defines a Fortran derived type, `fermion_sector2_type`, and associated procedures. It is designed to represent and manage fermionic Hilbert space sectors that are constructed as direct products of two individual `fermion_sector_type` components (from `fermion_hilbert_class.f90`). This structure is typically used for systems where spin-up and spin-down particle numbers (\(N_{up}\), \(N_{down}\)) are separately conserved quantum numbers.
*   **Role in Project:** Provides essential data structures and utility functions for working with many-body basis states within a product basis, such as \(|N_{up}, N_{down}\rangle = |N_{up}\rangle \otimes |N_{down}\rangle\). It facilitates the mapping between combined state indices (ranks) and the indices of states within the individual spin sectors. The module also handles the distribution of these product states across MPI processes for parallel computations.

## Key Components

*   **`TYPE fermion_sector2_type`**:
    Represents a Hilbert space sector formed by the direct product of a spin-up sector and a spin-down sector.
    *   `LOGICAL :: NAMBU = .false.`: A flag indicating if this product sector is to be interpreted in a Nambu (superconducting) context. Its direct usage within this module is minimal but is passed during construction.
    *   `LOGICAL :: not_physical = .false.`: Set to `.TRUE.` if the sector is unphysical (e.g., if \(N_{up} + N_{down} > 2 \times N_s\), or if constituent up/down sectors are unphysical).
    *   `TYPE(fermion_sector_type) :: up`: The `fermion_sector_type` object (from `fermion_hilbert_class.f90`) representing the spin-up subspace.
    *   `TYPE(fermion_sector_type) :: down`: The `fermion_sector_type` object representing the spin-down subspace.
    *   `CHARACTER(LEN=100) :: title = '\0'`: A descriptive title for the sector, typically formatted as "(nup, ndo) = (value1, value2)".
    *   `INTEGER :: norbs = 0`: Total number of single-particle orbitals in the product space (usually `up%norbs + down%norbs`, which is `2 * Ns` if `Ns` is the number of spatial sites).
    *   `INTEGER :: npart = 0`: Total number of particles in the sector (`up%npart + down%npart`).
    *   `INTEGER :: nstates = 0`: Total number of states in the full product space if there were no particle number constraints in sub-sectors (i.e., `up%nstates * down%nstates`).
    *   `INTEGER :: dimen = 0`: The actual dimension of this product sector (`up%dimen * down%dimen`).
    *   `INTEGER :: strideup = 0`, `stridedown = 0`: Strides used for indexing within the product basis. `strideup` is typically 1, and `stridedown` is `up%dimen`.
    *   `INTEGER, POINTER :: chunk(:) => NULL()`: An array storing the number of states in this product sector managed by each MPI process.
    *   `INTEGER, POINTER :: istatemin(:) => NULL()`: An array storing the starting global index of states for each MPI process's local chunk.
    *   `INTEGER, POINTER :: istatemax(:) => NULL()`: An array storing the ending global index of states for each MPI process's local chunk.

*   **`INTERFACE new_fermion_sector2`**: A generic interface for creating `fermion_sector2_type` objects.
    *   `MODULE PROCEDURE new_fermion_sector2_from_scratch(SEC, nup, ndo, Ns, NAMBU)`: The primary constructor. It creates a new `fermion_sector2_type` object `SEC`.
        *   It initializes `SEC%up` and `SEC%down` by calling `new_fermion_sector` (from `fermion_hilbert_class.f90`) with `nup` particles in `Ns` orbitals (for up) and `ndo` particles in `Ns` orbitals (for down).
        *   It then calculates the combined properties (`norbs`, `npart`, `dimen`, etc.) and calls `split_fermion_sector2` to set up MPI distribution.
    *   `MODULE PROCEDURE new_fermion_sector2_from_old(SECOUT, SECIN)`: Creates `SECOUT` as a deep copy of `SECIN`.

*   **`INTERFACE tabrankupdo`**: Generic interface for calculating arrays of ranks.
    *   `MODULE PROCEDURE rankupdo_strideup(rank, tabiup, ido, sector2)`: Calculates global ranks `rank(:)` for a fixed down-state index `ido` and an array of up-state indices `tabiup(:)`.
    *   `MODULE PROCEDURE rankupdo_stridedo(rank, iup, tabido, sector2)`: Calculates global ranks `rank(:)` for a fixed up-state index `iup` and an array of down-state indices `tabido(:)`.

*   **Public Management and Utility Subroutines & Functions:**
    *   **`SUBROUTINE copy_fermion_sector2(SECOUT, SECIN)`**: Performs a deep copy from `SECIN` to `SECOUT`.
    *   **`SUBROUTINE delete_fermion_sector2(sector2)`**: Destructor for `fermion_sector2_type` objects. It calls `delete_fermion_sector` for its `up` and `down` components and deallocates MPI-related arrays.
    *   **`FUNCTION is_in_fermion_sector2(state, sector2) AS LOGICAL`**: Checks if a combined integer `state` representation corresponds to a valid state within the product `sector2` (i.e., its up and down components are valid in their respective sub-sectors).
    *   **`FUNCTION rankupdo(iup, ido, sector2) AS INTEGER`**: Given `iup` (1-based rank in the up-sector `sector2%up%state` array) and `ido` (1-based rank in the down-sector `sector2%down%state` array), calculates the combined 1-based rank in the product basis using the formula: `iup + (ido - 1) * sector2%up%dimen`.
    *   **`FUNCTION rank2(state, sector2) AS INTEGER`**: Converts a combined integer `state` representation (where up and down parts are packed) into its 1-based rank within the `fermion_sector2_type` `sector2`. It first extracts the up and down integer states, finds their ranks in their respective sub-sectors, then uses `rankupdo`.
    *   **`FUNCTION stateupdo(iup, ido, sector2) AS INTEGER`**: Given `iup` (1-based index into `sector2%up%state`) and `ido` (1-based index into `sector2%down%state`), constructs the combined integer representation of the product state using the formula: `sector2%up%state(iup) + sector2%down%state(ido) * sector2%up%nstates`.
    *   **`FUNCTION state2(rank, sector2) AS INTEGER`**: The inverse of `rank2`. Converts a 1-based `rank` in the product sector `sector2` back into its combined integer state representation.
    *   **`SUBROUTINE read_raw_fermion_sector2(sector2, UNIT, NAMBU)` / `SUBROUTINE write_raw_fermion_sector2(sector2, UNIT)`**: For unformatted (binary) input/output of `fermion_sector2_type` definitions (primarily reads/writes the `npart` and `norbs` of constituent sectors).

*   **Internal Subroutines:**
    *   **`SUBROUTINE split_fermion_sector2(SEC)`**: Configures the MPI distribution for the product sector `SEC`. It sets the `chunk`, `istatemin`, and `istatemax` arrays for the product sector. The current implementation seems to make the `up` sector local to all processes (`SEC%up%chunk = SEC%up%dimen`) and distributes based on the `down` sector's distribution, then combines these for the product space.

## Important Variables/Constants

*   **Within `fermion_sector2_type`:**
    *   `up :: TYPE(fermion_sector_type)`: The `fermion_sector_type` object defining the spin-up subspace.
    *   `down :: TYPE(fermion_sector_type)`: The `fermion_sector_type` object defining the spin-down subspace.
    *   `dimen :: INTEGER`: The total dimension of the product Hilbert space sector, equal to `up%dimen * down%dimen`.
    *   `strideup :: INTEGER` (typically 1) and `stridedown :: INTEGER` (typically `up%dimen`): These are fundamental for mapping 2D (up_index, down_index) state indexing to a 1D global index in the product space.

## Usage Examples

```fortran
MODULE example_fermion_sector2_usage
  USE fermion_sector2_class
  ! Assuming fermion_hilbert_class is available for its types if we delve into up/down components
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: manage_product_basis_sector

CONTAINS
  SUBROUTINE manage_product_basis_sector
    TYPE(fermion_sector2_type) :: product_sector
    INTEGER, PARAMETER :: num_spatial_sites = 3
    INTEGER :: n_up = 1  ! Number of up-spin particles
    INTEGER :: n_down = 1 ! Number of down-spin particles
    INTEGER :: up_state_rank, down_state_rank, combined_rank_val
    INTEGER :: combined_int_state, recovered_product_rank

    ! Create a new product sector for (n_up=1, n_down=1) on Ns=3 spatial sites
    CALL new_fermion_sector2(product_sector, nup=n_up, ndo=n_down, Ns=num_spatial_sites)

    PRINT *, "Product Sector Initialized:"
    PRINT *, "  Title: ", TRIM(product_sector%title)
    PRINT *, "  Total dimension (dimen): ", product_sector%dimen
    PRINT *, "  Dimension of up-spin sector (product_sector%up%dimen): ", product_sector%up%dimen
    PRINT *, "  Dimension of down-spin sector (product_sector%down%dimen): ", product_sector%down%dimen

    ! Example: Get the rank of a state composed of the 1st up-state and 1st down-state
    up_state_rank = 1
    down_state_rank = 1
    combined_rank_val = rankupdo(up_state_rank, down_state_rank, product_sector)
    PRINT *, "Rank of product state (up_rank=", up_state_rank, ", down_rank=", down_state_rank, ") is: ", combined_rank_val

    ! Get the integer representation of this specific product state
    combined_int_state = stateupdo(up_state_rank, down_state_rank, product_sector)
    PRINT *, "Integer representation of this state: ", combined_int_state

    ! Convert the integer state back to its rank in the product sector
    recovered_product_rank = rank2(combined_int_state, product_sector)
    PRINT *, "Rank recovered from integer state ", combined_int_state, " is: ", recovered_product_rank

    ! Clean up
    CALL delete_fermion_sector2(product_sector)

  END SUBROUTINE manage_product_basis_sector
END MODULE example_fermion_sector2_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `fermion_hilbert_class`: This is a primary dependency. `fermion_sector2_type` is composed of two `fermion_sector_type` objects, and its constructor `new_fermion_sector2_from_scratch` calls `new_fermion_sector` from `fermion_hilbert_class` to initialize these components. Routines like `copy_fermion_sector` and `delete_fermion_sector` from `fermion_hilbert_class` are also used.
    *   `genvar`: For `log_unit` and `messages4` (used for debugging output).
    *   `common_def`: For string conversion utilities (`c2s`, `i2c`) and error handling (`create_seg_fault`).
    *   `quantum_algebra`: The comment mentions `fermion`, but it's not directly used in the public interface of this module.
    *   `mpi_mod`: The `split` function is used by the internal `split_fermion_sector2` routine, although the distribution logic in `split_fermion_sector2` appears to simplify how it uses the results from `split`.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules.

*   **Interactions with other components:**
    *   This module builds directly upon `fermion_hilbert_class` to define Hilbert space sectors suitable for systems where \(N_{up}\) and \(N_{down}\) are conserved (e.g., standard Hubbard models without spin-orbit coupling or pairing).
    *   **Sector Definition (`sector_class.f90`):** A higher-level `sector_type` object (from `sector_class.f90`) would likely contain an instance of `fermion_sector2_type` when the quantum numbers correspond to `(N_{up}, N_{down})`.
    *   **Operator Application (e.g., `apply_c.f90`):** Routines that apply operators in this product basis (like `apply_Cupdo` in `apply_c.f90`) heavily rely on the mapping functions provided by this module (`rankupdo`, `stateupdo`, `rank2`, `state2`) to correctly identify states and determine indices after an operator acts on either the up or down spin component of a state.
    *   **Hamiltonian Construction (`h_class.f90`):** When constructing the Hamiltonian matrix for a sector defined by `fermion_sector2_type`, these mapping functions are essential for indexing matrix elements correctly.
    *   **ED Solvers:** The `dimen` property and the MPI distribution information (`chunk`, `istatemin`, `istatemax`) are used by ED algorithms to manage memory and parallel operations.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
