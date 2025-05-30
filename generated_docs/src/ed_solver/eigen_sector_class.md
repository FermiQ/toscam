# `src/ed_solver/eigen_sector_class.f90`

## Overview

*   **Purpose:** This module defines Fortran derived types (`eigensector_type` and `eigensectorlist_type`) and a suite of associated procedures to manage and organize many-body eigenstates and their corresponding eigenvalues according to their quantum number sectors (e.g., fixed particle number, fixed total spin projection Sz).
*   **Role in Project:** It provides essential data structures for storing the results of Exact Diagonalization (ED) calculations, particularly when the Hamiltonian is block-diagonal with respect to certain conserved quantum numbers. An `eigensector_type` links a specific sector definition (from `sector_class.f90`) to the list of eigenpairs (from `eigen_class.f90`) found within that block. An `eigensectorlist_type` then manages a collection of these `eigensector_type` objects, often representing the low-lying energy spectrum of the quantum system across multiple relevant sectors.

## Key Components

*   **`TYPE eigensector_type`**:
    Represents the eigenpairs belonging to a single, specific quantum number sector.
    *   `type(sector_type) :: sector`: An object (from `sector_class.f90`) that defines the quantum numbers (e.g., particle number `npart`, total spin projection `Sz`, or specific up/down particle numbers `nup`, `ndo`) characterizing this sector.
    *   `type(eigenlist_type) :: lowest`: An object (from `eigen_class.f90`) that stores the list of all eigenpairs (eigenvalues and eigenvectors) computed for this particular `sector`.

*   **`TYPE eigensectorlist_type`**:
    Represents a collection of `eigensector_type` objects, typically covering all relevant sectors for a given problem.
    *   `INTEGER :: nsector = 0`: The number of distinct `eigensector_type` objects currently stored in the list.
    *   `type(eigensector_type), pointer :: es(:) => null()`: An allocatable array where each element is an `eigensector_type` object.

*   **`INTERFACE new_eigensector`**: A generic interface for creating `eigensector_type` objects.
    *   `MODULE PROCEDURE new_eigensector_from_scratch(es, sector_in)`: Creates a new `eigensector_type` object `es` associated with a given `sector_type` object `sector_in`. The `es%lowest` eigenlist is initialized (typically empty).
    *   `MODULE PROCEDURE copy_eigensector(esout, esin)`: Performs a deep copy of an `eigensector_type` object from `esin` to `esout`.

*   **Lifecycle and Management Subroutines:**
    *   **`SUBROUTINE new_eigensectorlist(list, nsector)`**: Initializes (or re-initializes) an `eigensectorlist_type` `list` to be able to hold `nsector` eigensectors. If `nsector > 0`, the `list%es` array is allocated.
    *   **`SUBROUTINE delete_eigensector(es)`**: Destructor for a single `eigensector_type` object `es`. It calls `delete_sector` for `es%sector` and `delete_eigenlist` for `es%lowest`.
    *   **`SUBROUTINE delete_eigensectorlist(list)`**: Destructor for an `eigensectorlist_type` `list`. It iterates through all stored `eigensector_type` objects, calling `delete_eigensector` for each, and then deallocates the `list%es` array.
    *   **`SUBROUTINE add_eigensector(es, list)`**: Adds a copy of an `eigensector_type` object `es` to the `eigensectorlist_type` `list`. The `list` is dynamically grown.
    *   **`SUBROUTINE remove_eigensector(sector, list)`**: Finds and removes the `eigensector_type` object from `list` whose `sector` component matches the provided `sector`.
    *   **`SUBROUTINE copy_eigensectorlist(listout, listin)`**: Performs a deep copy of an `eigensectorlist_type` from `listin` to `listout`.

*   **Analysis and Utility Subroutines & Functions:**
    *   **`SUBROUTINE filter_eigensector(list, window)`**: Iterates through each `eigensector_type` in the `list` and calls `filter_eigen` (from `eigen_class.f90`) on its `lowest` eigenlist to remove eigenpairs outside the specified energy `window`. If an eigensector becomes empty after filtering, it is removed from the `list`.
    *   **`FUNCTION gsenergy(list) AS REAL(DP)`**: Finds and returns the absolute lowest eigenvalue across all eigenpairs in all eigensectors within the `list`.
    *   **`SUBROUTINE nambu_energy_shift(list, back)`**: Applies (`back` is not present or false) or removes (`back` is true) global energy shifts (`energy_global_shift` and `energy_global_shift2` from `globalvar_ed_solver.f90`) to all eigenvalues in all eigensectors in the `list`. This is often used in Nambu formalism for superconductivity where the chemical potential reference might be handled via such shifts.
    *   **`FUNCTION partition(beta, list) AS REAL(DP)`**: Calculates the partition function \(Z = \sum_k \exp(-\beta (E_k - E_0))\) by summing the Boltzmann weights of all eigenpairs across all eigensectors in the `list`, relative to the global ground state energy `E0 = GSenergy(list)`.
    *   **`FUNCTION average_energy(beta, list) AS REAL(DP)`**: Calculates the thermal average energy \(\langle E \rangle = \frac{1}{Z} \sum_k (E_k - E_0) \exp(-\beta (E_k - E_0))\).
    *   **`FUNCTION nlowest(list) AS INTEGER`**: Returns the total count of individual eigenpairs stored across all eigensectors in the `list`.
    *   **`FUNCTION rank_sector_in_list(sector, list) AS INTEGER`**: Searches the `list` for an `eigensector_type` whose `sector` component matches the input `sector` (using `equal_sector` from `sector_class.f90`). Returns the 1-based index if found, otherwise 0.
    *   **`FUNCTION not_commensurate_sector(list, isector) AS LOGICAL` / `FUNCTION not_commensurate_sector_(sector) AS LOGICAL`**: Checks if a given sector's dimensions (specifically for `updo` basis type) are "commensurate," likely referring to whether its constituent up/down spin subspace dimensions are evenly divisible or align well with MPI process counts (`size2`). Returns true if not commensurate or if not an `updo` sector.

*   **I/O Subroutines:**
    *   **`SUBROUTINE print_eigensectorlist(list, UNIT)`**: Prints a summary of each eigensector in the `list` (sector title, dimension) followed by the eigenvalues within that sector (using `print_eigenlist`).
    *   **`SUBROUTINE write_raw_eigensectorlist(list, FILEOUT)` / `SUBROUTINE read_raw_eigensectorlist(list, FILEIN)`**: For unformatted (binary) input/output of an `eigensectorlist_type`, allowing for efficient saving and loading of the entire set of computed sectors and their eigenpairs.
    *   **Helper I/O (for reading specific vector formats directly into an eigenvector, not public):**
        *   `SUBROUTINE new_rcvec_from_file_r(unit, vec, sector)`: Reads real coefficients from `unit` and populates `vec`, mapping to states based on `sector`.
        *   `SUBROUTINE new_rcvec_from_file_c(unit, vec, sector)`: Reads complex coefficients similarly.

## Important Variables/Constants

*   **Within `eigensector_type`:**
    *   `sector :: type(sector_type)`: Defines the quantum numbers for this block of the Hilbert space.
    *   `lowest :: type(eigenlist_type)`: Contains all the eigenvalues and eigenvectors found for this specific `sector`.
*   **Within `eigensectorlist_type`:**
    *   `es(:) :: type(eigensector_type), pointer`: The array holding all the `eigensector_type` objects.
    *   `nsector :: INTEGER`: The number of `eigensector_type` objects in the `es` array.

## Usage Examples

```fortran
MODULE example_eigen_sector_usage
  USE eigen_sector_class
  USE sector_class, ONLY: sector_type, new_sector !, other_sector_initializers
  USE eigen_class, ONLY: eigen_type, add_eigen     !, other_eigen_initializers
  USE genvar, ONLY: DP

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: manage_ed_results_by_sector

CONTAINS
  SUBROUTINE manage_ed_results_by_sector
    TYPE(eigensectorlist_type) :: all_computed_sectors
    TYPE(eigensector_type)     :: current_sector_data
    TYPE(sector_type)          :: sector_def_A, sector_def_B
    REAL(DP)                   :: overall_ground_state_energy
    REAL(DP), PARAMETER        :: beta_temperature = 1.0_DP
    REAL(DP)                   :: system_partition_function
    REAL(DP)                   :: energy_window_filter(2)

    ! Initialize the list that will hold all eigensectors
    CALL new_eigensectorlist(all_computed_sectors, 0)

    ! --- Simulate obtaining results for Sector A ---
    ! (In reality, sector_def_A would be fully initialized based on N, Sz, etc.)
    CALL new_sector(sector_def_A) ! Simplified initialization
    CALL new_eigensector(current_sector_data, sector_def_A)
    ! (Assume diagonalization was done, and current_sector_data%lowest is populated)
    ! Example: CALL add_eigen(some_eigenpair_A, current_sector_data%lowest)
    CALL add_eigensector(current_sector_data, all_computed_sectors)
    CALL delete_eigensector(current_sector_data) ! Clean up temporary holder

    ! --- Simulate obtaining results for Sector B ---
    CALL new_sector(sector_def_B) ! Simplified initialization
    CALL new_eigensector(current_sector_data, sector_def_B)
    ! (Assume diagonalization for sector_def_B populates current_sector_data%lowest)
    CALL add_eigensector(current_sector_data, all_computed_sectors)
    CALL delete_eigensector(current_sector_data)

    ! Analyze the collected results
    IF (all_computed_sectors%nsector > 0) THEN
      overall_ground_state_energy = GSenergy(all_computed_sectors)
      system_partition_function = partition(beta_temperature, all_computed_sectors)

      PRINT *, "Overall Ground State Energy: ", overall_ground_state_energy
      PRINT *, "System Partition Function (Z): ", system_partition_function

      CALL print_eigensectorlist(all_computed_sectors)

      ! Filter results to an energy window relative to ground state
      energy_window_filter = [overall_ground_state_energy, overall_ground_state_energy + 5.0_DP]
      CALL filter_eigensector(all_computed_sectors, energy_window_filter)
      PRINT *, "After filtering, number of sectors: ", all_computed_sectors%nsector
      CALL print_eigensectorlist(all_computed_sectors)
    END IF

    ! Clean up
    CALL delete_eigensectorlist(all_computed_sectors)
    ! CALL delete_sector(sector_def_A) ! If it was dynamically allocated and not part of an eigensector
    ! CALL delete_sector(sector_def_B)

  END SUBROUTINE manage_ed_results_by_sector
END MODULE example_eigen_sector_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `eigen_class`: This is a fundamental dependency, as `eigensector_type` contains an `eigenlist_type` object. Many procedures in `eigen_sector_class` (e.g., `filter_eigensector`, `GSenergy`, `partition`) call corresponding routines in `eigen_class` to operate on the `lowest` component.
    *   `genvar`: For `DP` (double precision kind), `log_unit` (for logging/printing), `messages3`, `messages4` (for debug messages), `huge_` (large number constant), `size2` (number of MPI processes, used in `not_commensurate_sector`), `strongstop`.
    *   `sector_class`: Essential for the `sector_type` component of `eigensector_type`. Utilities from `sector_class` like `copy_sector`, `delete_sector`, `equal_sector`, `title_sector`, `dimen_func`, `rank_func` are used.
    *   `common_def`: For `utils_assert` (error checking), string utilities (`c2s`, `i2c`), messaging (`dump_message`), and file operations (`open_safe`, `close_safe`).
    *   `globalvar_ed_solver`: The `nambu_energy_shift` routine uses `energy_global_shift` and `energy_global_shift2` from this module.
    *   `stringmanip` (used in `filter_eigensector`): For `toString`.
    *   `string5` (used in `write_raw_eigensectorlist`): For `get_unit`.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules. No explicit third-party external libraries are directly used by this module.

*   **Interactions with other components:**
    *   **ED Solvers (e.g., `ed_arpack.f90`, `block_lanczos.f90`):** These routines diagonalize Hamiltonians for specific sectors. The results (eigenvalues/eigenvectors for a sector) would be stored by populating the `lowest` (`eigenlist_type`) component of an `eigensector_type` object.
    *   **Main Solver Logic:** Higher-level routines controlling the ED calculation would manage an `eigensectorlist_type` object, adding solutions from different quantum number sectors as they are computed.
    *   **Post-processing/Analysis Modules (e.g., `correlations.f90`, `density_matrix.f90`):** These modules would take an `eigensectorlist_type` object as input to access the complete set of relevant low-lying states and energies for calculating physical observables, Green's functions, density matrices, etc.
    *   The routines for calculating thermodynamic quantities like `partition` and `average_energy` directly use the stored eigenenergies.
    *   The `filter_eigensector` routine allows for focusing on states within a particular energy window relative to the ground state, which is often necessary for practical calculations.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
