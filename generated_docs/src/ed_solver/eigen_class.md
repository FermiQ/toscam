# `src/ed_solver/eigen_class.f90`

## Overview

*   **Purpose:** This module defines Fortran derived types (`eigen_type` and `eigenlist_type`) and a suite of associated procedures for robustly storing, managing, and manipulating individual eigenpairs (consisting of an eigenvalue and its corresponding eigenvector) as well as lists of these eigenpairs.
*   **Role in Project:** Provides the fundamental data structures for handling the results obtained from various eigenvalue calculation routines (such as those found in `ed_arpack.f90` or `block_lanczos.f90`). It enables the creation, deletion, copying, and manipulation (like filtering, sorting by rank, and I/O) of eigenvalues and their eigenvectors, which are essential for subsequent physical property calculations.

## Key Components

*   **`TYPE eigen_type`**:
    Represents a single eigenpair.
    *   `INTEGER :: rank = 0`: An identifier for the eigenpair, often corresponding to its order (e.g., 0 for ground state, 1 for first excited, etc., though usage might vary).
    *   `LOGICAL :: converged = .false.`: A flag indicating whether the iterative diagonalization procedure for this eigenpair has converged to the desired tolerance.
    *   `REAL(DP) :: val = 0.0_DP`: Stores the real-valued eigenvalue. (Assumes Hermitian matrices).
    *   `TYPE(rcvector_type) :: vec`: An object of `rcvector_type` (from `rcvector_class.f90`) that stores the eigenvector components (can be real or complex).
    *   `REAL(DP) :: lanczos_vecp(1000)`: Seems to be a fixed-size array potentially for storing Lanczos vectors during calculation or some projection, its specific use might be tied to a particular Lanczos implementation.
    *   `REAL(DP) :: rdist`: Possibly stores the residual norm or a distance metric related to the convergence of the eigenpair.
    *   `INTEGER :: lanczos_iter`: Number of Lanczos iterations taken to converge this eigenpair.
    *   `INTEGER :: dim_space`: The full dimension of the Hilbert space in which the eigenvector lives.
    *   `INTEGER :: dim_sector`: The dimension of the specific symmetry sector this eigenstate belongs to (if applicable).

*   **`TYPE eigenlist_type`**:
    Represents a list (collection) of eigenpairs.
    *   `INTEGER :: neigen = 0`: The current number of eigenpairs stored in the list.
    *   `TYPE(eigen_type), POINTER :: eigen(:) => NULL()`: An allocatable array of `eigen_type` objects.

*   **`INTERFACE new_eigen`**: A generic interface for creating `eigen_type` objects.
    *   `MODULE PROCEDURE allocate_new_eigen(eigen, n)`: Allocates the eigenvector `eigen%vec` to size `n`.
    *   `MODULE PROCEDURE new_eigen_from_scratch(eigen, val, vec, CONVERGED, RANK, no_vector)`: Initializes an `eigen_type` object with a given eigenvalue `val`, eigenvector `vec` (of `rcvector_type`), and optional convergence status and rank. `no_vector` can skip vector copy.
    *   `MODULE PROCEDURE copy_eigen(EIGENOUT, EIGENIN)`: Performs a deep copy of an `eigen_type` object from `EIGENIN` to `EIGENOUT`.

*   **`INTERFACE add_eigen`**: A generic interface for adding eigenpairs to an `eigenlist_type`.
    *   `MODULE PROCEDURE add_eigen_(eigen, list)`: Adds a single, fully prepared `eigen_type` object to the `list`. The list size is dynamically increased.
    *   `MODULE PROCEDURE add_eigen__(n, VECP, VALP, list)`: Adds `n` eigenpairs to the `list`. `VECP` is a 2D array where columns are eigenvectors, and `VALP` is a 1D array of corresponding eigenvalues. This is a convenience routine for bulk addition from raw diagonalization output.

*   **Lifecycle and Management Subroutines:**
    *   **`SUBROUTINE new_eigenlist(list, neigen)`**: Initializes or reinitializes an `eigenlist_type` to hold a specified `neigen` number of eigenpairs. If `neigen > 0`, the `list%eigen` array is allocated.
    *   **`SUBROUTINE delete_eigen(eigen)`**: Destructor for a single `eigen_type` object. It deallocates the `eigen%vec` component and resets other members.
    *   **`SUBROUTINE delete_eigenlist(list)`**: Destructor for an `eigenlist_type`. It iterates through all `eigen_type` objects in the list, calling `delete_eigen` for each, and then deallocates the `list%eigen` array itself.
    *   **`SUBROUTINE eigen_allocate_vec(eigen, vec, N)`**: Allocates `eigen%vec` either by copying an existing `rcvector_type` `vec` or by creating a new one of size `N`.
    *   **`SUBROUTINE remove_eigen(rank, list)`**: Removes the eigenpair identified by `rank` from the `list`.
    *   **`SUBROUTINE copy_eigenlist(listout, listin)`**: Performs a deep copy of an `eigenlist_type` from `listin` to `listout`.

*   **Utility Functions and Subroutines:**
    *   **`FUNCTION is_eigen_in_window(E_, window) AS LOGICAL`**: Returns `.TRUE.` if the eigenvalue `E_` falls within the energy range defined by `window(1)` and `window(2)`.
    *   **`SUBROUTINE filter_eigen(list, window)`**: Modifies an `eigenlist_type` `list` by removing all eigenpairs whose eigenvalues fall outside the specified energy `window`.
    *   **`FUNCTION min_eigen(list) AS REAL(DP)`**: Returns the minimum eigenvalue present in the `list`.
    *   **`FUNCTION sum_boltz(list, beta, E0) AS REAL(DP)`**: Calculates the sum of Boltzmann weights, \(\sum_i \exp(-\beta(E_i - E_0))\), for all eigenpairs in the `list`. Filters based on `dEmax0` if `FLAG_FULL_ED_GREEN` is true.
    *   **`FUNCTION sum_boltz_times_E(list, beta, E0) AS REAL(DP)`**: Calculates the sum \(\sum_i (E_i - E_0) \exp(-\beta(E_i - E_0))\). Also filters based on `dEmax0`.
    *   **`SUBROUTINE print_eigenlist(list, UNIT)`**: Prints the eigenvalues (and their ranks) from the `list` to a specified Fortran `UNIT` or the default log unit.
    *   **`SUBROUTINE write_raw_eigenlist(list, UNIT)` / `SUBROUTINE read_raw_eigenlist(list, UNIT)`**: For unformatted (binary) input/output of an `eigenlist_type`, allowing for efficient saving and loading of diagonalization results.
    *   **`FUNCTION rank_eigen_in_list(rank, list) AS INTEGER`**: Finds the array index (position) within `list%eigen(:)` that corresponds to the eigenpair with the identifier `rank`.

## Important Variables/Constants

*   **Within `eigen_type`:**
    *   `val :: REAL(DP)`: Stores the eigenvalue.
    *   `vec :: TYPE(rcvector_type)`: Stores the eigenvector. This relies on `rcvector_class` for its definition and can hold real or complex data.
    *   `converged :: LOGICAL`: Indicates if the eigenvalue calculation met the convergence criteria.
    *   `rank :: INTEGER`: A unique identifier or label for the eigenpair.
*   **Within `eigenlist_type`:**
    *   `eigen(:) :: TYPE(eigen_type), POINTER`: The allocatable array holding the collection of `eigen_type` objects.
    *   `neigen :: INTEGER`: The number of eigenpairs currently stored in the `eigen` array.

## Usage Examples

```fortran
MODULE example_eigen_class_usage
  USE eigen_class
  USE rcvector_class, ONLY: new_rcvector, delete_rcvector ! For eigenvector setup
  USE genvar, ONLY: DP
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: manage_eigen_results

CONTAINS
  SUBROUTINE manage_eigen_results
    TYPE(eigenlist_type) :: results_list
    TYPE(eigen_type)     :: temp_eigenpair
    TYPE(rcvector_type)  :: temp_vector
    REAL(DP)             :: example_eigenvalues(2)
    COMPLEX(DP)          :: example_eigenvectors(5, 2) ! Assuming dim=5, 2 vectors
    INTEGER              :: i, problem_dimension
    REAL(DP)             :: energy_window(2)
    REAL(DP)             :: partition_func, avg_energy_term

    problem_dimension = 5
    example_eigenvalues = [0.5_DP, 1.5_DP]
    ! ... (Code to fill example_eigenvectors) ...

    ! 1. Initialize an eigenlist
    CALL new_eigenlist(results_list, 0) ! Start empty

    ! 2. Add eigenpairs (e.g., from a fictional diagonalization)
    ! Using add_eigen__ for bulk addition
    CALL add_eigen__(2, example_eigenvectors, example_eigenvalues, results_list)

    ! Or add one by one
    CALL new_rcvector(temp_vector, problem_dimension)
    temp_vector%rc = [(0.1_DP*i, i=1,problem_dimension)] ! Dummy vector
    CALL new_eigen_from_scratch(temp_eigenpair, 2.5_DP, temp_vector, CONVERGED=.TRUE., RANK=3)
    CALL add_eigen(temp_eigenpair, results_list)
    CALL delete_eigen(temp_eigenpair) ! Clean up temp_eigenpair as it's copied

    PRINT *, "Number of eigenpairs in list: ", results_list%neigen
    CALL print_eigenlist(results_list)

    ! 3. Filter eigenpairs within an energy window
    energy_window = [0.0_DP, 2.0_DP]
    CALL filter_eigen(results_list, energy_window)
    PRINT *, "After filtering (0.0 to 2.0), number of eigenpairs: ", results_list%neigen
    CALL print_eigenlist(results_list)

    ! 4. Calculate Boltzmann sums (assuming E0 = 0.0, beta = 1.0)
    IF (results_list%neigen > 0) THEN
      partition_func = sum_boltz(results_list, 1.0_DP, 0.0_DP)
      avg_energy_term = sum_boltz_times_E(results_list, 1.0_DP, 0.0_DP)
      PRINT *, "Sum of Boltzmann weights (Z): ", partition_func
      PRINT *, "Sum of E*Boltzmann weights: ", avg_energy_term
    END IF

    ! 5. Clean up
    CALL delete_eigenlist(results_list)
    CALL delete_rcvector(temp_vector)

  END SUBROUTINE manage_eigen_results
END MODULE example_eigen_class_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `genvar`: Provides `DP` (double precision kind), `log_unit` (for printing/logging), and `huge_` (a large number constant used in `min_eigen`).
    *   `rcvector_class`: Essential for the `rcvector_type` used within `eigen_type` to store the actual eigenvector data. Procedures like `new_rcvector` and `delete_rcvector` from this class are used.
    *   `common_def` (used in `remove_eigen` and `filter_eigen`): For string conversion utilities `c2s` (char to string) and `i2c` (integer to char array).
    *   `globalvar_ed_solver` (used in `sum_boltz`, `sum_boltz_times_E`): For flags `dEmax0` and `FLAG_FULL_ED_GREEN` which control filtering in Boltzmann sums.
    *   `linalg` (used in `sum_boltz`, `sum_boltz_times_E`): For `dexpc` (complex exponential function, likely `exp` for real DP in this context).

*   **External Libraries:**
    *   Standard Fortran intrinsic modules. No explicit third-party external libraries are directly used by this module.

*   **Interactions with other components:**
    *   **Diagonalization Routines (e.g., `ed_arpack.f90`, `block_lanczos.f90`):** These modules are primary *producers* of `eigenlist_type` objects. They call `add_eigen` or `add_eigen__` to store the computed eigenvalues and eigenvectors.
    *   **Physics Calculation Modules (e.g., `correlations.f90`, `density_matrix.f90`):** These modules are primary *consumers* of `eigenlist_type` objects. They access the `eigen%val` (eigenvalues) and `eigen%vec%rc` (eigenvector components) to calculate various physical quantities, expectation values, Lehmann sums for dynamic correlation functions, etc.
    *   The `filter_eigen` routine allows higher-level modules to select a subset of eigenpairs based on an energy window, which can be useful for focusing on low-energy physics or specific spectral regions.
    *   The Boltzmann sum functions (`sum_boltz`, `sum_boltz_times_E`) are utilities for calculating thermodynamic averages or partition functions.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
