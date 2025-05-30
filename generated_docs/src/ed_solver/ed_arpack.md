# `src/ed_solver/ed_arpack.f90`

## Overview

*   **Purpose:** This module provides Fortran subroutines for matrix diagonalization, offering interfaces to both the ARPACK library (for finding a few eigenvalues of large, sparse matrices using iterative methods) and a full, direct diagonalization method (suitable for smaller matrices where all eigenpairs can be computed).
*   **Role in Project:** Serves as a numerical eigensolver module within an Exact Diagonalization (ED) framework. It allows the user or higher-level routines to choose the most appropriate diagonalization strategy based on the system size and the number of required eigenstates. It abstracts the details of interacting with ARPACK or constructing and diagonalizing the full matrix.

## Key Components

*   **`MODULE ed_arpack`**:
    The main module containing the diagonalization routines.

*   **Public Subroutines:**
    *   **`SUBROUTINE ARPACK_diago(lowest)`**:
        Performs diagonalization using the ARPACK library to find a specified number of the lowest eigenvalues and eigenvectors.
        *   **Argument:**
            *   `lowest (TYPE(eigenlist_type), INTENT(INOUT))`: An `eigenlist_type` object (from `eigen_class.f90`) that will be populated with the computed lowest eigenpairs.
        *   **Functionality:**
            1.  Retrieves the Hamiltonian dimension (`dimenvec`) from `dimen_H()` (from `h_class.f90`).
            2.  Determines the number of eigenpairs to compute (`Neigen_`), ensuring it's less than `dimenvec`. The number requested might differ slightly for real vs. complex arithmetic due to ARPACK's handling of complex conjugate pairs or specific modes.
            3.  Allocates arrays `VALP` for eigenvalues and `VECP` for eigenvectors.
            4.  If `Neigen_ > 0`, it calls an external ARPACK wrapper routine:
                *   `arpack_eigenvector_sym_matrix_` if `_complex` is defined (for complex Hamiltonian).
                *   `arpack_eigenvector_sym_matrix` if `_complex` is not defined (for real Hamiltonian).
                These wrappers interface with ARPACK's reverse communication routines. They are passed the problem dimension, tolerance, number of eigenvalues requested, and crucially, the matrix-vector product routines (`Hmultc` or `Hmultr` from `h_class.f90`). The `SR` (Smallest Real part) or `BE` (Both Ends, then filtered) options are typically used to find low-lying states.
            5.  After the ARPACK call returns `nconv` converged eigenpairs, it filters these results to keep only those within `dEmax` (an energy window from `globalvar_ed_solver`) of the lowest found eigenvalue `VALP(1)`.
            6.  The selected, converged eigenpairs (eigenvalues and corresponding eigenvectors) are then added to the `lowest` list using `add_eigen` (from `eigen_class.f90`).
            7.  Timers record the duration of the ARPACK computation and data storage.
            8.  Deallocates temporary arrays.
        *   **Note:** This routine depends on the `_arpack` preprocessor macro being defined at compile time to enable the calls to the ARPACK library.

    *   **`SUBROUTINE ED_diago(lowest)`**:
        Performs a full, direct diagonalization of the Hamiltonian matrix. This involves explicitly constructing the matrix in memory.
        *   **Argument:**
            *   `lowest (TYPE(eigenlist_type), INTENT(INOUT))`: An `eigenlist_type` object that will be populated with the computed eigenpairs.
        *   **Functionality:**
            1.  Retrieves the Hamiltonian dimension (`dimenvec`) from `dimen_H()`.
            2.  Allocates arrays `VECP(dimenvec, dimenvec)` to store the full Hamiltonian matrix and later its eigenvectors, and `VALP(dimenvec)` for eigenvalues.
            3.  **Constructs the full Hamiltonian matrix `VECP`**: It iterates `ii` from 1 to `dimenvec`. In each iteration `ii`, it prepares a unit vector `invec` (1 at position `ii`, 0 elsewhere), computes `outvec = H * invec` by calling `Hmult__` (from `h_class.f90`), and stores `outvec%rc` as the `ii`-th column of `VECP`.
            4.  Calls `eigenvector_matrix` (a routine from `matrix_class.f90` or similar, for dense matrix diagonalization like LAPACK's `dsyev` or `zheev`) to diagonalize the constructed matrix `VECP`. Eigenvalues are stored in `VALP`, and eigenvectors overwrite `VECP`.
            5.  Filters the resulting eigenpairs:
                *   If `FLAG_FULL_ED_GREEN` (from `globalvar_ed_solver`) is true, all `dimenvec` eigenpairs are kept.
                *   Otherwise, only eigenpairs within `dEmax` of the ground state energy `VALP(1)` are kept.
            6.  The selected eigenpairs are added to the `lowest` list.
            7.  Timers record the duration.
            8.  Deallocates temporary arrays.

## Important Variables/Constants

*   **`lowest :: TYPE(eigenlist_type)`**: The primary output argument for both diagonalization routines, used to store the list of computed eigenvalues and eigenvectors.
*   **Control Parameters (typically from `globalvar_ed_solver.f90`):**
    *   `neigen :: INTEGER`: The desired number of eigenpairs to compute (especially relevant for `ARPACK_diago`).
    *   `tolerance :: REAL(DP)`: The convergence tolerance for the eigenvalues in ARPACK.
    *   `dEmax :: REAL(DP)`: An energy window above the ground state. Eigenpairs outside this window might be discarded.
    *   `FLAG_FULL_ED_GREEN :: LOGICAL`: If true, `ED_diago` attempts to store all computed eigenpairs.
*   **External Routines for Hamiltonian Operations (from `h_class.f90`):**
    *   `dimen_H() :: INTEGER`: Function returning the dimension of the Hamiltonian matrix.
    *   `Hmultc` / `Hmultr` / `Hmult__`: Subroutines that perform the matrix-vector product \(y = Hx\).
*   **External ARPACK Wrappers (if `_arpack` is defined):**
    *   `arpack_eigenvector_sym_matrix_` (complex) / `arpack_eigenvector_sym_matrix` (real): These are assumed to be interfaces to the ARPACK library routines.
*   **External Dense Diagonalizer (for `ED_diago`):**
    *   `eigenvector_matrix` (from `matrix` module): A routine for full diagonalization of a dense matrix.

## Usage Examples

```fortran
MODULE example_ed_arpack_usage
  USE ed_arpack
  USE eigen_class, ONLY: eigenlist_type, new_eigen_list, delete_eigen_list, print_eigenlist
  USE h_class, ONLY: init_H ! To make Hmult__ and dimen_H() available
  USE globalvar_ed_solver, ONLY: Neigen, tolerance, dEmax, FLAG_FULL_ED_GREEN
  ! Other modules for setting up the Hamiltonian and global parameters

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: perform_diagonalization

CONTAINS
  SUBROUTINE perform_diagonalization(use_arpack)
    LOGICAL, INTENT(IN) :: use_arpack
    TYPE(eigenlist_type) :: calculated_states

    ! --- Setup Hamiltonian (e.g., CALL init_H(...)) ---
    ! --- Set global parameters like Neigen, tolerance, dEmax, FLAG_FULL_ED_GREEN ---
    ! Example:
    ! Neigen = 10
    ! tolerance = 1.0E-9
    ! dEmax = 1.0E10 ! Effectively keep all converged states by ARPACK or up to Neigen
    ! FLAG_FULL_ED_GREEN = .FALSE. ! For ED_diago, only keep states within dEmax of GS

    CALL new_eigen_list(calculated_states)

    IF (use_arpack) THEN
      PRINT *, "Diagonalizing using ARPACK..."
      CALL ARPACK_diago(calculated_states)
    ELSE
      PRINT *, "Performing full diagonalization (ED_diago)..."
      CALL ED_diago(calculated_states)
    END IF

    PRINT *, "Diagonalization complete. Number of eigenpairs found: ", calculated_states%neigen
    CALL print_eigenlist(calculated_states)

    CALL delete_eigen_list(calculated_states)
    ! --- Cleanup Hamiltonian ---

  END SUBROUTINE perform_diagonalization
END MODULE example_ed_arpack_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `genvar`: For `DP` (double precision kind) and `rank`.
    *   `common_def`: For timing utilities (`reset_timer`, `timer_fortran`).
    *   `eigen_class`: Provides `eigenlist_type` for storing results and `add_eigen`, `is_eigen_in_window` utilities.
    *   `globalvar_ed_solver`: Supplies global control parameters for diagonalization (`dEmax`, `neigen`, `tolerance`, `FLAG_FULL_ED_GREEN`).
    *   `h_class`: Essential for defining the Hamiltonian problem. It provides `dimen_H()` for matrix size and `Hmultc`/`Hmultr`/`Hmult__` for the matrix-vector product operation required by iterative solvers like Lanczos/ARPACK and for constructing the matrix in `ED_diago`.
    *   `rcvector_class` (used in `ED_diago`): For `rcvector_type` (real/complex vector) and its associated memory management routines.
    *   `matrix` (used in `ED_diago`): Provides the `eigenvector_matrix` routine for dense matrix diagonalization.
    *   `timer_mod` (used in `ED_diago`): For `start_timer`, `stop_timer`.

*   **External Libraries:**
    *   **ARPACK**: The `ARPACK_diago` subroutine is a wrapper around ARPACK library routines. This dependency is typically enabled by a preprocessor flag like `_arpack`. ARPACK is a software library for solving large-scale eigenvalue problems, often used when the matrix is sparse and only a subset of eigenvalues/eigenvectors are needed.
    *   **MPI (Message Passing Interface)**: Although not explicitly shown in `ARPACK_diago`'s direct MPI calls in this snippet, ARPACK itself can be built with MPI support, and the overall framework (`h_class` matrix-vector products) might be MPI parallel. `ED_diago` as shown is serial in its matrix construction and diagonalization phase, but the `Hmult__` it calls could be parallel.

*   **Interactions with other components:**
    *   This module provides the core numerical routines for solving the eigenvalue problem \(H|\psi\rangle = E|\psi\rangle\).
    *   `ARPACK_diago` is suited for large systems where constructing the full Hamiltonian matrix is infeasible. It relies on the "matrix-free" approach where only the action of \(H\) on a vector is needed (provided by `h_class`).
    *   `ED_diago` is suitable for smaller systems. It explicitly builds the Hamiltonian matrix in memory by repeatedly applying \(H\) to basis vectors (effectively \(H_{ij} = \langle i | H | j \rangle\)) and then calls a standard dense diagonalization routine.
    *   Both routines populate an `eigenlist_type` object with the computed eigenvalues and eigenvectors, which are then used by other modules for calculating physical quantities (e.g., by `correlations.f90` or `density_matrix.f90`).
    *   The choice between `ARPACK_diago` and `ED_diago` would typically be made by a higher-level control routine based on system size or user input.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
