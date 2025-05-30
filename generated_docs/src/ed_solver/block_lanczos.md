# `src/ed_solver/block_lanczos.f90`

## Overview

*   **Purpose:** This module implements the Block Lanczos algorithm, a numerical method used for finding several of the lowest (or highest) eigenvalues and their corresponding eigenvectors of a large, sparse, symmetric Hamiltonian matrix.
*   **Role in Project:** Serves as a parallel numerical diagonalization engine within a larger Exact Diagonalization (ED) solver framework. It interfaces with an external Block Lanczos library (BLZPACK, specifically the `BLZDRD` routine) and manages the iterative process. This includes orchestrating the matrix-vector multiplications (delegated to `Hmult__` from `h_class`) and handling MPI communication for distributed vectors in a parallel computing environment.

## Key Components

*   **`MODULE Block_Lanczos`**:
    The main module containing the Block Lanczos diagonalization subroutine.

*   **`SUBROUTINE Block_Lanczos_diagonalize(lowest)`**:
    The primary public subroutine that performs the Block Lanczos diagonalization.
    *   **Argument:**
        *   `lowest (TYPE(eigenlist_type), INTENT(INOUT))`: An `eigenlist_type` object (likely defined in `eigen_class.f90`) where the computed lowest eigenvalues and eigenvectors will be stored.
    *   **Functionality:**
        1.  **Initialization:**
            *   Retrieves problem and control parameters from other modules/global variables:
                *   Hamiltonian dimension from `dimen_H()` (from `h_class`).
                *   Number of desired eigenpairs (`Neigen_`), block size (`Block_size`), maximum Lanczos iterations (`Nitermax`), and convergence `tolerance` (from `globalvar_ed_solver`).
            *   Determines local data distribution for MPI using `split` (from `mpi_mod`).
            *   Sets up query parameters (`ISTOR_QUERY`, `RSTOR_QUERY`) for the `BLZDRD` routine.
        2.  **BLZPACK Workspace Query (Optional):**
            *   If `Block_size` is initially 0, it makes a query call to `BLZDRD` to determine optimal workspace array sizes (`LISTOR`, `LRSTOR`) and potentially an optimal `Block_size`.
        3.  **Memory Allocation:**
            *   Allocates memory for BLZPACK workspace arrays (`ISTOR`, `RSTOR`), Lanczos block vectors (`U`, `V`), and arrays to store resulting eigenvalues and residuals (`VALP`) and eigenvectors (`VECP`).
        4.  **Iterative Lanczos Process (Loop with `BLZDRD` calls):**
            *   The core of the algorithm involves repeated calls to the external `BLZDRD` routine.
            *   `LFLAG` (returned by `BLZDRD`) controls the loop:
                *   `LFLAG = 0`: Convergence achieved or maximum iterations reached; exit loop.
                *   `LFLAG < 0`: Error condition; stop.
                *   `LFLAG = 1`: `BLZDRD` requires the matrix-vector product \(V = H \cdot U\). The subroutine performs this:
                    *   Gathers the distributed block of Lanczos vectors `U` from all MPI processes into a full vector `vec_in` using `MPI_ALLGATHERV`.
                    *   Calls `Hmult__` (provided by `h_class`) to compute `vec_out = H \cdot vec_in`.
                    *   Scatters the relevant local portion of `vec_out` back to the `V` array on each MPI process.
                *   `LFLAG = 4`: `BLZDRD` requests initial starting vectors. These are provided by filling the `V` array with random numbers using `randomize_mat`.
        5.  **Result Extraction and Storage:**
            *   After the `BLZDRD` loop terminates successfully, retrieves the number of converged eigenpairs (`Nconverged`) and Lanczos iterations (`Niter`).
            *   Iterates through the converged eigenpairs stored in `VALP` (eigenvalues, residuals) and `VECP` (eigenvector components).
            *   Filters eigenpairs based on `dEmax` (maximum energy above the absolute minimum found).
            *   For each valid eigenpair:
                *   Creates a new `eigen_type` object.
                *   Stores the eigenvalue and convergence status.
                *   Gathers the full eigenvector from distributed `VECP` components across all MPI processes into `eigen%vec%rc` using `MPI_ALLGATHERV`.
                *   Adds the new `eigen_type` object to the output `lowest` list using `add_eigen`.
        6.  **Logging and Cleanup:**
            *   Prints summary information, warnings if convergence issues arise (e.g., fewer eigenpairs converged than requested).
            *   Deallocates all workspace arrays.
            *   Records timing information for the diagonalization process.

## Important Variables/Constants

*   **Input/Output:**
    *   `lowest :: TYPE(eigenlist_type)`: The primary output, an object that accumulates the computed low-lying eigenvalues and their corresponding eigenvectors.
*   **Key Control Parameters (typically from `globalvar_ed_solver` or derived):**
    *   `Block_size :: INTEGER`: The block size used in the Block Lanczos algorithm. A block size of 0 allows BLZPACK to choose one.
    *   `Nitermax :: INTEGER`: Maximum number of Block Lanczos iterations allowed.
    *   `Neigen_ :: INTEGER`: The target number of eigenvalues/eigenvectors to compute.
    *   `tolerance :: REAL(DP)`: The desired convergence tolerance for the eigenvalues.
    *   `dimen_H() :: INTEGER`: Function returning the total dimension of the Hamiltonian matrix (from `h_class`).
*   **BLZPACK Workspace and Result Arrays:**
    *   `ISTOR :: INTEGER, ALLOCATABLE(:)`: Integer workspace array for `BLZDRD`.
    *   `RSTOR :: REAL(DP), ALLOCATABLE(:)`: Real workspace array for `BLZDRD`.
    *   `U(:,:) :: REAL(DP), ALLOCATABLE`: Stores the block of Lanczos vectors.
    *   `V(:,:) :: REAL(DP), ALLOCATABLE`: Stores the result of \(H \cdot U\) or starting vectors.
    *   `VALP(:,:) :: REAL(DP), ALLOCATABLE`: Stores computed eigenvalues and residual norms.
    *   `VECP(:,:) :: REAL(DP), ALLOCATABLE`: Stores computed eigenvector components (distributed).

## Usage Examples

```fortran
MODULE example_block_lanczos_usage
  USE Block_Lanczos
  USE eigen_class, ONLY: eigenlist_type, new_eigen_list, delete_eigen_list, print_eigenlist
  USE h_class, ONLY: init_H ! Assuming H needs initialization for Hmult__ and dimen_H
  USE globalvar_ed_solver, ONLY: Neigen, Block_size, Nitermax, tolerance, demax ! Set these params
  USE mpi_mod, ONLY: init_mpi, finalize_mpi ! If running with MPI

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: run_diagonalization

CONTAINS
  SUBROUTINE run_diagonalization
    TYPE(eigenlist_type) :: ground_states_list

    CALL init_mpi() ! Initialize MPI environment if needed

    ! Initialize Hamiltonian parameters (specific to the problem being solved)
    ! This would make Hmult__ and dimen_H() from h_class ready for use.
    ! CALL init_H(...)

    ! Set global parameters for Block Lanczos (can be read from input or set directly)
    Neigen = 5      ! Request 5 lowest eigenpairs
    Block_size = 0  ! Let BLZPACK determine block size (or set a specific value > 0)
    Nitermax = 1000 ! Max iterations
    tolerance = 1.0E-8 ! Convergence tolerance
    demax = 1.0E10  ! Energy window above ground state to keep eigenpairs

    CALL new_eigen_list(ground_states_list) ! Initialize the list to store results

    PRINT *, "Starting Block Lanczos diagonalization..."
    CALL Block_Lanczos_diagonalize(ground_states_list)
    PRINT *, "Block Lanczos diagonalization finished."

    ! Process the results stored in ground_states_list
    CALL print_eigenlist(ground_states_list)

    CALL delete_eigen_list(ground_states_list)
    CALL finalize_mpi() ! Finalize MPI if needed
  END SUBROUTINE run_diagonalization

END MODULE example_block_lanczos_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `genvar`: For various global constants (`DP`, `log_unit`, etc.) and MPI-related variables (`iproc`, `rank`, `size2`, `no_mpi`).
    *   `common_def`: For utility functions like string conversion (`c2s`, `i2c`), messaging (`dump_message`), and timing (`reset_timer`, `timer_fortran`).
    *   `eigen_class`: For data structures `eigen_type` (single eigenpair) and `eigenlist_type` (list of eigenpairs), and associated procedures (`add_eigen`, `delete_eigen`, `new_eigen`).
    *   `globalvar_ed_solver`: Provides runtime control parameters for the Lanczos algorithm (e.g., `block_size`, `demax`, `neigen`, `nitermax`, `tolerance`).
    *   `h_class`: Crucially depends on this module for:
        *   `Hmult__`: The user-provided subroutine that computes the matrix-vector product \(H \cdot \text{vector}\).
        *   `dimen_H()`: A function that returns the dimension of the Hamiltonian matrix.
        *   `title_H_()`: A function that returns a descriptive title for the Hamiltonian (used in logging).
    *   `mpi_mod`: The `split` subroutine is used to determine data distribution among MPI processes.
    *   `random`: The `randomize_mat` subroutine is used to generate initial random vectors when requested by BLZPACK.

*   **External Libraries:**
    *   **BLZPACK**: This module is fundamentally a Fortran wrapper around the BLZPACK library's `BLZDRD` routine, which performs the core Block Lanczos iterations. The compilation of this code is conditional on the `_blocklanczos` preprocessor macro being defined.
    *   **MPI (Message Passing Interface)**: Used for parallel execution. Direct MPI calls (e.g., `MPI_ALLGATHERV`) are made for distributing and collecting vector data across processes. The standard Fortran MPI interface module (`mpi`) is used.

*   **Interactions with other components:**
    *   This module serves as the primary numerical engine for diagonalizing the Hamiltonian matrix in an ED solver.
    *   It takes the definition of the Hamiltonian implicitly through `h_class` (which provides `dimen_H` and `Hmult__`).
    *   It receives control parameters from `globalvar_ed_solver`.
    *   It returns the computed low-energy eigenvalues and eigenvectors in an `eigenlist_type` object, which can then be used by other parts of the solver for calculating physical observables, correlation functions, etc.
    *   The module effectively orchestrates the complex interplay between the user's Hamiltonian representation, the parallel processing environment (MPI), and the specialized algorithms within the external BLZPACK library.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
