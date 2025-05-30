# `src/ed_solver/bath_class_hybrid.f90`

## Overview

*   **Purpose:** This module provides the core numerical routines for (1) calculating the hybridization function \(\Delta(\omega)\) from a given set of bath parameters (encapsulated in `bath_type`), and (2) fitting these bath parameters to reproduce a target hybridization function. These processes are fundamental to self-consistency cycles in Dynamical Mean-Field Theory (DMFT) and other quantum embedding methods.
*   **Role in Project:** Implements the two-way transformation between a parameterized model of a non-interacting bath and its corresponding hybridization function. The `bath2hybrid` routine computes \(\Delta(\omega)\) from bath parameters. The `hybrid2bath` routine performs the more complex task of optimizing bath parameters to best match a given target \(\Delta(\omega)\), often obtained from the impurity solver in a previous DMFT iteration. The module can handle different frequency types (Matsubara for `FERMIONIC`, real for `RETARDED`) and can optionally include Cluster Perturbation Theory (CPT) parameters in the fitting process to model more complex environments.

## Key Components

*   **Module-Level Private Variables:**
    *   `hybrid2fit :: TYPE(correl_type)`: Stores the target hybridization function that the bath parameters are being fitted against during the `hybrid2bath` routine.
    *   `batht :: TYPE(bath_type)`: A temporary or working instance of `bath_type` used internally during the fitting procedure in `hybrid2bath` to hold the current set of parameters being evaluated.
    *   `icount_ :: INTEGER`: An internal counter, likely used for tracking iterations or calls within the fitting process, primarily for logging or conditional debug output.

*   **Public Subroutines:**
    *   **`FUNCTION fill_mat_from_vector(n, vec) RESULT(fill_mat_from_vector(n,n)) REAL(DP)`**:
        A utility function that constructs a symmetric `n x n` matrix from a flat vector `vec` of size `n*(n+1)/2`. The vector is assumed to contain the elements of the upper (or lower) triangle of the symmetric matrix.
    *   **`SUBROUTINE bath2hybrid(BATH, FREQTYPE, WRITE_HYBRID, short, cptvec, cpt_build_matrix, Vweight)`**:
        Calculates the hybridization function \(\Delta(\omega)\) for a given `BATH` object and stores it in `BATH%hybrid` (if `FREQTYPE` is `FERMIONIC`) or `BATH%hybridret` (if `FREQTYPE` is `RETARDED`).
        *   **Arguments:**
            *   `BATH (TYPE(bath_type))`: The input bath object with defined parameters.
            *   `FREQTYPE (CHARACTER(LEN=9), INTENT(IN))`: Specifies frequency type (`FERMIONIC` or `RETARDED`).
            *   `WRITE_HYBRID (LOGICAL, OPTIONAL, INTENT(IN))`: If true, prints/writes the computed hybridization.
            *   `short (LOGICAL, OPTIONAL, INTENT(IN))`: If true and `fit_nw > 0`, limits calculation to `fit_nw` frequencies.
            *   `cptvec (REAL(DP), TARGET, OPTIONAL)`: Array of CPT parameters, used if `ncpt_approx_tot > 0`.
            *   `cpt_build_matrix (LOGICAL, OPTIONAL)`: If true and CPT is active, reconstructs the full CPT hopping matrix `T_cpt`.
            *   `Vweight (REAL(DP), OPTIONAL)`: Output, sum of squares of CPT hopping parameters.
        *   **Functionality:**
            *   Transforms bath parameters (`Eb`, `Vbc`, `Pb`, `PVbc`) to Nambu representation using routines from `bath_class`.
            *   If CPT is active (`ncpt_approx_tot > 0`), incorporates `cptvec` into the effective bath Hamiltonian `Eb`.
            *   Computes \(\Delta(\omega) = V_{cb} (\omega \mathbf{1} - E_{bath}^{eff})^{-1} V_{bc}\) using matrix inversion or eigenvalue decomposition, depending on flags like `fast_fit` and `freeze_poles_delta`.
            *   Handles special cases like `fit_green_func_and_not_hybrid` where the Green's function \(G = (\omega - E_c - \Delta(\omega))^{-1}\) is fitted instead of \(\Delta(\omega)\).
    *   **`SUBROUTINE hybrid2bath(bath)`**:
        Optimizes the parameters of the input `bath` object (e.g., `bath%Eb`, `bath%Vbc`, and potentially CPT parameters) so that the hybridization function calculated from these parameters matches a target hybridization function. The target is initially provided in `bath%hybrid`.
        *   **Functionality:**
            *   Copies the target `bath%hybrid` to the internal `hybrid2fit`.
            *   If `fit_green_func_and_not_hybrid` is true, transforms `hybrid2fit` into a target Weiss field.
            *   Uses `batht` as a working copy for optimization.
            *   Converts `batht` parameters to a flat vector using `bath2vec` (from `bath_class_vec`).
            *   Calls `minimize_func_wrapper` (from `minimization_wrapping`), providing `distance_func` as the function to minimize.
            *   The optimized parameter vector is then converted back to `bath_type` structure using `vec2bath`.
            *   Finally, calls `bath2hybrid` to update `bath%hybrid` and `bath%hybridret` with the result from the optimized parameters.

*   **Internal Subroutines (used by `hybrid2bath`):**
    *   **`SUBROUTINE distance_func(dist, n, vec)`**: The objective function passed to the minimizer.
        *   Takes a flat vector `vec` of `n` parameters.
        *   Updates `batht` using `vec2bath(vec, batht)`.
        *   Calls `bath2hybrid(batht, FERMIONIC, ...)` to compute \(\Delta(\omega)\) for the current `batht`.
        *   Calculates the squared distance `dist` between `batht%hybrid` and the target `hybrid2fit` using `diff_correl`.
        *   If CPT is active, adds a penalty term `Vweight * cpt_lagrange` to `dist`.
    *   **`SUBROUTINE plot_some_results()`**: If `rank == 0`, attempts to plot comparisons of the fitted hybridization vs. the target. (Actual plotting commands are commented out).
    *   **`SUBROUTINE plot_it_()`**: If `verbose_graph` is true and `rank == 0`, plots intermediate fitting results. (Actual plotting commands are commented out).

## Important Variables/Constants

*   **Module-private `hybrid2fit :: TYPE(correl_type)`**: This variable is central to the `hybrid2bath` routine. It stores the target hybridization function (copied from the input `bath%hybrid`) that the fitting procedure aims to reproduce.
*   **Module-private `batht :: TYPE(bath_type)`**: This is a working copy of the `bath_type` object used within the `hybrid2bath` fitting loop. The minimizer modifies the parameters stored in `batht` (via its vector representation), and `distance_func` evaluates the quality of these parameters.
*   **`BATH :: TYPE(bath_type)` (Argument to `bath2hybrid` and `hybrid2bath`)**: The primary data structure holding the bath parameters. In `bath2hybrid`, it's an input whose parameters are used. In `hybrid2bath`, its `hybrid` component provides the target, and its parameters are updated with the optimized values.
*   **`cptvec :: REAL(DP), POINTER / INTENT(IN)` (Argument to `bath2hybrid`)**: An optional array containing CPT parameters, which are included in the effective bath Hamiltonian if the CPT approximation is active. These can also be part of the fitting variables in `hybrid2bath`.
*   **`vec :: REAL(DP), INTENT(IN)` (Argument to `distance_func`)**: A flat vector containing all current bath parameters (and possibly CPT parameters) being evaluated by the optimizer.

## Usage Examples

```fortran
MODULE example_hybrid_usage
  USE bath_class, ONLY: bath_type, read_bath, delete_bath
  USE bath_class_hybrid
  USE correl_class, ONLY: correl_type ! For my_bath%hybrid
  USE genvar, ONLY: FERMIONIC, RETARDED

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: test_hybrid_fitter

CONTAINS
  SUBROUTINE test_hybrid_fitter(bath_input_file_name, num_imp_sites)
    CHARACTER(LEN=*), INTENT(IN) :: bath_input_file_name
    INTEGER, INTENT(IN) :: num_imp_sites
    TYPE(bath_type) :: my_bath

    ! 1. Initialize bath, e.g., by reading from a file.
    !    Assume bath_input_file_name points to a file that also defines
    !    the target hybridization function in my_bath%hybrid (or it's set some other way).
    CALL read_bath(my_bath, bath_input_file_name, num_imp_sites)
    ! (Alternatively, my_bath%hybrid could be populated from an impurity solver's output)

    ! 2. Fit the bath parameters in 'my_bath' to match the target 'my_bath%hybrid'.
    PRINT *, "Starting hybrid2bath fit..."
    CALL hybrid2bath(my_bath)
    PRINT *, "Finished hybrid2bath fit."
    ! Now, my_bath contains the optimized parameters, and my_bath%hybrid/hybridret
    ! reflect the hybridization achieved with these optimized parameters.

    ! 3. Optionally, explicitly recalculate and view the hybridization from the fitted parameters.
    CALL bath2hybrid(my_bath, FERMIONIC, WRITE_HYBRID=.TRUE.)
    CALL bath2hybrid(my_bath, RETARDED, WRITE_HYBRID=.TRUE.)

    CALL delete_bath(my_bath)
  END SUBROUTINE test_hybrid_fitter

END MODULE example_hybrid_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `correl_class`: For `correl_type` and its associated functions (`glimpse_correl`, `average_correlations`, `copy_correl`, `new_correl`, `diff_correl`).
    *   `bath_class`: For `bath_type` definition and utility functions like `nambu_eb`, `nambu_vbc`.
    *   `genvar`: Provides global constants (`DP`, `FERMIONIC`, `RETARDED`) and MPI-related variables (`iproc`, `rank`, etc.).
    *   `globalvar_ed_solver`: Supplies numerous global flags and parameters that control the behavior of `bath2hybrid` and `hybrid2bath` (e.g., `fit_green_func_and_not_hybrid`, CPT flags like `ncpt_approx`, fitting control like `niter_search_max`).
    *   `impurity_class`: The `Eccc` variable (a constant part of impurity self-energy) is used if `fit_green_func_and_not_hybrid` is enabled.
    *   `masked_matrix_class`: For `delete_masked_matrix` and `masked_matrix_type`.
    *   `matrix`: For matrix operations like `bande_mat`, `eigenvector_matrix`, `invmat`, `matmul_x`, `write_array`.
    *   `bath_class_vec`: Provides `bath2vec` (converts `bath_type` to a flat parameter vector) and `vec2bath` (converts vector back to `bath_type`), which are essential for interfacing with the numerical minimizer.
    *   `minimization_wrapping`: The `minimize_func_wrapper` subroutine is called by `hybrid2bath` to perform the numerical optimization of bath parameters.
    *   `common_def`: For utilities like `dump_message`, `reset_timer`, `timer_fortran`.
    *   `mpi_mod`: For MPI communication (`mpibarrier`, `mpibcast`) used within `hybrid2bath` to synchronize parameters across processes if run in parallel.
    *   `random`: The `dran_tab` function is used for generating initial random CPT parameters in `hybrid2bath`.
    *   `stringmanip`: The `tostring` function is used in the internal plotting routines.
    *   `timer_mod`: For timing sections of the code like `start_timer`, `stop_timer`.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules.
    *   MPI (Message Passing Interface) if the parallel capabilities of the minimization routines (called by `hybrid2bath`) are utilized.

*   **Interactions with other components:**
    *   This module forms a bridge between the parameterized representation of the bath (`bath_type`) and its physical effect as a hybridization function (`correl_type`).
    *   `bath2hybrid` consumes a `bath_type` object and produces/updates its `hybrid` and `hybridret` components. This is used by the impurity solver to get the current effective bath.
    *   `hybrid2bath` is a core part of the DMFT self-consistency loop. It takes a `bath_type` object whose `hybrid` component is the *target* (e.g., from \(G_{latt}^{-1} - \Sigma_{imp}\) or similar DMFT update), and modifies the parameters within this `bath_type` object to make its calculated hybridization match this target.
    *   The fitting process in `hybrid2bath` relies on `minimization_wrapping` for the optimization algorithm and `bath_class_vec` for parameter vectorization.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
