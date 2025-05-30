# `src/ed_solver/dmft_solver_ed.f90`

## Overview

*   **Purpose:** Implements the core Exact Diagonalization (ED) solver functionality within the DMFT framework.
*   **Role in Project:** Performs the ED calculations, manages impurity problems, bath parameters, and computes Green's functions and self-energies. This module is central to solving the quantum impurity model that arises in Dynamical Mean-Field Theory.

## Key Components

*   `SUBROUTINE solver_ED_interface(nbathparam_in, iter_, mmu, g_out, self_out, hybrid_in, Eimp, sigw, rdens, impurity_, Himp, gw, frequ, flip_input_output, para_state_, retarded, restarted_, compute_all, tot_rep, spm_, corhop_, only_density, skip_fit, use_specific_set_parameters_, param_input_, param_output_, nbathparam_, freeze_poles_delta_, freeze_poles_delta_iter_, freeze_poles_delta_niter_, imp_causality, Jhund_, use_input_delta_instead_of_fit, Jhund_matrix)`: The main entry point for the ED solver within a DMFT loop. It takes numerous inputs including bath parameters, chemical potential, hybridization function, impurity definition, and control flags. It outputs the computed Green's function, self-energy, density, and other relevant quantities.
*   `SUBROUTINE stand_alone_ed()`: A routine for running the ED solver as a standalone program. It reads parameters from an input file named `PARAMS` (defining cluster size, interaction matrices, frequencies, etc.) and then calls the main ED solver routines. Output Green's functions and self-energies are written to files like `_sigma_output_full_1`, `_sigma_output_full_real_1`, `green_output_matsu_full1`, etc.
*   `SUBROUTINE initialize_ed_solver(mmu, bbeta, impurity_, Himp, rdelta_width, nmatsu_frequ, FLAGIMP, nstep, Nww, wwmin, wmax, para_state_, supersc_state_, bath_param_ed, rdelta_frequ_eta1_, rdelta_frequ_T_, rdelta_frequ_w0_, init, average_G_, UUmatrix_in_size, UUmatrix_in, fit_green_func_and_not_hybrid_, donot_compute_holepart_)`: Sets up the necessary data structures and parameters for the ED calculation. This includes defining the impurity model, initializing correlation functions, reading bath parameters (from `BATHfile` or input), and setting up the Anderson Impurity Model (AIM) object. It also handles Slater integrals if `Uc.dat.ED` is present.
*   `SUBROUTINE finalize_ed_solver(finalize_iter, finalize)`: Cleans up resources used by the ED solver, such as deallocating arrays associated with the impurity and AIM objects.
*   `FUNCTION flip_matrix(mat, flip_input_output) AS COMPLEX(KIND=DP)`: A utility function that conditionally flips elements of a complex matrix `mat`. If `flip_input_output` is true and indices `i, j` are in the second half of the matrix dimensions, the element `mat(i,j)` is replaced by `-conjg(mat(i,j))`. This is often used in particle-hole transformations or handling specific symmetry sectors.
*   `FUNCTION flip_matrix_Eimp(mat, flip_input_output) AS COMPLEX(KIND=DP)`: Similar to `flip_matrix`, but specifically tailored for the impurity energy matrix `Eimp`. If `flip_input_output` is true, it modifies the second block-diagonal part (e.g., spin-down block) by `Eimp(Nc+i, Nc+j) = -mat(Nc+j, Nc+i)`.
*   `SUBROUTINE cpt_extract_G(vec, sig, g_out, imaginary)`: Used when Cluster Perturbation Theory (CPT) is applied (`ncpt_approx /= 0`). It calculates the Green's function `g_out` from the self-energy `sig` using the CPT hopping matrix `T_cpt` via the Dyson-like equation `G = (omega * I - T_cpt - Sigma)^-1`.
*   **Included Files:**
    *   `"dmft_ed_solver_read_write.h"`: Likely contains subroutines for reading DMFT parameters (`read_DMFT_parameters`) and writing header/footer information for logs (`write_header`, `write_footer`).
    *   `"dmft_ed_solver_hubbardI.h"`: Suggests routines related to the Hubbard-I approximation, potentially including `hub1_run` which is called if `fmos_hub1` is true.
    *   `"dmft_ed_solver_fmos.h"`: Points to routines for an FMOS (Factorized Moment Solver or similar) approach, like `fmos_solver`, called if `fmos` is true.
    *   `"dmft_ed_solver_dump_fort.h"`: Likely for debugging, containing routines to dump data to Fortran binary files (e.g., `dump_to_fortb`).
    *   `"correlations_keldysh.h"`: Contains routines for calculations on the Keldysh contour, such as the `keldysh` subroutine.
    *   `"dmft_ed_solver_hubbardI_run.h"`: Possibly contains the implementation for `hub1_run` or related helper routines for Hubbard-I.

## Important Variables/Constants

*   `CHARACTER(LEN=100) :: BATHfile`: Stores the filename for reading bath parameters (e.g., "bath.dat").
*   `CHARACTER(LEN=100) :: CORRELfile`: Stores the filename for correlation data (e.g., "correl.dat").
*   `CHARACTER(LEN=100) :: EDfile`: Stores the filename for main ED solver settings/input (e.g., "ed.inp").
*   `LOGICAL :: start_from_old_gs`: A flag indicating whether the calculation should attempt to restart from a previously computed and saved ground state.
*   `CHARACTER(LEN=100) :: OLDGSFILE`: If `start_from_old_gs` is true, this variable holds the filename of the saved ground state (e.g., "gs.dat").
*   `TYPE(impurity_type), SAVE :: impurity`: A structure (derived type) that encapsulates all properties of the quantum impurity, such as the number of orbitals (`N`), interaction strengths (`U`, `J`), and local energy levels.
*   `TYPE(bath_type), SAVE :: bath`: A structure that holds the parameters describing the non-interacting bath coupled to the impurity. This includes bath energies, hybridization strengths, and the hybridization function itself (`hybrid%fctn`).
*   `TYPE(AIM_type), SAVE :: AIM`: A structure representing the complete Anderson Impurity Model, combining the `impurity` and `bath` objects.
*   `TYPE(eigensectorlist_type), SAVE :: GS`: A structure that stores the results of the exact diagonalization, primarily the ground state wavefunction and energy, but also potentially a list of excited states (eigen sectors).
*   `TYPE h1occup`: A derived type likely used within Hubbard-I calculations (`hub1_run`) to store Hamiltonian matrices (`Hn`), energies (`En`), and other sector-specific information for occupancy calculations.

## Usage Examples

Typically called by higher-level DMFT scripts (e.g., in `src/scripts/`) via the `solver_ED_interface` subroutine as part of an iterative self-consistency loop.
The `stand_alone_ed` subroutine allows for direct ED calculations. To use it:
1.  Prepare an input file named `PARAMS`.
2.  This file should contain:
    *   `cluster_problem_size` (integer)
    *   `Eimp_` matrices (spin-up, then spin-down blocks)
    *   `UUmatrix_loc` (interaction matrix)
    *   `nmatsu_frequ` (number of Matsubara frequencies for input Delta)
    *   `nmatsu_long` (number of Matsubara frequencies for internal calculations)
    *   `nreal_frequ` (number of real frequencies)
    *   `compute_all` (logical flag)
    *   `FLAG_ORDER_TYPE` (integer for para/ferro/AFM state)
    *   `retarded` (logical for retarded Green's function)
    *   `frequ_min`, `frequ_max` (real frequency range)
    *   `mmu` (chemical potential)
    *   `bbeta` (inverse temperature)
    *   `dU1` (impurity potential shift)
    *   `rdelta_width`, `rdelta_frequ_eta1_`, `rdelta_frequ_T_`, `rdelta_frequ_w0_` (broadening/delta function parameters)
    *   `my_iter_dmft` (iteration number, negative for first run)
    *   `JJhund` (Hund's coupling J)
    *   `use_input_delta` (logical)
    *   `no_cdw` (logical)
    *   `fit_green` (logical)
    *   `nohole` (logical)
    *   Optionally: `JJmatrix_loc` (full Hund's matrix), `cluster_full` (logical).
3.  Execute the program that calls `stand_alone_ed`.
Outputs will be generated in files like `_sigma_output_full_1`, `green_output_matsu_full_1`, etc.

Further examples detailing specific use cases are to be added by developers.

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   **Modules:**
        *   `AIM_class`: For `AIM_type` definition and related procedures.
        *   `bath_class`: For `bath_type` definition and bath manipulation.
        *   `impurity_class`: For `impurity_type`, `hamiltonian`, `web` types, and impurity setup.
        *   `eigen_sector_class`: For `eigensectorlist_type` to store eigenstates.
        *   `genvar`: Provides global constants like `DP` (double precision), `pi`, MPI rank/size.
        *   `mesh`: For frequency mesh generation (e.g., `build1dmesh`, `mirror_array`).
        *   `splines`: For resampling functions (e.g., `resampleit`).
        *   `globalvar_ed_solver`: Contains global parameters and flags specific to the ED solver.
        *   `matrix`: Utility functions for matrix operations (e.g., `diag`, `write_array`, `id`, `invmat`).
        *   `bath_class_hybrid`: For converting between bath parameters and hybridization function (`bath2hybrid`, `hybrid2bath`).
        *   `mpi_mod`: For MPI parallelization primitives (e.g., `mpibarrier`).
        *   `solver`: Generic solver framework, including `solve_aim` and `init_solver`.
        *   `correlations`: For defining and computing correlation functions (`G`, `GF`, `SNAMBU`, etc.).
        *   `correl_class`: For `average_correlations`.
        *   `linalg`: Provides linear algebra routines (e.g., `modi`).
        *   `common_def`: Common definitions and utilities (e.g., `utils_system_call`, `utils_assert`).
        *   `random`: For random number generation (`init_rantab`, `rand_init`).
        *   `stringmanip`: For string manipulation (e.g., `tostring`).
        *   `strings`: Additional string utilities (e.g., `remove`).
    *   **Include Files:**
        *   `"dmft_ed_solver_read_write.h"`: For I/O operations related to solver parameters.
        *   `"dmft_ed_solver_hubbardI.h"`: For Hubbard-I approximation specific code.
        *   `"dmft_ed_solver_fmos.h"`: For FMOS specific code.
        *   `"dmft_ed_solver_dump_fort.h"`: For dumping data for debugging.
        *   `"correlations_keldysh.h"`: For Keldysh formalism calculations.
        *   `"dmft_ed_solver_hubbardI_run.h"`: For running Hubbard-I calculations.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules.
    *   MPI (Message Passing Interface) library, used via `mpi_mod` for parallel computations if compiled with MPI support.
    *   Numerical libraries like LAPACK and BLAS might be implicitly used through other modules (e.g., `linalg` or during compilation if linked by the build system) for matrix diagonalizations and other intensive linear algebra operations.

*   **Interactions with other components:**
    *   This module serves as a core computational engine, typically invoked by higher-level DMFT scripts (e.g., potentially located in `src/scripts/`) or other driver programs that manage the self-consistency loop.
    *   It heavily relies on various `*_class.f90` modules (e.g., `AIM_class`, `bath_class`, `impurity_class`, `eigen_sector_class`, `correl_class`) for object-oriented data structures and associated methods that define and manipulate the quantum impurity problem.
    *   Utilizes a suite of utility modules (`genvar`, `mesh`, `splines`, `matrix`, `linalg`, `common_def`, `random`, `stringmanip`, `strings`) for tasks such as array manipulation, mesh generation, mathematical operations, and I/O.
    *   Shared global parameters for the ED solver are accessed from `globalvar_ed_solver`.
    *   Parallel execution is managed through `mpi_mod`.
    *   The generic `solver` module provides the framework for the `solve_aim` routine, which performs the core diagonalization.
    *   The `correlations` module is essential for calculating various types of Green's functions and self-energies once the AIM eigenstates are found.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
