# `src/ed_solver/globalvar_ed_solver.f90`

## Overview

*   **Purpose:** This module centralizes a wide array of global variables, parameters, and constants that are used to control the behavior, settings, and execution flow of the entire Exact Diagonalization (ED) solver suite.
*   **Role in Project:** Acts as a global configuration registry or a common block for the ED solver. Various modules within the solver access these variables to determine algorithmic choices, convergence criteria, physical model parameters, I/O settings, debugging levels, and specific feature flags. This allows for consistent behavior and easy tuning of the solver's operation.

## Key Components

*   **`MODULE globalvar_ed_solver`**:
    This module primarily consists of public variable declarations. It does not define any public subroutines or functions itself. The core components are the global variables it makes available to other parts of the ED solver system.

## Important Variables/Constants

This module declares a large number of variables. Below is a categorized selection of some of the key public variables that influence the ED solver's behavior:

*   **General Control & Debugging:**
    *   `ALL_FIRST_CALL :: LOGICAL, PARAMETER`: A compile-time parameter that likely forces "first call" initialization logic in various modules to run every time if true (perhaps for debugging or specific reset scenarios). Default is often `.FALSE.`.
    *   `verboseall :: LOGICAL, PARAMETER`: A global verbosity flag, often enabled by the `debug` preprocessor macro, to increase detailed output throughout the solver.
    *   `iterdmft :: INTEGER`: Current DMFT iteration number (if the ED solver is used within a DMFT loop).
    *   `niterdmft :: INTEGER`: Total number of DMFT iterations to perform.
    *   `istati :: INTEGER`: A global status variable, likely used to check the success of allocations/deallocations across different modules.
    *   `ed_num_eigenstates_print :: INTEGER`: Controls how many eigenstates are printed in detailed analyses (e.g., in `density_matrix.f90`).
    *   `print_qc :: LOGICAL`: A flag to enable printing of Quality Control (QC) information.

*   **Diagonalization Parameters (for Lanczos/ARPACK/Full ED):**
    *   `Nitermax :: INTEGER`: Maximum number of iterations for iterative solvers like Lanczos or ARPACK.
    *   `tolerance :: REAL(DP)`: Convergence tolerance for eigenvalues in iterative solvers.
    *   `dEmax, dEmax0 :: REAL(DP)`: Energy window above the ground state. Eigenstates outside this window might be discarded or handled differently.
    *   `Neigen :: INTEGER`: Target number of eigenvalues/eigenvectors to compute.
    *   `Block_size :: INTEGER`: Block size for the Block Lanczos algorithm. If 0, the algorithm might choose an optimal size.
    *   `which_lanczos :: CHARACTER(20)`: String to select a specific Lanczos implementation if multiple are available.
    *   `FLAG_FULL_ED_GREEN :: LOGICAL`: If `.TRUE.`, indicates that all Green's function components from a full ED should be stored/used.

*   **Bath Fitting & Hybridization Function Control:**
    *   `Niter_search_max, Niter_search_max_0 :: INTEGER`: Maximum iterations for the bath hybridization fitting procedure (initial and subsequent iterations).
    *   `search_step :: REAL(DP)`: Initial step size for the minimizer in bath fitting.
    *   `dist_max :: REAL(DP)`: Convergence tolerance for the distance metric in hybridization fitting.
    *   `fit_meth :: CHARACTER(20)`: Specifies the minimization method for bath fitting (e.g., "CONJUG_GRAD").
    *   `fit_nw :: INTEGER`: Number of frequency points to use when calculating the distance in fitting.
    *   `window_hybrid, window_hybrid2 :: INTEGER`: Define frequency windows for fitting.
    *   `window_weight :: INTEGER`: Relative weight for frequencies inside vs. outside the `window_hybrid`.
    *   `fit_green_func_and_not_hybrid :: LOGICAL`: If true, fits to match a Green's function rather than the hybridization function directly.
    *   `skip_fit_ :: LOGICAL`: If true, skips the numerical fitting of bath parameters.
    *   `use_specific_set_parameters :: LOGICAL`: If true, uses bath parameters directly from `param_input` without fitting.
    *   `param_input(:), param_output(:) :: REAL(DP), ALLOCATABLE`: Arrays to hold flat vectors of input and output bath parameters for fitting routines.
    *   `freeze_poles_delta :: LOGICAL`, `freeze_poles_delta_iter, freeze_poles_delta_niter :: INTEGER`: Parameters to control a pole-freezing scheme during bath fitting.
    *   `FROZEN_POLES(:) :: REAL(DP), ALLOCATABLE`: Array to store the positions of frozen poles.

*   **CPT (Cluster Perturbation Theory) Parameters:**
    *   `gen_cpt :: LOGICAL`: If true, enables CPT calculations/parameter generation.
    *   `ncpt_approx, ncpt_tot, ncpt_approx_tot, ncpt_para, ncpt_chain_coup :: INTEGER`: Integers defining the structure and size of the CPT cluster/approximation.
    *   `cpt_upper_bound, cpt_lagrange :: REAL(DP)`: Constraints and Lagrange multipliers for CPT fitting.
    *   `epsilon_cpt(:,:), T_cpt(:,:) :: COMPLEX(DP), ALLOCATABLE`: Matrices used in CPT calculations.

*   **Physical Model & Basis Parameters:**
    *   `beta_ED, beta_ED_ :: REAL(DP)`: Inverse temperature for ED calculations.
    *   `Jhund :: REAL(DP)`: Hund's coupling strength for multi-orbital impurities.
    *   `UUmatrix(:,:), JJmatrix(:,:) :: REAL(DP), ALLOCATABLE`: Matrices defining the Coulomb U and exchange J interactions for multi-orbital impurities.
    *   `flag_slater_int :: LOGICAL`: If true, indicates that Slater integrals are being used for interactions.
    *   `para_state, supersc_state, superconducting_state :: LOGICAL`: Flags indicating the physical state being modeled (paramagnetic, generic superconducting, or specific BCS-like superconductivity).
    *   `force_sz_basis, force_nupdn_basis :: LOGICAL`: Flags to enforce the use of Sz-adapted basis or (Nup, Ndown) product basis, respectively.
    *   `force_no_pairing, force_para_state, force_singlet_state :: LOGICAL`: Flags to enforce specific symmetries or turn off pairing terms.

*   **Correlation Function & Output Control:**
    *   `average_G :: INTEGER`: Flag controlling how Green's functions are averaged (e.g., over equivalent orbitals or sites).
    *   `MASK_AVERAGE(:,:), MASK_AVERAGE_SIGMA(:,:) :: INTEGER, ALLOCATABLE`: Integer masks defining which elements to include/exclude during averaging of Green's functions or self-energies.

*   **Keldysh Formalism Parameters:**
    *   `do_keldysh :: LOGICAL`: Enables calculations on the Keldysh contour (for non-equilibrium dynamics).
    *   `keldysh_t0, keldysh_tmax, keldysh_delta :: REAL(DP)`: Time parameters defining the Keldysh contour.
    *   `keldysh_n :: INTEGER`: Number of time steps on the contour.

*   **FMOS (Factorized Moment Solver) Parameters:**
    *   `fmos, fmos_fluc, fmos_hub1 :: LOGICAL`: Flags to enable FMOS or its variants (e.g., including fluctuations, Hubbard-I like).
    *   `fmos_iter :: INTEGER`: Number of iterations for the FMOS method.
    *   `fmos_mix :: REAL(DP)`: Mixing parameter for FMOS iterations.

*   **MPI & Performance Flags:**
    *   `OPEN_MP :: LOGICAL`: Flag related to OpenMP usage.
    *   `USE_TRANSPOSE_TRICK_MPI :: LOGICAL`: Enables an MPI optimization for matrix transpositions.
    *   `use_cuda_lanczos :: LOGICAL`: If true, attempts to use CUDA-accelerated Lanczos routines.

## Usage Examples

This module itself is not executed. Its variables are `USE`d by other modules within the ED solver. The typical workflow involves:
1.  An input-handling part of the program (e.g., a main driver routine or a dedicated input module) reads configuration from a file or command line.
2.  These values are used to set the public variables in `globalvar_ed_solver`.
3.  Other modules in the ED solver suite then `USE globalvar_ed_solver` and access these variables to guide their execution.

```fortran
! In some input module or main program:
USE globalvar_ed_solver
IMPLICIT NONE
! ... read input_tolerance, input_neigen from a file ...
tolerance = input_tolerance
Neigen = input_neigen
beta_ED = 100.0_DP ! Set temperature
! ... set other global variables ...

! In another module, e.g., a diagonalization routine:
USE globalvar_ed_solver
IMPLICIT NONE
! ...
IF (current_residual < tolerance) THEN
  ! Converged
END IF
! ... uses Neigen to know how many states to find ...
```

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   `genvar`: This module only uses `genvar` for the `DP` kind parameter definition for its real variables. It does not have procedural dependencies on other modules as it primarily serves as a container for global data.
*   **External Libraries:**
    *   None. This module consists only of variable declarations and parameter definitions.

*   **Interactions with other components:**
    *   **Central Configuration:** This module is "USEd" by nearly all other modules within the `src/ed_solver/` directory (and potentially by higher-level scripts or codes that use the ED solver). It provides a centralized way to access shared parameters and settings.
    *   **Customization of Behavior:** The values of these global variables dictate how other modules perform their tasks. For example:
        *   Diagonalization routines (`ed_arpack.f90`, `block_lanczos.f90`) use `Neigen`, `tolerance`, `Block_size`, etc.
        *   Bath fitting routines (`bath_class_hybrid.f90`) use `Niter_search_max`, `dist_max`, `fit_meth`, CPT parameters, etc.
        *   Correlation function calculations (`correlations.f90`) use flags like `average_G`, `do_keldysh`, and parameters like `beta_ED`.
        *   Basis and sector definition modules (`fermion_hilbert_class.f90`, `sector_class.f90`) might use flags like `force_sz_basis`.
    *   **Initialization:** The variables in this module are expected to be initialized by a higher-level part of the application, typically after reading user inputs or configuration files, before the main computational routines of the ED solver are called.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
