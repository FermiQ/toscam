# `src/ed_solver/correlations.f90`

## Overview

*   **Purpose:** This module manages the definition, computation, and output of various physical correlation functions essential for characterizing an Anderson Impurity Model (AIM) solved via Exact Diagonalization (ED). These include Green's functions, self-energies, and potentially higher-order correlation functions or dynamic susceptibilities.
*   **Role in Project:** Acts as a central orchestrator for post-diagonalization analysis. It initializes data structures for storing numerous correlation functions, then uses the eigenstates and eigenvalues obtained from the ED solver (`GS` of `eigensectorlist_type`) and the AIM definition (`AIM_type`) to compute these functions. It supports various representations (e.g., normal spin-resolved, Nambu for superconductivity, Matsubara frequency, retarded real frequency) and handles the output of these calculated quantities.

## Key Components

*   **`MODULE correlations`**:
    *   **Public Module Variables (Data Storage):** This module declares several public, `SAVE`d variables of `correl_type` and `green_type` that serve as global containers for the computed correlation functions. These include:
        *   `SNAMBU, SNAMBUret :: TYPE(correl_type)`: Nambu self-energy (Matsubara and retarded).
        *   `GNAMBU, GNAMBUN, GNAMBUret :: TYPE(green_type)`: Nambu Green's functions (Matsubara, a variant `GNAMBUN` possibly for non-interacting or specific components, and retarded). `GNAMBU` typically holds \(\langle T \Psi \Psi^\dagger \rangle\) where \(\Psi = (c_\uparrow, c^\dagger_\downarrow)^T\).
        *   `G(2), GN(2), GKret(2) :: TYPE(green_type)`: Arrays for spin-resolved Green's functions. `G(1)` for spin-up \(\langle T c_\uparrow c^\dagger_\uparrow \rangle\), `G(2)` for spin-down \(\langle T c_\downarrow c^\dagger_\downarrow \rangle\). `GN` likely for non-interacting versions, `GKret` for Keldysh retarded Green's functions.
        *   `GF(2), GNF(2) :: TYPE(green_type)`: Arrays for spin-resolved anomalous Green's functions (e.g., \(\langle T c_\uparrow c_\downarrow \rangle\), \(\langle T c^\dagger_\uparrow c^\dagger_\downarrow \rangle\)). `GF(1)` for \(\langle c c \rangle\) type, `GF(2)` for \(\langle c^\dagger c^\dagger \rangle\) type.
        *   `Spm :: TYPE(green_type)`: For spin-flip correlation functions, like \(\langle S^+ S^- \rangle\).
    *   **Private Module Variables (Internal State/Data):**
        *   `Gret(2), GFret(2)`: Internal storage for retarded normal and anomalous Green's functions.
        *   `Sz, N, Szret, Spmret, Nret`: For static and dynamic \(S_z\), number \(N\), and spin-flip \(S^{\pm}\) correlations.
        *   `P3, P4, P3ret, P4ret`: For correlations involving 3-site and 4-site permutation operators.
        *   `CHI, CHIret`: For susceptibilities.
        *   Various control flags and parameters like `SUPER` (superconducting flag), `Nc` (number of impurity sites), `beta`, `width`, `triplets` and `quadruplets` definitions (for permutation operators), `compute_Sz`, etc.

*   **Public Subroutines:**
    *   **`SUBROUTINE init_correlations(...)`**: (Implementation details are in the included `correlations_init.h`). This routine initializes all the public (and private) module-level `correl_type` and `green_type` variables. It sets up their matrix dimensions (`N` based on `AIM%Nc`), frequency meshes (`freq` component using input `beta`, `wmin`, `wmax`, `nmatsu_freq`, `nreal_freq`), titles, particle statistics (`FERMIONIC`/`BOSONIC`), and masks for storing matrix elements.
    *   **`SUBROUTINE compute_correlations(AIM, GS, retard, COMPUTE_ALL_CORREL, only_dens)`**: This is the main computational engine of the module.
        *   **Arguments:**
            *   `AIM (TYPE(AIM_type), INTENT(IN))`: The definition of the Anderson Impurity Model.
            *   `GS (TYPE(eigensectorlist_type), INTENT(IN))`: The solution from the ED solver, containing many-body eigenstates and eigenvalues.
            *   `retard (LOGICAL, INTENT(IN))`: If true, retarded (real-frequency) functions are computed in addition to Matsubara ones.
            *   `COMPUTE_ALL_CORREL (LOGICAL, OPTIONAL, INTENT(IN))`: If true, a comprehensive set of correlation functions is computed. Otherwise, a minimal set might be computed.
            *   `only_dens (LOGICAL, OPTIONAL, INTENT(IN))`: If true, computations might be restricted to those needed for density calculations.
        *   **Functionality:**
            *   Initializes operator application modules (`apply_C`, `apply_NS`, `apply_P`) on the first call.
            *   Orchestrates calls to internal routines (likely found in the included `correlations_dyn_stat.h`) to compute:
                *   Matsubara frequency Green's functions (normal and Nambu).
                *   Retarded Green's functions if requested.
                *   Bosonic correlation functions (dynamic or equal-time) like spin susceptibilities, charge susceptibilities, etc., if `COMPUTE_ALL_CORREL` is true.
                *   Reduced density matrix elements.
            *   Calls `write_static` to output equal-time quantities derived from the Green's functions or density matrix.
    *   **`SUBROUTINE compute_self_energy(S, G, AIM, retarded)`**: Calculates a self-energy `S` (of `correl_type`) from a given Green's function `G` (of `correl_type`) and the AIM definition. It uses the Dyson equation \(S = G_0^{-1} - G^{-1}\). \(G_0^{-1}\) is constructed from the non-interacting parts of the AIM: local impurity energies \(E_c\) (from `AIM%impurity%Ec`) and the hybridization function \(\Delta(\omega)\) (obtained from `AIM%bath` via `bath2hybrid`). Handles averaging based on `average_G` flags.
    *   **`SUBROUTINE Nambu_green(GNAMBU, G, GF)`**: Constructs the 2x2 Nambu block Green's function `GNAMBU%correl(ipm,mipm)%fctn` from spin-resolved normal Green's functions `G(1)%correl`, `G(2)%correl` and anomalous Green's functions `GF(1)%correl`, `GF(2)%correl`. The standard Nambu structure is:
        \[
        \hat{G}_{Nambu} = \begin{pmatrix} \langle c_\uparrow c^\dagger_\uparrow \rangle & \langle c_\uparrow c_\downarrow \rangle \\ \langle c^\dagger_\downarrow c^\dagger_\uparrow \rangle & \langle c^\dagger_\downarrow c_\downarrow \rangle \end{pmatrix}
        \]
        The routine correctly maps the components from `G` and `GF` into this structure, including necessary transformations (e.g., \(G_{22}(i\omega_n) = -\text{conj}(G_{\downarrow\downarrow}(-i\omega_n))\) for the \( (c^\dagger_\downarrow c_\downarrow) \) block). It also handles a case where only one component of `GNAMBU` might be computed directly, and others are filled by symmetry.
    *   **`SUBROUTINE write_info_correlations(...)`**: (Implementation in `correlations_write.h`). Writes summary information or metadata about the correlation functions that have been initialized and are available for computation.
    *   **`SUBROUTINE write_static(AIM, GS)`**: (Implementation in `correlations_write.h`). Computes and writes static (equal-time) quantities, such as local occupations, double occupancy, or local moments, derived from the computed Green's functions or by direct expectation values over the ground state `GS`.

*   **Included Files:**
    *   `"correlations_write.h"`: Expected to contain various subroutines for writing the computed correlation functions (stored in the module-level variables) to output files in appropriate formats.
    *   `"correlations_init.h"`: Expected to contain the main `init_correlations` subroutine and helper routines to set up the diverse array of `correl_type` and `green_type` objects.
    *   `"correlations_dyn_stat.h"`: Expected to house the core computational routines that calculate the dynamic (frequency-dependent) and static (equal-time) correlation functions. These routines would use the Lehmann representation, applying operators (via `apply_C`, `apply_NS`, `apply_P`) to the many-body eigenstates from `GS` to compute matrix elements and then sum them with appropriate energy denominators and Boltzmann factors.

## Important Variables/Constants

*   **Module-level `correl_type` and `green_type` variables** (e.g., `SNAMBU`, `GNAMBU`, `G(1)`, `G(2)`, `GF(1)`, `GF(2)`): These are the primary containers that store the results of the correlation function computations. Their structure and content are central to this module's purpose.
*   **Input arguments to `compute_correlations`:**
    *   `AIM :: TYPE(AIM_type)`: Provides the complete definition of the Anderson Impurity Model (energies, couplings, interactions).
    *   `GS :: TYPE(eigensectorlist_type)`: Contains the many-body eigenstates and eigenvalues obtained from the ED solver, which are essential for computing correlation functions via Lehmann sums.

## Usage Examples

```fortran
MODULE example_correlations_usage
  USE correlations
  USE aim_class, ONLY: AIM_type
  USE impurity_class, ONLY: impurity_type
  USE eigen_sector_class, ONLY: eigensectorlist_type
  ! ... other necessary modules for AIM and GS setup ...

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: run_correlation_computation

CONTAINS
  SUBROUTINE run_correlation_computation
    TYPE(AIM_type) :: aim_definition
    TYPE(impurity_type), TARGET :: impurity_model ! Part of aim_definition
    TYPE(eigensectorlist_type) :: ed_solution
    REAL(DP) :: target_beta, real_freq_min, real_freq_max, r_delta_width
    INTEGER :: num_matsubara_freq, num_real_freq

    ! --- Initialize aim_definition and ed_solution ---
    ! (This would involve setting up impurity_model, bath_model,
    !  then aim_definition, and running a diagonalization to get ed_solution)
    ! Example:
    ! impurity_model%Nc = 1
    ! target_beta = 50.0
    ! real_freq_min = -10.0; real_freq_max = 10.0; r_delta_width = 0.05
    ! num_matsubara_freq = 1024; num_real_freq = 512
    ! --- End of initialization ---


    ! 1. Initialize all correlation function structures
    CALL init_correlations(impurity_def=impurity_model, beta=target_beta, &
     &                     fileout_correl="my_correlations.dat", & ! A base name for some outputs
     &                     rdelta_width_in=r_delta_width, &
     &                     wmin_in=real_freq_min, wmax_in=real_freq_max, &
     &                     nmatsu_freq=num_matsubara_freq, nreal_freq=num_real_freq)
    PRINT *, "Correlation function structures initialized."

    ! 2. Compute the correlation functions using the ED solution
    PRINT *, "Starting computation of correlations..."
    CALL compute_correlations(AIM=aim_definition, GS=ed_solution, retard=.TRUE., COMPUTE_ALL_CORREL=.TRUE.)
    PRINT *, "Correlation computation finished."

    ! 3. Write out information (and potentially data via routines in correlations_write.h)
    CALL write_info_correlations()

    ! The computed data is now in the public module variables like SNAMBU%fctn, G(1)%correl(1,1)%fctn etc.
    ! For example, to access the (0,0) component of spin-up Matsubara Green's function:
    ! PRINT *, "G_up(0,0,iw=1): ", G(1)%correl(1,1)%fctn(1,1,1)


    ! (Cleanup of aim_definition, ed_solution would happen here)
  END SUBROUTINE run_correlation_computation

END MODULE example_correlations_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `correl_class`: Provides the `correl_type` definition.
    *   `green_class`: Provides the `green_type` definition (which is likely a container of multiple `correl_type` objects, e.g., for different particle/hole components or time orderings).
    *   `readable_vec_class`: For `readable_vec_type`, used for `vec_list` (likely for I/O of triplet/quadruplet definitions).
    *   `aim_class`: For `AIM_type`, which defines the physical model.
    *   `impurity_class`: Used by `compute_self_energy` (via `nambu_ec`) and `init_correlations`.
    *   `eigen_sector_class`: For `eigensectorlist_type`, which provides the eigenstates and eigenvalues needed for Lehmann sums.
    *   `apply_c`, `apply_ns`, `apply_p`: These modules are crucial as they provide the routines to apply \(c, c^\dagger, N, S_z, S^\pm, P_3, P_4\) operators to many-body states. These operations are fundamental for calculating the matrix elements \(\langle m |Operator| n \rangle\) required in the Lehmann sums for correlation functions. These are likely called extensively by routines within the included `correlations_dyn_stat.h`.
    *   `common_def`, `genvar`, `globalvar_ed_solver`, `timer_mod`: For utilities, global parameters, and timing.
    *   `masked_matrix_class`, `matrix`, `mask_class`: For matrix and mask operations.
    *   `bath_class_hybrid`: `bath2hybrid` is used in `compute_self_energy`.
    *   **Include Files:** `correlations_write.h`, `correlations_init.h`, `correlations_dyn_stat.h` contain significant portions of the module's logic.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules.

*   **Interactions with other components:**
    *   This module is a high-level consumer of the results from the core ED solver (which produces the `eigensectorlist_type` object `GS`).
    *   It uses the physical model definition (`AIM_type`) to correctly interpret states and operators.
    *   The core calculations (likely in `correlations_dyn_stat.h`) involve iterating over pairs of eigenstates `|m>`, `|n>`, applying relevant physical operators (using `apply_C`, `apply_NS`, `apply_P`) to get matrix elements like \(\langle m|A|n \rangle\), and then constructing the correlation functions using these matrix elements, energy differences \(E_m - E_n\), and Boltzmann factors \(e^{-\beta E_n}/Z\).
    *   The computed correlation functions (e.g., `SNAMBU`, `G(1)`) are then made available for further analysis, output (via routines in `correlations_write.h`), or use in self-consistency loops (e.g., the self-energy might be used to update the AIM parameters in a DMFT cycle).

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
