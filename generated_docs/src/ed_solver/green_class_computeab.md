# `src/ed_solver/green_class_computeab.f90`

## Overview

*   **Purpose:** This module provides a specialized subroutine, `compute_greenAB`, for computing specific components of a Green's function that involve correlations of two potentially different types of operators, generically denoted as \(A\) and \(B\). Examples include \(\langle T A(\tau) B(0) \rangle\) or \(\langle T A(\tau) B^\dagger(0) \rangle\). The calculation is based on the eigenstates and eigenvalues obtained from an Exact Diagonalization (ED) study of the system.
*   **Role in Project:** Serves as a dedicated computational routine within an ED solver framework, tailored for "AB-type" two-time correlation functions. It manages the application of operators \(A\) (or \(A^\dagger\)) and \(B\) (or \(B^\dagger\)) to the system's eigenstates, the computation of the necessary mixed matrix elements (e.g., \(\langle \Psi_k | A | \Psi_m \rangle \langle \Psi_m | B | \Psi_n \rangle\)), and the subsequent summation according to the Lehmann representation to yield the frequency-dependent Green's function components. The specific nature of operators \(A\) and \(B\) is abstracted through procedure arguments (interfaces).

## Key Components

*   **`MODULE green_class_computeAB`**:
    The main module containing the `compute_greenab` public interface, which maps to the `compute_greenAB` subroutine.

*   **Public Subroutines:**
    *   **`SUBROUTINE compute_greenAB(greenAB, greenBA, greenAA, greenBB, AIM, beta, GS, Asector, Bsector, applyA, applyB, COMPUTE_DYN)`**:
        This is the core subroutine for computing "AB-type" Green's function components.
        *   **Arguments:**
            *   `greenAB (TYPE(green_type), INTENT(INOUT))`: The `green_type` object (from `green_class.f90`) that will be populated with the computed \(\langle A...B...\rangle\) type correlation functions.
            *   `greenBA (TYPE(green_type), INTENT(INOUT))`: The `green_type` object for the corresponding \(\langle B...A...\rangle\) type correlations.
            *   `greenAA (TYPE(green_type), INTENT(IN))`: A precomputed `green_type` object for \(\langle A...A...\rangle\) correlations. Likely used by `reshuffle_vecAB` for normalization or symmetrization.
            *   `greenBB (TYPE(green_type), INTENT(IN))`: A precomputed `green_type` object for \(\langle B...B...\rangle\) correlations. Likely used by `reshuffle_vecAB`.
            *   `AIM (TYPE(AIM_type), INTENT(IN))`: The definition of the Anderson Impurity Model.
            *   `beta (REAL(DP), INTENT(IN))`: The inverse temperature \(1/(k_B T)\).
            *   `GS (TYPE(eigensectorlist_type), INTENT(IN))`: A list containing the eigenstates (and eigenvalues) of the AIM, from a prior diagonalization.
            *   `Asector (INTERFACE)`: An abstract interface for a user-provided subroutine `Asector(Asec, pm, sector)` that determines the sector `Asec` resulting from applying operator \(A\) (or \(A^\dagger\), via `pm`) to an initial `sector`.
            *   `Bsector (INTERFACE)`: Similar to `Asector`, but for operator \(B\).
            *   `applyA (INTERFACE)`: An abstract interface for `applyA(Aes, pm, MASK, es, rankr)` that applies operator \(A\) (or \(A^\dagger\)) to eigenstate `es%lowest%eigen(rankr)`, storing result in `Aes`.
            *   `applyB (INTERFACE)`: Similar to `applyA`, but for operator \(B\).
            *   `COMPUTE_DYN (LOGICAL, OPTIONAL, INTENT(IN))`: If `.TRUE.`, the dynamic (frequency-dependent) part of the Green's function is computed.
        *   **Functionality:**
            1.  Calls an internal `init_data()` routine (defined in the included `"green_class_computeAB.h"`).
            2.  Calculates the partition function `Zpart` and ground state energy `E0` from `GS`.
            3.  The primary computational logic is encapsulated within the included file `"green_class_computeAB.h"`. This hidden code is expected to:
                *   Iterate over initial eigenstates \(|\Psi_n\rangle\) from `GS` and intermediate states \(|\Psi_m\rangle\).
                *   Use the passed `applyA` and `applyB` procedures to generate states like \(A|\Psi_n\rangle\) and \(B|\Psi_m\rangle\) (or their adjoints).
                *   Calculate matrix elements of the form \(\langle \Psi_k | A | \Psi_m \rangle\) and \(\langle \Psi_m | B | \Psi_n \rangle\), or more complex multi-operator matrix elements if the Green's function involves more than two operators explicitly.
                *   Sum these contributions with appropriate energy denominators and Boltzmann factors according to the Lehmann representation to compute components for `greenAB` and `greenBA`.
            4.  The visible Fortran code shows loops over `isector`, `ipm`, `jpm`, `iph` that structure the calls to routines within the include file (e.g., `init_sector`, `print_label`, `do_ab`, `do_ba`). `do_ab` likely computes \(\langle A...B...\rangle\) terms and `do_ba` computes \(\langle B...A...\rangle\) terms.
            5.  Handles MPI parallelism if `FLAG_MPI_GREENS > 0`, distributing parts of the calculation and summing results using `mpisum`.
            6.  After the main computation, it calls `reshuffle_vecAB` (from `green_class_compute_symmetric.f90`). This routine likely uses the computed `greenAB`, `greenBA` along with the input `greenAA`, `greenBB` to form specific symmetric linear combinations or to ensure proper normalization/ordering of the vectorized Green's function components.
            7.  Finally, calls an internal `clean_everything()` routine.

*   **Included Files:**
    *   `"green_class_computeAB.h"`: This include file is central and contains the detailed implementation for calculating the "AB-type" Green's function components. It relies heavily on the `Asector`, `Bsector`, `applyA`, and `applyB` procedures passed as arguments.

## Important Variables/Constants

*   **Input/Output `greenAB`, `greenBA :: TYPE(green_type)`**: These `green_type` objects are populated with the computed correlation function data.
*   **Input `greenAA`, `greenBB :: TYPE(green_type)`**: These are likely used by the `reshuffle_vecAB` post-processing step, possibly for constructing specific linear combinations or for normalization purposes.
*   **Input `AIM :: TYPE(AIM_type)`**: Provides the definition of the Anderson Impurity Model.
*   **Input `GS :: TYPE(eigensectorlist_type)`**: Contains the many-body eigenstates and eigenvalues of the AIM, crucial for the Lehmann sums.
*   **Input Interfaces `Asector`, `Bsector`, `applyA`, `applyB`**: These are of paramount importance as they define the specific operators \(A\) and \(B\) whose correlation is being computed. By supplying different concrete Fortran subroutines that match these interface definitions, a wide variety of two-operator correlation functions can be calculated using this single framework.

## Usage Examples

The `compute_greenAB` subroutine is intended to be called by higher-level routines that manage the computation of various physical observables.

```fortran
MODULE example_compute_ab_usage
  USE green_class_computeAB
  USE green_class, ONLY: green_type, new_green, delete_green ! For setup/teardown
  USE aim_class, ONLY: AIM_type
  USE eigen_sector_class, ONLY: eigensectorlist_type
  USE sector_class, ONLY: sector_type ! For dummy interface procedures
  ! ... other necessary modules for AIM, GS, and operator definitions ...

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calculate_my_ab_greens_function

! Dummy procedures matching the Asector, Bsector, applyA, applyB interfaces
CONTAINS
  SUBROUTINE my_operator_A_sector(Asec, pm, sector)
    TYPE(sector_type), INTENT(INOUT) :: Asec; TYPE(sector_type), INTENT(IN) :: sector
    CHARACTER(LEN=1), INTENT(IN)   :: pm; Asec = sector ! Placeholder
  END SUBROUTINE
  SUBROUTINE my_operator_B_sector(Bsec, pm, sector)
    TYPE(sector_type), INTENT(INOUT) :: Bsec; TYPE(sector_type), INTENT(IN) :: sector
    CHARACTER(LEN=1), INTENT(IN)   :: pm; Bsec = sector ! Placeholder
  END SUBROUTINE
  SUBROUTINE my_apply_operator_A(Aes, pm, MASK, es, rankr)
    USE eigen_sector_class, ONLY: eigensector_type
    TYPE(eigensector_type), INTENT(INOUT) :: Aes; CHARACTER(LEN=1), INTENT(IN) :: pm
    TYPE(eigensector_type), INTENT(IN) :: es; INTEGER, INTENT(IN) :: rankr
    LOGICAL, INTENT(IN) :: MASK(:); Aes%sector = es%sector; Aes%lowest = es%lowest%eigen(rankr) ! Placeholder
  END SUBROUTINE
  SUBROUTINE my_apply_operator_B(Bes, pm, MASK, es, rankr)
    USE eigen_sector_class, ONLY: eigensector_type
    TYPE(eigensector_type), INTENT(INOUT) :: Bes; CHARACTER(LEN=1), INTENT(IN) :: pm
    TYPE(eigensector_type), INTENT(IN) :: es; INTEGER, INTENT(IN) :: rankr
    LOGICAL, INTENT(IN) :: MASK(:); Bes%sector = es%sector; Bes%lowest = es%lowest%eigen(rankr) ! Placeholder
  END SUBROUTINE

  SUBROUTINE calculate_my_ab_greens_function
    TYPE(green_type) :: gf_AB, gf_BA, gf_AA_ref, gf_BB_ref ! Assume these are initialized
    TYPE(AIM_type) :: model_def  ! Assume initialized
    REAL(DP) :: temp_beta        ! Assume initialized
    TYPE(eigensectorlist_type) :: ed_sol ! Assume populated

    ! --- Initialize all green_type objects, model_def, temp_beta, ed_sol ---
    ! (This is a complex setup step not detailed here)

    PRINT *, "Computing AB-type Green's function..."
    ! CALL compute_greenAB(greenAB=gf_AB, greenBA=gf_BA, greenAA=gf_AA_ref, greenBB=gf_BB_ref, &
    ! &                   AIM=model_def, beta=temp_beta, GS=ed_sol, &
    ! &                   Asector=my_operator_A_sector, Bsector=my_operator_B_sector, &
    ! &                   applyA=my_apply_operator_A, applyB=my_apply_operator_B, &
    ! &                   COMPUTE_DYN=.TRUE.)
    PRINT *, "AB-type Green's function computation finished."

    ! gf_AB and gf_BA would now contain the results.
    ! --- Cleanup of green_type objects and other types ---
  END SUBROUTINE calculate_my_ab_greens_function
END MODULE example_compute_ab_usage
```
*(Note: The example above is conceptual. Full initialization of all `green_type` arguments and realistic operator procedures are necessary for a working call.)*

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `genvar`: For `DP` (double precision kind), `messages3` (debugging messages), MPI variables (`rank`, `size2`).
    *   `aim_class`: For `AIM_type` (Anderson Impurity Model definition).
    *   `eigen_class`: For `eigen_type` (representing a single eigenstate or an operator applied to a state).
    *   `eigen_sector_class`: For `eigensector_type` (grouping eigenstates by sector) and `eigensectorlist_type` (list of all computed eigensectors `GS`), and utilities like `gsenergy`, `partition`.
    *   `globalvar_ed_solver`: For global flags controlling execution, such as `FLAG_MPI_GREENS`, `USE_TRANSPOSE_TRICK_MPI`.
    *   `green_class`: For `green_type`, the data structure this module populates.
    *   `H_class`: The `delete_H` routine is called (within the include file's `clean_everything`), suggesting that Hamiltonian context might be temporarily modified.
    *   `mask_class`: For `mask_type` (likely used within the included file for operator definitions or projections).
    *   `mpi_mod`: For `mpisum` (for summing results from different MPI processes).
    *   `sector_class`: For `sector_type` (used by the `Asector` and `Bsector` interfaces).
    *   `green_class_compute_symmetric`: The `reshuffle_vecAB` subroutine is called after the main computation, implying that results might be combined or symmetrized using `greenAA` and `greenBB` as inputs.
    *   **Include File:** Critically, `#include "green_class_computeAB.h"` contains a significant portion of the computational logic, including routines like `init_data`, `init_sector`, `print_label`, `do_ab`, `do_ba`, `scan_indices`, `fix_indices`, `loop_over_states`, and `clean_everything`.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules.
    *   MPI (Message Passing Interface) if `FLAG_MPI_GREENS > 0` and `size2 > 1`, for parallel summation of results.

*   **Interactions with other components:**
    *   This module is a specialized computational routine, likely called by higher-level master routines within the ED solver (e.g., from `correlations.f90` or its includes like `correlations_dyn_stat.h`).
    *   Its main purpose is to calculate Green's functions involving two, possibly different, operators \(A\) and \(B\). The generic nature is achieved through the `Asector`, `Bsector`, `applyA`, and `applyB` procedure arguments (interfaces). By supplying different concrete subroutines for these, various types of two-operator correlation functions can be computed.
    *   The results (\(\langle A...B...\rangle\) and \(\langle B...A...\rangle\)) are stored in the `greenAB` and `greenBA` `green_type` objects. These can then be used for further physical analysis or output.
    *   The core calculations (Lehmann sums, operator applications, matrix element computations) are primarily handled within the included `"green_class_computeAB.h"` file.
    *   The `reshuffle_vecAB` call at the end suggests that the raw computed data might be further processed, possibly to form specific linear combinations or ensure certain symmetry properties, using `greenAA` and `greenBB` as auxiliary inputs.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
