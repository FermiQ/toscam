# `src/ed_solver/green_class_computeaa.f90`

## Overview

*   **Purpose:** This module provides a specialized subroutine, `compute_greenAA`, for computing specific components of a Green's function. These components typically involve time-ordered correlations of the same type of operator, generically denoted as \(A\), such as \(\langle T A(\tau) A(0) \rangle\) or \(\langle T A(\tau) A^\dagger(0) \rangle\). The calculation uses results from an Exact Diagonalization (ED) study, namely the eigenstates and eigenvalues of the system.
*   **Role in Project:** Serves as a dedicated computational routine within the broader ED solver framework, specifically tailored for "AA-type" correlation functions. It orchestrates the process of applying the operator \(A\) (or its adjoint \(A^\dagger\)) to the system's eigenstates, calculating the necessary matrix elements, and performing the summations required by the Lehmann representation to obtain the frequency-dependent Green's function components. The specific nature of operator \(A\) is made generic through the use of procedure arguments (interfaces).

## Key Components

*   **`MODULE green_class_computeAA`**:
    The main module containing the `compute_greenaa` public interface, which maps to the `compute_greenAA` subroutine.

*   **Public Subroutines:**
    *   **`SUBROUTINE compute_greenAA(green, AIM, beta, GS, Asector, applyA, COMPUTE_DYN, keldysh_level, GS_out)`**:
        This is the core subroutine for computing "AA-type" Green's function components.
        *   **Arguments:**
            *   `green (TYPE(green_type), INTENT(INOUT))`: The `green_type` object (from `green_class.f90`) that will be populated with the computed Green's function components.
            *   `AIM (TYPE(AIM_type), INTENT(IN))`: The definition of the Anderson Impurity Model.
            *   `beta (REAL(DP), INTENT(IN))`: The inverse temperature \(1/(k_B T)\).
            *   `GS (TYPE(eigensectorlist_type), INTENT(IN))`: A list containing the eigenstates (and eigenvalues) of the AIM, obtained from a prior diagonalization.
            *   `Asector (INTERFACE)`: An abstract interface for a user-provided subroutine: `SUBROUTINE Asector(Asec, pm, sector)`. This subroutine must determine the resulting sector `Asec` (of `sector_type`) after an operator \(A\) (or \(A^\dagger\), distinguished by `pm = '+'` or `'-'`) acts on an initial `sector`.
            *   `applyA (INTERFACE)`: An abstract interface for a user-provided subroutine: `SUBROUTINE applyA(Aes, pm, MASK, es, rankr)`. This subroutine must apply the operator \(A\) (or \(A^\dagger\), based on `pm`) associated with a `MASK` to a specific eigenstate `es%lowest%eigen(rankr)` and store the resulting state(s) in `Aes` (of `eigensector_type`).
            *   `COMPUTE_DYN (LOGICAL, OPTIONAL, INTENT(IN))`: If `.TRUE.`, the dynamic (frequency-dependent) part of the Green's function is computed. Defaults to a value that likely enables dynamic computation if not present.
            *   `keldysh_level (INTEGER, OPTIONAL, INTENT(IN))`: An optional argument to trigger or control calculations on the Keldysh contour (for non-equilibrium Green's functions).
            *   `GS_out (TYPE(eigensectorlist_type), OPTIONAL)`: An optional argument, its specific use is not clear from the visible code but might be relevant for calculations within the included files that transform the ground state or relevant states.
        *   **Functionality:**
            1.  Calls an internal `init_data()` routine (likely defined in the include file).
            2.  The core computational logic is primarily encapsulated within the included file `"green_class_computeAA.h"`. This included code is expected to:
                *   Iterate over the relevant eigenstates \(|\Psi_n\rangle\) from `GS`.
                *   For each state, and for each component of the Green's function to be computed (e.g., \(\langle A A^\dagger \rangle\), \(\langle A^\dagger A \rangle\)), it uses the passed `applyA` procedure to generate states like \(A|\Psi_n\rangle\) or \(A^\dagger|\Psi_n\rangle\).
                *   It then computes matrix elements of the form \(\langle \Psi_m | A | \Psi_n \rangle\).
                *   These matrix elements, along with energy differences \(E_m - E_n\) and Boltzmann factors \(e^{-\beta E_n}/Z\), are used in Lehmann sums to calculate the frequency-dependent Green's function components, which are stored in `green%correl(ipm,jpm)%fctn` and static parts in `green%correlstat(ipm,jpm)%rc%mat`.
            3.  The visible Fortran structure shows loops iterating over `isector` (sectors in `GS`), `ipm` and `jpm` (indices for Green's function components like `(1,1), (1,2), (2,1), (2,2)`), `iph` (particle/hole channels), and `ieigen` (eigenstates within a sector). These loops likely set up contexts for calls into the routines within the include file.
            4.  If MPI is enabled (`FLAG_MPI_GREENS > 0`), it handles parallel computation by distributing loops (e.g., `DO kk = rank + 1, ktot, size2`) and then summing up results using `mpisum`.
            5.  After the primary computation, it calls `reshuffle_vecAA` (from `green_class_compute_symmetric.f90`), which likely applies symmetry operations or reorders the vectorized representation of the Green's function components.
            6.  Finally, it calls an internal `clean_everything()` routine.

*   **Included Files:**
    *   `"green_class_computeAA.h"`: This file is central and contains the detailed implementation for calculating the "AA-type" Green's function components. It makes extensive use of the `Asector` and `applyA` procedures passed as arguments.

## Important Variables/Constants

*   **Input `green :: TYPE(green_type)`**: The `green_type` object that will be populated with the computed correlation function data.
*   **Input `AIM :: TYPE(AIM_type)`**: The definition of the Anderson Impurity Model.
*   **Input `GS :: TYPE(eigensectorlist_type)`**: The list of many-body eigenstates and eigenvalues of the AIM.
*   **Input Interfaces `Asector` and `applyA`**: These are crucial. They abstract the specific nature of the operator \(A\). By providing different concrete Fortran subroutines that match these interface definitions, one can use `compute_greenAA` to calculate various types of correlation functions (e.g., single-particle Green's functions where \(A\) is a fermion annihilation/creation operator, spin-spin correlations if \(A\) is a spin operator, etc.).

## Usage Examples

The `compute_greenAA` subroutine would be called by a higher-level routine that orchestrates the calculation of different types of correlation functions.

```fortran
MODULE example_compute_aa_usage
  USE green_class_computeAA
  USE green_class, ONLY: green_type, new_green, delete_green ! For setup/teardown
  USE aim_class, ONLY: AIM_type
  USE eigen_sector_class, ONLY: eigensectorlist_type
  USE sector_class, ONLY: sector_type ! For dummy interface procedures
  ! ... other necessary modules for AIM, GS, and operator definitions ...

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calculate_my_aa_greens_function

! Dummy procedures matching the Asector and applyA interfaces for illustration
CONTAINS
  SUBROUTINE my_operator_sector_calculator(Asec, pm, sector)
    TYPE(sector_type), INTENT(INOUT) :: Asec
    TYPE(sector_type), INTENT(IN)    :: sector
    CHARACTER(LEN=1), INTENT(IN)   :: pm
    ! ... logic to determine Asec based on operator A and initial sector ...
    PRINT *, "MyAsector called with pm=", pm
    Asec = sector ! Placeholder
  END SUBROUTINE my_operator_sector_calculator

  SUBROUTINE my_operator_applicator(Aes, pm, MASK, es, rankr)
    USE eigen_sector_class, ONLY: eigensector_type
    TYPE(eigensector_type), INTENT(INOUT) :: Aes
    CHARACTER(LEN=1), INTENT(IN)        :: pm
    TYPE(eigensector_type), INTENT(IN)    :: es
    INTEGER, INTENT(IN)                   :: rankr
    LOGICAL, INTENT(IN)                   :: MASK(:)
    ! ... logic to apply operator A to es%lowest%eigen(rankr) and store in Aes ...
    PRINT *, "MyApplya called with pm=", pm, " for rank=", rankr
    Aes%sector = es%sector ! Placeholder
    Aes%lowest = es%lowest%eigen(rankr) ! Placeholder
  END SUBROUTINE my_operator_applicator

  SUBROUTINE calculate_my_aa_greens_function
    TYPE(green_type) :: specific_aa_gf
    TYPE(AIM_type) :: model_definition  ! Assume initialized
    REAL(DP) :: temperature_beta        ! Assume initialized
    TYPE(eigensectorlist_type) :: ed_results ! Assume populated by diagonalization

    LOGICAL :: compute_flags(2,2)

    ! --- Initialize model_definition, temperature_beta, ed_results ---
    ! Example:
    ! compute_flags = .FALSE.
    ! compute_flags(1,1) = .TRUE. ! Request <A A_dag> type component
    ! compute_flags(2,2) = .TRUE. ! Request <A_dag A> type component

    ! CALL new_green(specific_aa_gf, compute=compute_flags, title="My_AA_GF", N=..., Nw=..., etc.)
    ! (Full initialization of specific_aa_gf needed here)

    PRINT *, "Computing AA-type Green's function..."
    ! CALL compute_greenAA(green=specific_aa_gf, AIM=model_definition, beta=temperature_beta, &
    ! &                   GS=ed_results, Asector=my_operator_sector_calculator, &
    ! &                   applyA=my_operator_applicator, COMPUTE_DYN=.TRUE.)
    PRINT *, "AA-type Green's function computation finished."

    ! specific_aa_gf would now contain the results.
    ! CALL delete_green(specific_aa_gf)
    ! --- Cleanup other types ---
  END SUBROUTINE calculate_my_aa_greens_function

END MODULE example_compute_aa_usage
```
*(Note: The example above is conceptual as `specific_aa_gf` needs full initialization before the call, and the dummy operator procedures need actual logic.)*

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `genvar`: For `DP` (double precision kind), `log_unit`, `messages3` (debugging messages), MPI variables (`rank`, `size2`).
    *   `aim_class`: For `AIM_type` (Anderson Impurity Model definition).
    *   `common_def`: For `dump_message` (logging).
    *   `eigen_class`: For `eigen_type` (representing a single eigenstate or an operator applied to a state).
    *   `eigen_sector_class`: For `eigensector_type` (grouping eigenstates by sector) and `eigensectorlist_type` (list of all computed eigensectors `GS`).
    *   `globalvar_ed_solver`: For various global flags controlling execution, such as `donot_compute_holepart`, `FLAG_MPI_GREENS`, `USE_TRANSPOSE_TRICK_MPI`.
    *   `green_class`: For `green_type`, which is the data structure this module populates.
    *   `H_class`: The `delete_H` routine is called, suggesting that the Hamiltonian context might be modified or cleared by routines within the included file.
    *   `mask_class`: For `mask_type` (likely used within the included file for operator definitions or projections).
    *   `mpi_mod`: For `mpisum` (for summing results from different MPI processes).
    *   `sector_class`: For `sector_type` (used by the `Asector` interface).
    *   `green_class_compute_symmetric`: The `reshuffle_vecAA` subroutine is called after the main computation, implying that results might be symmetrized or reordered.
    *   **Include File:** Critically, `#include "green_class_computeAA.h"` contains a significant portion of the computational logic.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules.
    *   MPI (Message Passing Interface) if `FLAG_MPI_GREENS > 0` and `size2 > 1`, for parallel summation of results.

*   **Interactions with other components:**
    *   This module is typically invoked by higher-level routines within the ED solver system, likely from the main `correlations.f90` module (or its own included files like `correlations_dyn_stat.h`), which is responsible for calculating a variety of different correlation functions.
    *   The key interaction mechanism is through the `Asector` and `applyA` procedure arguments (interfaces). The calling routine provides specific subroutines that define the operator \(A\) (e.g., a fermion creation/annihilation operator for single-particle Green's functions, a spin operator for spin susceptibilities). This makes `compute_greenAA` generic for any "AA-type" correlation.
    *   The results of the computation (the Green's function components) are stored in the `green_type` object passed as an argument. This object can then be used for further analysis, written to disk, or used in self-consistency cycles (like in DMFT).
    *   The core calculations performed within the included `"green_class_computeAA.h"` file would involve:
        *   Iterating over the initial eigenstates \(|\Psi_n\rangle\) from `GS`.
        *   Applying the operator \(A\) (or \(A^\dagger\)) using the `applyA` procedure to get intermediate states \(A|\Psi_n\rangle\).
        *   Calculating matrix elements like \(\langle \Psi_m | A | \Psi_n \rangle\).
        *   Summing these contributions according to the Lehmann representation of the desired Green's function component, using eigenvalues from `GS` and the input `beta`.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
