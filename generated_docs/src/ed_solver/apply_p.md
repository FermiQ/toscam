# `src/ed_solver/apply_p.f90`

## Overview

*   **Purpose:** This module provides routines to apply cyclic permutation operators (specifically 3-cycles and 4-cycles of site indices) to many-body eigenstates. Such operations are typically used for symmetrizing wavefunctions with respect to site permutations or for calculating observables that involve these permutations, like scalar spin chirality.
*   **Role in Project:** Functions as a specialized computational kernel within an Exact Diagonalization (ED) solver. It takes an input eigenstate, definitions of site permutations (triplets or quadruplets), and applies these permutations to the eigenstate. The module supports different many-body basis representations (Sz-adapted and up/down product basis).

## Key Components

*   **`MODULE apply_P`**:
    *   **Module-Level Variables:**
        *   `INTEGER, ALLOCATABLE :: CYC3(:, :, :)`, `CYC4(:, :, :)`: Arrays to store the definitions of 3-cycles and 4-cycles. `CYC3(i,j,1)` stores the j-th site of the i-th 3-cycle (original order), and `CYC3(i,j,2)` stores a related permutation (e.g., inverse or a specific element first). Similar for `CYC4`.
        *   `INTEGER, ALLOCATABLE :: IMPiorbupdo(:, :)`, `IMPiorbsz(:, :)`: Pointers to arrays storing impurity orbital mappings, similar to those in `apply_C` and `apply_NS`. `IMPiorbupdo` is for up/down product basis, `IMPiorbsz` for Sz-adapted basis.
        *   `INTEGER :: nCYC3, nCYC4`: Number of defined 3-cycles and 4-cycles.
        *   `INTEGER :: Ns, Ns2`: Total number of sites in AIM and `Ns*2`.
        *   `LOGICAL, PARAMETER :: force_reset_list = .false.`: If true, clears the target eigenstate list before adding results.

*   **Public Subroutines:**
    *   **`SUBROUTINE init_apply_P(triplets, quadruplets, AIM)`**: Initializes the module.
        *   **Arguments:**
            *   `triplets(:, :) (INTEGER, INTENT(IN))`: Array defining 3-cycles (e.g., `triplets(i, 1:3)` lists the 3 site indices for the i-th cycle).
            *   `quadruplets(:, :) (INTEGER, INTENT(IN))`: Array defining 4-cycles.
            *   `AIM (TYPE(AIM_type), INTENT(IN))`: The Anderson Impurity Model definition for orbital mappings.
        *   **Functionality:** Stores the cycle definitions in `CYC3` and `CYC4` (including a second related permutation for each cycle). Initializes `IMPiorbupdo`, `IMPiorbsz`, `Ns`, `Ns2` from `AIM`.
    *   **`SUBROUTINE P_sector(Csec, pm, sector_in)`**: Determines the resulting sector `Csec` after applying a permutation. Since permutations of identical particles do not change particle number or total Sz, the output sector is simply a copy of `sector_in`. `pm` is a dummy argument.
    *   **`SUBROUTINE apply_P3(Ces, pm, MASK, es, esrank)`**: Applies selected 3-cycle permutations to the eigenstate `es%lowest%eigen(esrank)`.
        *   Calls the generic worker `apply_Pn` with `ncyc=3`.
        *   `pm` ('+' or '-') selects between the original cycle in `CYC3(:,:,1)` or the related one in `CYC3(:,:,2)`.
        *   `MASK` is a logical array indicating which of the defined 3-cycles to apply.
    *   **`SUBROUTINE apply_P4(Ces, pm, MASK, es, esrank)`**: Applies selected 4-cycle permutations, calling `apply_Pn` with `ncyc=4`.
    *   **`SUBROUTINE apply_CHI(Ces, pm, MASK, es, esrank)`**: Applies an operator seemingly related to scalar spin chirality or a similar three-site operator. It computes \(P_3^+(\text{state}) - P_3^-(\text{state})\), where \(P_3^+\) and \(P_3^-\) are results from `apply_Pn` with `ncyc=3` and `pm='+'` and `pm='-'` respectively. `MASK` selects which triplet to use.

*   **Internal Core Subroutines:**
    *   **`SUBROUTINE apply_Pn(Ces, ncyc, pm, MASK, es, esrank)`**: The main worker routine for applying permutations.
        *   Selects cycle definitions (`CYC3` or `CYC4`) based on `ncyc` and `pm`.
        *   Iterates through permutations selected by `MASK`.
        *   For each component of the input eigenstate `es%lowest%eigen(esrank)%vec%rc`:
            *   Converts the integer representation of the state to a `fermion_ket_type`.
            *   Calls `Psz` (if Sz-adapted basis) or `Pupdo` (if up/down product basis) to apply the permutation to the ket.
            *   The `Psz` and `Pupdo` routines return a list of kets (`Pket`) because a permutation can map a single basis state to a sum of basis states if intermediate operations are involved (though for simple permutations it should be one-to-one).
            *   Accumulates the resulting permuted kets (with their fermionic signs) into the output `eigen_out%vec%rc`.
        *   Adds the final `eigen_out` to `Ces%lowest`.
    *   **`SUBROUTINE Pupdo(Pket, multi, spin, ket)`**: Applies a permutation defined by `multi` (an array of site indices in the cycle) to a `ket` in the up/down product basis, for a specific `spin`. It builds up the permutation by applying a sequence of `hop` operations (\(c^\dagger_j c_i\)) between the sites in the cycle, storing intermediate results in the `Pket` array. The `Pket` array on output contains components of the permuted state.
    *   **`SUBROUTINE Psz(Pket, multi, ket)`**: Applies a permutation `multi` to a `ket` in the Sz-adapted (Nambu) basis. It applies a sequence of `destroy` and `create` operations for both spin-up and spin-down components (using `IMPiorbsz` for orbital indices derived from `multi`) to achieve the permutation of particles among the specified sites.

## Important Variables/Constants

*   **Module-level arrays for cycle definitions:**
    *   `CYC3(nCYC3, 3, 2)`: Stores definitions of 3-site permutations. `CYC3(i,:,1)` is one permutation, `CYC3(i,:,2)` is its related one (e.g., inverse).
    *   `CYC4(nCYC4, 4, 2)`: Stores definitions of 4-site permutations similarly.
*   **Module-level pointer arrays for orbital mappings (initialized by `init_apply_P`):**
    *   `IMPiorbupdo(Nc, 2)`: Orbital indices for up/down basis.
    *   `IMPiorbsz(Nc, 2)`: Orbital indices for Sz-adapted basis.
*   **`force_reset_list :: LOGICAL, PARAMETER = .false.`**: If true, clears the target eigenstate list in `Ces` before adding results.

## Usage Examples

```fortran
MODULE example_apply_p_usage
  USE AIM_class, ONLY: AIM_type
  USE apply_P
  USE eigen_sector_class, ONLY: eigensector_type
  USE sector_class, ONLY: sector_type, delete_sector_list
  USE eigen_class, ONLY: delete_eigenlist

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: test_apply_p_module

CONTAINS
  SUBROUTINE test_apply_p_module(aim_model, source_eigenstate, source_rank)
    TYPE(AIM_type), INTENT(IN) :: aim_model
    TYPE(eigensector_type), INTENT(IN) :: source_eigenstate ! Assume populated
    INTEGER, INTENT(IN) :: source_rank       ! Rank of specific eigenvec

    TYPE(eigensector_type) :: result_eigenstate
    LOGICAL, ALLOCATABLE :: active_triplets_mask(:)
    INTEGER, DIMENSION(1,3) :: triplet_defs_example ! Example for one triplet

    ! Define a 3-cycle, e.g., sites 1 -> 2 -> 3 -> 1 (using 1-based indexing for sites)
    triplet_defs_example(1,:) = [1, 2, 3]

    ! Initialize the apply_P module (assuming no 4-cycles for this example)
    CALL init_apply_P(triplets=triplet_defs_example, &
     &                quadruplets=RESHAPE([0,0,0,0],[1,4]), AIM=aim_model) ! Pass empty or correctly shaped empty array for quadruplets

    ! Mask to apply only the first defined triplet
    IF (nCYC3 > 0) THEN
      ALLOCATE(active_triplets_mask(nCYC3))
      active_triplets_mask = .FALSE.
      active_triplets_mask(1) = .TRUE.
    ELSE
      PRINT *, "No 3-cycles defined in init_apply_P."
      RETURN
    END IF

    ! Determine the resulting sector (does not change for permutations)
    CALL P_sector(result_eigenstate%sector, pm='+', sector_in=source_eigenstate%sector)
    ! (Initialize result_eigenstate%lowest if needed)

    ! Apply the P(1->2->3) permutation
    CALL apply_P3(Ces=result_eigenstate, pm='+', MASK=active_triplets_mask, &
     &            es=source_eigenstate, esrank=source_rank)
    PRINT *, "Applied 3-cycle permutation. Result in result_eigenstate."
    IF (ASSOCIATED(result_eigenstate%lowest)) CALL delete_eigenlist(result_eigenstate%lowest)


    ! Example for apply_CHI (scalar spin chirality like)
    CALL apply_CHI(Ces=result_eigenstate, pm=' ', MASK=active_triplets_mask, &
     &             es=source_eigenstate, esrank=source_rank)
    PRINT *, "Applied CHI operator. Result in result_eigenstate."
    IF (ASSOCIATED(result_eigenstate%lowest)) CALL delete_eigenlist(result_eigenstate%lowest)


    CALL delete_sector_list(result_eigenstate%sector)
    IF (ALLOCATED(active_triplets_mask)) DEALLOCATE(active_triplets_mask)
    ! (Deallocate CYC3, CYC4, IMPiorb... if init_apply_P is not called again for new defs)

  END SUBROUTINE test_apply_p_module
END MODULE example_apply_p_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `genvar`: For `DP` (double precision kind).
    *   `aim_class`: The `init_apply_P` subroutine requires an `AIM_type` object for orbital mappings.
    *   `sector_class`: For `sector_type` and its utility functions (`delete_sector`, `new_sector`, `dimen_func`).
    *   `eigen_sector_class`: For `eigensector_type` (container for eigenstates) and related utilities like `new_eigensector`, `delete_eigensector`.
    *   `eigen_class`: For `eigen_type` (representing a single eigenstate vector) and its manipulation routines (`add_eigen`, `delete_eigen`, `delete_eigenlist`, `new_eigen`).
    *   `fermion_sector2_class`: For `rankupdo` (state indexing in up/down product basis).
    *   `fermion_ket_class`: For `fermion_ket_type`, `new_ket`, and crucially for implementing permutations: `hop` (particle hops for `Pupdo`), `create`/`destroy` (for `Psz`), and `copy_ket`.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules.

*   **Interactions with other components:**
    *   This module provides highly specialized operator application routines used within an ED solver, primarily for calculations involving site permutations.
    *   It is initialized via `init_apply_P` with definitions of 3-site and 4-site cycles and with `AIM_type` data for orbital context.
    *   It operates on `eigensector_type` objects, transforming input eigenstates.
    *   The core logic in `Pupdo` (for up/down product basis) applies permutations by a sequence of effective hopping operations (\(c^\dagger_j c_i\)).
    *   The logic in `Psz` (for Sz-adapted basis) applies permutations by sequences of creation and destruction operators for both spin components to achieve the particle rearrangement.
    *   The `apply_CHI` subroutine is a higher-level function that combines results from `apply_P3` (e.g., \(P\) and \(P^{-1}\)) to compute quantities like scalar spin chirality.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
