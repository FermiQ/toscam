# `src/ed_solver/fermion_ket_class.f90`

## Overview

*   **Purpose:** This module defines a Fortran derived type, `fermion_ket_type`, to represent fermionic Fock states (often referred to as "kets" in quantum mechanics) using an integer bitwise representation. It provides a suite of core subroutines for creating and manipulating these states, such as applying creation (\(c^\dagger\)) and destruction (\(c\)) operators, performing particle hops, checking orbital occupation, and critically, calculating the fermionic sign changes associated with these operations.
*   **Role in Project:** Serves as a low-level utility crucial for state manipulation within an Exact Diagonalization (ED) solver or any framework dealing with many-fermion systems. It allows other modules (e.g., those that apply operators to states, like `apply_c.f90`, or those that construct Hamiltonian matrix elements) to work with fermionic states in an efficient integer representation while correctly handling the anti-commutation relations of fermion operators through meticulous sign calculations.

## Key Components

*   **`TYPE fermion_ket_type`**:
    Represents a single many-fermion Fock state.
    *   `INTEGER :: norbs`: The total number of single-particle orbitals available in the basis for this ket. This is essential for determining the range of bitwise operations and for fermionic sign calculations.
    *   `INTEGER :: state`: The integer representation of the Fock state. Typically, if bit `i-1` (0-indexed) of this integer is set (1), orbital `i` (1-indexed) is considered occupied; otherwise, it's empty.
    *   `INTEGER :: fermion_sign`: Stores the accumulated fermionic sign, which is typically +1 or -1. This sign is updated when creation or destruction operators are applied, reflecting the anti-commutation rules.
    *   `LOGICAL :: is_nil`: A flag that is `.TRUE.` if the ket represents a null state (e.g., after attempting to annihilate a particle from an already empty orbital or create a particle in an already occupied orbital according to Pauli exclusion).

*   **`TYPE fermion_ket_vec`**: (Not public from the module)
    Likely intended to represent a state vector as a linear combination of `fermion_ket_type` basis states, primarily for internal calculations or examples shown within the module (like `build_and_diag_hamiltonian_matrix`).
    *   `INTEGER :: dim`: Dimension of the vector (number of basis kets).
    *   `REAL(DP) :: energy, total_S, total_Sz, total_N`: Physical properties that could be associated with such a state vector.
    *   `REAL(DP), ALLOCATABLE :: vec(:)`: Array of coefficients for the linear combination of basis kets.

*   **`INTERFACE new_ket`**: A generic interface for creating `fermion_ket_type` objects.
    *   `MODULE PROCEDURE new_ket_from_state(ket, state, norbs)`: The primary constructor. It initializes a `fermion_ket_type` `ket` with a given integer `state` representation and `norbs`. The `fermion_sign` is set to 1, and `is_nil` to `.FALSE.`.
    *   `MODULE PROCEDURE new_ket_from_old(ket_out, ket_in)`: Creates `ket_out` as a copy of `ket_in` (effectively calls `copy_ket`).

*   **Public Manipulation Subroutines & Functions:**
    *   **`SUBROUTINE copy_ket(ket_out, ket_in)`**: Performs a direct copy of one `fermion_ket_type` (`ket_in`) to another (`ket_out`).
    *   **`SUBROUTINE create(ket_out, iorb, ket_in)`**: Applies the creation operator \(c^\dagger_{iorb}\) to `ket_in`, storing the result in `ket_out`. If orbital `iorb` is already occupied in `ket_in`, `ket_out%is_nil` is set to `.TRUE.`. Otherwise, the corresponding bit in `ket_out%state` is set, and `ket_out%fermion_sign` is updated by calling the internal `Fermion_sign` routine.
    *   **`SUBROUTINE destroy(ket_out, iorb, ket_in)`**: Applies the destruction operator \(c_{iorb}\) to `ket_in`, result in `ket_out`. If orbital `iorb` is empty in `ket_in`, `ket_out%is_nil` is set. Otherwise, the bit is cleared, and `fermion_sign` is updated.
    *   **`SUBROUTINE create_pair(ket_out, jorb, iorb, ket_in)`**: Applies \(c^\dagger_{jorb} c^\dagger_{iorb}\) to `ket_in` by sequential calls to `create`.
    *   **`SUBROUTINE destroy_pair(ket_out, jorb, iorb, ket_in)`**: Applies \(c_{jorb} c_{iorb}\) to `ket_in` by sequential calls to `destroy`.
    *   **`SUBROUTINE hop(ket_out, jorb, iorb, ket_in)`**: Applies the hopping operator \(c^\dagger_{jorb} c_{iorb}\) (moves a particle from `iorb` to `jorb`) by calling `destroy` then `create`.
    *   **`SUBROUTINE Fermion_sign(jsign, ket, iorb)`**: (This is an internal helper but fundamental to the public routines). Updates an input `jsign` based on the occupations of orbitals with index greater than `iorb` in the `ket%state`. `jsign` is multiplied by -1 for each occupied orbital `jorb > iorb`. *Note: The standard convention for the fermionic sign when an operator acts at `iorb` usually involves summing occupations of orbitals with index *less than* `iorb`. This implementation (summing `jorb > iorb`) calculates the sign associated with moving the operator at `iorb` to the "end" of the canonical ordering of operators.*
    *   **`FUNCTION is_occupied(iorb, ket) AS LOGICAL`**: Returns `.TRUE.` if orbital `iorb` (1-indexed) is occupied in `ket` (checks bit `iorb-1` of `ket%state` using `BTEST`).
    *   **`FUNCTION ni__(iorb, ket) AS INTEGER`**: Returns 1 if orbital `iorb` is occupied in `ket`, and 0 otherwise.
    *   **`FUNCTION noccupied(ket, nnn) AS INTEGER`**: Returns the total number of occupied orbitals in `ket`. If `nnn` is present, it counts up to orbital `nnn`.

*   **Private Helper/Analysis Functions (Illustrative, not part of public API):**
    The module also contains several private routines like `locate_fermion`, `build_and_diag_hamiltonian_matrix`, `build_basis_of_hilbert_space`, `apply_cdagger_dn_and_cdagger_up_on_state`, and functions to calculate total spin, particle number, Sz, and diagonal energy for a given `fermion_ket_type`. These appear to be for internal testing, examples, or perhaps for use in very small, specific sub-problems rather than general public use.

## Important Variables/Constants

*   **Within `fermion_ket_type`:**
    *   `state :: INTEGER`: The integer whose bit pattern represents the Fock state.
    *   `fermion_sign :: INTEGER`: Accumulates the sign changes due to fermionic anti-commutation. Crucial for correct matrix element calculations.
    *   `is_nil :: LOGICAL`: Indicates if an operation resulted in a null state (e.g., \(c_i |0\rangle_i\)).
    *   `norbs :: INTEGER`: Defines the number of orbitals, setting the context for bitwise operations and sign calculations.

## Usage Examples

```fortran
MODULE example_fermion_ket_usage
  USE fermion_ket_class
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: demonstrate_ket_operations

CONTAINS
  SUBROUTINE demonstrate_ket_operations
    TYPE(fermion_ket_type) :: ket_A, ket_B_annihilate, ket_C_create, ket_D_hop
    INTEGER, PARAMETER :: total_orbitals = 4
    INTEGER :: orbital_idx_1 = 1  ! 1-indexed orbital
    INTEGER :: orbital_idx_2 = 2
    INTEGER :: orbital_idx_3 = 3
    INTEGER :: initial_state_int ! e.g., |1010> (orbitals 1 and 3 occupied)

    ! For |1010> with norbs=4: (orbital 4 is MSB, orbital 1 is LSB for this example)
    ! Binary: 1*2^3 (orb4) + 0*2^2 (orb3) + 1*2^1 (orb2) + 0*2^0 (orb1) = 8 + 2 = 10
    ! If LSB is orbital 1: 0*2^3 (orb4) + 1*2^2 (orb3) + 0*2^1 (orb2) + 1*2^0 (orb1) = 4 + 1 = 5
    ! The code uses BTEST(state, iorb-1), so iorb=1 is bit 0 (LSB).
    ! So, |1010> means orbital 1 (bit 0) occupied, orbital 3 (bit 2) occupied. state = 1*2^0 + 1*2^2 = 1 + 4 = 5.
    initial_state_int = 5

    CALL new_ket(ket_A, state=initial_state_int, norbs=total_orbitals)
    PRINT *, "Initial Ket A: State=", ket_A%state, " (Occupation: ", &
     &        ni__(1, ket_A), ni__(2, ket_A), ni__(3, ket_A), ni__(4, ket_A), &
     &        "), Sign=", ket_A%fermion_sign

    ! 1. Annihilate a particle from orbital_idx_1 (which is occupied)
    CALL destroy(ket_B_annihilate, orbital_idx_1, ket_A)
    IF (ket_B_annihilate%is_nil) THEN
      PRINT *, "Annihilation on orbital ", orbital_idx_1, " resulted in nil state."
    ELSE
      PRINT *, "Ket B (c_1 |A>): State=", ket_B_annihilate%state, " (Occupation: ", &
       &        ni__(1, ket_B_annihilate), ni__(2, ket_B_annihilate), ni__(3, ket_B_annihilate), ni__(4, ket_B_annihilate), &
       &        "), Sign=", ket_B_annihilate%fermion_sign
    END IF

    ! 2. Create a particle on orbital_idx_2 (which is currently empty) in Ket A
    CALL create(ket_C_create, orbital_idx_2, ket_A)
    IF (ket_C_create%is_nil) THEN
      PRINT *, "Creation on orbital ", orbital_idx_2, " resulted in nil state."
    ELSE
      PRINT *, "Ket C (c^+_2 |A>): State=", ket_C_create%state, " (Occupation: ", &
       &        ni__(1, ket_C_create), ni__(2, ket_C_create), ni__(3, ket_C_create), ni__(4, ket_C_create), &
       &        "), Sign=", ket_C_create%fermion_sign
    END IF

    ! 3. Hop a particle from orbital_idx_3 to orbital_idx_4 in Ket C (state should be |1011> -> |1110>)
    ! Ket C state is |1110> (state=7), occupied: 1,2,3. Hopping from 3 to 4 (empty).
    ! destroy(orb3) on |1110> -> state |1010> (5), sign depends on orb1,2 count before 3.
    ! create(orb4) on |1010> -> state |1010> | (1<<3) = |1000|0101> = |1101> (13)
    IF (.NOT. ket_C_create%is_nil) THEN
      CALL hop(ket_D_hop, orbital_idx_4, orbital_idx_3, ket_C_create)
      IF (ket_D_hop%is_nil) THEN
        PRINT *, "Hop from ", orbital_idx_3, " to ", orbital_idx_4, " resulted in nil state."
      ELSE
        PRINT *, "Ket D (c^+_4 c_3 |C>): State=", ket_D_hop%state, " (Occupation: ", &
         &        ni__(1, ket_D_hop), ni__(2, ket_D_hop), ni__(3, ket_D_hop), ni__(4, ket_D_hop), &
         &        "), Sign=", ket_D_hop%fermion_sign
      END IF
    END IF

  END SUBROUTINE demonstrate_ket_operations
END MODULE example_fermion_ket_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `genvar`: For `DP` (double precision kind definition).
    *   `matrix` (used in private, non-public helper routines): For `eigenvector_matrix`, `write_array`.
    *   `linalg` (used in private, non-public helper routines): For `scalprod`.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules, particularly for bitwise operations:
        *   `BTEST(integer, bit_pos)`: Tests if a bit is set.
        *   `IBSET(integer, bit_pos)`: Sets a bit to 1.
        *   `IBCLR(integer, bit_pos)`: Clears a bit to 0.
        *   `SUM` with array constructor for `noccupied`.

*   **Interactions with other components:**
    *   This module provides the fundamental tools for representing and algebraically manipulating fermionic Fock states at a low level.
    *   **Operator Application Modules (e.g., `apply_c.f90`, `apply_ns.f90`, `apply_p.f90`):** These modules are primary consumers. They would use `fermion_ket_type` and its associated routines (`create`, `destroy`, `hop`, `ni__`, etc.) extensively. When these modules apply an operator to a many-body eigenstate (which is a superposition of Fock states), they iterate through the constituent kets, apply the operation using this class's methods to each ket, and then construct the resulting state.
    *   **Hamiltonian Construction Modules (e.g., `h_class.f90`):** When building the Hamiltonian matrix, routines would use this class to take a basis ket, apply Hamiltonian terms (which are composed of creation/destruction operators), and identify the resulting ket and matrix element (which includes the `fermion_sign`).
    *   **Hilbert Space Management (e.g., `fermion_hilbert_class.f90`):** While `fermion_hilbert_class` manages collections of states and their organization into sectors, it uses the integer representation of states that `fermion_ket_class` manipulates. For example, `fermion_hilbert_class%new_fermion_Hilbert_space` uses `noccupied` (from this module) to determine the particle number of a state.
    *   The integer representation of Fock states is a common and efficient technique in ED codes, allowing for compact storage and fast bitwise operations to simulate fermion creation and annihilation. The correct calculation and propagation of `fermion_sign` by this module is critical for ensuring that all calculations respect fermionic anti-commutation relations.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
