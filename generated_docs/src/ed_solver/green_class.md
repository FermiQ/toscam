# `src/ed_solver/green_class.f90`

## Overview

*   **Purpose:** This module defines a Fortran derived type, `green_type`, and associated procedures. The `green_type` is designed to encapsulate a complete Green's function, which may consist of several distinct components (e.g., particle-hole, particle-particle, hole-hole channels, or different time orderings).
*   **Role in Project:** It serves as a higher-level container that groups multiple `correl_type` objects (from `correl_class.f90`), each representing a specific component of a Green's function. This module provides functionalities for creating, managing, copying, deleting, and performing operations like I/O, padding (filling by symmetry), or slicing on the entire set of these related correlation components that constitute a full Green's function.

## Key Components

*   **`TYPE green_type`**:
    A Fortran derived type to represent a multi-component Green's function.
    *   `INTEGER :: N = 0`: The matrix dimension of the Green's function components (e.g., number of orbitals or sites).
    *   `INTEGER :: Nw = 0`: The number of frequency points at which the Green's function is defined.
    *   `CHARACTER(LEN=100) :: title = '\0'`: A descriptive title for the Green's function (e.g., "Impurity_GF", "Nambu_GF").
    *   `LOGICAL :: compute(2, 2) = RESHAPE((/.false., .false., .false., .false./), (/2, 2/))`: A 2x2 logical mask. `compute(ipm, jpm)` indicates whether the corresponding `correl(ipm, jpm)` component is actively computed or considered primary. The indices `ipm` and `jpm` (1 or 2) typically map to particle/hole or creation/annihilation operator combinations (e.g., (1,1) for \(\langle c c^\dagger \rangle\), (1,2) for \(\langle c c \rangle\), (2,1) for \(\langle c^\dagger c^\dagger \rangle\), (2,2) for \(\langle c^\dagger c \rangle\)). The default initialization in `new_green_*` routines often sets `compute(1,1)=.FALSE.` and `compute(1,2)=.FALSE.`, implying these components might be derived from `compute(2,2)` and `compute(2,1)` by symmetry.
    *   `TYPE(correl_type) :: correl(2, 2)`: A 2x2 array of `correl_type` objects. Each element `correl(ipm, jpm)` stores a specific component of the Green's function matrix as a function of frequency.
    *   `TYPE(masked_matrix_type) :: correlstat(2, 2)`: A 2x2 array of `masked_matrix_type` objects, intended to store the static (equal-time) values of the corresponding Green's function components.
    *   `COMPLEX(DP) or REAL(DP), POINTER :: Amean(:, :) => NULL(), Bmean(:, :) => NULL()`: Pointers that can store mean values of operators, potentially \(\langle A \rangle\) and \(\langle B \rangle\) if the Green's function represents a correlation like \(\langle T A(z) B(0) \rangle\). The type depends on whether `_complex` is defined.

*   **`INTERFACE new_green`**: A generic interface for creating `green_type` objects.
    *   `MODULE PROCEDURE new_green_rfreq(GREEN, compute, title, N, Nw, wmin, wmax, width, STAT, IMASK, AB, force_compute)`: Creates a `green_type` object for real frequencies. It initializes the four `correl(ipm,jpm)` components using `new_correl` from `correl_class` with real frequency parameters.
    *   `MODULE PROCEDURE new_green_ifreq(GREEN, compute, title, N, Nw, beta, STAT, IMASK, AB, force_compute)`: Creates a `green_type` object for Matsubara (imaginary) frequencies. It initializes the `correl(ipm,jpm)` components for Matsubara frequencies.
        *   Both constructors take `compute(2,2)` and optional `force_compute(2,2)` logical arrays to determine which of the four `correl` components are marked for direct computation. By default, they set `compute(1,1)=.FALSE.` and `compute(1,2)=.FALSE.`, implying these components might be derived from `(2,2)` and `(2,1)` using symmetries.

*   **Public Management and Utility Subroutines:**
    *   **`SUBROUTINE copy_green(GREENOUT, GREENIN)`**: Performs a deep copy of a `green_type` object from `GREENIN` to `GREENOUT`.
    *   **`SUBROUTINE delete_green(GREEN)`**: Destructor for `green_type` objects. It calls `delete_correl` for each active `correl(ipm,jpm)` component and `delete_masked_matrix` for each `correlstat(ipm,jpm)`.
    *   **`SUBROUTINE pad_green(GREENOUT, GREENIN)`**: "Pads" or fills in potentially uncomputed components of `GREENOUT` using symmetry relations based on computed components, possibly from `GREENIN` if provided. This is useful if only a minimal set of components were directly calculated. The symmetry relation mentioned is \(\langle a(z)b \rangle = (-1)^F \langle A(-z^*)B \rangle^*\), where F indicates fermionic/bosonic nature.
    *   **`SUBROUTINE write_green(GREEN)`**: Writes all active (where `GREEN%compute` or its symmetric counterpart is true) `correl` components of the `GREEN` object to appropriately named files by calling `write_correl` from `correl_class.f90`.
    *   **`SUBROUTINE glimpse_green(GREEN, UNIT)`**: Prints a "glimpse" (typically the matrix at a specific frequency `iwprint`) of all active `correl` components to a specified Fortran `UNIT`.
    *   **`SUBROUTINE slice_green(GREENOUT, GREENIN, rbounds, cbounds)`**: Extracts a sub-matrix (slice defined by `rbounds` and `cbounds`) from each active `correl` component of `GREENIN` and stores it in the corresponding component of `GREENOUT`.

## Important Variables/Constants

*   **Within `green_type`:**
    *   `correl(2,2) :: TYPE(correl_type)`: This 2x2 array is the core of `green_type`, holding the different components of the Green's function (e.g., \(G_{<c,c^\dagger>}\), \(G_{<c,c>}\), etc.), each as a `correl_type` object that includes frequency dependence and matrix structure.
    *   `correlstat(2,2) :: TYPE(masked_matrix_type)`: Stores the static (equal-time) counterparts of the `correl` components.
    *   `compute(2,2) :: LOGICAL`: A mask indicating which of the four `correl(ipm,jpm)` components are considered primary or are actively computed. This allows for efficiency by computing a minimal set and deriving others via symmetry (using `pad_green`).
    *   `N :: INTEGER`: The matrix dimension of each component.
    *   `Nw :: INTEGER`: The number of frequency points.
    *   `title :: CHARACTER(LEN=100)`: A name for this Green's function object.

## Usage Examples

```fortran
MODULE example_green_class_usage
  USE green_class
  USE correl_class, ONLY: correl_type ! For direct access to components if needed
  USE frequency_class, ONLY: freq_type ! For frequency details of components
  USE genvar, ONLY: FERMIONIC, DP      ! For statistics

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: manage_greens_function

CONTAINS
  SUBROUTINE manage_greens_function
    TYPE(green_type) :: electron_gf
    LOGICAL :: gf_component_flags(2,2)
    INTEGER :: num_orbitals_sys, num_freq_points_sys
    REAL(DP) :: beta_temperature_sys

    num_orbitals_sys = 3
    num_freq_points_sys = 512
    beta_temperature_sys = 20.0_DP

    ! Define which components of the Green's function are to be explicitly handled/computed.
    ! For example, compute all components explicitly for demonstration.
    gf_component_flags = .TRUE.

    ! Create a new Green's function object for Matsubara frequencies
    CALL new_green(electron_gf, compute=gf_component_flags, title="ElectronGF_Matsubara", &
     &             N=num_orbitals_sys, Nw=num_freq_points_sys, beta=beta_temperature_sys, &
     &             STAT=FERMIONIC)
    ! Note: The interface 'new_green' will dispatch to 'new_green_ifreq'.

    PRINT *, "Initialized Green's function: ", TRIM(electron_gf%title)
    PRINT *, "Matrix dimension (N): ", electron_gf%N
    PRINT *, "Number of frequency points (Nw): ", electron_gf%Nw
    PRINT *, "Active components (compute flags): ", electron_gf%compute

    ! At this point, electron_gf%correl(ipm,jpm)%fctn would be allocated and
    ! could be filled by solver routines.
    ! For example, electron_gf%correl(1,1)%fctn would be G_-- (<c c_dag>)
    ! electron_gf%correl(2,2)%fctn would be G_++ (<c_dag c>) (or G_<c_dag c> if following default)

    ! After computation and filling, one might write it out:
    ! CALL write_green(electron_gf)

    ! Clean up the Green's function object
    CALL delete_green(electron_gf)

  END SUBROUTINE manage_greens_function
END MODULE example_green_class_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `correl_class`: This is a fundamental dependency. `green_type` is essentially a structured container for four `correl_type` objects. Most operations on `green_type` (creation, I/O, padding, slicing) delegate to the corresponding routines in `correl_class` for each active component.
    *   `genvar`: For `DP` (double precision kind), `log_unit` (logging), `messages3` (debug message control), and `pm` (character array `['+','-']` used in naming components).
    *   `masked_matrix_class`: Used for the `correlstat` components (which are `masked_matrix_type`) and via `correl_class` for the `MM` component of `correl_type`.
    *   `frequency_class`: Used via `correl_class`, as each `correl_type` object has a `freq_type` component.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules.

*   **Interactions with other components:**
    *   **`correlations.f90`**: This higher-level module is a primary user of `green_class`. It declares public module variables of `green_type` (e.g., `G(2)`, `GF(2)`, `GNAMBU`, `GNAMBUN`, `GNAMBUret`, `Spm`) which are initialized using `new_green` and then populated by its computational routines.
    *   **Solver routines (within `correlations_dyn_stat.h` or similar):** These routines calculate the matrix elements and frequency dependence of different Green's function components and would fill the `fctn` arrays of the `correl_type` objects held within a `green_type` instance.
    *   The `compute(2,2)` logical mask within `green_type` allows for optimization. For example, only \(\langle T c c^\dagger \rangle\) might be directly computed, and \(\langle T c^\dagger c \rangle\) could be derived using symmetry relations via `pad_green`, reducing computational effort.
    *   The `correlstat(2,2)` components provide a place to store equal-time (static) values, which can be derived from the dynamic Green's functions or computed directly.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
