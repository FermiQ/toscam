# `src/utils/common_def.f90`

## Overview

*   **Purpose:** Provides common definitions and utility functions for the project.
*   **Role in Project:** Contains shared utility routines and data types used across various parts of the application.

## Key Components

*   `FUNCTION c2s(c) AS CHARACTER(LEN=SIZE(c))`: Converts an array of characters to a string.
*   `SUBROUTINE close_safe(unit_)`: Safely closes a file unit, checking if it was opened.
*   `SUBROUTINE create_seg_fault()`: Intentionally attempts to create a segmentation fault (likely for testing memory access issues).
*   `SUBROUTINE dump_message(UNIT, TEXT)`: Writes a message to a specified output unit (defaults to `log_unit`) or increments a counter if no text is provided.
*   `FUNCTION find_rank(intg, tabintg) AS INTEGER`: Finds the first index of an integer `intg` within an integer array `tabintg`. Returns 0 if not found or if the array is empty.
*   `INTERFACE get_address`: Provides functions to get the memory address of different data types.
    *   `FUNCTION get_addressr(d1) AS INTEGER(8)`: Gets address of a single precision real.
    *   `FUNCTION get_addressi(d1) AS INTEGER(8)`: Gets address of a 4-byte integer.
    *   `FUNCTION get_addressl(d1) AS INTEGER(8)`: Gets address of a logical variable.
    *   `FUNCTION get_addresschar(d1) AS INTEGER(8)`: Gets address of the first character of a string.
    *   `FUNCTION get_address_d(d1) AS INTEGER(8)`: Gets address of a double precision real.
    *   `FUNCTION get_address_dvec(d1) AS INTEGER(8)`: Gets address of the first element of a double precision real vector.
    *   `FUNCTION get_address_c(d1) AS INTEGER(8)`: Gets address of a double precision complex number.
    *   `FUNCTION get_address_cvec(d1) AS INTEGER(8)`: Gets address of the first element of a double precision complex vector.
    *   `FUNCTION get_address_ivec(d1) AS INTEGER(8)`: Gets address of the first element of a 4-byte integer vector.
*   `FUNCTION i2c(i) AS CHARACTER, POINTER :: i2c(:)`: Converts an integer to a dynamically allocated array of characters.
*   `SUBROUTINE open_safe(unit_, filename_, status_, action_, get_unit)`: Safely opens a file with specified status and action, optionally getting a free unit. Includes error checking.
*   `INTERFACE put_address`: Provides subroutines to put data at a given memory address for different data types.
    *   `SUBROUTINE put_addressr(j, d1)`: Puts a single precision real at address `j`.
    *   `SUBROUTINE put_addressi(j, d1)`: Puts a 4-byte integer at address `j`.
    *   `SUBROUTINE put_addressdvec(j, d1)`: Puts a vector of double precision reals starting at address `j`.
    *   `SUBROUTINE put_addresscvec(j, d1)`: Puts a vector of double precision complex numbers starting at address `j`.
    *   `SUBROUTINE put_addressivec(j, d1)`: Puts a vector of 4-byte integers starting at address `j`.
    *   `SUBROUTINE put_addresschar(j, d1)`: Puts characters of a string starting at address `j`.
    *   `SUBROUTINE put_addressd(j, d1)`: Puts a double precision real at address `j`.
    *   `SUBROUTINE put_addressc(j, d1)`: Puts a double precision complex number at address `j` (real and imaginary parts separately).
    *   `SUBROUTINE put_addressl(j, d1)`: Puts a logical variable at address `j`.
*   `SUBROUTINE reset_timer(clock_ref)`: Initializes or resets a timer by getting the current system clock count.
*   `SUBROUTINE skip_line(unit_, NLINES)`: Reads and discards a specified number of lines (default 1) from an input unit.
*   `SUBROUTINE stats_func(TEXT)`: Prints system memory usage statistics (from `vmstat`) and optional text messages. (Note: `vmstat` call is platform-dependent).
*   `SUBROUTINE system_call(cmd)`: Executes a system command, redirecting its output to the log file. (Note: platform-dependent, uses `system`, `cat`, `rm`, `fseek`).
*   `SUBROUTINE timer_fortran(clock_ref, TEXT, unit_)`: Calculates and prints the elapsed time since `clock_ref`, with an optional descriptive text, to a specified unit or `log_unit`.
*   `SUBROUTINE utils_abort(error_message, filename)`: Terminates the program using `error stop` and prints an error message to a specified file or standard output.
*   `SUBROUTINE utils_assert(assertion, error_message)`: Checks a logical assertion. If false, it calls `utils_abort` with the given error message.
*   `SUBROUTINE utils_system_call(command, abort, report)`: Executes a system command using `execute_command_line`. Optionally aborts on non-zero exit status and/or reports the command being executed.
*   `FUNCTION utils_unit() AS INTEGER`: Finds and returns an available (closed and existing) Fortran I/O unit number between 10 and 99.
*   `INTERFACE utils_qc_print`: Provides procedures for printing Quality Control (QC) information to standard output in a formatted way.
    *   `SUBROUTINE utils_qc_print_integer(tag, value)`: Prints an integer value with a tag.
    *   `SUBROUTINE utils_qc_print_real(tag, value)`: Prints a double precision real value with a tag.
    *   `SUBROUTINE utils_qc_print_string(tag, value)`: Prints a string value with a tag.

## Important Variables/Constants

*   The module `use genvar` which likely defines `DP` (double precision kind parameter, used in `utils_qc_print_real`, `create_seg_fault`, `put_addressdvec`, `put_addresscvec`, `put_addressd`, `put_addressc`, `get_address_d`, `get_address_dvec`, `get_address_c`, `get_address_cvec`, `timer_fortran`), `SP` (single precision kind parameter, used in `put_addressr`, `get_addressr`), `logfile` (likely default log unit, used in `utils_qc_print_*`), `log_unit` (used in `dump_message`, `timer_fortran`, `system_call`), `one_over_clock_rate` (used in `timer_fortran`), and `iproc` (used in `system_call`).

## Usage Examples

Refer to calling routines for usage. Examples to be added by developers.

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   `use genvar`: Imports various global variables and kind parameters such as `DP`, `SP`, `logfile`, `log_unit`, `one_over_clock_rate`, `iproc`.
*   **External Libraries:**
    *   Standard Fortran intrinsic modules (e.g., for `SYSTEM_CLOCK`, `execute_command_line`, `flush`, `ADJUSTL`, `LEN_TRIM`, `SIZE`, `REAL`, `AIMAG`, `ALLOCATE`).
    *   Relies on system calls available in the environment for `utils_system_call` (via `execute_command_line`) and `system_call` (via `system`, `cat`, `rm`, `fseek`), and `stats_func` (via `vmstat`). The availability and behavior of these system commands are platform-dependent.
*   **Interactions with other components:**
    *   Provides utility functions (error handling, file I/O, timing, system interaction, string manipulation, QC printing) that are likely called from many other modules in the project.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
