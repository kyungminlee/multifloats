# Proposed Naming Scheme for multifloats.fypp

Based on the current implementation in `@fsrc/multifloats.fypp`, the codebase already has a strong foundation for a naming scheme, primarily using a **Type-First (pseudo-Object-Oriented)** approach. However, there are a few inconsistencies—such as numeric arity suffixes (`_1`, `_2`), mixed use of Fortran intrinsic names vs. internal kind names (`dble` vs `dp`), and messy internal C-delegate names (`_full`, `c_dd_cdd_`).

Here is a proposed **stable, pretty, and consistent** naming scheme that standardizes the existing patterns and cleans up the rough edges.

## The Core Convention: Type Abbreviations
First, establish a strict vocabulary for types across all function names:
* **Custom types:** `dd` (real64x2), `cdd` (cmplx64x2)
* **Real kinds:** `dp`, `sp`, `qp`
* **Complex kinds:** `cdp`, `csp`, `cqp`
* **Integer kinds:** `int` (default), `i8`, `i16`, `i64`

---

## 1. Binary Operators (`+`, `-`, `*`, `/`, `**`, `==`, `<`)
**Pattern:** `{lhs_type}_{operation}_{rhs_type}`

The current scheme here is already perfect. It reads naturally from left to right.
* `dd_add_dp` (dd + dp)
* `cdd_mul_dd` (cdd * dd)
* `dd_eq_int` (dd == int)
* `dd_pow_i32` (dd ** int32)

## 2. Unary Operators & Math Functions
**Pattern:** `{arg_type}_{function}`

Also mostly perfect, but needs cleanup regarding the `_full` suffix currently used to avoid naming collisions with C-bindings.
* `dd_sin`, `cdd_exp`, `dd_abs`, `dd_neg`
* **Refinement:** Drop the `_full` suffix (e.g., `dd_sin_full`). Instead, rename the underlying C-binding functions to create a clean namespace boundary (see section 7) or inline the C-calls directly into the public routine.

## 3. Constructors (Building from scalars)
**Pattern:** `{result_type}_from_{arg1_type}[_{arg2_type}]`

**Current Flaw:** Using `_1` and `_2` to denote the number of arguments (e.g., `cdd_from_dp_2`, `cdd_from_dp_1`), which hides the argument types.
**Proposed Fix:** Explicitly state the argument types. If the function takes two arguments, list both abbreviations.
* `dd_from_dp` (1 arg)
* `cdd_from_dp_dp` (instead of `cdd_from_dp_2`)
* `cdd_from_dd_dp` (instead of `cdd_from_dd_dp_2`)
* `cdd_from_dp` (instead of `cdd_from_dp_1` — implicit zero imaginary part)

## 4. Type Casting (Overloading `dble`, `int`, `cmplx`, `real`)
**Pattern:** `{src_type}_to_{target_type}`

**Current Flaw:** Mixing intrinsic names with kind names (e.g., `dd_dble`, `dd_int`, `cdd_real`, `dd_cmplx_2`) and creating redundant functions for casting vs constructing.
**Proposed Fix:** Use the established type abbreviations for both the source and the target. Furthermore, **unify casts and constructors** where functionally identical to eliminate redundant code.
* `dd_to_dp` (unifies `dd_dble` and `dd_to_double`)
* `dd_to_int` (instead of `dd_int`)
* `cdd_to_dd` (instead of `cdd_real`)
* For complex creation (e.g. `interface cmplx`), do not create a separate `dd_dd_to_cdd` or `dd_cmplx_2`. Simply expose the explicitly typed constructor `cdd_from_dd_dd`.

## 5. Assignment (`=`)
**Pattern:** `{lhs_type}_assign_{rhs_type}`

This is consistent with the binary operators and reads as "LHS assigned from RHS".
* `dd_assign_dp` (dd = dp)
* `dp_assign_cdd` (dp = cdd)

## 6. Array Reductions & Transformational Functions
**Pattern:** `{type}_{function}_{rank}d[_dim]`

The current scheme is excellent and should be kept exactly as is.
* `dd_maxval_1d`
* `dd_sum_2d_dim`
* `cdd_findloc_3d`

## 7. Internal C-ABI Bindings (The cleanup)
**Pattern:** `c_{NAME}` where `NAME` is the exact C function name from the ABI.

**Current Flaw:** The C bindings are inconsistently named (`c_sindd`, `c_dd_cdd_sqrt`, `c_dd_jn`), which forced the use of the awkward `_full` suffix on the Fortran side.
**Proposed Fix:** Use `c_` as a universal "C interop" prefix, immediately followed by the exact C symbol name. This provides zero-friction traceability to the C++ source and perfectly leverages existing C naming conventions (where `c` often prefixes complex math functions).
* `c_sindd` (binds to C symbol `sindd`)
* `c_csindd` (binds to C symbol `csindd` for complex sine)
* `c_jndd` (binds to C symbol `jndd`)
* `c_matmuldd_mm` (binds to C symbol `matmuldd_mm`)

---

## Summary of Refactorings Required in `multifloats.fypp`:
1. Rename `cdd_from_*_1` to `cdd_from_*`.
2. Rename `cdd_from_*_2` to `cdd_from_*_*` (e.g., `cdd_from_dp_dp`).
3. Replace intrinsic names in casts with type abbreviations (`dd_dble` → `dd_to_dp`, `cdd_real` → `cdd_to_dd`), and consolidate redundant casting functions like `dd_to_double` into `dd_to_dp`.
4. Unify identical casts and constructors by directly exposing the explicit constructor (e.g., `cdd_from_dd_dd`) in standard interfaces like `interface cmplx`.
5. Rename `c_...` C-interfaces to `c_{C_SYMBOL}` exactly (e.g. `c_csindd`).
6. Remove `_full` from the elemental wrappers (e.g., `dd_sin_full` → `dd_sin` since the C binding is now safely out of the way as `c_sindd`).