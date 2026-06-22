# muDock Compilation & Portability Issues Report

This report documents two primary build-time issues encountered when compiling the `muDock` application on Ubuntu with `Clang` (versions 17.0.6 and 18.1.3). Both issues have been resolved using standard-compliant workarounds that preserve the original logic and performance characteristics.

---

## Executive Summary

| Issue ID | Description | Severity | Impacted Files | Resolution Status |
| :--- | :--- | :--- | :--- | :--- |
| **ERR-001** | Compiler Crash (SIGSEGV, Exit Code 139) on OpenMP capture of structured bindings. | **Critical** (Blocks build) | [autodock_protein.cpp](file:///home/olly/UNI/progetto_aca/muDock/mudock/src/chem/autodock_protein.cpp) | **Resolved** via reference aliases. |
| **ERR-002** | Invalid template arguments (`this->`) inside template specializations. | **Major** (Blocks build) | 7 source files (CPP, GH, and XSimd implementations) | **Resolved** via direct static member lookup. |

---

## Issue ERR-001: Compiler Segmentation Fault on OpenMP Capture

### Description
During compilation, the compiler fails with an **Internal Compiler Error (ICE)**:
```
clang++-17: error: clang frontend command failed with exit code 139
```
The backtrace indicates a crash inside Clang's semantic analysis phase during OpenMP parsing:
```
clang::Sema::isOpenMPPrivateDecl(clang::ValueDecl*, unsigned int, unsigned int) const
clang::Sema::tryCaptureVariable(...)
```

### Root Cause Analysis
The crash is triggered by C++17 **structured bindings** used inside a loop parallelized with OpenMP (`#pragma omp parallel for collapse(3)`):
```cpp
const auto [vector1, vector2, exp, disorder] = compute_hbon_geometries(...);
#pragma omp parallel for collapse(3) schedule(static)
for (std::size_t index_z = 0; index_z < size_z; ++index_z) {
    // ...
    cos_theta = -diff.product(vector1[i]).sum_components(); // <--- Capture of 'vector1' structured binding triggers crash
}
```
Clang 17 and 18 have a known bug/limitation where capturing structured bindings inside OpenMP loops fails to resolve correctly, leading to an unhandled null pointer dereference or stack corruption (SIGSEGV) in the compiler frontend.

### Fix / Workaround
We manually unpack the struct instead of using C++17 structured bindings. The fields are accessed via standard reference bindings (`const auto&`), which are fully compatible with Clang's OpenMP analyzer:

```diff
     // find out the geometries of HBonds from the protein
-    const auto [vector1, vector2, exp, disorder] =
+    // const auto [vector1, vector2, exp, disorder] =
+    //     compute_hbon_geometries(x, y, z, this->get_base_molecule().get_num_hbond(), get_elements(), graph);
+    const auto geometries =
         compute_hbon_geometries(x, y, z, this->get_base_molecule().get_num_hbond(), get_elements(), graph);
+    const auto& vector1  = geometries.vector1;
+    const auto& vector2  = geometries.vector2;
+    const auto& exp      = geometries.exp;
+    const auto& disorder = geometries.disorder;
```

---

## Issue ERR-002: Invalid Template Arguments (`this->`)

### Description
Compiler rejects template instantiations of `invoke_kernel` with the following message:
```
error: no matching member function for call to 'invoke_kernel'
    q->invoke_kernel<this->geom_region_name>(geom_transform, ...
```

### Root Cause Analysis
The member `geom_region_name` is defined as a static constexpr member of the kernel classes:
```cpp
static constexpr char geom_region_name[] = "geometric_trasformation";
```
However, the code references it inside template parameters as `this->geom_region_name`. In standard C++, template arguments must be compile-time constant expressions. Because `this` is a runtime pointer, any expression containing `this->` cannot be evaluated at compile time, prompting compliant compilers (like Clang 17 and GCC) to reject it.

### Fix / Workaround
We resolved this across all 7 affected files by accessing the static constant directly by name (which compiles as a constant expression):

```diff
-q->invoke_kernel<this->geom_region_name>(geom_transform, ...
+q->invoke_kernel<geom_region_name>(geom_transform, ...
```

### Impacted Files & Lines:
1. **Geometric Transform:**
   - [geom_transform_cpp.cpp:L69](file:///home/olly/UNI/progetto_aca/muDock/mudock/src/cpp_implementation/geom_transform_cpp.cpp#L69)
   - [geom_transform_gh.cpp:L70](file:///home/olly/UNI/progetto_aca/muDock/mudock/src/gh_implementation/geom_transform_gh.cpp#L70)
   - [geom_transform_xsimd.cpp:L70](file:///home/olly/UNI/progetto_aca/muDock/mudock/src/xsimd_implementation/geom_transform_xsimd.cpp#L70)
2. **ADT Scoring:**
   - [adt_score_cpp.cpp:L212](file:///home/olly/UNI/progetto_aca/muDock/mudock/src/cpp_implementation/adt_score_cpp.cpp#L212)
   - [adt_score_gh.cpp:L397](file:///home/olly/UNI/progetto_aca/muDock/mudock/src/gh_implementation/adt_score_gh.cpp#L397)
   - [adt_score_xsimd.cpp:L367](file:///home/olly/UNI/progetto_aca/muDock/mudock/src/xsimd_implementation/adt_score_xsimd.cpp#L367)
3. **Genetic Algorithm:**
   - [genetic_cpp.cpp:L179](file:///home/olly/UNI/progetto_aca/muDock/mudock/src/cpp_implementation/genetic_cpp.cpp#L179) (Finalize)
   - [genetic_cpp.cpp:L190](file:///home/olly/UNI/progetto_aca/muDock/mudock/src/cpp_implementation/genetic_cpp.cpp#L190) (Iterate)
   - [genetic_cpp.cpp:L202](file:///home/olly/UNI/progetto_aca/muDock/mudock/src/cpp_implementation/genetic_cpp.cpp#L202) (Initialize)

---

> [!NOTE]
> All changes preserve compatibility with older C++20 compilers and maintain original OpenMP parallel loop optimizations (including collapse semantics).
