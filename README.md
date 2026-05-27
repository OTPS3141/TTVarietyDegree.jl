# TTVarietyDegree

Exact computation of degrees of tensor train varieties
`V^{<=}_{D,d}` using the recursive Schur-Weingarten procedure.

The implementation is pure Julia and uses `Rational{BigInt}` throughout.
It does not hardcode the sample cases: the tests compute the recursive tail
polynomial, the final functional `f(P)`, and the degree from the input data.

## Installation

For local development, run this once from any Julia session:

```julia
using Pkg
Pkg.develop(path="/home/user/pCloudDrive/pCloud Backup/ENVY-OTPS/Documents/TENORS/Code/degree_TT_varieties/degree_random_alg_geom")
```

For a public GitHub repository, use the repository URL:

```julia
using Pkg
Pkg.add(url="https://github.com/OTPS3141/TTVarietyDegree.jl")
```

Replace the URL with the actual repository once it is pushed.

## Julia API

The main package API is:

```julia
degree_TT_variety(D, d; method=:schur_weingarten, reduce=true)
dimension_TT_variety(D, d)
compute_tail_polynomial(D, d; verbose=false)
compute_fP(P_coeffs, k, n)
f_schur(lam, k, n)
compress_vacuous_boundary_modes(D, d)
degree_subspace_variety_N3(D, d)
```

Example:

```julia
using TTVarietyDegree

D = [1, 2, 2, 1]
d = [3, 3, 3]

P = compute_tail_polynomial(D, d)
fP = compute_fP(P, D[end - 1], d[end])
deg = degree_TT_variety(D, d)
```

`compute_tail_polynomial` is quiet by default. To inspect the recursive
Schur-Weingarten steps, use:

```julia
P = compute_tail_polynomial(D, d; verbose=true)
```

By default, `degree_TT_variety(D, d)` first reduces vacuous boundary modes and
then applies the Schur-Weingarten recursion to the reduced signature.

For three-way reduced signatures, an optional subspace-variety formula is also
available:

```julia
deg = degree_TT_variety(D, d; method=:subspace)
```

Lower-level combinatorial helpers are intentionally not exported from the
package namespace. They are still available for inspection using qualified
names such as `TTVarietyDegree.partitions(5)` or
`TTVarietyDegree.compute_H_schur((3, 1), 2, 2, 2)`.

## Tests

Run the checks:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

or equivalently:

```bash
julia --project=. test/runtests.jl
```

## Command-Line Scripts

The runnable utilities live in `scripts/`. Compatibility wrappers remain at
the package root, so both command styles work.

For a direct computation, edit `D` and `d` at the top of
`scripts/compute_degree.jl` and run:

```bash
julia scripts/compute_degree.jl
```

or pass them as command-line inputs:

```bash
julia scripts/compute_degree.jl "[1,2,2,1]" "[3,3,3]"
```

By default this computes only the degree. To also print the recursive tail
polynomial and `f(P)`, use:

```bash
julia scripts/compute_degree.jl "[1,2,2,1]" "[3,3,3]" --tail
```

To print the intermediate recursive Schur-Weingarten steps while computing the
tail, add `--verbose`:

```bash
julia scripts/compute_degree.jl "[1,2,2,1]" "[3,3,3]" --tail --verbose
```

To explicitly use the optional subspace-variety formula when the reduced input
is three-way, pass:

```bash
julia scripts/compute_degree.jl "[1,3,3,3,3,1]" "[3,3,3,3,3]" --subspace
```

`compute_degree.jl` also starts a resident-memory watchdog. The default limit
is set near the top of the script as `MAX_RSS_GB = 8.0`. Override it on the
command line with:

```bash
julia scripts/compute_degree.jl "[1,2,2,1]" "[3,3,3]" --max-ram-gb=4
```

or disable the watchdog with:

```bash
julia scripts/compute_degree.jl "[1,2,2,1]" "[3,3,3]" --no-ram-limit
```

If the process exceeds the limit, it exits with code `99` and prints a message
to stderr. The watchdog reads `/proc/self/status`, so it is intended for Linux.

To regenerate the N=2,3,4 comparison tables with this implementation, run:

```bash
julia scripts/generate_degree_tables.jl
```

The output is written to `computed_tables/` as exact integer CSV files:

```text
computed_tables/TT_variety_degrees_schur_weingarten_N_2.csv
computed_tables/TT_variety_degrees_schur_weingarten_N_3.csv
computed_tables/TT_variety_degrees_schur_weingarten_N_4.csv
```

The row and column signatures are copied from the existing comparison-table
exports in `../intersection_degree_computations_and_parametrizations/degree_tables`.
Entries outside the TT admissibility range are written as `0`.

## Boundary Mode Reduction

Before computing a degree, `degree_TT_variety` removes vacuous full-rank
constraints at the two ends of the tensor train. This can make a large-looking
input much cheaper.

For a tensor train with

```text
D = [1, D1, D2, ..., D_{N-1}, 1]
d = [d1, d2, ..., dN]
```

the first TT rank condition is the rank condition on the first flattening

```text
C^d1  |  C^d2 tensor ... tensor C^dN.
```

This flattening has only `d1` rows, so the condition `rank <= D1` is vacuous
when `D1 = d1`. In that case the first mode can be merged with the second mode:

```text
D = [1, D1, D2, D3, ..., 1]
d = [d1, d2, d3, ..., dN]

becomes

D = [1, D2, D3, ..., 1]
d = [d1*d2, d3, ..., dN].
```

Geometrically, this is only a reshaping

```text
C^d1 tensor C^d2  =  C^(d1*d2).
```

The first rank condition was already automatic, and all remaining TT flattening
rank conditions are exactly the same conditions after this reshaping. Therefore
the projective variety and its degree are unchanged.

The same rule applies at the right end. The last flattening has only `dN`
columns, so the condition `rank <= D_{N-1}` is vacuous when
`D_{N-1} = dN`. Then the last two modes can be merged:

```text
D = [1, ..., D_{N-2}, D_{N-1}, 1]
d = [d1, ..., d_{N-2}, d_{N-1}, dN]

becomes

D = [1, ..., D_{N-2}, 1]
d = [d1, ..., d_{N-2}, d_{N-1}*dN].
```

The code repeats these two reductions until neither boundary rank is full.

Example:

```text
D = [1,3,3,3,3,1]
d = [3,3,3,3,3]
```

First, `D1 = d1 = 3`, so the left boundary mode is merged:

```text
D = [1,3,3,3,1]
d = [9,3,3,3]
```

Then the right boundary has `D3 = d4 = 3`, so the right boundary mode is merged:

```text
D = [1,3,3,1]
d = [9,3,9]
```

At this point no boundary rank is full: `3 < 9` on both sides. The reduced
problem is a three-way tensor format. By default, the code still applies the
same Schur-Weingarten mechanism to this reduced signature. For this example,
the computation is therefore done as

```text
D = [1,3,3,1]
d = [9,3,9]
```

using the Schur-Weingarten recursion.

There is also an optional shortcut for reduced three-way inputs

```text
D = [1,r,s,1]
d = [a,b,c].
```

Passing `method=:subspace` in Julia, or `--subspace` in
`scripts/compute_degree.jl`,
uses the standard subspace-variety coefficient formula. This is kept as a
separate verification route, not as the default.

## Notes

`compute_H_raw_monomial_by_weingarten` implements the Stiefel Weingarten
integral by summing over symmetric-group cycle classes. This is equivalent to
the entrywise expansion in the cheat sheet, but avoids materializing every
matrix-entry monomial before integration.

The main `compute_H_schur` routine goes one step further: it uses character
orthogonality to compute the Schur coefficients of the Weingarten integral
directly. This is still the same Schur-Weingarten calculation, but it avoids
constructing the intermediate raw monomial or power-sum polynomial.

`scripts/compute_degree.jl` does not compute `P(CC*)` or `f(P)` unless
`--tail` is passed. With `--tail`, it still skips printing `P(CC*)` when the
generic tail expansion would require a large Weingarten degree.
