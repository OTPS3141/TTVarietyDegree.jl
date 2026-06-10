# TTVarietyDegree

Exact computation of degrees of tensor train varieties
`V_{D,d}` using the recursive Schur-Weingarten procedure as in the paper 

The implementation is pure Julia and uses `Rational{BigInt}` throughout.

## Installation

Using the registered package

```julia
using Pkg
Pkg.add("TTVarietyDegree")
```

or via public repository

```julia
using Pkg
Pkg.add(url="https://github.com/OTPS3141/TTVarietyDegree.jl")
```


## Julia API

The main package API is:

```julia
degree_TT_variety(D, d; method=:schur_weingarten, reduce=true)
dimension_TT_variety(D, d)
compute_tail_polynomial(D, d; verbose=false)
compute_fP(P_coeffs, k, n)
f_schur(lam, k, n)
compress_vacuous_boundary_modes(D, d)
```

Example use of API functions (assuming existing TTVarietyDegree.jl installation):

```julia
using TTVarietyDegree

degree_TT_variety([1,2,2,1], [3,3,3])
```

`compute_tail_polynomial` is quiet by default. To inspect the recursive
Schur-Weingarten steps, use:

```julia
P = compute_tail_polynomial(D, d; verbose=true)
```

By default, `degree_TT_variety(D, d)` first reduces vacuous boundary modes and
then applies the Schur-Weingarten recursion to the reduced signature.



## Command-Line Scripts

The runnable utilities live in `scripts/`. Compatibility wrappers remain at
the package root, so both command styles work.

Pass `D` and `d` (first and second argument) as command-line inputs:

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




