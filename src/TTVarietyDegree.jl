"""
Exact degree computations for tensor train varieties.

The public API computes projective dimensions, recursive Schur tail
polynomials, the final functional ``f(P)``, and degrees of tensor train
varieties from a rank signature `D = [D_0, ..., D_N]` and mode dimensions
`d = [d_1, ..., d_N]`. The default algorithm is the exact Schur-Weingarten
recursion with `Rational{BigInt}` coefficients.
"""
module TTVarietyDegree

export compute_tail_polynomial,
       f_schur,
       compute_fP,
       degree_TT_variety,
       dimension_TT_variety,
       compress_vacuous_boundary_modes

# Partition keys are stored as tuples, for example (6, 4).
const Partition = Tuple{Vararg{Int}}

# All rational coefficients use exact BigInt numerator/denominator arithmetic.
const Rat = Rational{BigInt}

# Memorization tables for the combinatorics and Schur-Weingarten subroutines.
const _partition_cache = Dict{Tuple{Int, Int, Int}, Vector{Partition}}()
const _permutation_cache = Dict{Int, Vector{Vector{Int}}}()
const _class_permutation_cache = Dict{Int, Dict{Partition, Vector{Vector{Int}}}}()
const _rim_hook_cache = Dict{Tuple{Partition, Int}, Vector{Tuple{Partition, Int}}}()
const _character_cache = Dict{Tuple{Partition, Partition}, BigInt}()
const _weingarten_cache = Dict{Tuple{Partition, Int, Int}, Rat}()
const _weingarten_class_sum_cache = Dict{Tuple{Int, Int}, Dict{Tuple{Partition, Partition}, Rat}}()
const _product_cycle_distribution_cache = Dict{Tuple{Int, Partition}, Dict{Tuple{Partition, Partition}, BigInt}}()
const _complete_h_cache = Dict{Tuple{Int, Int}, Dict{Partition, BigInt}}()
const _schur_raw_cache = Dict{Tuple{Partition, Int}, Dict{Partition, BigInt}}()
const _schur_monomial_cache = Dict{Tuple{Partition, Int}, Dict{Partition, BigInt}}()
const _transition_cache = Dict{Tuple{Int, Int}, Tuple{Vector{Partition}, Vector{Vector{Rat}}}}()
const _H_cache = Dict{Tuple{Partition, Int, Int, Int}, Dict{Partition, Rat}}()


# Input: nonnegative integer n.
# Output: n! as a BigInt.
# Purpose: avoid overflow in hook-length, Grassmannian, and Weingarten formulas.
_bigfactorial(n::Integer) = factorial(big(n))

# Input: sparse coefficient dictionary dict.
# Output: the same dictionary, mutated so entries with zero coefficients are gone.
# Purpose: keep sparse polynomial and expansion dictionaries compact.
function _cleanup!(dict::Dict{K, V}) where {K, V}
    for key in collect(keys(dict))
        if iszero(dict[key])
            delete!(dict, key)
        end
    end
    return dict
end

# Input: partition-like iterable lam, possibly with trailing zeroes.
# Output: normalized partition tuple with trailing zeroes removed.
# Throws: ArgumentError if any part is negative.
function trim_trailing_zeroes(lam)
    values = Int[x for x in lam]
    while !isempty(values) && values[end] == 0
        pop!(values)
    end
    if any(x -> x < 0, values)
        throw(ArgumentError("partitions cannot contain negative parts"))
    end
    return Tuple(values)
end

# Input: partition lam and target number of parts length.
# Output: tuple of length length obtained by appending zeroes.
# Throws: ArgumentError if lam has more than length nonzero parts.
function pad_partition(lam, length::Integer)
    length < 0 && throw(ArgumentError("padding length must be nonnegative"))
    trimmed = trim_trailing_zeroes(lam)
    if Base.length(trimmed) > length
        throw(ArgumentError("partition $trimmed has length greater than $length"))
    end
    return Tuple(vcat(collect(trimmed), zeros(Int, Int(length) - Base.length(trimmed))))
end

# Input: partition lam.
# Output: size |lam|, with the empty partition () treated as size 0.
partition_size(lam) = isempty(lam) ? 0 : sum(lam)

# Input: n, maximum allowed part max_part, and maximum length max_length.
# Output: vector of partitions of n satisfying those bounds.
# Purpose: internal cached recursion behind partitions(...).
function _partitions_cached(n::Int, max_part::Int, max_length::Int)
    key = (n, max_part, max_length)
    cached = get(_partition_cache, key, nothing)
    cached !== nothing && return cached

    out = Partition[]
    if n == 0
        push!(out, ())
    elseif max_length > 0 && max_part > 0
        for first in min(n, max_part):-1:1
            for rest in _partitions_cached(n - first, min(first, n - first), max_length - 1)
                push!(out, (first, rest...))
            end
        end
    end

    _partition_cache[key] = out
    return out
end

# Input: nonnegative integer n and optional max_length.
# Output: vector of all partitions of n, optionally restricted to len <= max_length.
# Order: lexicographic by first parts descending, e.g. (4), (3,1), ...
function partitions(n::Integer; max_length=nothing)
    n < 0 && throw(ArgumentError("n must be nonnegative"))
    nn = Int(n)
    ml = max_length === nothing ? nn : Int(max_length)
    return _partitions_cached(nn, nn, ml)
end

# Input: partition lam.
# Output: BigInt vector containing the hook length h_ij for every box (i,j).
# Purpose: used by hook-content and hook-length formulas.
function hook_lengths(lam)
    lam = trim_trailing_zeroes(lam)
    hooks = BigInt[]
    for i in eachindex(lam)
        for j in 1:lam[i]
            arm = lam[i] - j
            leg = count(r -> r >= j, lam[(i + 1):end])
            push!(hooks, big(arm + leg + 1))
        end
    end
    return hooks
end

# Input: partition lam.
# Output: BigInt product of all hook lengths of lam.
# Purpose: denominator in the hook-length formula for f^lam.
function _hook_product(lam)
    prod = big(1)
    for hook in hook_lengths(lam)
        prod *= hook
    end
    return prod
end

# Input: partition lam and number of variables k.
# Output: exact rational s_lam(1^k); returns 0 if length(lam) > k.
# Formula: product over boxes (k+j-i)/h_ij.
function schur_at_ones(lam, k::Integer)
    k < 0 && throw(ArgumentError("number of variables must be nonnegative"))
    lam = trim_trailing_zeroes(lam)
    Base.length(lam) > k && return zero(Rat)

    value = one(Rat)
    for i in eachindex(lam)
        for j in 1:lam[i]
            content = Int(k) + j - i
            content == 0 && return zero(Rat)
            arm = lam[i] - j
            leg = count(r -> r >= j, lam[(i + 1):end])
            hook = arm + leg + 1
            value *= big(content) // big(hook)
        end
    end
    return value
end

# Input: integers k,n specifying Gr(k,n), with 0 <= k <= n.
# Output: integer projective degree of the Grassmannian.
# Formula: (k(n-k))! divided by the hook product of the k x (n-k) rectangle.
function degree_grassmannian(k::Integer, n::Integer)
    k = Int(k)
    n = Int(n)
    (k < 0 || n < 0 || k > n) && throw(ArgumentError("expected 0 <= k <= n"))
    m = n - k
    denom = big(1)
    for i in 1:k, j in 1:m
        denom *= big(k - i + m - j + 1)
    end
    value = _bigfactorial(k * m) // denom
    denominator(value) == 1 || error("Grassmannian degree was not integral")
    return numerator(value)
end


# Input: integer vector values.
# Output: true iff all entries are nonnegative and weakly decreasing.
# Purpose: validates candidate partitions after rim-hook removals.
function _is_partition_vector(values::Vector{Int})
    all(x -> x >= 0, values) || return false
    for i in 1:(Base.length(values) - 1)
        values[i] >= values[i + 1] || return false
    end
    return true
end

# Input: set of Young-diagram boxes, each stored as (row, column).
# Output: true iff the boxes form one edge-connected component.
# Purpose: rim hooks must be connected skew strips.
function _is_connected_strip(strip::Set{Tuple{Int, Int}})
    isempty(strip) && return false
    start = first(strip)
    stack = [start]
    seen = Set{Tuple{Int, Int}}()
    while !isempty(stack)
        box = pop!(stack)
        box in seen && continue
        push!(seen, box)
        i, j = box
        for nb in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1))
            if nb in strip && !(nb in seen)
                push!(stack, nb)
            end
        end
    end
    return Base.length(seen) == Base.length(strip)
end

# Input: set of Young-diagram boxes.
# Output: true iff no four boxes form a 2-by-2 square.
# Purpose: rim hooks are connected skew strips without 2-by-2 blocks.
function _has_no_two_by_two(strip::Set{Tuple{Int, Int}})
    for (i, j) in strip
        if ((i + 1, j) in strip) && ((i, j + 1) in strip) && ((i + 1, j + 1) in strip)
            return false
        end
    end
    return true
end

# Input: partition lam and positive integer q.
# Output: vector of (mu, height), where mu is lam after removing a rim hook
# of q boxes and height is the number of rows touched by that hook.
# Purpose: Murnaghan-Nakayama character recursion.
function _rim_hook_removals(lam::Partition, q::Int)
    key = (lam, q)
    cached = get(_rim_hook_cache, key, nothing)
    cached !== nothing && return cached

    lam = trim_trailing_zeroes(lam)
    rows = Base.length(lam)
    counts = zeros(Int, rows)
    out = Tuple{Partition, Int}[]

    # Input: current row-wise removal counts in counts.
    # Output: appends a valid (mu,height) to out, or does nothing.
    # Purpose: validate partition shape, connectedness, and no-2-by-2 condition.
    function check_candidate()
        mu_vec = [lam[i] - counts[i] for i in 1:rows]
        _is_partition_vector(mu_vec) || return

        strip = Set{Tuple{Int, Int}}()
        for i in 1:rows
            for j in (mu_vec[i] + 1):lam[i]
                push!(strip, (i, j))
            end
        end
        Base.length(strip) == q || return
        _has_no_two_by_two(strip) || return
        _is_connected_strip(strip) || return

        height = count(c -> c > 0, counts)
        push!(out, (trim_trailing_zeroes(mu_vec), height))
    end

    # Input: row index and remaining number of boxes to remove.
    # Output: enumerates all count vectors summing to q via recursion.
    # Purpose: brute-force row-wise rim-hook candidates.
    function rec(row::Int, remaining::Int)
        if row > rows
            remaining == 0 && check_candidate()
            return
        end
        for c in 0:min(lam[row], remaining)
            counts[row] = c
            rec(row + 1, remaining - c)
        end
        counts[row] = 0
    end

    q > 0 && q <= partition_size(lam) && rec(1, q)
    _rim_hook_cache[key] = out
    return out
end

# Input: representation partition lam and conjugacy class cycle_type of S_|lam|.
# Output: BigInt value chi^lam(cycle_type).
# Method: Murnaghan-Nakayama recursion using cached rim-hook removals.
function symmetric_group_character(lam, cycle_type)
    lam = trim_trailing_zeroes(lam)
    cycle_type = trim_trailing_zeroes(sort(collect(cycle_type), rev=true))
    partition_size(lam) == partition_size(cycle_type) || return big(0)

    key = (lam, cycle_type)
    cached = get(_character_cache, key, nothing)
    cached !== nothing && return cached

    value = big(0)
    if isempty(cycle_type)
        value = isempty(lam) ? big(1) : big(0)
    else
        q = cycle_type[1]
        rest = Tuple(cycle_type[2:end])
        for (mu, height) in _rim_hook_removals(lam, q)
            sign = iseven(height - 1) ? big(1) : big(-1)
            value += sign * symmetric_group_character(mu, rest)
        end
    end

    _character_cache[key] = value
    return value
end

# Input: integer partition, usually a cycle type mu.
# Output: BigInt z_mu = prod_i i^(m_i(mu)) * m_i(mu)!.
# Purpose: centralizer size and denominator in power-sum/character formulas.
function z_value(partition)
    counts = Dict{Int, Int}()
    for part in partition
        part > 0 || continue
        counts[part] = get(counts, part, 0) + 1
    end
    value = big(1)
    for (part, mult) in counts
        value *= big(part)^mult * _bigfactorial(mult)
    end
    return value
end

# Input: partition lam and integer n.
# Output: BigInt C_lam(n) = prod over boxes (i,j) of n+j-i.
# Purpose: content denominator in the unitary Weingarten formula.
function _content_product(lam, n::Int)
    prod = big(1)
    for i in eachindex(lam)
        for j in 1:lam[i]
            prod *= big(n + j - i)
        end
    end
    return prod
end

# Input: partition lam of p.
# Output: BigInt dimension f^lam of the irreducible S_p representation.
# Formula: p! divided by the hook product of lam.
function _symmetric_group_irrep_dimension(lam)
    p = partition_size(lam)
    value = _bigfactorial(p) // _hook_product(lam)
    denominator(value) == 1 || error("hook-length formula did not produce an integer")
    return numerator(value)
end

# Input: cycle_type partition of p.
# Output: BigInt size of the conjugacy class in S_p, equal to p!/z_cycle_type.
function _class_size(cycle_type::Partition)
    p = partition_size(cycle_type)
    return _bigfactorial(p) ÷ z_value(cycle_type)
end



# Input: Schur label lam from the previous recursion level, previous rank
# D_prev, current rank D_cur, and physical dimension d_r.
# Output: dictionary mu => coefficient of s_mu in H_lam^{(r)}.
# Method: direct Schur-basis Stiefel-Weingarten integral via character
# orthogonality, avoiding raw monomial expansion in the main path.
function compute_H_schur(lam, D_prev::Integer, D_cur::Integer, d_r::Integer)
    lam = trim_trailing_zeroes(lam)
    key = (lam, Int(D_prev), Int(D_cur), Int(d_r))
    cached = get(_H_cache, key, nothing)
    cached !== nothing && return cached

    p = partition_size(lam)
    if p == 0
        out = Dict{Partition, Rat}(() => one(Rat))
        _H_cache[key] = out
        return out
    end

    n_r = Int(D_prev) * Int(d_r)
    prefactor = schur_at_ones(lam, Int(D_prev)) / _symmetric_group_irrep_dimension(lam)
    out = Dict{Partition, Rat}()

    for mu in partitions(p; max_length=Int(D_cur))
        weighted_inner_product = big(0)
        for kappa in partitions(p)
            chi_lam = symmetric_group_character(lam, kappa)
            chi_lam == 0 && continue
            chi_mu = symmetric_group_character(mu, kappa)
            chi_mu == 0 && continue
            weighted_inner_product +=
                _class_size(kappa) * big(d_r)^Base.length(kappa) * chi_lam * chi_mu
        end

        if !iszero(weighted_inner_product)
            coeff = prefactor * weighted_inner_product / _content_product(mu, n_r)
            !iszero(coeff) && (out[mu] = coeff)
        end
    end

    _H_cache[key] = out
    return out
end

# Input: Schur expansion T_coeffs, determinant power m_r, and variable count D_r.
# Output: Schur expansion of det(Sigma)^m_r * T_coeffs.
# Effect: pads each partition to D_r parts and adds m_r to every part.
function determinant_shift(T_coeffs, m_r::Integer, D_r::Integer)
    m_r = Int(m_r)
    D_r = Int(D_r)
    out = Dict{Partition, Rat}()
    for (mu, coeff) in T_coeffs
        padded = pad_partition(mu, D_r)
        lam = trim_trailing_zeroes(x + m_r for x in padded)
        out[lam] = get(out, lam, zero(Rat)) + coeff
    end
    return _cleanup!(out)
end

# Input: TT signature D=[D_0,...,D_N] and d=[d_1,...,d_N].
# Output: N = length(d) if the basic conditions pass.
# Throws: ArgumentError for wrong lengths, nonpositive entries, bad endpoints,
# or local Stiefel infeasibility D_r > D_{r-1} d_r.
function _validate_TT_input(D, d)
    N = Base.length(d)
    Base.length(D) == N + 1 || throw(ArgumentError("D must have length length(d)+1"))
    D[1] == 1 || throw(ArgumentError("expected D_0 = 1"))
    D[end] == 1 || throw(ArgumentError("expected D_N = 1"))
    all(x -> x > 0, D) || throw(ArgumentError("all D entries must be positive"))
    all(x -> x > 0, d) || throw(ArgumentError("all d entries must be positive"))
    for r in 1:(N - 1)
        n_r = D[r] * d[r]
        n_r >= D[r + 1] || throw(ArgumentError("expected D_$r <= D_$(r - 1)d_$r"))
    end
    return N
end

"""
    compute_tail_polynomial(D, d; verbose=false, io=stdout)

Compute the recursive tail polynomial ``F_{N-1}`` for the tensor train
signature `D = [D_0, ..., D_N]` and `d = [d_1, ..., d_N]`.

The return value is a dictionary `lam => coeff` representing

```math
F_{N-1}(CC^*) = \\sum_\\lambda c_\\lambda s_\\lambda(CC^*).
```

Partitions are stored as tuples, for example `(6, 4)`, and coefficients are
exact `Rational{BigInt}` values. Set `verbose=true` to print each local
Schur-Weingarten averaging step to `io`.
"""
function compute_tail_polynomial(D, d; verbose::Bool=false, io::IO=stdout)
    D = Int[x for x in D]
    d = Int[x for x in d]
    N = _validate_TT_input(D, d)

    F = Dict{Partition, Rat}(() => one(Rat))
    verbose && println(io, "For r = 0 we have F = $F")
    for r in 1:(N - 1)
        verbose && println(io, "---------r = $r -------------")
        D_prev = D[r]
        D_cur = D[r + 1]
        d_r = d[r]
        m_r = D_prev * d_r - D_cur

        T = Dict{Partition, Rat}()
        for (lam, coeff) in F
            H_lam = compute_H_schur(lam, D_prev, D_cur, d_r)
            verbose && println(io, "--H-function for lam = $lam: H_lam = $H_lam with D_prev = $D_prev, D_cur = $D_cur, d_r = $d_r")
            for (mu, hcoeff) in H_lam
                T[mu] = get(T, mu, zero(Rat)) + coeff * hcoeff
                verbose && println(io, "--Adding contribution from lam = $lam, mu = $mu: coeff=$coeff, hcoeff=$hcoeff, total for mu = $mu is now $(T[mu])")
            end
        end
        F = determinant_shift(T, m_r, D_cur)
        verbose && println(io, "For r = $r we have F = $F")
    end

    return F
end

"""
    f_schur(lam, k, n)

Evaluate the final-tail functional on one Schur polynomial:

```math
f(s_\\lambda(CC^*))
=
s_\\lambda(1^k)
\\prod_{i=1}^k \\frac{(\\lambda_i+n-i)!}{(n-i)!}.
```

Here `C` is a `k x n` complex matrix and missing parts of `lam` are treated as
zero. The result is returned exactly as a `Rational{BigInt}`.
"""
function f_schur(lam, k::Integer, n::Integer)
    k = Int(k)
    n = Int(n)
    lam = trim_trailing_zeroes(lam)
    Base.length(lam) > k && return zero(Rat)
    n >= k || throw(ArgumentError("the final tail formula expects n >= k"))

    padded = pad_partition(lam, k)
    value = schur_at_ones(padded, k)
    for i in 1:k
        value *= _bigfactorial(padded[i] + n - i) // _bigfactorial(n - i)
    end
    return value
end

"""
    compute_fP(P_coeffs, k, n)

Apply the final-tail functional to a Schur expansion `P_coeffs`.

`P_coeffs` should be a dictionary `lam => coeff`, such as the output of
`compute_tail_polynomial`. The returned value is the exact integer `f(P)` as a
`BigInt`.
"""
function compute_fP(P_coeffs, k::Integer, n::Integer)
    total = zero(Rat)
    for (lam, coeff) in P_coeffs
        total += coeff * f_schur(lam, k, n)
    end
    denominator(total) == 1 || error("f(P) was not integral: $total")
    return numerator(total)
end

"""
    compress_vacuous_boundary_modes(D, d)

Merge full-rank boundary modes without changing the tensor train variety.

If `D_1 == d_1`, the first TT rank condition is vacuous and the first two mode
dimensions are replaced by `d_1*d_2`. The analogous right-boundary rule is
applied when `D_{N-1} == d_N`. The process repeats until neither boundary rank
is full.

Returns `(D_reduced, d_reduced)`.
"""
function compress_vacuous_boundary_modes(D, d)
    D = Int[x for x in D]
    d = Int[x for x in d]
    _validate_TT_input(D, d)

    changed = true
    while changed && Base.length(d) > 1
        changed = false

        while Base.length(d) > 1 && D[2] == d[1]
            d = vcat([d[1] * d[2]], d[3:end])
            D = vcat([D[1]], D[3:end])
            changed = true
        end

        while Base.length(d) > 1 && D[end - 1] == d[end]
            d = vcat(d[1:(end - 2)], [d[end - 1] * d[end]])
            D = vcat(D[1:(end - 2)], [D[end]])
            changed = true
        end
    end

    return D, d
end



"""
    degree_TT_variety(D, d; method=:schur_weingarten, reduce=true)

Compute the projective degree of the tensor train variety with rank signature
`D = [D_0, ..., D_N]` and mode dimensions `d = [d_1, ..., d_N]`.

By default, vacuous full-rank boundary modes are first compressed and the exact
Schur-Weingarten recursion is applied to the reduced signature. 

Returns the degree as a `BigInt`.
"""
function degree_TT_variety(D, d; method::Symbol=:schur_weingarten, reduce::Bool=true, verbose::Bool=false)
    D = Int[x for x in D]
    d = Int[x for x in d]
    _validate_TT_input(D, d)

    if reduce
        reduced_D, reduced_d = compress_vacuous_boundary_modes(D, d)
        if (reduced_D, reduced_d) != (D, d)
            return degree_TT_variety(reduced_D, reduced_d; method=method, reduce=false, verbose=verbose)
        end
    end

    N = Base.length(d)
    P_coeffs = compute_tail_polynomial(D, d; verbose=verbose)
    fP = compute_fP(P_coeffs, D[N], d[N])

    degree = fP // big(1)
    for r in 1:(N - 1)
        n_r = D[r] * d[r]
        m_r = n_r - D[r + 1]
        degree *= degree_grassmannian(D[r + 1], n_r)
        degree /= _bigfactorial(D[r + 1] * m_r)
    end

    denominator(degree) == 1 || error("degree was not integral: $degree")
    return numerator(degree)
end

"""
    dimension_TT_variety(D, d)

Compute the projective dimension

```math
sum_{r=1}^{N-1} D_r(D_{r-1}d_r-D_r) + D_{N-1}d_N - 1
```

for the tensor train signature `D,d`.
"""
function dimension_TT_variety(D, d)
    D = Int[x for x in D]
    d = Int[x for x in d]
    N = _validate_TT_input(D, d)

    running = 0
    for r in 1:(N - 1)
        n_r = D[r] * d[r]
        m_r = n_r - D[r + 1]
        running += D[r + 1] * m_r
    end
    return running + D[N] * d[N] - 1
end

function degree_segre_variety(d::Vector{Int}, N::Int)
    numerator = _bigfactorial(sum(d))
    denominator = big(1)
    for i in 1:(N)
        denominator *= _bigfactorial(d[i])
    end
    return numerator ÷ denominator
end

end # module TTVarietyDegree
