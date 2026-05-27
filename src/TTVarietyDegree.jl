"""
Exact symbolic degree computations for tensor train varieties.

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
       compress_vacuous_boundary_modes,
       degree_subspace_variety_N3

# Partition keys are stored as tuples, for example (6, 4).
const Partition = Tuple{Vararg{Int}}

# All rational coefficients use exact BigInt numerator/denominator arithmetic.
const Rat = Rational{BigInt}

# Memoization tables for the combinatorics and Schur-Weingarten subroutines.
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

# Input: integer n.
# Output: exact rational n//1 as Rational{BigInt}.
# Purpose: small convenience for constructing exact coefficients.
_rat(n::Integer) = big(n) // big(1)

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

# Input: cycle_type, a partition of p.
# Output: one permutation vector of 1:p with exactly that cycle type.
# Convention: permutation vectors use 1-based images, so perm[i] is pi(i).
function _representative_permutation(cycle_type::Partition)
    p = partition_size(cycle_type)
    perm = collect(1:p)
    start = 1
    for len in cycle_type
        if len > 1
            for i in start:(start + len - 2)
                perm[i] = i + 1
            end
            perm[start + len - 1] = start
        end
        start += len
    end
    return perm
end

# Input: nonnegative integer p.
# Output: vector of all permutations of 1:p, each stored as a 1-based image vector.
# Note: factorial-size; only used by diagnostic/legacy routines at small p.
function _permutations_list(p::Int)
    cached = get(_permutation_cache, p, nothing)
    cached !== nothing && return cached

    out = Vector{Int}[]
    values = collect(1:p)

    # Input: current position pos in the permutation under construction.
    # Output: appends all completions to out.
    # Purpose: standard in-place backtracking over all permutations.
    function backtrack(pos::Int)
        if pos > p
            push!(out, copy(values))
            return
        end
        for i in pos:p
            values[pos], values[i] = values[i], values[pos]
            backtrack(pos + 1)
            values[pos], values[i] = values[i], values[pos]
        end
    end
    if p == 0
        push!(out, Int[])
    else
        backtrack(1)
    end
    _permutation_cache[p] = out
    return out
end

# Input: permutations p and q of the same size.
# Output: composition p after q, i.e. r[i] = p[q[i]].
# Purpose: cycle-type computations for diagnostic product distributions.
function _compose(p::Vector{Int}, q::Vector{Int})
    n = Base.length(q)
    out = Vector{Int}(undef, n)
    @inbounds for i in 1:n
        out[i] = p[q[i]]
    end
    return out
end

# Input: permutation vector p.
# Output: inverse permutation vector inv satisfying inv[p[i]] = i.
function _inverse_perm(p::Vector{Int})
    inv = similar(p)
    @inbounds for i in eachindex(p)
        inv[p[i]] = i
    end
    return inv
end

# Input: permutation vector perm.
# Output: cycle type as a partition, sorted in descending order.
function _cycle_type(perm::Vector{Int})
    n = Base.length(perm)
    seen = falses(n)
    lengths = Int[]
    for i in 1:n
        if !seen[i]
            cur = i
            len = 0
            while !seen[cur]
                seen[cur] = true
                len += 1
                cur = perm[cur]
            end
            push!(lengths, len)
        end
    end
    sort!(lengths, rev=true)
    return Tuple(lengths)
end

# Input: permutation vector perm.
# Output: +1 for even permutations and -1 for odd permutations.
# Formula: sign = (-1)^(p - number_of_cycles).
function _perm_sign(perm::Vector{Int})
    return iseven(Base.length(perm) - Base.length(_cycle_type(perm))) ? 1 : -1
end

# Input: nonnegative integer p.
# Output: dictionary mapping cycle_type => all permutations in that conjugacy class.
# Note: factorial-size; retained for diagnostic/legacy computations.
function _class_permutations(p::Int)
    cached = get(_class_permutation_cache, p, nothing)
    cached !== nothing && return cached

    classes = Dict{Partition, Vector{Vector{Int}}}()
    for perm in _permutations_list(p)
        ctype = _cycle_type(perm)
        push!(get!(classes, ctype, Vector{Int}[]), perm)
    end
    _class_permutation_cache[p] = classes
    return classes
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

# Input: cycle type of pi in S_p, the degree p, and ambient dimension n.
# Output: exact rational unitary Weingarten value Wg_n(pi).
# Formula: 1/p! * sum_rho f^rho chi^rho(pi)/C_rho(n).
function unitary_weingarten(cycle_type_pi, p::Integer, n::Integer)
    p = Int(p)
    n = Int(n)
    cycle_type_pi = trim_trailing_zeroes(sort(collect(cycle_type_pi), rev=true))
    partition_size(cycle_type_pi) == p || throw(ArgumentError("cycle type $cycle_type_pi does not have size $p"))
    p == 0 && return one(Rat)

    key = (cycle_type_pi, p, n)
    cached = get(_weingarten_cache, key, nothing)
    cached !== nothing && return cached

    total = zero(Rat)
    for rho in partitions(p)
        Base.length(rho) <= n || continue
        f_rho = _symmetric_group_irrep_dimension(rho)
        chi = symmetric_group_character(rho, cycle_type_pi)
        chi == 0 && continue
        total += (f_rho * chi) // _content_product(rho, n)
    end
    value = total / _bigfactorial(p)
    _weingarten_cache[key] = value
    return value
end

# Input: Weingarten degree p and ambient dimension n.
# Output: dictionary (kappa, eta) => class-summed Weingarten coefficient.
# Purpose: replaces explicit sums over permutations by character orthogonality.
function _weingarten_class_sum_matrix(p::Int, n::Int)
    key = (p, n)
    cached = get(_weingarten_class_sum_cache, key, nothing)
    cached !== nothing && return cached

    parts = partitions(p)
    matrix = Dict{Tuple{Partition, Partition}, Rat}()
    for kappa in parts
        for eta in parts
            total = zero(Rat)
            for rho in parts
                Base.length(rho) <= n || continue
                chi_eta = symmetric_group_character(rho, eta)
                chi_eta == 0 && continue
                chi_kappa = symmetric_group_character(rho, kappa)
                chi_kappa == 0 && continue
                total += (chi_eta * chi_kappa) // _content_product(rho, n)
            end
            matrix[(kappa, eta)] = total / z_value(eta)
        end
    end

    _weingarten_class_sum_cache[key] = matrix
    return matrix
end

# Input: degree p and a cycle type nu used to choose a representative phi.
# Output: dictionary (cycle_type(sigma), cycle_type(phi*sigma)) => count.
# Note: diagnostic/legacy helper; factorial-size because it enumerates S_p.
function _product_cycle_count_distribution(p::Int, nu::Partition)
    key = (p, nu)
    cached = get(_product_cycle_distribution_cache, key, nothing)
    cached !== nothing && return cached

    phi = _representative_permutation(nu)
    dist = Dict{Tuple{Partition, Partition}, BigInt}()
    for sigma in _permutations_list(p)
        kappa = _cycle_type(sigma)
        zeta = _cycle_type(_compose(phi, sigma))
        dkey = (kappa, zeta)
        dist[dkey] = get(dist, dkey, big(0)) + 1
    end

    _product_cycle_distribution_cache[key] = dist
    return dist
end

# Input: power-sum type nu, previous rank D_prev, physical dimension d_r,
# and Stiefel ambient dimension n_r = D_prev*d_r.
# Output: dictionary eta => coefficient of p_eta(Sigma_t) after integration.
# Method: Stiefel-Weingarten expansion compressed by character sums.
function _expected_power_sum_coefficients(nu::Partition, D_prev::Int, d_r::Int, n_r::Int)
    p = partition_size(nu)
    p == 0 && return Dict{Partition, Rat}(() => one(Rat))

    # This is the entrywise Stiefel-Weingarten expansion after summing the
    # row/block/color Kronecker constraints by character orthogonality.
    parts = partitions(p)
    wg_sums = _weingarten_class_sum_matrix(p, n_r)
    out = Dict{Partition, Rat}()

    for kappa in parts
        block_factor = big(d_r)^Base.length(kappa)
        block_factor == 0 && continue

        row_sum = zero(Rat)
        for rho in parts
            schur_dim = schur_at_ones(rho, D_prev)
            iszero(schur_dim) && continue
            chi_nu = symmetric_group_character(rho, nu)
            chi_nu == 0 && continue
            chi_kappa = symmetric_group_character(rho, kappa)
            chi_kappa == 0 && continue
            f_rho = _symmetric_group_irrep_dimension(rho)
            row_sum += schur_dim * chi_nu * chi_kappa / f_rho
        end
        iszero(row_sum) && continue

        factor = block_factor * _class_size(kappa) * row_sum
        for eta in parts
            contribution = factor * wg_sums[(kappa, eta)]
            if !iszero(contribution)
                out[eta] = get(out, eta, zero(Rat)) + contribution
            end
        end
    end

    return _cleanup!(out)
end

# Input: power_coeffs mapping eta => coefficient of p_eta, total degree,
# and number of variables.
# Output: dictionary mu => coefficient of s_mu in the same symmetric polynomial.
# Formula: p_eta = sum_mu chi^mu(eta) s_mu, with length(mu) bounded.
function _power_sum_to_schur(power_coeffs::Dict{Partition, Rat}, degree::Int, num_variables::Int)
    degree == 0 && return Dict{Partition, Rat}(() => get(power_coeffs, (), zero(Rat)))

    out = Dict{Partition, Rat}()
    for mu in partitions(degree; max_length=num_variables)
        coeff = zero(Rat)
        for (eta, power_coeff) in power_coeffs
            chi = symmetric_group_character(mu, eta)
            chi == 0 && continue
            coeff += power_coeff * chi
        end
        !iszero(coeff) && (out[mu] = coeff)
    end
    return out
end

# Input: Schur label lam and local TT data D_prev, D_cur, d_r, n_r.
# Output: dictionary eta => coefficient of p_eta in H_lam.
# Role: slower explicit Frobenius-character route used by raw diagnostic output.
function _H_power_sum_coefficients(lam::Partition, D_prev::Int, D_cur::Int, d_r::Int, n_r::Int)
    degree = partition_size(lam)
    degree == 0 && return Dict{Partition, Rat}(() => one(Rat))

    total = Dict{Partition, Rat}()
    for nu in partitions(degree)
        char = symmetric_group_character(lam, nu)
        char == 0 && continue
        coeff_nu = char // z_value(nu)
        moment = _expected_power_sum_coefficients(nu, D_prev, d_r, n_r)
        for (eta, coeff) in moment
            total[eta] = get(total, eta, zero(Rat)) + coeff_nu * coeff
        end
    end

    return _cleanup!(total)
end

# Input: power-sum expansion in num_variables diagonal variables.
# Output: raw monomial dictionary gamma => coefficient, where gamma is an
# ordered exponent tuple of length num_variables.
function _power_sum_to_raw(power_coeffs::Dict{Partition, Rat}, num_variables::Int)
    num_variables < 0 && throw(ArgumentError("number of variables must be nonnegative"))
    raw = Dict{Partition, Rat}()
    zero_exp = ntuple(_ -> 0, num_variables)

    for (eta, coeff) in power_coeffs
        poly = Dict{Partition, Rat}(zero_exp => one(Rat))
        for q in eta
            next = Dict{Partition, Rat}()
            for (gamma, gamma_coeff) in poly
                for i in 1:num_variables
                    values = collect(gamma)
                    values[i] += q
                    key = Tuple(values)
                    next[key] = get(next, key, zero(Rat)) + gamma_coeff
                end
            end
            poly = next
        end

        for (gamma, gamma_coeff) in poly
            raw[gamma] = get(raw, gamma, zero(Rat)) + coeff * gamma_coeff
        end
    end

    return _cleanup!(raw)
end

# Input: Schur label lam and local TT data D_prev, D_cur, d_r, n_r.
# Output: raw monomial expansion of H_lam(Sigma_t) in diagonal variables.
# Use: diagnostic compatibility with the entrywise cheat-sheet algorithm.
function compute_H_raw_monomial_by_weingarten(lam, D_prev::Integer, D_cur::Integer, d_r::Integer, n_r::Integer)
    lam = trim_trailing_zeroes(lam)
    power_coeffs = _H_power_sum_coefficients(lam, Int(D_prev), Int(D_cur), Int(d_r), Int(n_r))
    return _power_sum_to_raw(power_coeffs, Int(D_cur))
end

# Input: raw monomial dictionary gamma => coefficient.
# Output: monomial-symmetric dictionary eta => common coefficient of all
# distinct permutations of eta.
# Throws: error if the raw polynomial is not symmetric.
function raw_monomial_to_monomial_symmetric(raw)
    out = Dict{Partition, Rat}()
    for (gamma, coeff) in raw
        iszero(coeff) && continue
        eta = trim_trailing_zeroes(sort(collect(gamma), rev=true))
        if haskey(out, eta)
            out[eta] == coeff || error("raw polynomial is not symmetric at monomial type $eta")
        else
            out[eta] = coeff
        end
    end
    return out
end

# Input: nonnegative integer n and number of parts k.
# Output: vector of all weak k-part compositions of n.
# Purpose: raw expansion of complete homogeneous symmetric polynomials.
function _compositions(n::Int, k::Int)
    if k == 0
        return n == 0 ? [()] : Partition[]
    elseif k == 1
        return [(n,)]
    end

    out = Partition[]
    for first in 0:n
        for rest in _compositions(n - first, k - 1)
            push!(out, (first, rest...))
        end
    end
    return out
end

# Input: degree q and number of variables k.
# Output: raw monomial dictionary for h_q(x_1,...,x_k) with BigInt coefficients.
# Convention: q < 0 returns the zero polynomial.
function _complete_homogeneous_raw(q::Int, k::Int)
    q < 0 && return Dict{Partition, BigInt}()
    key = (q, k)
    cached = get(_complete_h_cache, key, nothing)
    cached !== nothing && return cached

    out = Dict{Partition, BigInt}()
    for comp in _compositions(q, k)
        out[comp] = big(1)
    end
    _complete_h_cache[key] = out
    return out
end

# Input: number of variables k.
# Output: sparse raw polynomial dictionary representing constant 1.
function _poly_one(k::Int)
    return Dict{Partition, BigInt}(ntuple(_ -> 0, k) => big(1))
end

# Input: output polynomial out, input polynomial poly, integer scale.
# Output: mutates and returns out after adding scale*poly.
# Coefficients: BigInt sparse raw polynomial coefficients.
function _poly_add_scaled!(out::Dict{Partition, BigInt}, poly::Dict{Partition, BigInt}, scale::Integer)
    for (exp, coeff) in poly
        out[exp] = get(out, exp, big(0)) + big(scale) * coeff
    end
    return _cleanup!(out)
end

# Input: sparse raw polynomial dictionaries a and b.
# Output: sparse raw polynomial dictionary for a*b.
# Keys: ordered exponent tuples of equal variable count.
function _poly_mul(a::Dict{Partition, BigInt}, b::Dict{Partition, BigInt})
    out = Dict{Partition, BigInt}()
    for (ea, ca) in a
        for (eb, cb) in b
            exp = map(+, ea, eb)
            out[exp] = get(out, exp, big(0)) + ca * cb
        end
    end
    return _cleanup!(out)
end

# Input: partition mu and number of variables.
# Output: raw monomial dictionary for the Schur polynomial s_mu.
# Method: Jacobi-Trudi determinant expanded through h_q polynomials.
function _schur_raw_polynomial(mu::Partition, num_variables::Int)
    mu = trim_trailing_zeroes(mu)
    key = (mu, num_variables)
    cached = get(_schur_raw_cache, key, nothing)
    cached !== nothing && return cached

    if Base.length(mu) > num_variables
        out = Dict{Partition, BigInt}()
        _schur_raw_cache[key] = out
        return out
    end

    ell = Base.length(mu)
    if ell == 0
        out = _poly_one(num_variables)
        _schur_raw_cache[key] = out
        return out
    end

    out = Dict{Partition, BigInt}()
    for perm in _permutations_list(ell)
        sign = _perm_sign(perm)
        prod_poly = _poly_one(num_variables)
        valid = true
        for i in 1:ell
            q = mu[i] - i + perm[i]
            hq = _complete_homogeneous_raw(q, num_variables)
            if isempty(hq)
                valid = false
                break
            end
            prod_poly = _poly_mul(prod_poly, hq)
        end
        valid && _poly_add_scaled!(out, prod_poly, sign)
    end

    _cleanup!(out)
    _schur_raw_cache[key] = out
    return out
end

# Input: partition mu and number of variables.
# Output: dictionary eta => coefficient of m_eta in s_mu.
# Purpose: gives Kostka transition data from Schur to monomial basis.
function _schur_monomial_coeffs(mu::Partition, num_variables::Int)
    mu = trim_trailing_zeroes(mu)
    key = (mu, num_variables)
    cached = get(_schur_monomial_cache, key, nothing)
    cached !== nothing && return cached

    raw = _schur_raw_polynomial(mu, num_variables)
    out = Dict{Partition, BigInt}()
    for (gamma, coeff) in raw
        eta = trim_trailing_zeroes(sort(collect(gamma), rev=true))
        if haskey(out, eta)
            out[eta] == coeff || error("Schur polynomial expansion is not symmetric at $eta")
        else
            out[eta] = coeff
        end
    end

    _schur_monomial_cache[key] = out
    return out
end

# Input: homogeneous degree and number of variables.
# Output: (basis, rows), where basis is the partition list and rows is the
# exact matrix for s_mu = sum_eta K_{mu,eta} m_eta.
function _transition_matrix(degree::Int, num_variables::Int)
    key = (degree, num_variables)
    cached = get(_transition_cache, key, nothing)
    cached !== nothing && return cached

    basis = partitions(degree; max_length=num_variables)
    n = Base.length(basis)
    matrix = [zero(Rat) for _ in 1:n, _ in 1:n]
    for (j, mu) in enumerate(basis)
        coeffs = _schur_monomial_coeffs(mu, num_variables)
        for (i, eta) in enumerate(basis)
            matrix[i, j] = get(coeffs, eta, big(0)) // big(1)
        end
    end

    rows = [Vector{Rat}(matrix[i, :]) for i in 1:n]
    result = (basis, rows)
    _transition_cache[key] = result
    return result
end

# Input: square rational matrix A and rational vector b.
# Output: rational vector x solving A*x=b.
# Throws: error if A is singular.
function _solve_rational_system(A::Vector{Vector{Rat}}, b::Vector{Rat})
    n = Base.length(b)
    M = [Rat[A[i]...; b[i]] for i in 1:n]

    for col in 1:n
        pivot = findfirst(r -> !iszero(M[r][col]), col:n)
        pivot === nothing && error("singular transition matrix")
        pivot_row = pivot + col - 1
        if pivot_row != col
            M[col], M[pivot_row] = M[pivot_row], M[col]
        end

        pivot_value = M[col][col]
        for j in col:(n + 1)
            M[col][j] /= pivot_value
        end

        for r in 1:n
            r == col && continue
            factor = M[r][col]
            iszero(factor) && continue
            for j in col:(n + 1)
                M[r][j] -= factor * M[col][j]
            end
        end
    end

    return Rat[M[i][n + 1] for i in 1:n]
end

# Input: monomial-symmetric expansion Gamma, total degree, and variable count.
# Output: Schur coefficient dictionary h_mu satisfying Gamma_eta =
# sum_mu h_mu K_{mu,eta}.
function monomial_symmetric_to_schur(Gamma, degree::Integer, num_variables::Integer)
    degree = Int(degree)
    num_variables = Int(num_variables)
    degree == 0 && return Dict{Partition, Rat}(() => get(Gamma, (), zero(Rat)))

    basis, transition = _transition_matrix(degree, num_variables)
    b = Rat[get(Gamma, eta, zero(Rat)) for eta in basis]
    coeffs = _solve_rational_system(transition, b)

    out = Dict{Partition, Rat}()
    for (mu, coeff) in zip(basis, coeffs)
        !iszero(coeff) && (out[mu] = coeff)
    end
    return out
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
function compute_tail_polynomial(D, d; verbose::Bool=true, io::IO=stdout)
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

# Input: sparse polynomial poly, variable indices i,j, and target exponent tuple.
# Output: sparse polynomial for poly*(x_i-x_j), discarding monomials above target.
# Role: helper for the optional three-way subspace coefficient extraction.
function _multiply_linear_difference(poly, i::Int, j::Int, target::Partition)
    out = Dict{Partition, BigInt}()
    for (exp, coeff) in poly
        exp_i = collect(exp)
        exp_i[i] += 1
        if exp_i[i] <= target[i]
            key = Tuple(exp_i)
            out[key] = get(out, key, big(0)) + coeff
        end

        exp_j = collect(exp)
        exp_j[j] += 1
        if exp_j[j] <= target[j]
            key = Tuple(exp_j)
            out[key] = get(out, key, big(0)) - coeff
        end
    end
    return _cleanup!(out)
end

# Input: sparse polynomial poly, variable indices xi,zj, multiplicity, target.
# Output: poly multiplied by the truncated expansion of
# (1 - x_xi - z_zj)^(-multiplicity).
# Role: helper for the optional three-way subspace formula.
function _multiply_pair_segre_factor(poly, xi::Int, zj::Int, multiplicity::Int, target::Partition)
    out = Dict{Partition, BigInt}()
    for (exp, coeff) in poly
        max_u = target[xi] - exp[xi]
        max_v = target[zj] - exp[zj]
        for u in 0:max_u, v in 0:max_v
            q = u + v
            series_coeff = binomial(big(multiplicity + q - 1), q) * binomial(big(q), u)
            values = collect(exp)
            values[xi] += u
            values[zj] += v
            key = Tuple(values)
            out[key] = get(out, key, big(0)) + coeff * series_coeff
        end
    end
    return _cleanup!(out)
end

"""
    degree_subspace_variety_N3(D, d)

Compute the degree of a three-way tensor train variety using the standard
subspace-variety coefficient formula.

The input must have the form

```julia
D = [1, r, s, 1]
d = [a, b, c]
```

The variety is the subspace variety obtained by choosing an `r`-dimensional
subspace in the first tensor factor and an `s`-dimensional subspace in the
third tensor factor. The degree is computed as the coefficient integral over
`Gr(r,a) x Gr(s,c)`.

This is an optional verification route. The default degree computation uses
the Schur-Weingarten recursion instead.
"""
function degree_subspace_variety_N3(D, d)
    D = Int[x for x in D]
    d = Int[x for x in d]
    N = _validate_TT_input(D, d)
    N == 3 || throw(ArgumentError("degree_subspace_variety_N3 expects length(d) == 3"))

    a, b, c = d
    r = D[2]
    s = D[3]
    r <= a || throw(ArgumentError("expected D_1 <= d_1"))
    s <= c || throw(ArgumentError("expected D_2 <= d_3"))

    target = Tuple(vcat([a - i for i in 1:r], [c - j for j in 1:s]))
    poly = Dict{Partition, BigInt}(ntuple(_ -> 0, r + s) => big(1))

    if r > 1
        for i in 1:(r - 1), j in (i + 1):r
            poly = _multiply_linear_difference(poly, i, j, target)
        end
    end

    if s > 1
        for i in 1:(s - 1), j in (i + 1):s
            poly = _multiply_linear_difference(poly, r + i, r + j, target)
        end
    end

    for i in 1:r, j in 1:s
        poly = _multiply_pair_segre_factor(poly, i, r + j, b, target)
    end

    return get(poly, target, big(0))
end


"""
    degree_TT_variety(D, d; method=:schur_weingarten, reduce=true)

Compute the projective degree of the tensor train variety with rank signature
`D = [D_0, ..., D_N]` and mode dimensions `d = [d_1, ..., d_N]`.

By default, vacuous full-rank boundary modes are first compressed and the exact
Schur-Weingarten recursion is applied to the reduced signature. Passing
`method=:subspace` uses `degree_subspace_variety_N3` after reduction, and
therefore requires the reduced input to be three-way.

Returns the degree as a `BigInt`.
"""
function degree_TT_variety(D, d; method::Symbol=:schur_weingarten, reduce::Bool=true, verbose::Bool=false)
    D = Int[x for x in D]
    d = Int[x for x in d]
    _validate_TT_input(D, d)
    method in (:schur_weingarten, :subspace) ||
        throw(ArgumentError("method must be :schur_weingarten or :subspace"))

    if reduce
        reduced_D, reduced_d = compress_vacuous_boundary_modes(D, d)
        if (reduced_D, reduced_d) != (D, d)
            return degree_TT_variety(reduced_D, reduced_d; method=method, reduce=false, verbose=verbose)
        end
    end

    N = Base.length(d)

    if method == :subspace
        N == 3 || throw(ArgumentError("method=:subspace requires a three-way signature after reduction"))
        return degree_subspace_variety_N3(D, d)
    end

    P_coeffs = compute_tail_polynomial(D, d, verbose=verbose)
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
\\sum_{r=1}^{N-1} D_r(D_{r-1}d_r-D_r) + D_{N-1}d_N - 1
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

end
