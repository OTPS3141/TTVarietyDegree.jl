#!/usr/bin/env julia

using Pkg
const PACKAGE_ROOT = dirname(@__DIR__)
Pkg.activate(PACKAGE_ROOT; io=devnull)

using TTVarietyDegree

# Edit these two lines for an interactive/file-based computation.
D = [1, 2, 2, 1]
d = [3, 3, 3]

TAIL_WEINGARTEN_LIMIT = 20
SHOW_TAIL_DATA = false
METHOD = :schur_weingarten
REDUCE = true
VERBOSE = false
MAX_RSS_GB = 8.0
MEMORY_POLL_SECONDS = 0.25
MEMORY_ABORT_EXIT_CODE = 99

function print_usage(io=stdout)
    println(io, "Usage:")
    println(io, "  julia --project=. compute_degree.jl \"[1,2,2,1]\" \"[3,3,3]\" [options]")
    println(io)
    println(io, "Options:")
    println(io, "  --tail, --with-tail              also print P(CC*) and f(P)")
    println(io, "  --verbose                        print recursive Schur-Weingarten steps with --tail")
    println(io, "  --subspace, --method=subspace    use the optional three-way subspace formula")
    println(io, "  --non-reduced                    disable vacuous boundary-mode reduction")
    println(io, "  --max-ram-gb=<GB>                set resident-memory watchdog limit")
    println(io, "  --no-ram-limit                   disable the resident-memory watchdog")
end

function parse_int_vector(text::AbstractString)
    cleaned = replace(strip(text), r"^\[" => "", r"\]$" => "")
    isempty(strip(cleaned)) && return Int[]
    return [parse(Int, strip(x)) for x in split(cleaned, ",")]
end

function parse_float_flag(arg::AbstractString, prefix::AbstractString)
    startswith(arg, prefix) || return nothing
    value = strip(arg[(lastindex(prefix) + 1):end])
    isempty(value) && error("Missing numeric value for ", prefix)
    return parse(Float64, value)
end

function current_rss_bytes()
    status_path = "/proc/self/status"
    isfile(status_path) || return nothing

    for line in eachline(status_path)
        if startswith(line, "VmRSS:")
            fields = split(line)
            length(fields) >= 2 || return nothing
            return parse(Int, fields[2]) * 1024
        end
    end

    return nothing
end

function format_bytes(bytes::Integer)
    gib = bytes / 1024.0^3
    return string(round(gib; digits=2), " GiB")
end

function start_memory_watchdog(max_rss_gb)
    max_rss_gb === nothing && return nothing
    max_rss_gb <= 0 && return nothing
    isfile("/proc/self/status") || begin
        println(stderr, "RAM limit disabled: /proc/self/status is unavailable on this system.")
        return nothing
    end

    max_bytes = floor(Int, max_rss_gb * 1024.0^3)

    rss = current_rss_bytes()
    if rss !== nothing && rss > max_bytes
        println(stderr, "Aborting: resident memory is already ", format_bytes(rss),
                ", above limit ", format_bytes(max_bytes), ".")
        println(stderr, "Increase the limit with --max-ram-gb=<GB> or disable it with --no-ram-limit.")
        flush(stderr)
        exit(MEMORY_ABORT_EXIT_CODE)
    end

    return Timer(0; interval=MEMORY_POLL_SECONDS) do timer
        rss = current_rss_bytes()
        rss === nothing && return

        if rss > max_bytes
            close(timer)
            println(stderr)
            println(stderr, "Aborting: resident memory reached ", format_bytes(rss),
                    ", above limit ", format_bytes(max_bytes), ".")
            println(stderr, "Increase the limit with --max-ram-gb=<GB> or disable it with --no-ram-limit.")
            flush(stderr)
            exit(MEMORY_ABORT_EXIT_CODE)
        end
    end
end

function format_schur_polynomial(P)
    terms = String[]
    for lam in sort(collect(keys(P)); rev=true)
        coeff = P[lam]
        push!(terms, string(coeff, " * s_", lam))
    end
    return isempty(terms) ? "0" : join(terms, " + ")
end

function max_tail_weingarten_degree(D, d)
    running = 0
    max_degree = 0
    for r in 1:(length(d) - 1)
        max_degree = max(max_degree, running)
        m_r = D[r] * d[r] - D[r + 1]
        running += D[r + 1] * m_r
    end
    return max_degree
end

flags = Set(arg for arg in ARGS if startswith(arg, "-"))
data_args = [arg for arg in ARGS if !startswith(arg, "-")]

allowed_flags = Set([
    "--tail",
    "--with-tail",
    "--subspace",
    "--method=subspace",
    "--schur-weingarten",
    "--method=schur_weingarten",
    "--method=schur-weingarten",
    "--no-ram-limit",
    "--non-reduced",
    "--verbose",
    "--help",
    "-h"
])
known_flag(arg) = arg in allowed_flags ||
                  startswith(arg, "--max-ram-gb=") ||
                  startswith(arg, "--max-rss-gb=")
unknown_flags = Set(arg for arg in flags if !known_flag(arg))
if !isempty(unknown_flags)
    error("Unknown flag(s): ", join(sort(collect(unknown_flags)), ", "))
end

if "--help" in flags || "-h" in flags
    print_usage()
    exit()
end

if !isempty(intersect(flags, allowed_flags))
    if "--tail" in flags || "--with-tail" in flags
        global SHOW_TAIL_DATA = true
    end
    if "--verbose" in flags
        global VERBOSE = true
    end
    if "--subspace" in flags || "--method=subspace" in flags
        global METHOD = :subspace
    end
    if "--non-reduced" in flags
        global REDUCE = false
    end
    if "--schur-weingarten" in flags ||
       "--method=schur_weingarten" in flags ||
       "--method=schur-weingarten" in flags
        global METHOD = :schur_weingarten
    end
end

if "--no-ram-limit" in flags
    global MAX_RSS_GB = nothing
else
    for arg in flags
        value = parse_float_flag(arg, "--max-ram-gb=")
        if value !== nothing
            global MAX_RSS_GB = value
            continue
        end

        value = parse_float_flag(arg, "--max-rss-gb=")
        if value !== nothing
            global MAX_RSS_GB = value
            continue
        end
    end
end

if length(data_args) == 2
    global D = parse_int_vector(data_args[1])
    global d = parse_int_vector(data_args[2])
elseif length(data_args) != 0
    print_usage(stderr)
    error("expected either zero data arguments or exactly D and d")
end

memory_watchdog = start_memory_watchdog(MAX_RSS_GB)

reduced_D, reduced_d = compress_vacuous_boundary_modes(D, d)
tail_D, tail_d = REDUCE ? (reduced_D, reduced_d) : (D, d)
dim = dimension_TT_variety(D, d)
deg = degree_TT_variety(D, d; method=METHOD, reduce=REDUCE, verbose=VERBOSE)

println("D = ", D)
println("d = ", d)
println("method = ", METHOD)
println("RAM limit = ", MAX_RSS_GB === nothing ? "disabled" : string(MAX_RSS_GB, " GiB RSS"))
if (reduced_D, reduced_d) != (D, d) && REDUCE
    println("reduced D = ", reduced_D)
    println("reduced d = ", reduced_d)
end
println("dimension = ", dim)
println("degree = ", deg)

if SHOW_TAIL_DATA
    tail_degree = max_tail_weingarten_degree(tail_D, tail_d)
    if tail_degree > TAIL_WEINGARTEN_LIMIT
        println("P(CC*) skipped: generic tail expansion would need Weingarten degree ", tail_degree)
        exit()
    end

    P = compute_tail_polynomial(tail_D, tail_d; verbose=VERBOSE)
    k = tail_D[end - 1]
    n = tail_d[end]
    fP = compute_fP(P, k, n)
    println("P(CC*) = ", format_schur_polynomial(P))
    println("f(P) = ", fP)
end
