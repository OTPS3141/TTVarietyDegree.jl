#!/usr/bin/env julia

using Pkg
const PACKAGE_ROOT = dirname(@__DIR__)
Pkg.activate(PACKAGE_ROOT; io=devnull)

using TTVarietyDegree

const SOURCE_TABLE_DIR = joinpath(
    dirname(PACKAGE_ROOT),
    "intersection_degree_computations_and_parametrizations",
    "degree_tables",
)

const OUTPUT_TABLE_DIR = joinpath(PACKAGE_ROOT, "computed_tables")

const TEMPLATE_FILES = Dict(
    2 => "TT_variety_degrees_N_2.csv",
    3 => "TT_variety_degrees_metric_N_3.csv",
    4 => "TT_variety_degrees_metric_N_4.csv",
)

function parse_csv_line(line::AbstractString)
    cells = String[]
    buf = IOBuffer()
    in_quotes = false
    i = firstindex(line)

    while i <= lastindex(line)
        c = line[i]
        if c == '"'
            if in_quotes && i < lastindex(line) && line[nextind(line, i)] == '"'
                print(buf, '"')
                i = nextind(line, i)
            else
                in_quotes = !in_quotes
            end
        elseif c == ',' && !in_quotes
            push!(cells, String(take!(buf)))
        else
            print(buf, c)
        end
        i = nextind(line, i)
    end

    push!(cells, String(take!(buf)))
    return cells
end

function csv_quote(value)
    text = string(value)
    return "\"" * replace(text, "\"" => "\"\"") * "\""
end

function parse_int_vector(text::AbstractString)
    return [parse(Int, m.match) for m in eachmatch(r"-?\d+", text)]
end

function read_table_signature_template(path::AbstractString)
    rows = [parse_csv_line(line) for line in readlines(path) if !isempty(strip(line))]
    isempty(rows) && error("empty template table: $path")

    d_values = [parse_int_vector(cell) for cell in rows[1][2:end]]
    D_values = [parse_int_vector(row[1]) for row in rows[2:end]]
    return D_values, d_values
end

function is_admissible_TT_signature(D, d)
    length(D) == length(d) + 1 || return false
    D[1] == 1 && D[end] == 1 || return false
    all(>(0), D) && all(>(0), d) || return false

    N = length(d)
    for r in 1:(N - 1)
        left_dim = prod(d[1:r])
        right_dim = prod(d[(r + 1):N])
        D[r + 1] <= min(left_dim, right_dim) || return false
        D[r + 1] <= D[r] * d[r] || return false
        D[r + 1] <= d[r + 1] * D[r + 2] || return false
    end

    return true
end

function table_degree(D, d; method::Symbol)
    is_admissible_TT_signature(D, d) || return big(0)
    return degree_TT_variety(D, d; method=method)
end

function output_name(N::Int, method::Symbol)
    method_suffix = method == :subspace ? "subspace" : "schur_weingarten"
    return "TT_variety_degrees_$(method_suffix)_N_$(N).csv"
end

function write_degree_table(N::Int; method::Symbol=:schur_weingarten)
    haskey(TEMPLATE_FILES, N) || error("no template configured for N = $N")
    template_path = joinpath(SOURCE_TABLE_DIR, TEMPLATE_FILES[N])
    D_values, d_values = read_table_signature_template(template_path)

    mkpath(OUTPUT_TABLE_DIR)
    out_path = joinpath(OUTPUT_TABLE_DIR, output_name(N, method))

    open(out_path, "w") do io
        println(io, join(vcat(["D"], [csv_quote(d) for d in d_values]), ","))

        for D in D_values
            entries = String[csv_quote(D)]
            for d in d_values
                push!(entries, string(table_degree(D, d; method=method)))
            end
            println(io, join(entries, ","))
        end
    end

    return out_path
end

function parse_method(args)
    if "--subspace" in args || "--method=subspace" in args
        return :subspace
    end
    return :schur_weingarten
end

function parse_N_values(args)
    values = Int[]
    for arg in args
        startswith(arg, "--") && continue
        push!(values, parse(Int, arg))
    end
    return isempty(values) ? [2, 3, 4] : values
end

method = parse_method(ARGS)
N_values = parse_N_values(ARGS)

for N in N_values
    path = write_degree_table(N; method=method)
    println(path)
end
