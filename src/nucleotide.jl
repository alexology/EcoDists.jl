const _DNA_ALPHABET = codeunits("ACGTRYSWKMBDHVN-")
const _DNA_CANONICAL = codeunits("ACGT")

function _dna_bytes(sequence::AbstractString)
    isascii(sequence) || throw(ArgumentError("aligned DNA must use ASCII nucleotide characters"))
    bytes = Vector{UInt8}(codeunits(uppercase(sequence)))
    all(b -> b in _DNA_ALPHABET, bytes) ||
        throw(ArgumentError("aligned DNA must contain only A/C/G/T, IUPAC ambiguity codes, or '-'"))
    return bytes
end

function _alignment_columns(mask, length)
    isnothing(mask) && return collect(1:length)
    values = collect(mask)
    if eltype(values) <: Bool
        Base.length(values) == length || throw(DimensionMismatch("alignment mask must match sequence length"))
        columns = findall(values)
    else
        all(x -> x isa Integer && !(x isa Bool), values) ||
            throw(ArgumentError("alignment mask must contain booleans or integer column indices"))
        columns = Int.(values)
        issorted(columns) && allunique(columns) || throw(ArgumentError("mask indices must be sorted and unique"))
        all(i -> 1 <= i <= length, columns) || throw(ArgumentError("alignment mask index out of bounds"))
    end
    isempty(columns) && throw(ArgumentError("alignment mask must retain at least one column"))
    return columns
end

"""
    NucleotidePDistance(sequences; mask=nothing, gaps=:error, ambiguities=:error)

Fit one fixed mask to a nonempty collection of aligned homologous DNA strings.
The distance is the fraction of mismatched retained columns. Alignment lengths
must agree and be positive. Lowercase is accepted. Invalid characters (including
`?`, `.`, `U`, and whitespace) are rejected even outside the mask.

`mask` is a boolean vector or sorted unique column indices. At retained columns,
`gaps=:error` rejects `-`, `:mask` removes its column across ALL training sequences,
and `:state` treats it as a fifth state. `ambiguities=:error` rejects IUPAC
ambiguities; `:mask` removes their columns across ALL training sequences.
No pairwise deletion occurs. At least one column must remain.

The fitted `Distances.Metric` accepts aligned strings with this same mask.
New sequences with gaps/ambiguities still present at retained columns are rejected
(except gaps with `:state`); comparisons never change the mask. Removing columns
can make distinct original strings have distance zero (a pseudometric).
"""
struct NucleotidePDistance <: Distances.Metric
    original_length::Int
    columns::Vector{Int}
    gap_state::Bool
end

function NucleotidePDistance(sequences; mask=nothing, gaps=:error, ambiguities=:error)
    gaps in (:error, :mask, :state) || throw(ArgumentError("gaps must be :error, :mask, or :state"))
    ambiguities in (:error, :mask) || throw(ArgumentError("ambiguities must be :error or :mask"))
    sequences isa AbstractString && throw(ArgumentError("provide a collection of DNA strings"))
    seqs = collect(sequences)
    isempty(seqs) && throw(ArgumentError("at least one aligned sequence is required"))
    all(x -> x isa AbstractString, seqs) || throw(ArgumentError("sequences must be nonmissing DNA strings"))
    encoded = [_dna_bytes(s) for s in seqs]
    length0 = length(first(encoded))
    length0 > 0 || throw(ArgumentError("aligned sequences cannot be empty"))
    all(s -> length(s) == length0, encoded) || throw(DimensionMismatch("aligned sequence lengths differ"))
    columns = _alignment_columns(mask, length0)
    retained = Int[]
    for j in columns
        has_gap = any(s -> s[j] == UInt8('-'), encoded)
        has_ambiguity = any(s -> !(s[j] in _DNA_CANONICAL) && s[j] != UInt8('-'), encoded)
        has_gap && gaps == :error && throw(ArgumentError("gap in retained alignment column $j"))
        has_ambiguity && ambiguities == :error && throw(ArgumentError("ambiguous base in retained alignment column $j"))
        (has_gap && gaps == :mask || has_ambiguity && ambiguities == :mask) && continue
        push!(retained, j)
    end
    isempty(retained) && throw(ArgumentError("no comparable columns remain in the fixed alignment mask"))
    return NucleotidePDistance(length0, retained, gaps == :state)
end

"""
    alignment_mask(d::NucleotidePDistance)

Return a copy of the fixed retained alignment column indices.
"""
alignment_mask(d::NucleotidePDistance) = copy(d.columns)

function _retained_dna(d::NucleotidePDistance, sequence::AbstractString)
    bytes = _dna_bytes(sequence)
    length(bytes) == d.original_length || throw(DimensionMismatch("sequence length differs from fitted alignment"))
    out = bytes[d.columns]
    all(b -> b in _DNA_CANONICAL || (d.gap_state && b == UInt8('-')), out) ||
        throw(ArgumentError("new gap or ambiguity at a retained column; the fitted mask cannot change"))
    return out
end

function (d::NucleotidePDistance)(a::AbstractString, b::AbstractString)
    x, y = _retained_dna(d,a), _retained_dna(d,b)
    return count(i -> x[i] != y[i], eachindex(x)) / length(x)
end
Distances.result_type(::NucleotidePDistance, ::Type, ::Type) = Float64

# Cache the exact p-distance matrix as integer numerators and one common
# denominator. This avoids floating-point costs inside optimal transport.
struct _PDistanceMatrix <: AbstractMatrix{Float64}
    mismatches::Matrix{Int}
    sites::Int
end
Base.size(m::_PDistanceMatrix) = size(m.mismatches)
Base.getindex(m::_PDistanceMatrix, i::Int, j::Int) = m.mismatches[i,j] / m.sites

function _p_distance_matrix(d::NucleotidePDistance, sequences)
    encoded = [_retained_dna(d,s) for s in sequences]
    n = length(encoded)
    counts = zeros(Int,n,n)
    for j in 1:n, i in (j+1):n
        counts[i,j] = counts[j,i] = count(k -> encoded[i][k] != encoded[j][k], eachindex(encoded[i]))
    end
    return _PDistanceMatrix(counts, length(d.columns))
end

# Validate even a singleton's characters (the generic metric implementation
# would fill its diagonal with zero without evaluating the input). Encode each
# sequence once, rather than once for every pair in the distance matrix.
Distances.pairwise(d::NucleotidePDistance,sequences::AbstractVector{<:AbstractString}) =
    Matrix(_p_distance_matrix(d,sequences))
