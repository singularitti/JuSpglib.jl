"""
# module CAPI



# Examples

```jldoctest
julia>
```
"""
module CAPI

using CoordinateTransformations
using DataStructures: counter

export symmetry_operations, international_symbol, schoenflies_symbol

# Load library after build
if isfile(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
    include(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
else
    error("Spglib not properly installed. Please run Pkg.build(\"JuSpglib\")")
end

"""
version number of the underlying C-API
"""
const version = VersionNumber(ccall((:spg_get_major_version, spglib), Cint, ()),
  ccall((:spg_get_minor_version, spglib), Cint, ()),
  ccall((:spg_get_micro_version, spglib), Cint, ()),
)

function symmetry_operations(lattice::AbstractMatrix, positions::AbstractMatrix, types::AbstractVector; symprec::Real = 1e-8)
    if size(positions, 2) != length(types)
        error("Number of positions and types do not match")
    end
    if size(positions, 1) != 3
        error("Operating in 3D here")
    end
    maxsize::Integer = 52
    rotations = Array{Cint}(undef, 3, 3, maxsize)
    translations = Array{Cdouble}(undef, 3, maxsize)

    type_indices = convert(Vector{Cint}, (collect ∘ values ∘ counter)(types))
    clattice = convert(Matrix{Cdouble}, lattice)
    cpositions = convert(Matrix{Cdouble}, positions)
    numops = ccall((:spg_get_symmetry, spglib), Cint,
        (Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        rotations, translations, maxsize, clattice, cpositions, type_indices, length(type_indices), symprec)
    if numops == 0
        error("Could not determine symmetries")
    end
    [AffineMap(transpose(rotations[:, :, i]), translations[:, i]) for i in 1:numops]
end

function international_symbol(lattice::AbstractMatrix, positions::AbstractMatrix, types::AbstractVector; symprec::Real = 1e-8)
    result = zeros(Cchar, 11)

    type_indices = convert(Vector{Cint}, (collect ∘ values ∘ counter)(types))
    clattice = convert(Matrix{Cdouble}, lattice)
    cpositions = convert(Matrix{Cdouble}, positions)
    numops = ccall((:spg_get_international, spglib), Cint,
        (Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        result, clattice, cpositions, type_indices, length(type_indices), symprec)

    if numops == 0
        error("Could not determine internation symbol")
    end
    join(convert(Array{Char}, result[1:findfirst(iszero, result) - 1]))
end

function schoenflies_symbol(lattice::AbstractMatrix, positions::AbstractMatrix, types::AbstractVector; symprec::Real = 1e-8)
    result = zeros(Cchar, 11)

    type_indices = convert(Vector{Cint}, (collect ∘ values ∘ counter)(types))
    clattice = convert(Matrix{Cdouble}, lattice)
    cpositions = convert(Matrix{Cdouble}, positions)
    numops = ccall((:spg_get_schoenflies, spglib), Cint,
        (Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        result, clattice, cpositions, type_indices, length(type_indices), symprec)

    if numops == 0
        error("Could not determine Schoenflies symbol")
    end
    join(convert(Array{Char}, result[1:findfirst(iszero, result) - 1]))
end

end