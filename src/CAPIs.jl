"""
# module CAPIs



# Examples

```jldoctest
julia>
```
"""
module CAPIs

using CoordinateTransformations
using DataStructures: counter

export spglib_version, get_symmetry, get_international, get_schoenflies

# Load library after build
if isfile(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
    include(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
else
    error("Spglib not properly installed. Please run Pkg.build(\"JuSpglib\")")
end

macro lazy_version(cfuncname)
    return :(ccall(($cfuncname, spglib), Cint, ()))
end

"""
This returns version number of spglib.
"""
const spglib_version = VersionNumber(@lazy_version(:spg_get_major_version), @lazy_version(:spg_get_minor_version), @lazy_version(:spg_get_micro_version))

function get_symmetry(lattice::AbstractMatrix, positions::AbstractMatrix, types::AbstractVector; symprec::Real = 1e-8)
    size(positions, 2) != length(types) && throw(DimensionMismatch("The number of positions and atomic types do not match!"))
    size(positions, 1) != 3 && error("Operations in 3D space is supported here!")

    maxsize::Integer = 52
    rotations = Array{Cint}(undef, 3, 3, maxsize)
    translations = Array{Cdouble}(undef, 3, maxsize)

    type_indices = convert(Vector{Cint}, (collect ∘ values ∘ counter)(types))
    clattice = convert(Matrix{Cdouble}, lattice)
    cpositions = convert(Matrix{Cdouble}, positions)
    numops = ccall((:spg_get_symmetry, spglib), Cint,
        (Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        rotations, translations, maxsize, clattice, cpositions, type_indices, length(type_indices), symprec)

    numops == 0 && error("Could not determine symmetries!")

    [AffineMap(transpose(rotations[:, :, i]), translations[:, i]) for i in 1:numops]
end

function get_international(lattice::AbstractMatrix, positions::AbstractMatrix, types::AbstractVector; symprec::Real = 1e-8)
    result = zeros(Cchar, 11)

    type_indices = convert(Vector{Cint}, (collect ∘ values ∘ counter)(types))
    clattice = convert(Matrix{Cdouble}, lattice)
    cpositions = convert(Matrix{Cdouble}, positions)
    numops = ccall((:spg_get_international, spglib), Cint,
        (Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        result, clattice, cpositions, type_indices, length(type_indices), symprec)

    numops == 0 && error("Could not determine the international symbol!")

    join(convert(Array{Char}, result[1:findfirst(iszero, result) - 1]))
end

function get_schoenflies(lattice::AbstractMatrix, positions::AbstractMatrix, types::AbstractVector; symprec::Real = 1e-8)
    result = zeros(Cchar, 11)

    type_indices = convert(Vector{Cint}, (collect ∘ values ∘ counter)(types))
    clattice = convert(Matrix{Cdouble}, lattice)
    cpositions = convert(Matrix{Cdouble}, positions)
    numops = ccall((:spg_get_schoenflies, spglib), Cint,
        (Ptr{Cchar}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cdouble),
        result, clattice, cpositions, type_indices, length(type_indices), symprec)

    numops == 0 && error("Could not determine the Schoenflies symbol!")

    join(convert(Array{Char}, result[1:findfirst(iszero, result) - 1]))
end

end