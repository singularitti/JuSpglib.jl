module JuSpglib

using Reexport

export spglib_version

include("../deps/build.jl")

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

include("DataStructure.jl")
@reexport using .DataStructure

include("CAPIs.jl")
@reexport using .CAPIs

end # module
