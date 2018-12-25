module JuSpglib

using Reexport

include("../deps/build.jl")
include("CAPIs.jl")
@reexport using .CAPIs

end # module
