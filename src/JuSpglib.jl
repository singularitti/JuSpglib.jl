module JuSpglib

using Reexport

include("../deps/build.jl")
include("capi.jl")
@reexport using .CAPI

end # module
