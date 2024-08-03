module DiscreteDiff

# Write your package code here.
export d1x

include("boundary.jl")

include("center.jl")

include("scheme.jl")

include("compact.jl")

end
