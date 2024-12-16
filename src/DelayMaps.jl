module DelayMaps

using RadiiPolynomial

include("DelayStructures/delay_structs.jl")
    export DelayMap, ExplicitDelayMap

include("rotation.jl")
    export f

end
