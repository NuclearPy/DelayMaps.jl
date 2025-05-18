module DelayMaps

using RadiiPolynomial

include("DelayStructures/delay_structs.jl")
    export DelayMap, ExplicitDelayMap, Jacobian

include("DelayStructures/ExplicitDelayMap.jl")
    export ExplicitDelayMap

include("rotation.jl")
    export f

end
