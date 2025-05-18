"""
    DelayMap{T <: Function, S <: Function}

Abstract type for all delay maps.
"""
abstract type DelayMap{T<: Function, S <: Function} end

# Delay maps are functions from C[-1,1] -> C[-1,1]
const PHASE_SPACE_INTERVAL_START = -1.0
const PHASE_SPACE_INTERVAL_END = 1.0
const PHASE_SPACE_INTERVAL_LENGTH = PHASE_SPACE_INTERVAL_END - PHASE_SPACE_INTERVAL_START