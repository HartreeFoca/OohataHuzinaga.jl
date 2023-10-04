module OohataHuzinaga
using LinearAlgebra
using SpecialFunctions
using TimerOutputs
using BasisSets

include("auxiliary.jl")
include("overlap.jl")
include("kinetic.jl")
include("boys.jl")
include("attraction.jl")
include("repulsion.jl")
include("hartreefock.jl")

export doublefactorial
export gaussianproduct
export normalization

export overlap
export kinetic
export boys
export attraction
export attraction_2
export repulsion

export computeenergy

end
