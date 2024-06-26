module OohataHuzinaga
using LinearAlgebra
using SpecialFunctions
using TimerOutputs
using BasisSets
using StaticArrays

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
export Operator
export oei
export kinetic
export boys
export attraction
export repulsion
export Results

export rhf

end
