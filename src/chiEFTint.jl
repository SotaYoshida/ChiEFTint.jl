module chiEFTint

using Base.Threads
using Combinatorics
using LaTeXStrings
using LinearAlgebra
using Printf
using PyCall
using SpecialFunctions
using TimerOutputs
using ThreadPools
@pyimport matplotlib.pyplot as plt
using WignerSymbols
include("plt_misc.jl")
include("op_func.jl")
include("main.jl")
export make_chiEFTint
export hw_formula
include("../parameters.jl")
end




