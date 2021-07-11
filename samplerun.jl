include("./src/chiEFTint.jl")
using .chiEFTint

#make_chiEFTint(;is_plot=true)
make_chiEFTint()

## for profiling
#using StatProfilerHTML
#@profilehtml make_chiEFTint()
