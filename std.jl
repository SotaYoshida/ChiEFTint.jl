include("./src/chiEFTint.jl")
using .chiEFTint

sps,dictsps,dictTBMEs = readsnt("usdb.snt",18)

sps,dictsps,dictTBMEs = readsnt("gxpf1a.snt",42)
@time jj_std(sps,dictsps,dictTBMEs;fname="gxpf1a")

sps,dictsps,dictTBMEs = readsnt("nn_em500n3lo_A28hw20_pf.snt",42)
@time jj_std(sps,dictsps,dictTBMEs;fname="nn")
sps,dictsps,dictTBMEs = readsnt("nn3nf_em500n3lo_A28hw20_pf.snt",42)
@time jj_std(sps,dictsps,dictTBMEs;fname="nn3nf")
