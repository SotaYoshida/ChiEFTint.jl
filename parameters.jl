# global valiables are specified in this file
# See https://docs.julialang.org/en/v1/manual/variables-and-scoping/#Constants
# for more details about "const"

const n_mesh = 40
const emax = 2
const Nnmax= 20
const pmax_fm = 4
const chi_order = 3 #0:LO 1:NLO 2:NNLO 3:N3LO
const calc_NN = true
const calc_3N = false #density-dependent 3NF

"""
hbar-omega can be specify by hand
const hw = 20.0 
or using formula
const hw = hw_formula(Anum,fnum)
#fnum=1: 41A^(-1/3) MeV  
#funm=2: 45A^(-1/3)-25A^(-2/3) MeV (used for sd-shell)
"""
const Anum = 4
const hw = 20.0
#Anum = 28; const hw = hw_formula(Anum,2)
#const hw = hw_formula(Anum,1)

## SRG evolution (srg_lambda is in fm^{-1}
#const srg = false
const srg = true
const srg_lambda = 2.0

## file name and format for TBME 
tx = "bare";if srg; tx ="srg"*string(srg_lambda);end
#const tbme_fmt = "snt"
const tbme_fmt = "myg"
const fn_tbme = "tbme_em500n3lo_"*tx*"hw"*string(round(Int64,hw))*"emax"*string(emax)*"."*tbme_fmt


const calc_monopole = false
const calc_std = false

## target {[n,l,j]} for TBMEs
const target_nlj=[] # no-core
#const target_nlj=[[1,0,1],[0,2,3],[0,2,5]] # sd-shell part
#const target_nlj=[[1,1,1],[1,1,3],[0,3,5],[0,3,7]] # pf-shell part

##for valence space operators
#-(v_chi_order=1): vsNLO
#-(v_chi_order=3): vsN3LO (not implemnted)
const v_chi_order = 0 # 0: free-space only
const n_mesh_P = 10
const Pmax_fm = 3


## constants (no need to modify basically...)
const jmax = 6
const lmax = jmax+1
const lcmax = jmax+1 
const Mp = 938.272
const Mn = 939.5653
const Mm = (Mp+Mn)/2
const Ms = [Mp,Mm,Mn]
const Lambchi = 500 # cutoff Lambda
const itts = [-2,0,2]
const hc = 197.327053
const gA = -1.29
const Fpi = 92.4
const hc2 = hc^2
const hc3 = hc^3
const hc4 = hc^4
const gA2 = gA^2
const gA4 = gA^4
const mpis = [139.5702,134.9766,139.5702]
### Fermi momentum for density-dependent NN from 3NFs
const kF = 1.35
