# parameters
const n_mesh = 40
const emax = 3
const jmax = 6
const lmax = jmax+1
const Nnmax= 20
const pmax_fm = 4
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
const gA2 = gA^2
const gA4 = gA^4
const mpis = [139.5702,134.9766,139.5702]
const chi_order = 3 #0:LO 1:NLO 2:NNLO 3:N3LO

const calc_NN = true
const calc_3N = true

# hbar-omega can be specify by hand
const hw = 20.0 
# or using formula: hw_formula(Anum,fnum)
#fnum=1: 41A^(-1/3) MeV  
#funm=2: 45A^(-1/3)-25A^(-2/3) MeV (used for sd-shell)
const Anum = 28
#const hw = hw_formula(Anum,2)
#println("hw $hw")

# file name and format for TBME
const fn_tbme = "tmp_em500n3lo_A"*string(Anum)*"hw"*string(round(Int64,hw))*".snt"
const tbme_fmt = "snt" # "myg"/"snt"
#const target_nlj=[] # all
#const target_nlj=[[1,0,1],[0,2,3],[0,2,5]] # when you only need sd-shell part
const target_nlj=[[1,1,1],[1,1,3],[0,3,5],[0,3,7]] # when you only need pf-shell part emax must be >= 3

#const fn_tbme = "tbme_em500NLOvs_A"*string(Anum)*"hw"*string(round(Int64,hw))*".dat"
#const fn_tbme = "tbme_em500_2n3n_A"*string(Anum)*"hw"*string(round(Int64,hw))*".dat"
#const tbme_fmt = "myg" # "myg"/"snt"

const calc_monopole = true
const calc_std = true

#for valence space operators
#-(v_chi_order=1): vsNLO
#-(v_chi_order=3): vsN3LO not implemnted now
const v_chi_order = 0 # 0: free-space only
const n_mesh_P = 10
const Pmax_fm = 3
const lcmax = jmax+1 

### density-dependent NN from 3NFs
const kF = 1.35
const pnm = false
