# parameters
const n_mesh = 30
const jmax = 6
const lmax = jmax+1
const Nnmax= 20
const emax = 2
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
const hc3 = hc^3
const gA2 = gA^2
const gA4 = gA^4
const mpis = [139.5702,134.9766,139.5702]
const chi_order = 3 #0:LO 1:NLO 2:NNLO 3:N3LO

#const target_nlj=[] # all
const target_nlj=[[1,0,1],[0,2,3],[0,2,5]] # when you only need sd-shell part

# hbar-omega can be specify by hand
#const hw = 20.0 
# or using formula: hw_formula(Anum,fnum)
#fnum=1: 41A^(-1/3) MeV  
#funm=2: 45A^(-1/3)-25A^(-2/3) MeV (used for sd-shell)
const Anum = 28
const hw = hw_formula(Anum,2)
println("hw $hw")

# file name and format for TBME
const fn_tbme = "tmp_em500NLOvs_A"*string(Anum)*"hw"*string(round(Int64,hw))*".snt"
const tbme_fmt = "snt" # "myg"/"snt"

#for valence space operators
#-(v_chi_order=1): vsNLO
#-(v_chi_order=3): vsN3LO not implemnted now
const v_chi_order = 1
const n_mesh_P = 25
const Pmax_fm = 3
const lcmax = jmax+1 

