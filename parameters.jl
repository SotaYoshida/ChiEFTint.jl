# parameters
const n_mesh = 30
const jmax = 6
const lmax = jmax+1
const Nnmax= 20
const emax = 1
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

# hbar-omega can be specify by hand
#const hw = 20.0 
# or using formula: hw_formula(Anum,fnum)
#fnum=1: 41A^(-1/3) MeV  
#funm=2: 45A^(-1/3)-25A^(-2/3) MeV (used for sd-shell)
const Anum = 28
const hw = hw_formula(Anum,2)

# file name and format for TBME
const fn_tbme = "tbme_em500_NLOvs.dat" 
const tbme_fmt = "snt" # "myg"/"snt"

## LO 10^4 GeV^-2
const C0_1S0  = -0.147167
const C0_3S1  = -0.11897839697675265
const C_CSB = 0.00049950
const C_CIB = 0.00069075
## NLO 10^4 GeV^-4
const C2_3S1  = 0.76
const C2_3P0  = 1.487
const C2_1P1  = 0.656
const C2_3P1  = -0.630
const C2_1S0  = 2.380
const C2_3SD1 = 0.826
const C2_3P2  = -0.538
## NNLO GeV GeV^-1
const c1_NNLO = -0.81
const c2_NNLO =  2.80
const c3_NNLO = -3.20
const c4_NNLO =  5.40
## N3LO 10^4 GeV^-6
const hD_1S0 = -2.545
const D_1S0 = -16.0
const D_1P1 = 5.25
const D_3P2 = 2.295
const D_3P1 = 2.35
const D_3P0 = 0.245
const hD_3S1 = 7.0
const D_3S1 = 6.55
const hD_3SD1 = 2.25
const D_3SD1 = 6.61
const D_3D1 = -2.80
const D_1D2 = -1.77
const D_3D2 = -1.46
const D_3PF2 = -0.465
const D_3D3 = 5.66
## GeV^-2
const d12 = 3.06
const d3 = -3.27
const d5 = 0.45
const d145 = -5.65

#for valence space operators
#v_chi_order=0: free-space
#v_chi_order=1: vsNLO
#v_chi_order=3: vsN3LO not implemnted now
const v_chi_order = 0
const n_mesh_P = 25
const Pmax_fm = 3
const lcmax = jmax+1 
const c_vs_1 =  0.19162
const c_vs_2 = -0.28374
const c_vs_3 =  0.02680
const c_vs_4 = -0.34499
const c_vs_5 = -0.16235
