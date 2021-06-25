# betaChiralEFT.jl

**Julia code for Chiral EFT interactions**


[README in Japanese (日本語はこちら)](https://github.com/SotaYoshida/betaChiralEFT.jl/blob/main/README_JP.md)

## Quick Start

1. Prepare Julia environment

2. Download betaChiralEFT.jl  
```
$git clone https://github.com/SotaYoshida/betaChiralEFT.jl
```

3. Install other Julia Packages
```
$cd betaChiralEFT.jl
$julia src/package_install.jl
```

4. Run the sample script
```
$julia -t 12 sample_run.jl
```
One can specify the number of execution threads like this.  
(The -t/--threads command line argument requires at least Julia >= 1.5.)  

All the parameters are given in "parameters.jl".
Edit it by your hand and play with it!


## What is supported now?  
 
* Free-space NN interaction (Entem&Machleidt type, up to N3LO) => Ref[1]  
    -The author confirmed that results (up to N3LO) agree with those by widely used Fortran codes by R. Machleidt.
* Valence space NN interaction (up to NLOvs, not optimized yet) => Ref[2]

* One can get NN interaction in (relative) momentum space or TBMEs in HO basis (by Talmi-Moshinsky trans. [3]).

## To be implemented:  

* Higher order terms (N4LO, N3LOvs, etc.)

* Density-dependent effective NN interaction from 3NFs

* SRG

* and many more

Any contributions, suggestions, feedbacks, and "Issues&Pull requests" are welcomed.

## Licence  

MIT License, Copyright (c) 2021 Sota Yoshida

## References

[1] [R. Machleidt, D. R. Entem, Phys.Rept.503:1-75,2011](https://www.sciencedirect.com/science/article/pii/S0370157311000457)  
[2] [L. Huth et al., Phys. Rev. C 98, 044301 (2018)](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.98.044301)  
[3] [G.P.Kamuntavičius et al.,](https://www.sciencedirect.com/science/article/pii/S0375947401011010)  


## How to cite  

in prep.
