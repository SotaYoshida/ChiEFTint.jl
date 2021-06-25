# betaChiralEFT.jl

**Julia code to generate Chiral EFT interactions (up to N3LO for now)**


[README in Japanese (日本語はこちら)](https://github.com/SotaYoshida/betaChiralEFT.jl/blob/main/README_JP.md)


For free-space interactions, this code reproduces the interactions with Fortran codes, which are written by Ruprecht Machleidt and widely used in this community.


## What is supported now?  
 
* Free-space NN interaction (Entem&Machleidt type, up to N3LO) => Ref[1]

* Valence space NN interaction (up to NLOvs, not optimized yet) => Ref[2]

## To be implemented:  

* Higher order term for N4LO, N3LOvs

* Density-dependent effective NN interaction from 3NFs

* SRG

* and many more

Any contributions, suggestions, feedbacks, and "Issues&Pull requests" are welcomed.

## Licence  

MIT License, Copyright (c) 2021 Sota Yoshida

## References

[1] [R. Machleidt, D. R. Entem, Phys.Rept.503:1-75,2011](https://www.sciencedirect.com/science/article/pii/S0370157311000457)
[2] [L. Huth et al., Phys. Rev. C 98, 044301 (2018)](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.98.044301)

## How to cite  

in prep.
