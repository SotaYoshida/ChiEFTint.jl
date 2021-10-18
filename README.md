# ChiEFTint.jl (Beta) 

**Julia code for Chiral EFT interactions**

Note: The source code is currently undergoing final checks and is only available by request.  
Please feel free to contact me!

[README in Japanese (日本語はこちら)](https://github.com/SotaYoshida/ChiEFTint.jl/blob/main/README_JP.md)

## Quick Start

1. Prepare Julia environment

2. Clone chiEFTint.jl  
```
$git clone https://github.com/SotaYoshida/ChiEFTint.jl
```

3. Install other Julia Packages
```
$cd chiEFTint.jl
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

If you are new to Julia and have some troubles regarding "PyCall" or "pyimport",
plase comment out the ```using PyCall``` in ```src/chiEFTint.jl``` and ```src/bayesopt.jl```.

## What are supported now?  
 
* Free-space NN interaction (Entem&Machleidt type, up to N3LO) [1]  
    - The author confirmed that results (up to N3LO) agree with those by widely used Fortran codes by R. Machleidt.
    - One can get NN interaction in (relative) momentum space or TBMEs in HO basis (by Talmi-Moshinsky trans. [2]).

* Similarity Renormalization Group (SRG) to soften the NN interactions 
    - Only the two-body part (induced forces are not taken into account)
    - Slow (due to naive implemantation, could be resolved by introducing efficient ODE solvers)

* Valence space NN interaction [3]
    - up to NLOvs, not optimized yet
    - different from the so-called "effective interaction" for a valence space

* Density-dependent effective 2NF from three-nucleon force (3NF)  [4]
  

## To be implemented:  

* NNN (+SRG in three-body space)

* Higher order terms (N4LO, N3LOvs, etc.)

* and many more

Any contributions, suggestions, feedbacks, and "Issues&Pull requests" are welcomed.

## Licence  

MIT License, Copyright (c) 2021 Sota Yoshida

## References

[1] [R. Machleidt, D. R. Entem, Phys.Rept.503 (2011) 1-75](https://www.sciencedirect.com/science/article/pii/S0370157311000457)  
[2] [G.P.Kamuntavičius et al., Nucl. Phys. A 695 (2001) 191-201](https://www.sciencedirect.com/science/article/pii/S0375947401011010)  
[3] [L. Huth et al., Phys. Rev. C 98, 044301 (2018)](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.98.044301)  
[4] [M.Kohno, Phys. Rev. C 88, 064005 (2013)](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.88.064005), [Erratum: Phys. Rev. C 96, 059903 (2017)](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.96.059903) 

## How to cite  

in prep.
