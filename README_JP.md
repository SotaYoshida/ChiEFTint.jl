# betaChiralEFT.jl

**Julia code to generate Chiral EFT interactions (up to N3LO for now)**


## What is supported now?  
 
* Entem&Machleidt型のNN相互作用 => Ref[1]  
    -N3LOまで、流通するR.Machleidt氏のFortranコードと一致。
* バレンス系の演算子を含むNN相互作用 ("NLOvs"まで, まだ全然最適化されていません) => Ref[2]

* (relative)なmomentum spaceでのポテンシャルと、殻模型で言うところのTBMEs(Talmi-Moshinsky変換[3])が生成できます

## To be implemented:  

* 高次の項 (N4LO, N3LOvs, etc.)

* 3体力の足を畳んだ有効2体力

* SRG

* その他多数

あらゆる寄与,改良,ご提案,ご意見,"Issues&Pull requests"を歓迎します

## Licence  

MIT License, Copyright (c) 2021 Sota Yoshida

## References

[1] [R. Machleidt, D. R. Entem, Phys.Rept.503:1-75,2011](https://www.sciencedirect.com/science/article/pii/S0370157311000457)  
[2] [L. Huth et al., Phys. Rev. C 98, 044301 (2018)](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.98.044301)  
[3] [G.P.Kamuntavičius et al.,](https://www.sciencedirect.com/science/article/pii/S0375947401011010)  


## How to cite  

準備中
