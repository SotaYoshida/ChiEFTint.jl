# ChiralEFTint.jl (Beta)

**Chiral EFT(カイラル摂動論)による核子間相互作用を計算するJuliaコード**

現在公開に向けた確認作業を行っています。  
そのため、ソースコードは当面リクエストによる公開としていますが、必要な方は気軽にコンタクトしてください。

## What are supported now?  
 
* Entem&Machleidt型のNN相互作用 => Ref[1]  
    - N3LOまで、流通するR.Machleidt氏のFortranコードと一致
    - 運動量空間(部分波展開)のポテンシャルもしくはHO基底のTBME(Talmi-Moshinsky 変換[2])が得られます

*  Similarity Renormalization Group (SRG)  
    - NN-only (induced 3NFは未実装)
    - 遅い(ODE solverを導入すれば高速化の余地あり)
    
* バレンス系の演算子を含むNN相互作用  => Ref[3]  
    - "NLOvs"まで (チェックや最適化はこれから)
    - いわゆる"有効相互作用"とは異なります (詳しくはRef[3]

* 有効2体化3体力 => Ref[4]


## To be implemented:  

* 高次の項 (N4LO, N3LOvs, etc.)

* その他多数

あらゆる寄与,改良,ご提案,ご意見,"Issues&Pull requests"を歓迎します

## Licence  

MIT License, Copyright (c) 2021 Sota Yoshida

## References

[1] [R. Machleidt, D. R. Entem, Phys.Rept.503 (2011) 1-75](https://www.sciencedirect.com/science/article/pii/S0370157311000457)  
[2] [G.P.Kamuntavičius et al., Nucl. Phys. A 695 (2001) 191-201](https://www.sciencedirect.com/science/article/pii/S0375947401011010)  
[3] [L. Huth et al., Phys. Rev. C 98, 044301 (2018)](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.98.044301)  
[4] [M.Kohno, Phys. Rev. C 88, 064005 (2013)](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.88.064005), [Erratum: Phys. Rev. C 96, 059903 (2017)](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.96.059903)  


## How to cite  

準備中
