function wigner9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
    s= 0.0
    xmin = max(abs(j1-j9),abs(j2-j6),abs(j4-j8))
    xmax = min(abs(j1+j9),abs(j2+j6),abs(j4+j8))
    @inbounds for x =xmin:xmax
        t  = (-1)^(2*x) * (2*x+1)
        t *= wigner6j(Float64,j1,j4,j7,j8,j9,x)
        t *= wigner6j(Float64,j2,j5,j8,j4,x,j6)
        t *= wigner6j(Float64,j3,j6,j9,x,j1,j2)
        s += t 
    end
    return s
end

#  9j-symbols for specific cases j2=j5=1/2
function s_wigner9j(j1,j3,j4,j6,j7,j8,j9) 
    s= 0.0
    xmin = max(abs(j1-j9),abs(Int(1//2-j6)),abs(j4-j8))
    xmax = min(abs(j1+j9),abs(Int(1//2+j6)),abs(j4+j8))
    @inbounds for x =xmin:xmax
        t  = (-1.0)^(2*x) * (2.0*x+1.0)
        t *= wigner6j(Float64,j1,j4,j7,j8,j9,x)
        t *= wigner6j(Float64,1//2,1//2,j8,j4,x,j6)
        t *= wigner6j(Float64,j3,j6,j9,x,j1,1//2)
        s += t
    end
    return s
end

function TMtrans(dLECs,xr,wr,xrP,wrP,Rnl,RNL,
                 nTBME,infos,izs_ab,
                 Numpn,V12ab,arr_numst,to
                 ;calc_relcm=false)
    Nrmax = 2*emax
    ichan = length(infos)
    nume = zeros(Int64,Nrmax+1,Nrmax+1)
    numknn = [ [ [ zeros(Int64,div(K1,2)+1,div(KK-K1,2)+1) for K1=0:KK]
                 for Lam=0:KK] for KK=0:Nrmax]
    Ndim = [ zeros(Int64,i) for i=1:Nrmax+1]
    Transbk= Float64[]
    
    sp_P5_9j = [ wigner9j(1,1,0,1//2,1//2,0,1//2,1//2,0) wigner9j(1,1,0,1//2,1//2,1,1//2,1//2,1);
                 wigner9j(1,1,2,1//2,1//2,0,1//2,1//2,0) wigner9j(1,1,2,1//2,1//2,1,1//2,1//2,1)]
    cg1s = [clebschgordan(Float64,1,0,1,0,0,0),
            clebschgordan(Float64,1,0,1,0,2,0)]
    nofst = 0
    nljsnt = [ [0,0] ]; deleteat!(nljsnt,1)
    for temax = 0:emax  #2n +l = temax
        for l = temax%2:2:temax
            n = div(temax - l,2)
            jmin = 2*l-1
            if jmin < 1; jmin=1;end
            for j=jmin:2:2*l+1
                push!(nljsnt,[n,l,j])
                nofst += 1
            end
        end
    end
    
    X9,U6 = prepareX9U6(Nrmax,to)
    ##START: Transall in nn_int.f90
    f_mb,g_mb,w_mb = def_fgw()    
    num=0
    for KK=0:Nrmax
        for Lam=0:KK
            nume[KK+1,Lam+1]=num #nume(KK,Lam)=num
            for K1=0:KK
                K2=KK-K1
                if K2 < K1; continue;end
                for N1=0:div(K1,2) 
                    L1=K1-2*N1
                    for N2=0:div(K2,2)
                        L2=K2-2*N2
                        for K3=0:KK
                            K4=KK-K3
                            if K4 < K3;continue;end
                            for N3=0:div(K3,2)#N3=0,K3/2
                                L3=K3-2*N3
                                for N4=0:div(K4,2)
                                    L4=K4-2*N4
                                    if abs(L1-L2) > Lam || (L1+L2) < Lam;continue;end
                                    if 2*N1+L1+2*N2+L2 != KK;continue;end
                                    if abs(L3-L4) > Lam || (L3+L4) < Lam;continue;end
                                    if 2*N3+L3+2*N4+L4 != KK;continue;end
                                    num=num+1
                                    if KK <= Nrmax                                        
                                        fmosh = gmosh(N1,L1,N2,L2,N3,L3,N4,L4,Lam,1.0,f_mb,g_mb,w_mb)
                                        push!(Transbk,fmosh)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    println("dim. of HO brackets = $num")
    for KK=0:Nrmax
        tarr = numknn[KK+1]        
        for Lam=0:KK
            num = 0
            tar2 = tarr[Lam+1]
            for K1=0:KK
                K2=KK-K1
                tar3 = tar2[K1+1]
                if K2 < K1; continue;end
                for N1=0:div(K1,2) # do N2=0,K2/2
                    L1=K1-2*N1
                    for N2=0:div(K2,2)
                        L2=K2-2*N2
                        if abs(L1-L2) >Lam || L1+L2 < Lam;continue;end
                        if 2*N1+L1+2*N2+L2 != KK; continue;end
                        num=num+1
                        tar3[N1+1,N2+1]=num
                    end
                end
            end
            Ndim[KK+1][Lam+1]=num
        end
    end
    ##END: Transall in nn_int.f90
    ##START: Trans_TM in nn_int.f90
    nljdict = Dict(0=>0); delete!(nljdict,0)
    io=open(fn_tbme,"w")
    if tbme_fmt  == "myg"
        println(io,"tza,  a, tzb,  b, tzc,  c, tzd,  d, J12,   <ab|v|cd>,   <ab|p_ip_j|cd>/mhw,    <ab|r_ir_j|cd>mw/h")
    elseif tbme_fmt=="snt"
        if target_nlj!=[]
            nljdict = def_sps_snt(emax,target_nlj)[2]
        end
    end
    @inbounds for ich=1:length(infos)
        izz,ip,Jtot,ndim=infos[ich]
        pnrank = Int(div(izz,2))+2
        izs = izs_ab[ich]        
        vv = zeros(Float64,ndim,ndim)
        p_pro = zeros(Float64,ndim,ndim)
        r_pro = zeros(Float64,ndim,ndim)
        @timeit to "vtrans" @inbounds @qthreads for i = 1:ndim
            t2v=zeros(Int64,2)
            t5v=zeros(Int64,5)
            iza,ia,izb,ib = izs[i]
            @inbounds for j = 1:i
                izc,ic,izd,id= izs[j]                
                v12 = vtrans(pnrank,izz,ip,Jtot,
                             iza,ia,izb,ib,izc,ic,izd,id,
                             nljsnt,V12ab,X9,U6,Nrmax,
                             nume,numknn,Ndim,Transbk,
                             arr_numst,t2v,t5v,to)
                if v_chi_order >= 1
                    vvs = zeros(Float64,5) 
                    NLOvs(dLECs,vvs,xr,wr,xrP,wrP,Rnl,RNL,
                          cg1s,sp_P5_9j,nljsnt,pnrank,ip,
                          X9,U6,t5v,
                          Jtot,iza,ia,izb,ib,
                          izc,ic,izd,id,to)
                    for i=1:5
                        v12 += vvs[i]
                    end
                end
                vv[i, j] = v12;vv[j, i] = v12                
                if calc_relcm # rel-cm part (e.g. NCSM)
                    # not yet implemented
                end   
            end
        end
        if tbme_fmt == "myg"
            write_tbme(io,ndim,izs,Jtot,vv,nljsnt,nljdict;ofst=nofst)
        else
            write_tbme(io,ndim,izs,Jtot,vv,nljsnt,nljdict;ofst=nofst)
        end
    end
    close(io)
    return nothing
end


function Nfac_jj(na,la,ja,nb,lb,jb)
    s =1.0/ sqrt(1.0+delta(na,nb)*delta(la,lb)*delta(ja,jb))
    return s
end

"""
    HObracket(n1,l1,n2,l2,n,l,N,L,Lam;d=1)

To calc. harmonic oscillator brackets.
This function has large overlap to other functions (gmosh, etc.),
and it will be merged in near future...

Reference:
[1] B.Buck& A.C.merchant, Nucl.Phys. A600 (1996) 387-402
[2] G.P.Kamuntavicius et al., Nucl.Phys. A695 (2001) 191-201
<EL,el:Lam|e1l1,e2l2:Lam>d
d=1: two-body

phase is needed to reproduce Tab.1 of [1]
"""
function HObracket(n1,l1,n2,l2,n,l,N,L,Lam;d=1)
    e1 = 2*n1+l1
    e2 = 2*n2+l2
    e  = 2*n+l
    E  = 2*N+L
    t1 = d^((e1-e)/2) * (1+d)^(-(e1+e2)/2)
    phase = (-1)^(N+n+n1+n2)
    sum = 0.0
    for ed = 0:min(e2,e)
        ec = e2-ed
        ea = E-ec
        eb = e1-ea
        for la = ea:-2:0
            na = div(ea-la,2)
            for lb = eb:-2:0
                nb = div(eb-lb,2)
                for lc = ec:-2:0
                    nc = div(ec-l,2)
                    for ld = ed:-2:0
                        nd = div(ed-ld,2)
                        t2 = (-d)^ed * wigner9j(la,lb,l1,lc,ld,l2,L,l,Lam)
                        if t2 == 0.0;continue;end
                        t2 *= Ghob(e1,l1,ea,la,eb,lb) * Ghob(e2,l2,ec,lc,ed,ld)
                        if t2 == 0.0;continue;end
                        t2 *= Ghob(E,L,ea,la,ec,lc) *  Ghob(e,l,eb,lb,ed,ld)
                        if t2 == 0.0;continue;end
                        sum += t2
                    end
                end
            end
        end
    end
    return Float64(t1*sum)*phase
end

function trinomial(n1,na,nb)
    return doublefactorial(n1) /( doublefactorial(na) * doublefactorial(nb))
end

function Ghob(n1,l1,na,la,nb,lb)
    t = hat(la) *hat(lb) * clebschgordan(Float64,la,0,lb,0,l1,0)
    if t ==0.0
        return 0.0        
    else
        if n1-l1<0 || na-la<0|| nb-lb<0;return 0.0;end
        c = trinomial(n1-l1,na-la,nb-lb) * trinomial(n1+l1+1,na+la+1,nb+lb+1)
        c = sqrt(c)
        t *= c
    end
    return t
end

function tri_check(ja,jb,jc)
    TF= true
    if ja+jb < jc ; TF=false;return TF ;end
    if jb+jc < ja ; TF=false;return TF ;end
    if jc+ja < jb ; TF=false;return TF ;end
    if abs(ja-jb) > jc ; TF=false;return TF ;end
    if abs(jb-jc) > ja ; TF=false;return TF ;end
    if abs(jc-ja) > jb ; TF=false;return TF ;end
    return TF
end

# Naive translation of gmosh in Fortran
# gmosh is to prepare Mochinsky brakets
function gmosh(n,l,nc,lc,n1,l1,n2,l2,lr,d,
                f_mb,g_mb,w_mb)
    ret = 0.0
    if n+n+nc+nc+l+lc-n1-n1-n2-n2-l1-l2 != 0 ;return ret;end
    if l+lc-lr < 0 || l1+l2-lr < 0 ;return ret;end
    if abs(l-lc)-lr > 0 || abs(l1-l2)-lr > 0  ;return ret;end
    dl=log(d)
    d1l=log(d+1.0)

    bb =  f_mb[n1+1]+f_mb[n2+1]+f_mb[n+1]-f_mb[nc+1]
    bb += g_mb[n1+l1+1]+g_mb[n2+l2+1]-g_mb[n+l+1]-g_mb[nc+lc+1]
    ba =  w_mb[l1+1]+w_mb[l2+1]+w_mb[lc+1]+w_mb[l+1]
    ba += f_mb[l1+l2-lr+1]+f_mb[l+lc+lr+2]+f_mb[l+lc-lr+1]+f_mb[lc+lr-l+1]
    ba += f_mb[lr+l-lc+1]-f_mb[l1+l2+lr+2]-f_mb[l1+lr-l2+1]-f_mb[l2+lr-l1+1]-l*d1l

    ip=lr+n+n1+n2
    p=1+2*(div(ip,2)*2-ip)

    anorm=p*exp(0.5*(bb+ba))
    y=0.0
    j1f=l+1

    for j1=1:j1f
        j2=l+2-j1
        k1i=abs(l1-j1+1)+1
        k1f=l1+j1
        for k1=k1i:2:k1f
            m1f=n1-div(j1+k1-l1,2)+2
            if m1f-1 < 0 ;continue;end
            k2i=max(abs(l2-j2+1),abs(lc-k1+1))+1
            k2f=min(l2+j2,lc+k1)
            if k2i-k2f > 0 ;continue;end
            for k2=k2i:2:k2f
                m2f=n2-div(j2+k2-l2,2)+2
                if m2f-1 < 0 ;continue;end
                ip=j2-1+div(l1+k1+j1+l2+k2+j2,2)
                p=1+2*(div(ip,2)*2-ip)
                bc = 0.5*((k1+j2-2)*dl-(k1+k2-2)*d1l)
                bc += f_mb[k1+l1-j1+1]+f_mb[k1+k2-lc-1]+f_mb[k2+l2-j2+1]-f_mb[k1+l1+j1]-f_mb[k1+k2+lc]
                bc += -f_mb[k2+l2+j2]+w_mb[k1]+w_mb[k2]+f_mb[div(k1+l1+j1,2)]
                bc += f_mb[div(k1+k2+lc,2)]+f_mb[div(k2+l2+j2,2)]              
                bc += -f_mb[div(k1+l1-j1,2)+1]-f_mb[div(l1+j1-k1,2)+1]
                bc += -f_mb[div(j1+k1-l1,2)]-f_mb[div(k1+k2-lc,2)]
                bc += -f_mb[div(k2+lc-k1,2)+1]-f_mb[div(lc+k1-k2,2)+1] 
                bc += -f_mb[div(k2+l2-j2,2)+1]-f_mb[div(l2+j2-k2,2)+1]
                bc += -f_mb[div(j2+k2-l2,2)]
                #println("bc $bc")

                cfac=p*exp(bc)
                sxy=0.0
                ixf=min(k1+k1,k1+k2-lc)-1
                for ix=1:ixf
                    iyi=max(1,ix+j1+l2-k1-lr)
                    iyf=min(l2+l2+1,l1+l2-lr+1,l2+lc+ix-k1-j2+2)
                    if iyi-iyf > 0;continue;end
                    for iy=iyi:iyf
                        ip=ix+iy
                        p=1+2*(div(ip,2)*2-ip)
                        bxy =  f_mb[k1+k1-ix]+f_mb[l2+l2-iy+2]
                        bxy += f_mb[k2+lc-k1+ix]+f_mb[l1+lr-l2+iy] 
                        bxy += -f_mb[ix]-f_mb[iy]-f_mb[k1+k2-lc-ix]
                        bxy += -f_mb[l1+l2-lr-iy+2]-f_mb[k1-l2+lr-j1+iy-ix+1]
                        bxy += -f_mb[l2-k1+lc-j2+ix-iy+3]
                        sxy += p*exp(bxy)
                    end
                end
                s=cfac*sxy
                sm=0.0
                for m1=1:m1f
                    m2i=max(1,nc-m1-div(k1+k2-lc,2)+3)
                    if m2i-m2f > 0;continue;end
                    for m2=m2i:m2f
                        ip=m1+m2
                        p=1+2*(div(ip,2)*2-ip)
                        bm =  (m1-1)*dl-(m1+m2-2)*d1l+g_mb[1]
                        bm += g_mb[m1+m2+div(k1+k2+lc,2)-2]-g_mb[k1+m1-1]
                        bm += -g_mb[k2+m2-1]+f_mb[m1+m2+div(k1+k2-lc,2)-2]
                        bm += -f_mb[m1]-f_mb[m2]-f_mb[n1-m1-div(j1+k1-l1,2)+3]
                        bm += -f_mb[n2-m2-div(j2+k2-l2,2)+3]-f_mb[m1+m2-nc+div(k1+k2-lc,2)-2]
                        sm=sm+p*exp(bm)
                    end
                end
                y += s*sm
            end
        end
    end
    ret = anorm*y
    return ret
end

function def_fgw(;maxjj=200)
    f_mb = zeros(Float64,maxjj)
    g_mb = zeros(Float64,maxjj)
    w_mb = zeros(Float64,maxjj)
    f_mb[1]=0.0
    g_mb[1]=log(0.5)
    w_mb[1]=0.0
    for i=2:maxjj
       a=i-1
       f_mb[i]=f_mb[i-1]+log(a)
       g_mb[i]=g_mb[i-1]+log(a+0.5)
       w_mb[i]=log(a+a+1.0)
    end
    return f_mb,g_mb,w_mb
end

function vtrans(pnrank,izz,ip,Jtot,iza,ia,izb,ib,izc,ic,izd,id,
                   nljsnt,V12ab,X9,U6,Nrmax,
                   nume,numknn,Ndim,Transbk,
                   arr_numst,t2v,t5v,to)
    ret = 0.0
    na,la,jda = nljsnt[ia]; nb,lb,jdb = nljsnt[ib]
    nc,lc,jdc = nljsnt[ic]; nd,ld,jdd = nljsnt[id]
    lrmax = jmax + 1
    Nrmax = 2*emax   
    mab=2*na+la+2*nb+lb; mcd=2*nc+lc+2*nd+ld
    TF = false
    if izz != iza+izb || izz != izc+izd; TF;end
    if Jtot > div(jda+jdb,2) || Jtot < abs(div(jda-jdb,2)); TF=true;end
    if Jtot > div(jdc+jdd,2) || Jtot < abs(div(jdc-jdd,2)); TF=true;end
    if (-1)^(la+lb) != ip || (-1)^(lc+ld) != ip; TF=true;end
    if (izz==2 || izz==-2) && ia==ib && (-1)^(Jtot)==-1; TF=true;end
    if (izz==2 || izz==-2) && ic==id && (-1)^(Jtot)==-1; TF=true;end
    if TF; return ret;end   
    U6_j = U6[Jtot+1]    
    for S=0:1
        tX9 = X9[S+1]
        U6_s = U6_j[S+1]
        tarr_numst = arr_numst[pnrank][S+1]
        lmax1=min(Jtot+S,la+lb)
        lmin1=max(abs(Jtot-S),abs(la-lb))
        if lmin1 > lmax1;continue;end
        lmax2=min(Jtot+S,lc+ld)
        lmin2=max(abs(Jtot-S),abs(lc-ld))
        if lmin2 > lmax2;continue;end
        @inbounds for Lam=lmin1:lmax1
            Ja = jda-2*la; Jb = jdb-2*lb
            t5v[1] = la;t5v[2] = Ja;t5v[3] = lb;
            t5v[4] = Jb;t5v[5] = Lam
            ttX9 = tX9[Jtot+1]            
            x1= get(ttX9,t5v,0.0) * (-1)^Lam
            U6_lam1 = U6_s[Lam+1]
            @inbounds for Lamp=lmin2:lmax2
                Jc = jdc-2*lc; Jd = jdd-2*ld
                t5v[1] = lc;t5v[2] = Jc;t5v[3] = ld
                t5v[4] = Jd;t5v[5] = Lamp
                x2=get(ttX9,t5v,0.0) * (-1)^Lamp
                kncu=div(mab,2)
                U6_lam2 = U6_s[Lamp+1]                
                @inbounds for Ncm=0:kncu
                    klcmax=min((mab-2*Ncm),(mcd-2*Ncm))
                    if klcmax < 0;continue;end
                    @inbounds for Lcm=0:klcmax
                        klrmi1=abs(Lcm-Lam)
                        klrma1=min(lrmax,Lcm+Lam)
                        U6_Lcm1 = U6_lam1[Lcm+1]
                        U6_Lcm2 = U6_lam2[Lcm+1]
                        @inbounds for lr1=klrmi1:klrma1
                            nx1=mab-2*Ncm-(lr1+Lcm)
                            nr1=div(nx1,2)
                            if nr1 < 0 || nx1!=2*nr1;continue;end
                            y1 = 0.0
                            if 2*nr1+lr1+2*Ncm+Lcm <= Nrmax && mab <= Nrmax
                                KK = 2*nr1+lr1 + 2*Ncm+Lcm 
                                s_numknn = numknn[KK+1][Lam+1]
                                ndim = Ndim[KK+1][Lam+1]
                                y1=trbknum(nr1,lr1,Ncm,Lcm,
                                           na,la,nb,lb,Lam,
                                           nume,s_numknn,ndim,Transbk,t2v)
                            end                            
                            if abs(y1) < 1.e-10;continue;end
                            U6_lr1 = U6_Lcm1[lr1+1]
                            klrmi2=abs(Lcm-Lamp); klrma2=Lcm+Lamp
                            @inbounds for lr2=klrmi2:klrma2
                                if lr1%2 != lr2%2;continue;end
                                if lr2 > lrmax;continue;end
                                nx2=mcd-2*Ncm-(lr2+Lcm); nr2=div(nx2,2)
                                if nr2 < 0 || (nx2!=2*nr2);continue;end
                                y2 = 0.0
                                if 2*nr2+lr2+2*Ncm+Lcm <= Nrmax && mcd <= Nrmax
                                    KK = 2*nr2+lr2 + 2*Ncm+Lcm 
                                    s_numknn = numknn[KK+1][Lamp+1]
                                    ndim = Ndim[KK+1][Lamp+1]
                                    y2=trbknum(nr2,lr2,Ncm,Lcm,
                                               nc,lc,nd,ld,Lamp,
                                               nume,s_numknn,ndim,Transbk,t2v)
                                end
                                if abs(y2) < 1.e-10;continue;end
                                mj1=abs(lr1-S); mj2=abs(lr2-S); mj3=abs(Jtot-Lcm)
                                kjmin=max(mj1,mj2,mj3)
                                kjmax=min(lr1+S,lr2+S,Jtot+Lcm,6)
                                if kjmin > kjmax;continue;end
                                U6_lr2 = U6_Lcm2[lr2+1]
                                sumv=0.0
                                @inbounds for Jrel=kjmin:kjmax
                                    zu1=U6_lr1[Jrel-abs(lr1-S)+1]
                                    if abs(zu1) < 1.e-10;continue;end
                                    zu2=U6_lr2[Jrel-abs(lr2-S)+1]
                                    if abs(zu2) < 1.e-10;continue;end
                                    izfac=0
                                    if izz==-2 || izz==2
                                        izfac=1+(-1)^(lr1+S)
                                    end
                                    if izz==0;izfac=1;end
                                    rv12=0.0
                                    if izfac!=0
                                        num= tarr_numst[Jrel+1][lr1-abs(Jrel-S)+1][lr2-abs(Jrel-S)+1]
                                        rv12=V12ab[num][nr1+1,nr2+1]
                                    end
                                    sumv += zu1*zu2*izfac*rv12
                                end
                                zxy=x1*x2*y1*y2
                                ret += sumv*zxy
                            end
                        end
                    end
                end
            end
        end
    end    
    if pnrank!= 2
        Nab = Nfac_jj(na,la,jda,nb,lb,jdb)
        Ncd = Nfac_jj(nc,lc,jdc,nd,ld,jdd)
        ret *= Nab*Ncd
    end
    return ret
end

function prepareX9U6(Nrmax,to)
    jrange = max(Nrmax+1,2*jmax+2)
    X9 = [ [ Dict( [0,0] => 0.0) for J=0:jrange ] for S=0:1]
    for S=0:1;for J=0:jrange; delete!(X9[S+1][J+1], [0,0]);end;end

    jrmax = jmax
    lrmax = jrmax+1
    hit6 = 0; hit9=0
    U6 = [[[[[ zeros(Float64,lr+iss-abs(lr-iss)+1) for lr=0:lrmax ] for lc=0:Nrmax] for lam=0:Nrmax] for iss=0:1] for jj=0:Nrmax+1]
    for lc =0:Nrmax
        for lr =0:lrmax
            for lam=0:Nrmax
                if lam < abs(lc-lr) || lam > lc+lr; continue;end
                for iss=0:1
                    for jj=abs(lam-iss):lam+iss
                        for jr=abs(lr-iss):lr+iss
                            if jr > jrmax;continue;end
                            if jr < abs(lc-jj) || jr > lc+jj;continue;end
                            sfac = sqrt((2.0*lam+1.0)*(2.0*jr+1.0))*(-1.0)^(lc+lr+iss+jj)
                            tmp =  sfac * wigner6j(Float64,lc,lr,lam,iss,jj,jr)
                            U6[jj+1][iss+1][lam+1][lc+1][lr+1][jr-abs(lr-iss)+1] = tmp
                            hit6 += 1
                        end
                    end
                end
            end
        end
    end
    for la=0:Nrmax
        for lb=0:Nrmax
            for nja=-1:2:1
                jda = 2*la + nja 
                if jda < 0;continue;end
                for njb=-1:2:1
                    jdb = 2*lb + njb 
                    if jdb < 0;continue;end
                    for lam=abs(la-lb):la+lb
                        if lam > Nrmax;continue;end
                        for iss=0:1
                            for jj=abs(lam-iss):lam+iss
                                sfac = sqrt((jda+1.0)*(jdb+1.0)*(2.0*lam+1.0)*(2.0*iss+1.0))
                                X9[iss+1][jj+1][
                                     [la,nja,lb,njb,lam]
                                 ] = sfac .* s_wigner9j(la,jda//2,
                                                        lb,jdb//2,
                                                        lam,iss,jj)
                                hit9 += 1
                            end
                        end
                    end
                end
            end
        end
    end
    return X9,U6
end   

function trbknum(Nr,Lr,Nc,Lc,Na,La,Nb,Lb,Lam,nume,s_numknn,ndim,Transbk,t2v)
    ret = 0.0    
    L1=0;L2=0; L3=0;L4=0; N1=0; N2=0;N3=0; N4=0
    phase=0.0
    Kr=2*Nr+Lr
    Kc=2*Nc+Lc
    Ka=2*Na+La
    Kb=2*Nb+Lb
    if (Kr+Kc != Ka+Kb) || abs(Lr-Lc) > Lam || Lr+Lc < Lam || abs(La-Lb) > Lam || La+Lb < Lam
        return ret
    end
    if Kr <= Kc && Ka <= Kb
        N1=Nr; L1=Lr; N2=Nc;L2=Lc; N3=Na; L3=La
        N4=Nb; L4=Lb; phase=1.0
    elseif Kr > Kc && Ka <= Kb
        N1=Nc; L1=Lc; N2=Nr;L2=Lr; N3=Na; L3=La
        N4=Nb; L4=Lb; phase=(-1.0)^(La+Lam)
    elseif Kr <= Kc && Ka > Kb
        N1=Nr; L1=Lr; N2=Nc;L2=Lc; N3=Nb; L3=Lb
        N4=Na; L4=La; phase=(-1.0)^(Lc+Lam)
    elseif Kr > Kc && Ka > Kb
        N1=Nc; L1=Lc; N2=Nr;L2=Lr; N3=Nb; L3=Lb
        N4=Na; L4=La; phase=(-1.0)^(Lr+La)
    end
    K1=2*N1+L1; K3=2*N3+L3; KK=Kr+Kc
    numknn_1 = s_numknn[K1+1][N1+1,N2+1] 
    numknn_2 = s_numknn[K3+1][N3+1,N4+1]
    num=nume[KK+1,Lam+1]+ndim*(numknn_1-1) +numknn_2
    ret = Transbk[num]*phase
    #println("num $num Ndim $ndim n1 $numknn_1 n2 $numknn_2 \n",
    #        "KK $KK Lam $Lam for1 $K1 $N1 $N2 for2 $K3 $N3 $N4 \n")
    return ret
end 

const l2l = [ wigner3j(Float64,l,2,l,0,0,0) for l=0:8]
const l2lnd =[[ wigner3j(Float64,l1,2,l2,0,0,0) for l2=0:8] for l1=0:8]
"""

cg1s = (1,0,l,0|l',0)
cg2s = (2,0,l,0|l',0)

"""
function prep_wsyms()
    lmax = 7
    dim = lmax+3
    cg1s = zeros(Float64,dim,dim)
    cg2s = zeros(Float64,dim,dim)
    d6_121 = zeros(Float64,dim,dim,dim)
    d6_222 = zeros(Float64,dim,dim,dim)
    d6_21 = zeros(Float64,dim,dim,dim,dim)
    d9_12 = zeros(Float64,dim,dim,dim,dim)
    for li = 0:dim-1
        for lj =0:dim-1
            cg1s[li+1,lj+1] = clebschgordan(Float64,1,0,li,0,lj,0)
            cg2s[li+1,lj+1] = clebschgordan(Float64,2,0,li,0,lj,0)
            for j = 0:dim-1
                d6_121[li+1,lj+1,j+1] = wigner6j(Float64,1,2,1,li,lj,j)
                d6_222[li+1,lj+1,j+1] = wigner6j(Float64,2,li,lj,j,2,2)
                for jp = 0:dim-1
                    d6_21[li+1,lj+1,j+1,jp+1] = wigner6j(Float64,li,lj,2,j,jp,1)
                    d9_12[li+1,lj+1,j+1,jp+1] = wigner9j(1,li,lj,1,j,jp,2,2,2)
                end
            end
        end
    end
    return cg1s,cg2s,d6_121,d6_21,d6_222,d9_12
end

""" 
slightly different from readsnt() in ShellModel.jl
"""
function readsnt(sntf,Anum) 
    f = open(sntf,"r");tlines = readlines(f);close(f)
    lines = rm_comment(tlines)
    line = lines[1]
    lp,ln,cp,cn = map(x->parse(Int,x),rm_nan(split(line," ")))
    p_sps = [[0,0]];deleteat!(p_sps,1)
    n_sps = [[0,0]];deleteat!(n_sps,1)
    dictsps = Dict([0,0,0,0]=>0);delete!(dictsps,[0,0,0,0])
    nls = []
    nlhit=0
    @inbounds for i = 1:lp
        ith,n,l,j,tz = map(x->parse(Int,x),rm_nan(split(lines[1+i]," "))[1:5])
        push!(p_sps,[n,l,j,tz])
        if ([n,l,tz] in nls)==false
            nlhit +=1
            push!(nls,[n,l,tz])
        end
        dictsps[[n,l,j,tz]] = i

    end
    @inbounds for i = 1:ln
        ith, n,l,j,tz = map(x->parse(Int,x),rm_nan(split(lines[1+i+ln]," "))[1:5])
        if ([n,l,tz] in nls)==false
            nlhit +=1
            push!(nls,[n,l,tz])
        end
        push!(n_sps,[n,l,j,tz])
        dictsps[[n,l,j,tz]] = i + lp
    end
    sps = vcat(p_sps,n_sps)    
    nsp,zero = map(x->parse(Int,x),rm_nan(split(lines[1+ln+lp+1]," "))[1:2])
    SPEs = [ [0.0 for i=1:lp],[0.0 for i=1:ln]]
    @inbounds for i = 1:nsp
        idx=0; j=i
        if i<=lp;idx=1;else;idx=2;j-=lp;end
        SPEs[idx][j] =parse(Float64,rm_nan(split(lines[1+ln+lp+1+i]," "))[3])
    end
    ntbme,massop,Aref,p = rm_nan(split(lines[1+ln+lp+1+nsp+1]," "))
    ntbme = parse(Int,ntbme); massop=parse(Int,massop)

    dictTBMEs=[ Dict([0,0,0,0,0]=>0.0) for pnrank=1:3]
    for i=1:3
        delete!(dictTBMEs[i],[0,0,0,0,0])
    end
    @inbounds for ith = 1:ntbme
        i,j,k,l,totJ,TBME= rm_nan(split(lines[1+ln+lp+1+nsp+1+ith], " "))
        i = parse(Int,i);j = parse(Int,j);k = parse(Int,k);l = parse(Int,l);
        totJ = parse(Int,totJ)
        nth = 0
        if i<=lp && j<=lp
            nth = 1
        elseif i>lp && j > lp
            nth = 3
        elseif i<=lp && j>lp
            nth = 2
        else
            println("err");exit()
        end
        TBME = parse(Float64,TBME)
        ## snt file must be "ordered"; a<=b & c=d & a<=c
        dictTBMEs[nth][[i,j,k,l,totJ]] = TBME        
    end
    return sps,dictsps,dictTBMEs
end

function jj_std(sps,dictsps,dictTBMEs;fname="")
    if fname=="";fname = "monopole"*fname*".dat"
    else;fname = "monopole_"*fname*".dat";end
    ln = length(sps)
    nmax = maximum([sps[i][1] for i=1:ln])
    nmin = minimum([sps[i][1] for i=1:ln])
    lmax = maximum([sps[i][2] for i=1:ln])
    lmin = minimum([sps[i][2] for i=1:ln])
    jmax = maximum([sps[i][3] for i=1:ln])
    jmin = minimum([sps[i][3] for i=1:ln])
    emax = maximum( [ 2*sps[i][1]+sps[i][2] for i=1:length(sps)])
    emin = minimum( [ 2*sps[i][1]+sps[i][2] for i=1:length(sps)])
    Jmax = div(2*jmax+1,2)
    Lmin = 0; Lmax = 2*lmax
    monodict = Dict(["",""]=> [[0.0,[0.0]]]);delete!(monodict,["",""])
    cpnrank=["pp","pn","nn"]

    keya = [0,0,0,0]; keyb = [0,0,0,0];keyc=[0,0,0,0]; keyd=[0,0,0,0]   # nljtz
    key_abcdJ = [0,0,0,0,0]
    tvjj = zeros(Float64,5)
    for pnrank =1:3
        tdict = dictTBMEs[pnrank]
        for tkey in keys(tdict)
            a,b,c,d,J = tkey
            na,la,ja2,iza = sps[a]; nb,lb,jb2,izb = sps[b]
            nc,lc,jc2,izc = sps[c]; nd,ld,jd2,izd = sps[d]
            keya[1]=na; keya[2]=la; keya[3]=ja2; keya[4]=iza
            keyb[1]=nb; keyb[2]=lb; keyb[3]=jb2; keyb[4]=izb
            keyc[1]=nc; keyc[2]=lc; keyc[3]=jc2; keyc[4]=izc
            keyd[1]=nd; keyd[2]=ld; keyd[3]=jd2; keyd[4]=izd
            v_target = tdict[tkey]            
            tvjj .= 0.0
            hats = (-1)^J * sqrt((ja2+1.0)*(jb2+1.0)*(jc2+1.0)*(jd2+1.0))
            for L = abs(la-lb):la+lb
                for S=0:1
                    for Lp = abs(lc-ld):lc+ld
                        for Sp = 0:1
                            Jfac  = wigner9j(la,1//2,ja2//2,lb,1//2,jb2//2,L,S,J)
                            Jfac *= wigner9j(lc,1//2,jc2//2,ld,1//2,jd2//2,Lp,Sp,J)
                            for jpa2 = 2*la-1:2:2*la+1
                                if jpa2 <=0;continue;end
                                for jpb2 = 2*lb-1:2:2*lb+1
                                    if jpb2 <=0;continue;end
                                    for jpc2 = 2*lc-1:2:2*lc+1
                                        if jpc2 <=0;continue;end
                                        for jpd2 = 2*ld-1:2:2*ld+1
                                            if jpd2 <=0;continue;end
                                            hatps = sqrt((jpa2+1.0)*(jpb2+1.0)*(jpc2+1.0)*(jpd2+1.0))
                                            LLSS = (2*L+1)*(2*Lp+1)*(2*S+1)*(2*Sp+1)
                                            for Jp = max(abs(Lp-Sp),abs(L-S)):min(Lp+Sp,L+S)
                                                d9j_ab = wigner9j(la,1//2,jpa2//2,lb,1//2,jpb2//2,L,S,Jp)
                                                d9j_cd = wigner9j(lc,1//2,jpc2//2,ld,1//2,jpd2//2,Lp,Sp,Jp)
                                                keya[3] = jpa2; keyb[3]=jpb2;keyc[3] = jpc2; keyd[3]=jpd2                                                
                                                ap = dictsps[keya];bp = dictsps[keyb];cp = dictsps[keyc];dp = dictsps[keyd]
                                                key_abcdJ[1] = ap;key_abcdJ[2] = bp
                                                key_abcdJ[3] = cp;key_abcdJ[4] = dp;key_abcdJ[5] = Jp
                                                ta=ap; tb=bp;tc=cp;td=dp
                                                phase = 1.0
                                                if ap > bp
                                                    ta = bp; tb = ap
                                                    phase *= (-1)^(div(sps[ap][3]+sps[bp][3],2)+Jp+1)
                                                end
                                                if cp > dp
                                                    tc = dp; td = cp
                                                    phase *= (-1)^(div(sps[cp][3]+sps[dp][3],2)+Jp+1)
                                                end
                                                fa=ta; fb=tb;fc=tc;fd=td
                                                if fa > fc || (fa==fc && fb > fd)
                                                    fa=tc;fb=td;fc=ta;fd=tb
                                                end
                                                key_abcdJ[1] = fa;key_abcdJ[2] = fb
                                                key_abcdJ[3] = fc;key_abcdJ[4] = fd;key_abcdJ[5] = Jp
                                                vjj = get(tdict,key_abcdJ,false)                                                
                                                if vjj==false;continue; end
                                                #println("abcd' $ap $bp $cp $dp J' $Jp $nvjj $vjj")
                                                Jpfac = (-1)^Jp * (2*Jp+1) * d9j_ab * d9j_cd
                                                for k = 0:2 
                                                    tidx = k+1; if S==Sp && k==1;tidx=4;end                                                
                                                    kfac = (2*k+1) * wigner6j(Float64,L,S,J,Sp,Lp,k) * wigner6j(Float64,L,S,Jp,Sp,Lp,k)
                                                    tvjj[tidx] += hats*hatps*kfac*Jfac*Jpfac * LLSS * vjj *phase
                                                end
                                            end
                                        end
                                    end
                                end                                
                            end
                        end
                    end
                end
            end
            tvjj[5] = sum( @views tvjj[1:4])
            if abs(v_target-tvjj[5]) > 1.e-9
                println("Error TBME(in jj) mismatch!: $tkey  v ",
                        @sprintf("%12.4e",v_target),"; sum ", @sprintf("%12.4e",tvjj[5]))
            end
            if a==c && b==d && ( a <= b )
                c = cpnrank[pnrank]
                ca = string(na)*chara_l[la+1]*string(ja2)
                cb = string(nb)*chara_l[lb+1]*string(jb2)
                if get(monodict,[c,ca,cb],false) == false
                    monodict[[c,ca,cb]] = [ [J*1.0,copy(tvjj)] ]
                else
                    push!(monodict[[c,ca,cb]],[J*1.0,copy(tvjj)])
                end
            end
        end
        monopole(monodict,fname)
    end
    return nothing
end


function monopole(monodict,fname)
    io = open(fname,"w")
    nume = zeros(Float64,5)
    deno = 0.0
    for key in keys(monodict)
        tarr = monodict[key]
        nume .= 0.0
        deno = 0.0
        for tmp in tarr
            J, V =tmp
            nume .+= (2*J+1.0) .* V
            deno += 2*J +1.0
        end
        vm = nume ./ deno
        tx = key[1]*":"*key[2]*"-"*key[3]*@sprintf("%20.12f", vm[1])*@sprintf("%20.12f", vm[2])*@sprintf("%20.12f", vm[3])*@sprintf("%20.12f", vm[4])*@sprintf("%20.12f", vm[5])
        println(io,tx)
    end
    close(io)
    return nothing
end

# function sample_std()
#     sps,dictsps,dictTBMEs = readsnt("usdb.snt",18)
#     spin_tensor_decomposition(sps,dictsps,dictTBMEs)
# end

function nlj_from_str(str)
    lstr = match(r"[a-z]",str).match
    n,j = split(str,lstr)
    l = -1
    for k=1:length(chara_L)
        if lstr==chara_l[k];l=k;break;end
    end    
    l = l-1
    n = parse(Int64,n)
    j = parse(Int64,j)
    return n,l,j   
end
