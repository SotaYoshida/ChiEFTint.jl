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

# for specific cases
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

const chara_L = ["S","P","D","G","H","I","J","K","L","M","N"]
chara_SLJ(S,L,J) = "{}^"*ifelse(S==0,"1","3")*chara_L[L+1]*"_{"*string(J)*"}"
delta(a,b) = ifelse(a==b,1.0,0.0)
hat(a) = sqrt(2.0*a+1.0)

function genLaguerre(n,alpha,x)
    if n==0
        return 1
    elseif n==1
        return -x + (alpha+1)
    else
        s = 0.0
        for i=0:n
            bfac = gamma(n+alpha+1) / ( gamma(n-i+1) * gamma(alpha+i+1))
            s += (-1)^i * x^i / factorial(i) * bfac
        end        
        return s
    end
end

function Gauss_Legendre(xmin,xmax,n;eps=3.e-16) 
    m = div(n+1,2)
    x = zeros(Float64,n); w = zeros(Float64,n)
    xm = 0.5 * (xmin+xmax); xl = 0.5 * (xmax-xmin)
    p1 = 0.0;p2=0.0;p3 =0.0; pp = 0.0; z = 0.0; z1=0.0
    for i =1:m
        z1 = 0.0
        z = cos((pi*(i-0.25))/(n+0.5))
        hit = 0        
        while true
            hit += 1
            p1 = 1.0; p2 = 0.0
            for j=1:n
                p3=p2
                p2=p1
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j
            end
            pp = n*(z*p1-p2)/(z^2 -1.0)
            z1 = z
            z -= p1/pp
            if abs(z-z1) < eps; break;end
            if hit > 100; println("warn! in Gausse_Legendre");exit();end
        end
        x[i]=xm-xl*z
        x[n+1-i]=xm+xl*z
        w[i]=2.0*xl/( (1.0-z^2) * pp^2)
        w[n+1-i]=w[i]
    end
    return x,w
end

function make_chiEFTint(;is_plot=false)
    to = TimerOutput()
    @timeit to "prep." begin
        # prep. momentum mesh
        xr_fm,wr = Gauss_Legendre(0.0,pmax_fm,n_mesh)
        xr = xr_fm .* hc
        numst2,dict_numst,arr_numst = bstate()
        V12mom = [ zeros(Float64,n_mesh,n_mesh) for i=1:length(numst2)]
        ## prep. radial functions
        rmass = Mp*Mn/(Mp+Mn)
        br = sqrt(hc^2 /(rmass*hw))
        Rnl = Rnl_all_ab(lmax,br,n_mesh,xr_fm)

        ## prep. for valence space oparators
        xrP=[0.0]; wrP=[0.0]
        xrP_fm,wrP = Gauss_Legendre(0.0,Pmax_fm,n_mesh_P)
        xrP = xrP_fm .* hc   
        RNL = Rnl_all_ab(lcmax,br,n_mesh_P,xrP_fm)
        
        ## prep. for partial-wave decompositon
        lsjs = [ [ [J,J,0,J],[J,J,1,J],[J+1,J+1,1,J],
                   [J-1,J-1,1,J],[J+1,J-1,1,J],[J-1,J+1,1,J]]
                 for J = 0:jmax]
        llpSJ_s = [ [0,0,1,1],[1,1,1,0],[1,1,0,1],[1,1,1,1],
                    [0,0,0,0],[0,2,1,1],[1,1,1,2]]
        tllsj = [0,0,0,0,0]
        opfs = [ zeros(Float64,11) for i=1:5]#T,SS,C,LS,SL
        f_ss!(opfs[2]);f_c!(opfs[3])
        ## prep. Gauss point for integrals
        ts, ws = Gauss_Legendre(-1,1,20)
    end

    tpnrank = 2 # 1:pp 2:pn 3:nn pnrank for momplot 
    @timeit to "calcV" begin
        # ***Leading Order (LO)***
        ### OPEP
        @timeit to "OPEP" begin
            OPEP(ts,ws,xr,V12mom,dict_numst,to,lsjs,llpSJ_s,tllsj,opfs)
        end        
        ### contact term
        LECs_LO_contact = [C0_1S0,C0_3S1]
        LO(xr,LECs_LO_contact,V12mom,dict_numst,to)


        if chi_order >= 1
            # ***Next-to-Leading Order(NLO)***
            LECs_NLO_contact = [C2_3S1,C2_3P0,C2_1P1,C2_3P1,
                                C2_1S0,C2_3SD1,C2_3SD1,C2_3P2]
            ### contact term Q^2
            NLO(xr,LECs_NLO_contact,V12mom,dict_numst,to)

            ### TPE terms (NLO,NNLO,N3LO)
            @timeit to "TPE" begin
                tpe(ts,ws,xr,V12mom,dict_numst,to,llpSJ_s,lsjs,tllsj,opfs)
            end
            # *** N3LO ***
            if chi_order >= 3
                ### contact Q^4
                LECs_N3LO_contact = [[hD_1S0,D_1S0],[D_3P0],[D_1P1],[D_3P1],
                                     [hD_3S1,D_3S1],[D_3D1],
                                     [hD_3SD1,D_3SD1],[hD_3SD1,D_3SD1],
                                     [D_1D2],[D_3D2],[D_3P2],[D_3PF2],
                                     [D_3PF2],[D_3D3]]                 
                N3LO(xr,LECs_N3LO_contact,V12mom,dict_numst,to)
            end
        end
    end
    #write_vmom(xr,V12mom,dict_numst[tpnrank],tpnrank,llpSJ_s;label=tl)
    if is_plot; momplot(xr,V12mom,dict_numst[tpnrank],tpnrank,llpSJ_s);end
    
    #transform mom. int. to HO matrix element
    @timeit to "Vrel" begin
        V12ab = Vrel(V12mom,numst2,xr_fm,wr,n_mesh,Rnl,to)
    end
    @timeit to "HOtrans" begin
        jab_max = 4 * emax + 2
        Numpn= Dict( [0,0,0,0] => 0 ) ;delete!(Numpn,[0,0,0,0])   
        infos,izs_ab,nTBME = make_sp_state(jab_max,Numpn)
        println("# of channels 2bstate ",length(infos)," #TBME = $nTBME")
        TMtrans(xr,wr,xrP,wrP,Rnl,RNL,
                nTBME,infos,izs_ab,
                Numpn,V12ab,arr_numst,to)
    end

    show(to, allocations = true,compact = false);println("")
    
    return nothing    
end

#############
##basis_states.jl
########

#In this file, most functions are naive translation from the Fortran codes,
#which are used in the references below, and its modification by Takayuki Miyagi (@TRIUMF)
#Refs:
#D.R. Entem and R. Machleidt, PRC 68, 041001 (2003)
#R. Machleidt and D.R. Entem, Physics Reports 503, 1 (2011).
#

function make_sp_state(jab_max,Numpn)
    kh = Dict( [0,0] => 0 ) ;delete!(kh,[0,0])
    kn = Dict( [0,0] => 0 ) ;delete!(kn,[0,0])
    kl = Dict( [0,0] => 0 ) ;delete!(kl,[0,0])
    kj = Dict( [0,0] => 0 ) ;delete!(kj,[0,0])    
    maxsps = Int((emax+1) * (emax+2) / 2)
    println("# of sp states $maxsps")    
    n = 0
    for NL=0:emax
        for L=0:NL
            if (NL-L) % 2 == 1;continue;end
            Nn= div(NL-L,2)
            for IS=-1:2:1
                jd=2*L+IS
                if jd < 0; continue;end
                n=n+1
                kh[[-1,n]]=3;  kh[[1,n]]=3
                kn[[-1,n]]=Nn; kn[[1,n]]=Nn
                kl[[-1,n]]=L;  kl[[1,n]]=L
                kj[[-1,n]]=jd; kj[[1,n]]=jd
                Numpn[[-1,Nn,L,jd]]=n
                Numpn[[ 1,Nn,L,jd]]=n
            end
        end
    end 
    infos,izs_ab,nTBME = get_twq_2b(kh,kn,kl,kj,maxsps,jab_max)
    return infos,izs_ab,nTBME
end

function get_twq_2b(kh,kn,kl,kj,maxsps,jab_max)
    infos = [ [0,0,0,0] ];deleteat!(infos,1)
    izs_ab = [ [[0,0,0,0]] ];deleteat!(izs_ab,1)
    ichan = 0
    nTBME = 0
    for izz = -2:2:2
        for ipp = -1:2:1
            for jjx = 0:2*emax+1
                if jjx  > div(jab_max,2);continue;end
                ndim = 0
                tl = [ [0,0,0,0]];deleteat!(tl,1)
                for iza = -1:2:1
                    for izb = iza:2:1
                        if iza + izb != izz;continue;end
                        for ia = 1:maxsps
                            la = kl[[iza, ia]]; ja = kj[[iza, ia]]
                            ibmin = 1
                            if iza == izb;ibmin = ia; end                            
                            for ib = ibmin:maxsps
                                lb = kl[[izb, ib]]; jb = kj[[izb, ib]]                                
                                if (-1)^(la+lb) != ipp;continue;end
                                if tri_check(ja, jb, 2*jjx) == false; continue;end
                                if iza==izb && ia==ib && jjx % 2 == 1;continue;end
                                ndim = ndim + 1
                                push!(tl,[iza,ia,izb,ib])
                            end                            
                        end
                    end
                end
                if ndim != 0
                    ichan += 1
                    #println("ichan $ichan ndim $ndim")
                    push!(infos,[izz,ipp,jjx,ndim])
                    nTBME += Int(ndim*(ndim+1)/2)
                    push!(izs_ab,tl)
                end
            end
        end
    end
    return infos,izs_ab,nTBME
end
#infos => iz12,ipty,jtb,ndim
function TMtrans(xr,wr,xrP,wrP,Rnl,RNL,
                 nTBME,infos,izs_ab,
                 Numpn,V12ab,arr_numst,to
                 ;CM=false)
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
    if tbme_fmt == "snt"
        #2n +l = emax
        for temax = 0:emax
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
    end
    
    @timeit to "9j6j" X9,U6 = prepareX9U6(Nrmax,to)
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
    io=open(fn_tbme,"w")
    if tbme_fmt  == "myg"
        println(io,"tza,  a, tzb,  b, tzc,  c, tzd,  d, J12,   <ab|v|cd>,   <ab|p_ip_j|cd>/mhw,    <ab|r_ir_j|cd>mw/h")
    else
        nothing
    end
    @inbounds for ich=1:length(infos)
        izz,ip,Jtot,ndim=infos[ich]
        pnrank = Int(div(izz,2))+2
        izs = izs_ab[ich]        
        vv = zeros(Float64,ndim,ndim)
        p_pro = zeros(Float64,ndim,ndim)
        r_pro = zeros(Float64,ndim,ndim)
        #@timeit to "vtrans" @inbounds @qthreads for i = 1:ndim
        @timeit to "vtrans" @inbounds for i = 1:ndim
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
                    NLOvs(vvs,xr,wr,xrP,wrP,Rnl,RNL,
                          cg1s,sp_P5_9j,nljsnt,pnrank,ip,
                          X9,U6,t5v,
                          Jtot,iza,ia,izb,ib,
                          izc,ic,izd,id,to)
                    for i=1:5
                        v12 += vvs[i]
                    end
                end
                vv[i, j] = v12;vv[j, i] = v12                
                if CM # needed for CM motion
                    #call X12cmmat(izz, ip, Jtot, iza, ia, izb, ib, izc, ic, izd, id, rcm, pcm)
                    #call X12relmat(izz, ip, Jtot, iza, ia, izb, ib, izc, ic, izd, id, rrel, prel)
                    #p_pro[i, j] = pcm - prel; p_pro[j, i] = pcm - prel
                    #r_pro[i, j] = rcm - rrel; r_pro[j, i] = rcm - rrel
                end   
            end
        end
        write_tbme(io,ndim,izs,Jtot,vv,nljsnt;ofst=nofst)
    end
    close(io)
    ##END: Trans_TM in nn_int.f90
    return nothing
end

function write_tbme(io,ndim,izs,Jtot,vv,nljsnt;ofst=0)
    if tbme_fmt == "myg"
        @inbounds for i = 1:ndim
            iza,ia,izb,ib = izs[i]
            for j = 1:i
                izc,ic,izd,id= izs[j]
                print(io,@sprintf("%4i", iza))
                print(io,@sprintf("%4i", ia))
                print(io,@sprintf("%4i", izb))
                print(io,@sprintf("%4i", ib))
                print(io,@sprintf("%4i", izc))
                print(io,@sprintf("%4i", ic))
                print(io,@sprintf("%4i", izd))
                print(io,@sprintf("%4i", id))
                print(io,@sprintf("%4i", Jtot))
                println(io,@sprintf("%20.10e", vv[i,j]))
            end
        end
    elseif tbme_fmt == "snt"
        @inbounds for i = 1:ndim
            iza,ia,izb,ib = izs[i]
            for j = 1:i
                izc,ic,izd,id= izs[j]
                a = ifelse(iza==-1,ia,ia+ofst)
                b = ifelse(izb==-1,ib,ib+ofst)
                c = ifelse(izc==-1,ic,ic+ofst)
                d = ifelse(izd==-1,id,id+ofst)
                ra = a; rb = b; rc=c; rd=d
                fac = 1.0
                if a > b
                    ra = b; rb = a
                    fac *= (-1.0)^(div(nljsnt[a][3]+nljsnt[b][3],2)+Jtot+1)
                end
                if c > d 
                    ra = d; rd = c
                    fac *= (-1.0)^(div(nljsnt[c][3]+nljsnt[d][3],2)+Jtot+1)
                end
                fa = ra; fb =rb; fc=rc; fd=rd
                if ra > rc 
                    fa=rc;fb=rd;fc=ra;fd=rb
                else
                    if rb > rd
                        fa=rc;fb=rd;fc=ra;fd=rb
                    end
                end
                
                print(io,@sprintf("%5i", fa))
                print(io,@sprintf("%5i", fb))
                print(io,@sprintf("%5i", fc))
                print(io,@sprintf("%5i", fd))
                print(io,@sprintf("%6i", Jtot))
                println(io,@sprintf("%18.10f", vv[i,j]))
            end
        end
        
    end
    return nothing
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

function print_vec(s,v;ine=false)
    s *= " "
    for i = 1:length(v)
        if ine
            s *= @sprintf "%9.1e" v[i]
        else
            #s *= @sprintf "%10.4f" v[i]
            s *= @sprintf "%12.4e" v[i]
            #s *= @sprintf "%25.15f" v[i] 
	end
    end
    println(s)
end

###########
### potentials.f
###########

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

function freg(p,pp,n)
    Lambchi
    if n==0;return 1.0;end
    if n==-1 # for 3D3 pw
        return 0.5 *( exp( - (p/Lambchi)^4 - (pp/Lambchi)^4 )
                      + exp( - (p/Lambchi)^6 - (pp/Lambchi)^6 ))
    end
    return exp( - (p/Lambchi)^(2*n) - (pp/Lambchi)^(2*n) )
end

### function-forms for pw
f1(x,y,c,cdum) = c
fx(x,y,c,cdum) = c * x
fy(x,y,c,cdum) = c * y
fxxpyy(x,y,c,cdum) = c * (x^2 + y^2)
fxy(x,y,c,cdum) = c *x*y
fx2(x,y,c,cdum) = c * x^2
fy2(x,y,c,cdum) = c * y^2
f_442(x,y,c1,c2) = c1 * (x^4 + y^4) + c2 *(x^2 * y^2)
f_x42(x,y,c1,c2) = c1 * x^4  + c2 *(x^2 * y^2)
f_y42(x,y,c1,c2) = c1 * y^4  + c2 *(x^2 * y^2)
f_x2y2(x,y,c,cdum) = c * x^2 * y^2
f_x3y(x,y,c,cdum) = c * x^3 * y
f_xy3(x,y,c,cdum) = c * x * y^3
f_31(x,y,c,cdum) = c * (x^3 * y + x * y^3)

fp_P2(p,ell,pp,ellp,P) = P^2
fp_ddP(p,ell,pp,ellp,P) = P * (delta(ell,1)*delta(ellp,0)*p
                                 + delta(ell,0)*delta(ellp,1)*pp)
##

# ** Leading Order (LO)** Q^0 contact term
function LO(xr,LECs,V12mom,dict_numst,to)
    C_1S0,C_3S1 = LECs
    fac_q = hc^3 * 1.e-2 / (2*pi)^3
    pfunc = f1
    n_reg = 3    
    ress = [2*C_CIB+C_CSB,0.0,2*C_CIB-C_CSB]
    for pnrank = 1:3
        TF = true
        res = ress[pnrank]
        ##1S0 (LO)
        l=0;lp=0;S=0;J=0
        LEC = (C_1S0+res) * fac_q
        calc_Vmom!(pnrank,V12mom,dict_numst[pnrank],
                   xr,LEC,LEC,l,lp,S,J,pfunc,n_reg,to)
        
        ## 3S1 (LO)
        l=0;lp=0;S=1;J=1
        if pnrank%2==1;continue;end
        LEC = (C_3S1) * fac_q
        calc_Vmom!(pnrank,V12mom,dict_numst[pnrank],
                   xr,LEC,LEC,l,lp,S,J,pfunc,n_reg,to)
    end
    return nothing
end

# ** Next to Leading Order (NLO) ** Q^2 contact term
# LECs=>C2_3S1,C2_3P0,C2_1P1,C2_3P1,C2_1S0,C2_3SD1,C2_3DS1,C2_3P2
function NLO(xr,LECs,V12mom,dict_numst,to;
             n_reg=2)
    fac = hc^3  *  1.e-8  / (2*pi)^3  * 2.0 ###factor x2.0 needed to reproduce NLOonly

    funcs = [fxxpyy,fxy,fxy,fxy,fxxpyy,fx2,fy2,fxy] #fx2
    llpSJ_s = [ [0,0,1,1],[1,1,1,0],[1,1,0,1],[1,1,1,1],
                [0,0,0,0],[0,2,1,1],[2,0,1,1],[1,1,1,2]] #
    for pnrank = 1:3
        @inbounds for n=1:8            
            LEC = LECs[n] * fac
            pfunc = funcs[n]
            l,lp,S,J = llpSJ_s[n]
            if pnrank%2 == 1 && (l+S+1) % 2 != 1;continue;end
            calc_Vmom!(pnrank,V12mom,dict_numst[pnrank],xr,LEC,LEC,
                       l,lp,S,J,pfunc,n_reg,to)
        end
    end
    return nothing
end


function N3LO(xr,LECs_N3LO,V12mom,dict_numst,to)
    fac = hc^3  *  1.e-14  / (2*pi)^3  
    funcs = [ f_442,f_31,f_31,f_31,
              f_442,f_x2y2,
              f_x42,f_y42,
              f_x2y2,f_x2y2,f_31,f_x3y,f_xy3,f_x2y2]
    llpSJ_s = [ [0,0,0,0],[1,1,1,0],[1,1,0,1],[1,1,1,1],
                [0,0,1,1],[2,2,1,1],
                [0,2,1,1],[2,0,1,1],
                [2,2,0,2],[2,2,1,2],[1,1,1,2],[1,3,1,2],[3,1,1,2],[2,2,1,3]]
    nregs = [2,3,2,4,2,2,2,2,4,2,2,4,4,-1]
    for pnrank = 1:3
        for (n,LECs) in enumerate(LECs_N3LO)
            l,lp,S,J = llpSJ_s[n]            
            if pnrank%2 == 1 && (l+S+1) % 2 != 1;continue;end
            LEC = 0.0; LEC2 = 0.0
            if length(LECs)== 1;
                LEC =LEC2 = LECs[1]
            else
                LEC = LECs[1]; LEC2 = LECs[2]
            end                
            pfunc = funcs[n]
            n_reg = nregs[n]
            LEC *= fac; LEC2 *= fac
            calc_Vmom!(pnrank,V12mom,dict_numst[pnrank],xr,LEC,LEC2, 
                       l,lp,S,J,pfunc,n_reg,to)
        end
    end
    return nothing
end

function OPEP(ts,ws,xr,V12mom,dict_numst,to,lsjs,llpSJ_s,tllsj,opfs)
    hc3 = hc^3
    tVs = zeros(Float64,6)
    tVs_ch = zeros(Float64,6)
    v1d = zeros(Float64,6)
    opfs = zeros(Float64,8) # sq
    mpi0 = mpis[2]; mpi02 = mpi0^2
    mpipm = mpis[1]; mpipm2 = mpipm^2

    for pnrank =1:3
        tdict = dict_numst[pnrank]
        MN = Ms[pnrank];dwn = 1.0/MN;sq_dwn=dwn^2
        fff = pi / ((2*pi)^3 * MN^2)
        coeff = -(MN*gA/(2*Fpi))^2
        itt = itts[pnrank]
        tllsj[1] = itt
        CIB=false;hit = 0
        TF = true
        while TF 
            @inbounds for J=0:jmax
                lsj = lsjs[J+1]
                f_idx = 6
                if J==0; f_idx = 3;end                
                @inbounds for i= 1:n_mesh
                    x = xr[i];xdwn = x * dwn;sq_xdwn= xdwn^2
                    ex = sqrt(1.0+sq_xdwn)
                    @inbounds for j = 1:n_mesh
                        y = xr[j]; ydwn = y*dwn;sq_ydwn= ydwn^2
                        ey = sqrt(1.0+sq_ydwn)
                        nfac = 1.0/(x* y* sq_dwn)
                        ree = 1.0/sqrt(ex*ey) * freg(x,y,4) 
                        f_sq!(opfs,xdwn,ydwn)
                        if pnrank != 2
                            cib_lsj_opep(opfs,x,y,mpi02,1,J,pnrank,nfac,ts,ws,tVs)
                        else
                            cib_lsj_opep(opfs,x,y,mpi02,1,J,pnrank,nfac,ts,ws,tVs)
                            cib_lsj_opep(opfs,x,y,mpipm2,2,J,pnrank,nfac,ts,ws,tVs)
                        end
                        fc = fff * coeff * ree
                        t_fc= fc * hc3 
                        @inbounds for idx = 1:f_idx
                            @views tllsj[2:5] .= lsj[idx]
                            tl,tlp,tS,tJ = lsj[idx] 
                            if pnrank%2 == 1 && (tl+tS+1)%2 != 1;continue;end
                            V12idx = get(tdict,tllsj,-1)
                            if V12idx == -1;continue;end
                            tfac = tVs[idx] * t_fc
                            V = V12mom[V12idx]                            
                            V[i,j] += tfac
                        end
                    end
                end
            end
            break
            # if pnrank%2==1;TF=false;end
            # if CIB;TF=false;end
            # CIB = true
            # hit += 1           
            # mpi = mpis[1+hit] # for static (LO) OPEP
            #else #pnrank=2 (pn)
            #     if pigammacorrection #mpi+
            #         break
            #     else # gA=-1.29 is done up to here
            #         cgA = -0.06217
            #         coeff = -(Mm*cgA/(2.0*Fpi))^2
            #         pigammacorrection = true
            #     end
            # end
        end
    end
    return nothing
end

function set_pjs!(J,pjs,ts)
    for i=1:length(pjs);pjs[i] .= 0.0;end
    if J ==0
        pjs[1] .= 1.0; pjs[3] .= 0.0
    else
        pjs[3] .= 1.0
        for (i,t) in enumerate(ts)        
            pjs[1][i] = t
        end
        if J>1
            for (i,t) in enumerate(ts)
                pj = pjs[1][i]
                pjm1 = pjs[3][i]            
                for tJ = 2:J
                    a = t * pj
                    b = a-pjm1
                    pjm1 = pj
                    pj = -b/tJ + b+a
                end
                pjs[1][i] = pj
                pjs[3][i] = pjm1
            end
        end
    end
    for (i,t) in enumerate(ts)
        pjs[2][i] = pjs[1][i] * t
        pjs[4][i] = pjs[2][i] * t
        pjs[6][i] = pjs[4][i] * t
        pjs[5][i] = pjs[3][i] * t
        pjs[7][i] = pjs[5][i] * t
    end
    return nothing
end

function tpe(ts,ws,xr,V12mom,dict_numst,to,llpSJ_s,lsjs,tllsj,opfs)
    hc3 = hc^3
    tVs = zeros(Float64,6)
    v1d = zeros(Float64,6)
    mmpi = sum(mpis)/3.0
    mpi0 = mpis[2]; mpi02 = mpi0^2

    gis = [ zeros(Float64,7) for i=1:9]#Vt/Wt/Vs/Ws/Vc/Wc/Vls/Wls/Vsl
    pjs = [zeros(Float64,length(ts)) for i=1:7]
    
    for pnrank =1:3
        tdict = dict_numst[pnrank]
        MN = Ms[pnrank];dwn = 1.0/MN
        fff = pi / ((2*pi)^3 * MN^2)
        
        nd_mpi = mmpi/MN
        nd_mpi2 = nd_mpi^2
        nd_mpi4 = nd_mpi2^2
        nd_mpi6 = nd_mpi2^3
        nd_mpi8 = nd_mpi2^4
        Fpi2 = (Fpi/MN)^2
        Fpi4 = (Fpi/MN)^4
        Fpi6 = (Fpi/MN)^6
        c1 = c1_NNLO * MN * 1.e-3
        c2 = c2_NNLO * MN * 1.e-3
        c3 = c3_NNLO * MN * 1.e-3
        c4 = c4_NNLO * MN * 1.e-3
        r_d12 = d12 * MN^2 * 1.e-6
        r_d3 = d3 * MN^2 * 1.e-6
        r_d5 = d5 * MN^2 * 1.e-6
        r_d145 = d145 * MN^2 * 1.e-6        
        itt = itts[pnrank]
        tllsj[1] = itt
        @inbounds for J=0:jmax
            lsj = lsjs[J+1]
            f_idx = 6
            if J==0; f_idx = 3;end
            set_pjs!(J,pjs,ts)            
            @inbounds for i= 1:n_mesh
                x = xr[i];xdwn = x*dwn
                xdwn2 = xdwn^2       
                ex = sqrt(1.0+xdwn2)
                @inbounds for j = 1:n_mesh
                    y = xr[j];ydwn = y*dwn;ydwn2=ydwn^2
                    ey = sqrt(1.0+ydwn2)
                    k2=0.5*(xdwn2 + ydwn2)
                    ree = 1.0/sqrt(ex*ey)
                    fc = fff * ree * hc3 * freg(x,y,2)
                    single_tpe(nd_mpi,nd_mpi2,nd_mpi4,nd_mpi6,nd_mpi8,Fpi2,Fpi4,Fpi6,
                               c1,c2,c3,c4,r_d12,r_d3,r_d5,r_d145,
                               J,pnrank,ts,ws,xdwn,ydwn,xdwn2,ydwn2,k2,pjs,
                               gis,opfs,fc,f_idx,tVs,lsj,tllsj,
                               tdict,V12mom,i,j,to)
                end
            end
        end
    end
    return nothing
end

function Nfac_jj(na,la,ja,nb,lb,jb)
    s =1.0/ sqrt(1.0+delta(na,nb)*delta(la,lb)*delta(ja,jb))
    return s
end

function tm_hats(Lam,Lamp,S,Sp,Jrel,Jrelp,ja,jb,jc,jd)
    t = hat(Lam)^2 * hat(Lamp)^2 * hat(S) * hat(Sp)
    t *= hat(Jrel) *hat(Jrelp) * hat(ja/2) * hat(jb/2) * hat(jc/2) * hat(jd/2)
    return t
end

# pfunc(p,p',P)
# 1: P^2 2: P^2 3: (ddp +ddp') 4: dd p-dd p' 5: P^2
function calc_pPfac(xr,wr,xrP,wrP,
                    n,ell,np,ellp,
                    N,L,Np,Lp,
                    Rnl,RNL,
                    pnrank,dwn,pfunc,to)
    Rnlx = @views Rnl[n+1,ell+1,:]
    Rnly = @views Rnl[np+1,ellp+1,:]
    RNL1 = @views RNL[N+1,L+1,:]
    RNL2 = @views RNL[Np+1,Lp+1,:]    
    r = 0.0
    @inbounds for (i,x) in enumerate(xr)
        xdwn2= (x*dwn)^2
        ex = sqrt(1.0+xdwn2)
        xRw = Rnlx[i]* wr[i] * xdwn2
        @inbounds for (j,y) in enumerate(xr)
            ydwn2= (y*dwn)^2
            ey = sqrt(1.0+ydwn2)
            ree = 1.0/sqrt(ex*ey)
            yRw = Rnly[j]* wr[j] *ydwn2
            @inbounds for (k,P) in enumerate(xrP)
                tRNL = RNL1[k] * RNL2[k] * (P *dwn)^2 * wrP[k] 
                PRw = pfunc(x,ell,y,ellp,P) * tRNL                 
                r += xRw * yRw* ree * PRw
            end
        end
    end
    return r 
end

"""

not yet optimized
"""
function NLOvs(vs,xr,wr,xrP,wrP,Rnl,RNL,
               cg1s,sp_P5_9j,nljsnt,pnrank,ip,
               X9,U6,t5v,
               Jtot,iza,ia,izb,ib,izc,ic,izd,id,to)
    U6_j = U6[Jtot+1]
    na,la,jda = nljsnt[ia]; nb,lb,jdb = nljsnt[ib]
    nc,lc,jdc = nljsnt[ic]; nd,ld,jdd = nljsnt[id]
    mab = 2*na+la+2*nb+lb; mcd = 2*nc+lc+2*nd+ld

    Nab = Nfac_jj(na,la,jda,nb,lb,jdb)
    Ncd = Nfac_jj(nc,lc,jdc,nd,ld,jdd)
    lrmax = jmax + 1
    elab=2*na+la+2*nb+lb; elabp=2*nc+lc+2*nd+ld
    Nrmax = 2*emax
    MN = Ms[pnrank];dwn = 1.0/MN
    mfac = 1.0    
    rc_vs_1 = c_vs_1 * mfac
    rc_vs_2 = c_vs_2 * mfac
    rc_vs_3 = c_vs_3 * mfac
    rc_vs_4 = c_vs_4 * mfac
    rc_vs_5 = c_vs_5 * mfac
    for S=0:1
        tX9 = X9[S+1]
        U6_S = U6_j[S+1]
        for Sp=0:1
            U6_Sp = U6_j[Sp+1]
            tX9p = X9[Sp+1]
            if S ==Sp # term 1,2,5
                ell = ellp = 0
                pfunc = fp_P2
                minLam = max(abs(Jtot-S),abs(la-lb))
                maxLam = min(Jtot+S,la+lb)
                if minLam > maxLam;continue;end
                for Lam = minLam:maxLam
                    L = Lam                        
                    minLamp = max(abs(Jtot-Sp),abs(lc-ld))
                    maxLamp = min(Jtot+Sp,lc+ld)
                    if minLamp > maxLamp;continue;end
                    for Lamp = minLamp:maxLamp
                        Lp = Lamp
                        nja = jda-2*la; njb = jdb-2*lb
                        t5v[1]=la;t5v[2]=nja;t5v[3]=lb
                        t5v[4]=njb;t5v[5]=Lam
                        x1 = get(tX9[Jtot+1],t5v,0.0)
                        
                        njc = jdc-2*lc; njd = jdd-2*ld
                        t5v[1]=lc;t5v[2]=njc;t5v[3]=ld
                        t5v[4]=njd;t5v[5]=Lamp
                        x2 = get(tX9p[Jtot+1],t5v,0.0)
                        if x1==0.0 || x2 == 0.0;continue;end
                        t9j = x1*x2* ifelse((Lam+Lamp)%2==0,1.0,-1.0)
                        t6j = 1.0 # due to ell=ellp=0
                        for N=0:div(elab,2)
                            if (elab-L) % 2 !=0;continue;end
                            n = div(elab-(2*N+L),2);if n<0;continue;end
                            for Np=0:div(elabp,2)
                                if (elabp-Lp) % 2 !=0;continue;end
                                np = div(elabp-(2*Np+Lp),2); if np<0;continue;end
                                HOB1=HObracket(na,la,nb,lb,n,ell,N,L,Lam)
                                HOB2=HObracket(nc,lc,nd,ld,np,ellp,Np,Lp,Lamp)
                                wsyms = t9j * t6j
                                coeff = HOB1*HOB2 * wsyms
                                if coeff == 0.0;continue;end
                                pPfac = calc_pPfac(xr,wr,xrP,wrP,
                                                   n,ell,np,ellp,
                                                   N,L,Np,Lp,
                                                   Rnl,RNL,
                                                   pnrank,dwn,pfunc,to)
                                for Jrel =0:jmax
                                    Jrelp = Jrel
                                    if Jrel == S
                                        rv1 = rc_vs_1 * coeff * 4*pi * pPfac
                                        rv2 = (2*S*(S+1)-3) * rc_vs_2 * coeff * 4*pi * pPfac 
                                        vs[1] += rv1
                                        vs[2] += rv2
                                    end                                    
                                    sum5 = 0.0
                                    for a=0:2:2
                                        t = hat(a)
                                        t *= cg1s[div(a,2)+1]
                                        t *= clebschgordan(Float64,Lp,0,a,0,L,0)
                                        t *= wigner6j(Float64,Jrel,L,Jtot,Lp,Jrelp,a)
                                        t *= sp_P5_9j[div(a,2)+1,S+1]
                                        t *= (-1.0)^(L+Jtot+Jrel)
                                        sum5 += t
                                    end
                                    fac5 = 24*pi *hat(Lp) * hat(S)^2 * sum5
                                    vs[5] += rc_vs_5 * fac5 * pPfac                                    
                                end
                            end
                        end
                    end
                end
            else
                ## term 3,4 S+Sp=1
                pfunc = fp_P2
                for ellidx=1:2
                    ell = ellp = 0
                    if ellidx == 1
                        ell = 1 ; ellp = 0 # p' ~0
                    else
                        ell = 0 ; ellp = 1 # p ~ 0
                    end

                    minLam = max(abs(Jtot-S),abs(la-lb))
                    maxLam = min(Jtot+S,la+lb)
                    if minLam > maxLam;continue;end
                    for Lam = minLam:maxLam
                        minLamp = max(abs(Jtot-Sp),abs(lc-ld))
                        maxLamp = min(Jtot+Sp,lc+ld)
                        if minLamp > maxLamp;continue;end
                        for L = abs(Lam-ell):Lam+ell
                            for Lamp = minLamp:maxLamp
                                for Lp = abs(Lamp-ellp):Lamp+ellp
                                    nja = jda-2*la; njb = jdb-2*lb
                                    njc = jdc-2*lc; njd = jdd-2*ld
                                    t5v[1]=la;t5v[2]=nja;t5v[3]=lb;
                                    t5v[4]=njb;t5v[5]=Lam
                                    x1 = get(tX9[Jtot+1],t5v,0.0)
                                    t5v[1]=lc;t5v[2]=njc;t5v[3]=ld;
                                    t5v[4]=njd;t5v[5]=Lamp
                                    x2 = get(tX9p[Jtot+1],t5v,0.0)
                                    if x1*x2 == 0.0;continue;end
                                    t9j= x1*x2*ifelse((Lam+Lamp)%2==0,1.0,-1.0)
                                    for N=0:div(elab,2)
                                        if (elab-L) % 2 !=0;continue;end
                                        n = div(elab-(2*N+L),2);if n<0;continue;end
                                        for Np=0:div(elabp,2)
                                            if (elabp-Lp) % 2 !=0;continue;end
                                            np = div(elabp-(2*Np+Lp),2); if np<0;continue;end
                                            HOB1=HObracket(na,la,nb,lb,n,ell,N,L,Lam)
                                            HOB2=HObracket(nc,lc,nd,ld,np,ellp,Np,Lp,Lamp)
                                            if HOB1 * HOB2 == 0.0; continue;end
                                            @timeit to "pP" pPfac = calc_pPfac(xr,wr,xrP,wrP,
                                                                               n,ell,np,ellp,
                                                                               N,L,Np,Lp,
                                                                               Rnl,RNL,
                                                                               pnrank,dwn,pfunc,to)
                                            for Jrel = abs(ell-S):ell+S
                                                for Jrelp = abs(ellp-Sp):ellp+Sp
                                                    u6j =U6_S[Lam+1][L+1][ell+1][Jrel-abs(ell-S)+1]
                                                    u6jp=U6_Sp[Lamp+1][Lp+1][ellp+1][Jrelp-abs(ellp-Sp)+1]
                                                    t6j = u6j * u6jp
                                                    coeff = HOB1*HOB2 * t6j * t9j
                                                    if coeff == 0.0;continue;end                                                    
                                                    fac34 = 12.0*pi *hat(Jrel) * hat(Jrelp) * hat(Lp)
                                                    fac34 *= clebschgordan(Float64,Lp,0,1,0,L,0)
                                                    fac34 *= wigner6j(Float64,Jrel,L,Jtot,Lp,Jrelp,1)
                                                    fac34 *= wigner9j(Float64,1,1,1,ell,S,Jrel,ellp,Sp,Jrelp)
                                                    fac4 = fac34 
                                                    fac3 = 2.0*sqrt(2.0) * fac34
                                                    fac3 *= (-1)^(L+Jtot+Jrelp+Sp+1)
                                                    fac4 *= (-1)^(L+Jtot+Jrelp)
                                                    vs[3] += rc_vs_3 * fac3 * pPfac                                    
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
        end
    end
    if pnrank != 2
        vs .*= Nab * Ncd
    end
    return nothing
end

function f_CM_NLO(iterm,x,y,P,s,sp,l,lp,j,jp,L,Lp,J)
    ret = 0.0
    if iterm==1
        ds = delta(s,sp) * delta(l,lp) * delta(l,0) * delta(j,jp) * delta(j,s)* delta(L,Lp)
        ret = 4.0 * pi * P^2 * ds
    elseif iterm == 2
        ds = delta(s,sp) * delta(l,lp) * delta(l,0) * delta(j,jp) * delta(j,s)* delta(L,Lp)
        ret = 4.0 * pi * P^2 * ds * (2.0*s*(s+1)-3)
    elseif iterm == 3
        ds = delta(s+sp,1) *(delta(l,1)*delta(lp,0)*x
                             + delta(l,0)*delta(lp,1)*y) *P
        if ds == 0.0
            ret=0.0
        else
            t1 = ds * 24 * pi *sqrt(2.0) * hat(j) * hat(jp) * hat(Lp)
            t2 = (-1)^(L+J+jp+sp+1)
            t3 = clebschgordan(Float64,Lp,0,1,0,L,0)
            t4 = wigner6j(Float64,j,L,J,Lp,jp,1)
            t5 = wigner_9j(Float64,1,1,1,l,s,j,lp,sp,jp)
            ret = t1 * t2 * t3 * t4 * t5
        end
    elseif iterm == 4
        ds = delta(s+sp,1) *(delta(l,1)*delta(lp,0)*x
                             - delta(l,0)*delta(lp,1)*y) *P
        if ds == 0.0
            ret=0.0
        else
            t1 = ds * 12 * pi *sqrt(2.0) * hat(j) * hat(jp) * hat(Lp)
            t2 = (-1)^(L+J+jp)
            t3 = clebschgordan(Float64,Lp,0,1,0,L,0)
            t4 = wigner6j(Float64,j,L,J,Lp,jp,1)
            t5 = wigner_9j(Float64,1,1,1,l,s,j,lp,sp,jp)
            ret = t1 * t2 * t3 * t4 * t5
        end
    elseif iterm == 5
        ds = delta(s+sp,1) * delta(l,lp) * delta(l,0) * delta(j,jp)
        if ds == 0; return 0.0;end
        t1 = 24 * pi * ds * hat(Lp) * hat(s)^2 * P^2 * (-1)^(L+J+j)
        sum = 0.0
        for a=0:2:2
            tmp = hat(a) *  clebschgordan(Float64,1,0,1,0,a,0) *clebschgordan(Float64,L,0,Lp,0,a,0)
            tmp *= wigner6j(j,L,J,Lp,j,a) * wigner9j(1,1,a,1//2,1//2,s,1//2,1//2,s)
            sum += tmp
        end
        ret = sum * t1
    end
    return ret
end

function write_vmom(xr,V12mom,tdict,pnrank,llpSJ_s;label="")
    tx1d=""
    itt = itts[pnrank]    
    for i= 1:n_mesh
        x = xr[i]
        for j = 1:n_mesh
            y = xr[j]
            if y != x; continue;end
            tx = @sprintf("%18.8e", x/hc)
            for tmp in llpSJ_s
                l,lp,S,J = tmp
                V12idx = get(tdict,[itt,l,lp,S,J],-1)
                if V12idx==-1
                    tx *= @sprintf("%18.8e",0.0)
                else
                    v = V12mom[V12idx][i,j]
                    tx *= @sprintf("%18.8e",v)
                end
            end
            tx *= "\n"
            tx1d *= tx
        end
    end        
    io = open("vmom_1d_"*label*".dat","w")
    println(io,rstrip(tx1d))
    close(io)
    return nothing
end

function cib_lsj_opep(opfs,x,y,mpi2,nterm,J,pnrank,nfac,
                      ts,ws,tVs)
    z = (mpi2+x^2+y^2) / (2.0*x*y)                        
    QJ   = Q_L2(z,J,ts,ws)
    QJm1 = 0.0
    if J>0;QJm1=Q_L2(z,J-1,ts,ws);end
    IJ0 = nfac * QJ
    IJ1 = nfac * (z * QJ -delta(J,0))
    IJ2 = nfac * (J*z* QJ + QJm1) /(J+1) 
    IJ3 = nfac * sqrt(J/(J+1)) * (z* QJ - QJm1)    
    v1  = opfs[1] * IJ0 + opfs[2] *IJ1
    v2  = opfs[3] * IJ0 + opfs[4] *IJ2
    v3 = opfs[5] * IJ0 + opfs[6] *IJ1
    v4 = opfs[4] * IJ0 + opfs[3] *IJ2
    v5 = opfs[7] * IJ3
    v6 = -v5                        
    if J==0; v2=v4=v5=v6=0.0;end
    if J%2==1 && pnrank!=2; v1 =0.0;end
    if J%2==0 && pnrank!=2; v2 =0.0;end
    if J%2!=0 && pnrank!=2; v3=v4=v5=v6=0.0;end
    v34 = -sqrt(J*(J+1)) *(v3-v4)
    v56 = sqrt(J*(J+1)) * (v5+v6)
    d2j1 = 1.0/(2*J+1)
    if nterm == 1
        tVs[1] = v1 
        tVs[2] = v2 
        tVs[3] = d2j1 * ((J+1)* v3 + J*v4-v56)
        tVs[4] = d2j1 * ( J*v3 + (J+1)*v4 +v56) 
        tVs[5] = -d2j1 * (v34-(J+1)*v5+J*v6)
        tVs[6] = -d2j1 * (v34+J*v5-(J+1)*v6)
        tVs .*= ifelse(pnrank==2,-1.0,1.0)
    else
        is = J%2 + 1
        it = is%2 +1
        ttis = ifelse(is==2,-2.0,2.0)
        ttit = ifelse(it==2,-2.0,2.0)
        tVs[1] += ttis * v1
        tVs[2] += ttit * v2
        tVs[3] += d2j1 * ((J+1)* (ttis*v3) + J*(ttis*v4)-(ttis*v56))
        tVs[4] += d2j1 * ( J*(v3*ttis) + (J+1)*(ttis*v4) +(ttis*v56)) 
        tVs[5] += -d2j1 * ((ttis*v34)-(J+1)*(ttis*v5)+J*(ttis*v6))
        tVs[6] += -d2j1 * ((ttis*v34)+J*(ttis*v5)-(J+1)*(ttis*v6))
    end
    return nothing 
end

""" 
To calculate Legendre functions of second kind 
by Gauss-Legendre quadrature
"""
function Q_L2(z,J,ts,ws)
    s = 0.0
    @inbounds for (i,t) in enumerate(ts)
        s += ws[i] * (1.0-t^2)^J / (z-t)^(J+1) 
    end
    return s / 2.0^(J+1)
end
  
function calc_Vmom!(pnrank,V12mom,tdict,xr,LEC,LEC2,
                    l,lp,S,J,pfunc,n_reg,to;
                    pnranks=[1,2,3],co="",is_n3lo=false)
    itt = itts[pnrank]; MN = Ms[pnrank]; dwn = 1.0/MN
    V12idx = get(tdict,[itt,lp,l,S,J],-1)
    if V12idx == -1;return nothing;end
    V = V12mom[V12idx]
    @inbounds for i= 1:n_mesh
        x = xr[i]; ex = sqrt(1.0+(dwn.*x)^2)
        @inbounds for j = 1:n_mesh
            y = xr[j]; ey = sqrt(1.0+(dwn.*y)^2)
            ree = 1.0/sqrt(ex*ey)
            fac = pfunc(x,y,LEC,LEC2) * freg(x,y,n_reg) *ree
            V[i,j] += fac
        end
    end
    return nothing
end

"""
    Rnl_all_ab(lmax,br,n_mesh,xr_fm)

Returns array for radiul functions 
(prop to generalized Laguerre polynomials)
HO w.f. in momentum space.
Rnlk(l,n,k)=sqrt(br) * R(n,L,Z) *Z with Z=br*k (k=momentum in fm^-1)
"""
function Rnl_all_ab(lmax,br,n_mesh,xr_fm)
    Rnl = zeros(Float64,Nnmax+1,lmax+1,n_mesh)
    for l=0:lmax
        for kidx=1:n_mesh
            pb = br * xr_fm[kidx]
            pb2 = pb^2
            fexp = exp(-0.5*pb2)
            fpow = pb^(l+1)
            for n=0:Nnmax
                fac = sqrt(2.0*br*factorial(n) / gamma(n+l+3/2)) * fexp * fpow
                fac *= genLaguerre(n,l+1//2,pb2)
                Rnl[n+1,l+1,kidx]= fac
            end
        end
    end
    return Rnl
end

function bstate()
    #iz,lz1,lz2,isz,jz
    numst2 = [ [0,0,0,0,0] ];deleteat!(numst2,1)
    dict_numst = [ Dict([0,0]=>1) for i=1:3]
    for i=1:3; delete!(dict_numst[i],[0,0]); end
    arr_numst = [[[[ zeros(Int64,j+iss-abs(j-iss)+1) for ll1=abs(j-iss):j+iss ] for j=0:jmax ] for iss=0:1 ] for pnrank=1:3]
    
    num=0
    ## pp iz = -2
    itt = 1
    for iss=0:1
        for j=0:jmax
            for ll1=abs(j-iss):j+iss
                for ll2=abs(j-iss):j+iss
                    if (-1)^ll1 != (-1)^ll2;continue;end
                    if itt != Int( (1+(-1)^(ll1+iss))/2);continue;end
                    if (-1)^(ll1+iss+itt) != -1;continue;end
                    num=num+1
                    push!(numst2,[-2,ll1,ll2,iss,j])
                    dict_numst[1][[-2,ll1,ll2,iss,j]] = num
                    arr_numst[1][iss+1][j+1][ll1-abs(j-iss)+1][ll2-abs(j-iss)+1] = num
                end
            end
        end
    end
    ##  nn iz = 2
    pnrank=3
    for iss=0:1
        for j=0:jmax
            for ll1=abs(j-iss):j+iss
                for ll2=abs(j-iss):j+iss
                    if (-1)^ll1 != (-1)^ll2;continue;end
                    if itt != Int( (1+(-1)^(ll1+iss))/2);continue;end
                    if (-1)^(ll1+iss+itt) != -1;continue;end
                    num=num+1
                    push!(numst2,[2,ll1,ll2,iss,j])
                    dict_numst[pnrank][[2,ll1,ll2,iss,j]] = num
                    arr_numst[pnrank][iss+1][j+1][ll1-abs(j-iss)+1][ll2-abs(j-iss)+1] = num
                end
            end
        end
    end
    ## pn iz = 0
    pnrank=2
    for iss=0:1
        for ktt=1:2
            if (iss==0 && ktt==1) || (iss==1 && ktt==2);itt=1;end
            if (iss==0 && ktt==2) || (iss==1 && ktt==1);itt=0;end
            for j=0:jmax
                for ll1=abs(j-iss):j+iss
                    for ll2=abs(j-iss):j+iss
                        if (-1)^ll1 != (-1)^ll2;continue;end
                        if itt != Int((1+(-1)^(ll1+iss))/2) ;continue;end
                        if (-1)^(ll1+iss+itt) != -1;continue;end
                        num=num+1
                        push!(numst2,[0,ll1,ll2,iss,j])
                        dict_numst[pnrank][[0,ll1,ll2,iss,j]] = num
                        arr_numst[pnrank][iss+1][j+1][ll1-abs(j-iss)+1][ll2-abs(j-iss)+1] = num
                    end
                end
            end
        end
    end
    println("# of two-body states $num")
    return numst2,dict_numst,arr_numst
end

function momplot(xr,V12mom,tdict,pnrank,llpSJ_s;ctext="",fpath="")
    tfdat = []
    if fpath != ""; xf,yfs = compfdat(fpath); end
    
    itt = itts[pnrank]
    tv = zeros(Float64,n_mesh)
    for vidx = 1:7
        l,lp,S,J= llpSJ_s[vidx]
        V12idx = get(tdict,[itt,l,lp,S,J],-1)
        if V12idx==-1;println("V12idx==$V12idx");continue;end
        V = V12mom[V12idx]
        tx = ""
        cS = ifelse(S==0,"1","3")
        if l==lp
            tx=ctext*cS*chara_L[l+1]*string(J)
        else
            tx =ctext*cS*chara_L[l+1]*string(J)
            tx *= "_"*cS*chara_L[lp+1]*string(J)
        end
        for i =1:n_mesh;  tv[i] = V[i,i]; end
        if fpath != ""; tfdat = [xf,yfs[vidx]];end
        pw_plt(tx,xr,V,tv,pnrank;fdat=tfdat)
    end
    return nothing
end


function Vrel(V12mom,numst2,xr_fm,wr,n_mesh,Rnl,to)
    nstmax = length(numst2)
    V12ab = [zeros(Float64,Nnmax+1,Nnmax+1) for i=1:nstmax]
    x = zeros(Float64,Nnmax+1,n_mesh)
    tl = [0,0,0]
    for num = 1:nstmax
        Vtmp = V12mom[num]
        Vab = V12ab[num]
        iz,l1,l2,isz,jz = numst2[num]
        @inbounds for n1 = 0:Nnmax
            tR = @views Rnl[n1+1,l1+1,:]
            tx = @views x[n1+1,:]
            @inbounds for k2=1:n_mesh
                sum=0.0
                @inbounds for k1=1:n_mesh
                    vk1k2=Vtmp[k1,k2]
                    if vk1k2 == 0.0; continue;end
                    pkx1=tR[k1]
                    sum += pkx1*vk1k2*xr_fm[k1]*wr[k1]
                end
                #x[n1+1,k2]=sum*wr[k2]*xr_fm[k2]
                tx[k2]=sum*wr[k2]*xr_fm[k2]
            end
        end
        @inbounds for n1=0:Nnmax
            @inbounds for n2=0:Nnmax
                tR = @views Rnl[n2+1,l2+1,:]
                vsum=0.0
                @inbounds for k2=1:n_mesh
                    pkx2=tR[k2]
                    vsum += x[n1+1,k2]*pkx2
                end
                phase=(-1.0)^(n1+n2) 
                Vab[n1+1,n2+1]=phase*vsum #+ vcoul[num][n1,n2] # coulomb is not implemented now                
            end
        end
    end
    return V12ab
end

function hw_formula(A,fnum)
    hw = 0.0
    if fnum == 2
        #J. Blomqvist and A. Molinari, Nucl. Phys. A106, 545 (1968).
        hw = 45.0 * (A^(-1.0/3.0)) -25.0 * (A^(-2.0/3.0))
    else
        hw = 41.0 * (A^(-1.0/3.0))
    end
    return hw
end
