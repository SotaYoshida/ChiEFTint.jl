
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

    LECs =zeros(Float64,39) 
    idxLECs=Dict{String,Int64}()
    dLECs=Dict{String,Float64}()
    read_LECs!(LECs,idxLECs,dLECs;initialize=true)

    @timeit to "calcV" begin
        # ***Leading Order (LO)***
        ### OPEP
        @timeit to "OPEP" begin
            OPEP(ts,ws,xr,V12mom,dict_numst,to,lsjs,llpSJ_s,tllsj,opfs)
        end
        ### contact term
        LO(xr,dLECs,V12mom,dict_numst,to)
        if chi_order >= 1
            # ***Next-to-Leading Order(NLO)***
            ### contact term Q^2
            NLO(xr,dLECs,V12mom,dict_numst,to)

            ### TPE terms (NLO,NNLO,N3LO)
            @timeit to "TPE" begin
                tpe(dLECs,ts,ws,xr,V12mom,dict_numst,to,llpSJ_s,lsjs,tllsj,opfs)
            end
            # *** N3LO ***
            ### contact Q^4
            if chi_order >= 3
                N3LO(xr,dLECs,V12mom,dict_numst,to)
            end
        end
    end

    if is_plot
        @timeit to "momplot" begin
            tpnrank = 2 # 1:pp 2:pn 3:nn pnrank for momplot 
            #write_vmom(xr,V12mom,dict_numst[tpnrank],tpnrank,llpSJ_s)
            momplot(xr,V12mom,dict_numst[tpnrank],tpnrank,llpSJ_s)
        end
    end
    
    #transform mom. int. to HO matrix element
    @timeit to "Vtrans" begin
        V12ab = Vrel(V12mom,numst2,xr_fm,wr,n_mesh,Rnl,to)
        jab_max = 4 * emax + 2
        Numpn= Dict( [0,0,0,0] => 0 ) ;delete!(Numpn,[0,0,0,0])   
        infos,izs_ab,nTBME = make_sp_state(jab_max,Numpn)
        println("# of channels 2bstate ",length(infos)," #TBME = $nTBME")
        @timeit to "TM" TMtrans(dLECs,xr,wr,xrP,wrP,Rnl,RNL,
                                nTBME,infos,izs_ab,
                                Numpn,V12ab,arr_numst,to)
    end

    show(to, allocations = true,compact = false);println("")
    
    return nothing    
end


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


#############
##basis_states.jl
########

#In this file, most functions are naive translation from the Fortran codes,
#which are used in the references below,
#and its modification by Takayuki Miyagi (@TRIUMF)
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


function def_sps_snt(emax,target_nlj)
    nljsnt = [ [0,0] ]; deleteat!(nljsnt,1)    
    tzs = [-1,1]
    dict = Dict(0=>0); delete!(dict,0)
    idx = 0; hit = 0
    for pn = 1:2
        tz = tzs[pn]
        for temax = 0:emax
            for l = temax%2:2:temax
                n = div(temax - l,2)
                jmin = 2*l-1
                if jmin < 1; jmin=1;end
                for j=jmin:2:2*l+1
                    push!(nljsnt,[n,l,j,tz])
                    idx += 1
                    for (k,tmp) in enumerate(target_nlj)
                        if tmp[1]==n && tmp[2]==l && tmp[3]==j
                            dict[idx] = k + ifelse(pn==2,length(target_nlj),0)
                        end
                    end
                end
            end
        end
    end
    return nljsnt,dict
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


function read_LECs!(LECs,idxLECs,dLECs;initialize=false,inpf="LECs.jl")
    if initialize
        f=open(inpf,"r");lines=readlines(f);close(f)
        hit = 0
        for line in lines
            tl = split(rstrip(line))
            if length(tl) == 0;continue;end
            if tl[1] != "const"; continue;end
            hit += 1
            LEC = parse(Float64,tl[4])
            LECs[hit] = LEC
            idxLECs[tl[2]] = hit
            dLECs[tl[2]] = LEC
        end
    else
        for tkey in keys(idxLECs)
            idx = idxLECs[tkey]
            dLECs[tkey] = LECs[idx]
        end
    end
    return nothing
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


