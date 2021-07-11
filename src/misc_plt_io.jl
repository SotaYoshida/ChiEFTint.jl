const chara_L = ["S","P","D","F","G","H","I","J","K","L","M","N"]
const chara_l = ["s","p","d","f","g","h","i","j","k","l","m","n"]
chara_SLJ(S,L,J) = "{}^"*ifelse(S==0,"1","3")*chara_L[L+1]*"_{"*string(J)*"}"
delta(a,b) = ifelse(a==b,1.0,0.0)
hat(a) = sqrt(2.0*a+1.0)

function delta_arr(a,b)
    hit = 0
    for i=1:length(a)
        hit += ifelse(a[i]==b[i],1,0)
    end
    return ifelse(length(a)==hit,1.0,0.0)
end

function c_orbit(tarr)
    n,l,j=tarr
    tx = string(n)*lowercase(chara_L[l+1])*string(j)
    return tx
end

const Mat_pw_NLO = 4*pi*[ 1.0  0.25 -3.0 -3/4  0.0 -1.0 -0.25;
                          -2/3  1/6 -2/3  1/6 -2/3  2.0 -0.5;
                          -2/3  1/6  2.0 -1/2  0.0  2/3 -1/6;
                          -2/3  1/6 -2/3  1/6 -1/3 -4/3  1/3;
                          1.0   1/4  1.0  1/4  0.0  1/3  1/12;
                          0.0   0.0  0.0  0.0  0.0  -2*sqrt(2)/3 -sqrt(2)/6;
                          -2/3  1/6 -2/3  1/6  1/3  0.0  0.0]

const InvMat_pw_NLO = 1.0/(4*pi) .* [0.125 -0.0625 -3/16 -0.1875 0.375 0.0 -0.3125;
                                     0.5 0.25 0.75 0.75 1.5 0.0 1.25;
                                     -0.125 -0.0625 0.1875 0.0 0.125 1/(4*sqrt(2)) -0.125;
                                     -0.5 0.25 -0.75 0.0 0.5 sqrt(2)/2 0.5;
                                     0.0 -0.5 0.0 -0.75 0.0 0.0 1.25;
                                     0.0 0.125 0.0 -0.1875 0.0 -3/(4*sqrt(2)) 0.0625;
                                     0.0 -0.5 0.0 0.75 0.0 -3/sqrt(2) -0.25]
function rm_comment(lines)
    nlines = []
    for line in lines
        line = strip(line)
        if length(line)>1
            if startswith(line,"!")
                continue
            end
        end
        push!(nlines,line)
    end
    return nlines
end


function rm_nan(array)
    na = []
    for tmp in array
        if tmp != "";push!(na,tmp); end
    end
    return na
end

function write_tbme(io,ndim,izs,Jtot,vv,nljsnt,nljdict;ofst=0)
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
                tv = vv[i,j]
                izc,ic,izd,id= izs[j]
                cpnrank = ""
                if iza==izb==-1;cpnrank="pp";end
                if iza==izb==1;cpnrank="nn";end
                if iza != izb ==1;cpnrank="pn";end
                a = ifelse(iza==-1,ia,ia+ofst)
                b = ifelse(izb==-1,ib,ib+ofst)
                c = ifelse(izc==-1,ic,ic+ofst)
                d = ifelse(izd==-1,id,id+ofst)
                
                if target_nlj != []
                    if get(nljdict,a,0) == 0; continue;end
                    if get(nljdict,b,0) == 0; continue;end
                    if get(nljdict,c,0) == 0; continue;end
                    if get(nljdict,d,0) == 0; continue;end
                end
                
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
                ca = c_orbit(nljsnt[fa-ifelse(iza==1,ofst,0)])
                cb = c_orbit(nljsnt[fb-ifelse(izb==1,ofst,0)])
                cc = c_orbit(nljsnt[fc-ifelse(izc==1,ofst,0)])
                cd = c_orbit(nljsnt[fd-ifelse(izd==1,ofst,0)])
                if target_nlj != []
                    fa = nljdict[fa]
                    fb = nljdict[fb]
                    fc = nljdict[fc]
                    fd = nljdict[fd]
                end
                print(io,@sprintf("%5i", fa))
                print(io,@sprintf("%5i", fb))
                print(io,@sprintf("%5i", fc))
                print(io,@sprintf("%5i", fd))
                print(io,@sprintf("%6i", Jtot))
                println(io,@sprintf("%18.10f", tv))
            end
        end
    end
    return nothing
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

"""
fdat: fortran data for closs-check 
(SY confirmed that both match up to N3LO (w/o valence op.))
"""
function pw_plt(tx,xr,z,zds,pnrank;fdat=[])
    tls = ["pp","pn","nn"]
    cpnrank = tls[pnrank]    
    xr *= 1.0/hc; yr = copy(xr)
    fig = plt.figure(figsize=(10,4))
    axs = [fig.add_subplot(121),fig.add_subplot(122)]
    axs[1].set_xlabel(latexstring("p ")*" [fm"*latexstring("^{-1}")*"]")
    axs[1].set_ylabel(latexstring("p'")*" [fm"*latexstring("^{-1}")*"]")
    axs[1].contourf(xr, yr, z)
    axs[2].set_xlabel(latexstring("p=p' ")*" [fm"*latexstring("^{-1}")*"]")
    axs[2].set_ylabel(latexstring("V(p,p)")*" [MeV fm"*latexstring("^3")*"]")
    axs[2].plot(xr,zds,marker="o",markersize=2)
    if fdat != []
        axs[2].plot(fdat[1],fdat[2],marker="x",markersize=2,alpha=0.4,label="Fortran")
        axs[2].legend()
    end
    axs[2].grid(color="gray",linestyle="dotted")
    plt.savefig("pic/chiEFT_"*tx*"_"*cpnrank*".pdf",pad_inches=0)
    plt.close()
end

function compfdat(inpf)
    f = open(inpf,"r");lines = readlines(f);close(f)
    xs = Float64[]
    ys = [ Float64[] for i=1:7]
    for line in lines
        tl = split(line)
        tmp = map(x->parse(Float64,x),tl)
        push!(xs,tmp[1])
        for i=2:length(tmp); push!(ys[i-1],tmp[i]);end
    end
    return xs,ys
end
