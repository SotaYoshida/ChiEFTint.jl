"""
fdat: fortran data for closs-check 
(SY confirmed that both match up to N3LO (w/o valence op.))
"""
function pw_plt(tx,xr,z,zds,pnrank;fdat=[])
    tls = ["_pp","_pn","_nn"]
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
