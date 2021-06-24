function f_sq!(opf,xdwn,ydwn)
    opf .= 0.0
    opf[1] = -2.0 * (xdwn^2 + ydwn^2)
    opf[2] = 4.0 * xdwn * ydwn
    opf[3] = -opf[1]; opf[4] = -opf[2]
    opf[5] =  opf[2]; opf[6] =  opf[1]
    opf[7] = 2.0*( xdwn^2 - ydwn^2)
    opf[8] = -opf[7]
    return nothing
end

function f_ss!(opf)
    opf .= 0.0
    opf[1] = -6.0
    opf[2] = 0.0
    opf[3] = 2.0
    opf[4] = 0.0
    opf[5] = 0.0
    opf[6] = opf[3]
    opf[7] = -opf[3]
    opf[8] = opf[7]
    return nothing
end

function f_c!(opf)
    opf .= 0.0
    opf[1] = 2.0
    opf[2] = 0.0
    opf[3] = opf[1]
    opf[4] = opf[2]
    opf[5] = opf[2]
    opf[6] = opf[1]
    opf[7] = -opf[1]
    opf[8] = opf[7]
    return nothing
end

function f_ls!(opf,xdwn,ydwn)
    xy2 = 2.0 * xdwn * ydwn
    opf .= 0.0
    opf[4] = -xy2
    opf[5] = -xy2
    opf[10] = xy2
    opf[11] = -xy2
    return nothing
end

function f_sl!(opf,xdwn,ydwn)
    xxyy2 = 2.0 * xdwn^2 * ydwn^2
    opf .= 0.0
    opf[1] = -xxyy2
    opf[3] = -xxyy2
    opf[6] =  xxyy2
    opf[7] = -xxyy2
    opf[8] = -xxyy2
    opf[9] = 2.0 * xxyy2
    return nothing
end

function single_tpe(nd_mpi,nd_mpi2,nd_mpi4,nd_mpi6,nd_mpi8,Fpi2,Fpi4,Fpi6,
                    c1,c2,c3,c4,r_d12,r_d3,r_d5,r_d145,
                    J,pnrank,ts,ws,xdwn,ydwn,xdwn2,ydwn2,k2,pjs,
                    gis,opfs,fc,f_idx,tVs,lsj,tllsj,
                    tdict,V12mom,V_i,V_j,to)
    f_T,f_SS,f_C,f_LS,f_SL = opfs
    f_sq!(f_T,xdwn,ydwn);f_ls!(f_LS,xdwn,ydwn);f_sl!(f_SL,xdwn,ydwn)
    for i=1:length(gis); gis[i] .= 0.0; end
    @inbounds for (n,t) in enumerate(ts) 
        q2 = xdwn2 + ydwn2 -2.0*xdwn*ydwn*t
        q = sqrt(q2)
        w2 = 4.0*nd_mpi2 + q2
        w = sqrt(w2)
        tw2 = 2.0*nd_mpi2 + q2
        Lq = w/q * log((w+q)/(2.0*nd_mpi) )
        Aq = atan(q/(2.0*nd_mpi))  / (2.0*q)
        w2Aq = w2 * Aq
        tw2Aq = tw2 * Aq

        ## Tensor term: Vt
        gi = gis[1]
        f_NLO_Vt =  -3.0* gA4  /  (64.0 * pi^2 * Fpi4)
        f_NNLO_Vt = 9.0 * gA4 / (512.0 * pi * Fpi4)
        f_N3LO_Vt = gA4 / (32.0 * pi^2 * Fpi4)
        f_N3LO_2l_Vt = -gA2 * r_d145 /(32.0*pi^2 *Fpi4)
        # NLO
        tmp_s = Lq * f_NLO_Vt
        # NNLO
        if chi_order >= 2
            it_pi = (nd_mpi + w^2 * Aq) / 3.0
            tmp_s +=  (tw2Aq +it_pi) * f_NNLO_Vt
        end
        # N3LO
        if chi_order >= 3
           tmp_s += Lq * (k2 +3.0/8.0 *q2 + nd_mpi4 /w2) *f_N3LO_Vt
           tmp_s += w2 * Lq * f_N3LO_2l_Vt
        end
        for nth=1:7; gi[nth] += pjs[nth][n] *tmp_s *ws[n] ; end

        ##Tensor term: Wt
        if chi_order >=2
            gi = gis[2]       
            fac_NNLO = - gA2 / (32.0 * pi * Fpi4)
            f_N3LO_a_Wt = c4^2 / (96.0*pi^2 *Fpi4)
            f_N3LO_b_Wt = - c4 / (192.0*pi^2 *Fpi4)
            f_N3LO_c_Wt = 1.0/ (1536.0* pi^2 * Fpi4)        
            f_N3LO_2l_Wt= gA4 / (2048.0 * pi^2 * Fpi6)
            term1 = (c4+0.25) * w2
            term2 = -0.125*gA2 * (10.0 *nd_mpi2 + 3.0 *q2)
            it_pi = gA2 /8.0 * (nd_mpi + w2Aq)
            tmp_s = ((term1+term2)*Aq +it_pi) * fac_NNLO            
            if chi_order >= 3                
                tmp_s += Lq * w2 * f_N3LO_a_Wt
                tmp_s += Lq * ( gA2 *( 16.0 * nd_mpi2 + 7.0 * q2) -w2) * f_N3LO_b_Wt
                tmp_s += Lq * ( 4.0*gA4 *(7.0*nd_mpi2 + 17.0/4.0 * q2
                                          + 4.0 * nd_mpi4 /w2)
                                -32.0*gA2 *(nd_mpi2 + 7.0/16.0 * q2) + w2) * f_N3LO_c_Wt
                tmp_s += f_N3LO_2l_Wt * w2Aq *( w2*Aq + 2.0*nd_mpi * (1.0+2.0*gA2))
            end
            for nth=1:7; gi[nth] += pjs[nth][n] * tmp_s * ws[n];end
        end

        # sigma-sigma term: Vs
        gi = gis[3]
        f_NLO_Vs = -3.0 * gA4  / (64.0 * pi^2 * Fpi4)
        f_NNLO_Vs = 9.0 * gA4 / (512.0 * pi * Fpi4)
        f_N3LO_Vs = gA4 / (32.0 * pi^2 * Fpi4)
        f_N3LO_2l_Vs = -gA2 * r_d145 /(32.0*pi^2 *Fpi4)
        tmp_s = Lq * f_NLO_Vs
        if chi_order >= 2 # NNLO
            it_pi =  (nd_mpi + w2Aq)/3.0
            tmp_s += (tw2Aq+it_pi) * f_NNLO_Vs
        end
        if chi_order >= 3
           tmp_s +=  Lq * (k2 +3.0/8.0 *q2 + nd_mpi4 /w2) *f_N3LO_Vs
           tmp_s +=  w2 * Lq * f_N3LO_2l_Vs
        end
        tmp_s *= -q2
        for nth=1:7; gi[nth] += pjs[nth][n] * tmp_s * ws[n];end

        # sigma-sigma term: Ws
        if chi_order >= 2
            gi = gis[4]
            f_NNLO_Ws = - gA2 / (32.0 * pi * Fpi4)
            f_N3LO_a_Ws = c4^2 / (96.0*pi^2 *Fpi4)
            f_N3LO_b_Ws = - c4 / (192.0*pi^2 *Fpi4)
            f_N3LO_c_Ws = 1.0/ (1536.0* pi^2 * Fpi4)
            f_N3LO_2l_Ws= gA4 / (2048.0 * pi^2 * Fpi6)
            term1 = (c4+0.25) * w2Aq
            term2 = -0.125*gA2 * (10.0 *nd_mpi2 + 3.0 *q2) *Aq
            it_pi = gA2 /8.0 * (nd_mpi + w2Aq)
            tmp_s =  (term1+term2+it_pi) * f_NNLO_Ws
            if chi_order >= 3
                tmp_s += Lq * w2 * f_N3LO_a_Ws
                tmp_s += Lq * ( gA2 *(16.0*nd_mpi2 + 7.0 * q2) -w2) * f_N3LO_b_Ws
                tmp_s += Lq * ( 4.0*gA4 *(7.0*nd_mpi2 +
                                          17.0/4.0 * q2 + 4.0 * nd_mpi4 /w2)
                                -32.0 *gA2 * (nd_mpi2 + 7.0/16.0 * q2)
                                + w2) * f_N3LO_c_Ws
                tmp_s += f_N3LO_2l_Ws * w2Aq *( w2Aq + 2.0*nd_mpi * (1.0+2.0*gA2))
            end        
            tmp_s *= -q2        
            for nth=1:7; gi[nth] += pjs[nth][n] * tmp_s * ws[n];end
        end

        # Central term: Vc
        if chi_order >= 2
            gi = gis[5]        
            f_NNLO_Vc = 3.0 * gA2 / (16.0 * pi * Fpi4)
            f_N3LO_f1_Vc = 3.0/ (16.0 * pi^2 * Fpi4)
            f_N3LO_f6_Vc = - gA2 / (32.0 * pi^2 * Fpi4)
            f_N3LO_2l_Vc = 3.0 * gA4 /(1024.0 *pi^2 * Fpi6)                             
            term1 = gA2 * nd_mpi^5 / (16.0*w^2)
            term2 = - ( 2.0 * nd_mpi2 * (2.0*c1 -c3) -q2*(c3+3.0/16.0 * gA2)) * tw2Aq
            it_pi = -gA2 /16.0 * (nd_mpi * w2 + tw2 * tw2Aq )
            tmp_s =  (term1+term2+it_pi) * f_NNLO_Vc
            if chi_order >= 3
                brak1  = ((c2*w2)/6.0 + c3*tw2 -4.0*c1*nd_mpi2)^2 + (c2^2 *w2^2)/45.0
                brak6  = (c2-6.0*c3) * q2^2 + 4.0*(6.0*c1+c2-3.0*c3) *q2 *nd_mpi2
                brak6 += 6.0*(c2-2.0*c3) * nd_mpi2^2
                brak6 += 24.0 *(2.0*c1+c3) * nd_mpi2^3 / w2 
                Mm2cor = Lq *(2.0*nd_mpi2^4 / w^4 +
                          8.0*nd_mpi2^3 / w2 -q2^2
                              -2.0*nd_mpi2^2 ) + nd_mpi2^3 /(2.0 *w2)
                term3s  = Lq *( brak1 * f_N3LO_f1_Vc +  brak6 * f_N3LO_f6_Vc)
                term3s += (f_N3LO_f6_Vc *gA2) *Mm2cor
                tmp_s += term3s
                tmp_s += f_N3LO_2l_Vc * tw2Aq * (
                    (nd_mpi2 + 2.0 * q2)*(2.0*nd_mpi + tw2Aq) +4.0*gA2 * nd_mpi *tw2)
            end
            for nth=1:7; gi[nth] += pjs[nth][n] * ws[n] * tmp_s;end        
        end

        # Central term: Wc
        gi = gis[6]
        f_NLO_Wc = -1.0 / (384.0 * pi^2 * Fpi4)
        f_NNLO_Wc = gA2 / (128.0 * pi * Fpi4)
        f_N3LO_a_Wc = -c4 / (192.0 * pi^2 * Fpi4)
        f_N3LO_b_Wc = -1.0 / (768.0 * pi^2 * Fpi4)
        f_N3LO_2l_Wc= 1.0/(18432.0 * pi^4 * Fpi6)

        brak = 4.0*nd_mpi2*(5.0 * gA4 - 4.0 * gA2 -1.0)
        brak +=  q2 * (23.0*gA4 -10.0*gA2 -1.0) +48.0* gA4 *nd_mpi4 / w2
        tmp_s = Lq * brak * f_NLO_Wc
         
        if chi_order >= 2 # NNLO
            term1 = 3.0* gA2 * nd_mpi^5 / w^2
            term2 = - (4.0*nd_mpi^2 + 2.0*q2 -gA2*(4.0*nd_mpi2 +3.0*q2)) * tw2Aq
            it_pi = gA2 * ( nd_mpi * w2 + tw2 * tw2Aq)
            tmp_s += (term1+term2+it_pi) * f_NNLO_Wc
        end
        if chi_order >= 3 # N3LO
            terma = (gA2 *(8.0*nd_mpi2+5.0*q2) +w2) * q2 * Lq
            brak1 = Lq * ( 8.0*gA2 * (11.0/4.0 *q2^2 +5.0*nd_mpi2*q2
                                      + 3.0*nd_mpi2^2 -6.0*nd_mpi2^3 /w2
                                      -k2*(8.0*nd_mpi2 + 5.0*q2))
                           + 4.0*gA4 *(k2*(20.0*nd_mpi2 +7.0*q2 -16.0*nd_mpi2^2 / w2)
                                       +16.0 * nd_mpi2^4 /w2^2 + 12.0 *nd_mpi2^3 /w2
                                       -27.0/4.0* q2^2-11.0*nd_mpi2*q2-6.0*nd_mpi2^2)
                           +(q2-4.0*k2)*w2)
            brak2 =  16.0*gA4 *nd_mpi2^3 / w2
            termb = brak1 + brak2
            tmp_s += (terma * f_N3LO_a_Wc + termb * f_N3LO_b_Wc)
            # two-loop correction (EM:D.20)
            term1 = 192.0* pi^2 *Fpi2*w2*r_d3* (2.0*gA2*tw2 -3.0/5.0 *(gA2-1.0) *w2)
            brak = 384.0 * pi^2 *Fpi2 * (tw2 * r_d12 + 4.0*nd_mpi2 * r_d5)
            brak += Lq * (4.0*nd_mpi2 *(1+2.0*gA2) + q2 * (1.0+5.0*gA2))
            brak += - (q2 * (5.0+13.0*gA2)/3.0 + 8.0 * nd_mpi2 * (1.0+2.0*gA2))
            term2 = (6.0*gA2 *tw2 - (gA2 -1.0) *w2)  * brak
            tmp_s += (term1 + term2) * Lq * f_N3LO_2l_Wc
        end
        for nth=1:7; gi[nth] += pjs[nth][n] * ws[n] * tmp_s;end
        
        # LS term: Vls
        if chi_order >= 2
            gi = gis[7]
            f_NNLO_Vls = 3.0 * gA4 / (32.0 * pi * Fpi4)
            f_N3LO_a_Vls = c2 *gA2 /(8.0* pi^2 * Fpi4)
            f_N3LO_b_Vls = gA4 /(4.0* pi^2 * Fpi4)
            tmp_s =  tw2Aq * f_NNLO_Vls
            if chi_order >= 3
               tmp_s += f_N3LO_a_Vls * w2 * Lq
               tmp_s += f_N3LO_b_Vls * Lq * (11.0/32.0 * q2 + nd_mpi4 / w2)
            end
            for nth=1:7; gi[nth] += pjs[nth][n] * ws[n] *tmp_s;end
        end
        # LS term: Wls   
        if chi_order >= 2
            gi =gis[8]
            f_NNLO_Wls = gA2 * (1.0- gA2) / (32.0 *pi * Fpi4)
            f_N3LO_a_Wls = -c4 /(48.0*pi^2 * Fpi4)
            f_N3LO_b_Wls = 1.0 /(256.0*pi^2 * Fpi4) 
            tmp_s = w2Aq * f_NNLO_Wls
            if chi_order >= 3
                tmp_s += f_N3LO_a_Wls*Lq * (gA2 *(8.0*nd_mpi2 + 5.0*q2) + w2)
                tmp_s += f_N3LO_b_Wls*Lq * (16.0*gA2 *(nd_mpi2 + 3.0/8.0 *q2)
                                        + 4.0/3.0 * gA4 *(4.0*nd_mpi4 /w2
                                                          -11.0/4.0 * q2
                                                          -9.0*nd_mpi2)
                                        -w2)
            end
           for nth=1:7; gi[nth] += pjs[nth][n] * ws[n] * tmp_s;end
        end
        # sigma-L term: Vsl        
        if chi_order >=3
            gi = gis[9]
            f_N3LO = gA4 /(32.0* pi^2 * Fpi4)
            tmp_s = f_N3LO * Lq
            for nth=1:7; gi[nth] += pjs[nth][n] * ws[n] *tmp_s;end
        end        
    end
    
    #Vt
    calc_IJ_V(J,pnrank,gis[1],f_T,fc,f_idx,tVs,lsj,tllsj,
              tdict,V12mom,V_i,V_j,to)
    #Wt
    calc_IJ_V(J,pnrank,gis[2],f_T,fc,f_idx,tVs,lsj,tllsj,
              tdict,V12mom,V_i,V_j,to;isodep=true)
    #Vs
    calc_IJ_V(J,pnrank,gis[3],f_SS,fc,f_idx,tVs,lsj,tllsj,
              tdict,V12mom,V_i,V_j,to)
    #Ws
    calc_IJ_V(J,pnrank,gis[4],f_SS,fc,f_idx,tVs,lsj,tllsj,
              tdict,V12mom,V_i,V_j,to;isodep=true)     
    #Vc
    calc_IJ_V(J,pnrank,gis[5],f_C,fc,f_idx,tVs,
              lsj,tllsj,tdict,V12mom,V_i,V_j,to)
    #Wc
    calc_IJ_V(J,pnrank,gis[6],f_C,fc,f_idx,tVs,
              lsj,tllsj,tdict,V12mom,V_i,V_j,to;isodep =true)
    #Vls
    calc_IJ_V(J,pnrank,gis[7],f_LS,fc,f_idx,tVs,lsj,tllsj,
              tdict,V12mom,V_i,V_j,to;addtype="ls")    
    #Wls
    calc_IJ_V(J,pnrank,gis[8],f_LS,fc,f_idx,tVs,lsj,tllsj,
              tdict,V12mom,V_i,V_j,to;isodep=true,addtype="ls")
    #Vsl
    calc_IJ_V(J,pnrank,gis[9],f_SL,fc,f_idx,tVs,lsj,tllsj,
              tdict,V12mom,V_i,V_j,to;addtype="sl")    
end

function calc_IJ_V(J,pnrank,gi,opf,fc,f_idx,tVs,lsj,tllsj,
                   tdict,V12mom,V_i,V_j,to;
                   isodep=false,addtype="")
    if J==0;gi[3]=gi[5]=gi[7]=0.0;end
    IJ0 = gi[1]; IJ1 = gi[2]
    IJ2 = (J* gi[2] + gi[3]) /(J+1) 
    IJ3 = sqrt(J/(J+1)) * (gi[2]-gi[3])

    V0  = opf[1] * IJ0 + opf[2] *IJ1
    V1  = opf[3] * IJ0 + opf[4] *IJ2
    V12 = opf[5] * IJ0 + opf[6] *IJ1
    V34 = opf[4] * IJ0 + opf[3] *IJ2
    V55 = opf[7] * IJ3
    V66 = opf[8] * IJ3

    if addtype == "ls"
        IJ4 = gi[4] 
        IJ5 = (J* gi[4] + gi[5]) /(J+1) 
        IJ6 = sqrt(J/(J+1)) * (gi[4]-gi[5])
        V0  += opf[9] * IJ4 
        V1  += opf[10] * IJ1 + opf[9] * IJ5 
        V12 += opf[10] * IJ4
        V34 += opf[10] * IJ5 +opf[9] * IJ1
        e1 = opf[11] * IJ6
        V55 += e1
        V66 += e1
    elseif addtype == "sl"
        IJ4 = gi[4] 
        IJ5 = (J* gi[4] + gi[5]) /(J+1) 
        IJ6 = sqrt(J/(J+1)) * (gi[4]-gi[5])
        IJ10 = gi[6]
        IJ11 = (J * gi[6] + gi[7]) / (J+1)
        IJ12 = sqrt(J/(J+1)) * (gi[6]-gi[7])
        V0  += opf[6] * IJ4 
        V1  += opf[1] * IJ4 + opf[9] * IJ5 
        V12 += opf[1] * IJ10
        V34 += opf[9] * IJ1 +opf[1] * IJ11
        e1 = opf[6] * IJ12
        V55 += e1
        V66 += e1
    end        
    transV_into_lsj(J,pnrank,tVs,
                    V0,V1,V12,V34,V55,V66;isodep=isodep)
    @inbounds for idx = 1:f_idx
        @views tllsj[2:5] .= lsj[idx]
        V12idx = get(tdict,tllsj,-1)
        if V12idx == -1;continue;end
        tfac = tVs[idx] * fc
        V12mom[V12idx][V_i,V_j] += tfac
    end
    return nothing
end

function transV_into_lsj(J,pnrank,Vs,v1,v2,v3,v4,v5,v6;isodep=false)
    ttis=1.0;ttit=1.0
    d2j1 = 1.0/(2*J+1)
    if J==0; v2=v4=v5=v6=0.0;end
    if J%2==1 && pnrank!=2; v1 =0.0;end
    if J%2==0 && pnrank!=2; v2 =0.0;end
    if J%2!=0 && pnrank!=2; v3=v4=v5=v6=0.0;end
    v34 = -sqrt(J*(J+1)) *(v3-v4)
    v56 = sqrt(J*(J+1)) * (v5+v6)
    if isodep
        is = J%2 + 1
        it = is%2 +1
        ttis = ifelse(is==2,-3.0,1.0)
        ttit = ifelse(it==2,-3.0,1.0)
    end
    Vs[1] = v1 * ttis
    Vs[2] = v2 * ttit
    Vs[3] = d2j1 * ((J+1)*v3 + J*v4-v56) * ttis
    Vs[4] = d2j1 * ( J*v3+(J+1)*v4 +v56) * ttis
    Vs[5] = -d2j1 * (v34-(J+1)*v5+J*v6) * ttis
    Vs[6] = -d2j1 * (v34+J*v5-(J+1)*v6) * ttis
    return nothing
end
