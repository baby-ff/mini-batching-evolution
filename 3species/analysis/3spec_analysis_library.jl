function hist2d(vecx::Vector{Float64},vecy::Vector{Float64},nbins::Int64)
    x_edges = range(minimum(vecx), maximum(vecx), length=nbins+1)
    y_edges = range(minimum(vecy), maximum(vecy), length=nbins+1)
    h = fit(Histogram, (vecx, vecy), (x_edges,y_edges))
    counts = h.weights./sum(h.weights)
    counts[findall(x->x==0,counts)] .= NaN
    midpoints(edges) = [0.5 * (edges[i] + edges[i+1]) for i in 1:length(edges)-1]
    x_mids = midpoints(x_edges)
    y_mids = midpoints(y_edges)
    
    return copy(x_mids), copy(y_mids), copy(counts)
end;

function Fit(x::Float64,y::Float64,avg_vec::Matrix{Float64})
    fa = avg_vec[1,1]
    fb = avg_vec[2,1]
    fc = avg_vec[3,1]
    return fc + (fa - fc)*x + (fb - fc)*y
end;

function Var(x::Float64,y::Float64,var_vec::Matrix{Float64})
    Va = var_vec[1,1]
    Vb = var_vec[2,1]
    Vc = var_vec[3,1]
    return Vc + (Va - Vc)*x + (Vb - Vc)*y
end;

function ff(x::Float64,y::Float64,avg_vec::Matrix{Float64})
    fa = avg_vec[1,1]
    fb = avg_vec[2,1]
    fc = avg_vec[3,1]
    return fc^2 + (fa^2 - fc^2)*x + (fb^2 - fc^2)*y
end;

function VV(x::Float64,y::Float64,var_vec::Matrix{Float64})
    Va = var_vec[1,1]
    Vb = var_vec[2,1]
    Vc = var_vec[3,1]
    return Vc^2 + (Va^2 - Vc^2)*x + (Vb^2 - Vc^2)*y
end; 

function fV(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    fa = avg_vec[1,1]
    fb = avg_vec[2,1]
    fc = avg_vec[3,1]
    Va = var_vec[1,1]
    Vb = var_vec[2,1]
    Vc = var_vec[3,1]
    return fc*Vc + (fa*Va - fc*Vc)*x + (fb*Vb - fc*Vc)*y
end;

function VVV(x::Float64,y::Float64,var_vec::Matrix{Float64})
    Va = var_vec[1,1]
    Vb = var_vec[2,1]
    Vc = var_vec[3,1]
    return Vc^3 + (Va^3 - Vc^3)*x + (Vb^3 - Vc^3)*y
end;

function ffV(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    fa = avg_vec[1,1]
    fb = avg_vec[2,1]
    fc = avg_vec[3,1]
    Va = var_vec[1,1]
    Vb = var_vec[2,1]
    Vc = var_vec[3,1]
    return fc^2*Va + (fa^2*Va - fc^2*Va)*x + (fb^2*Vb - fc^2*Va)*y
end;

function fVV(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    fa  = avg_vec[1,1]
    fb  = avg_vec[2,1]
    fc  = avg_vec[3,1]
    Va = var_vec[1,1]
    Vb = var_vec[2,1]
    Vc = var_vec[3,1]
    return fc*Vc*Vc + (fa*Va*Va - fc*Vc*Vc)*x + (fb*Vb*Vb - fc*Vc*Vc)*y
end;

function flz(x::Float64,y::Float64,avg_vec::Matrix{Float64})
    fa = avg_vec[1,1]
    fb = avg_vec[2,1]
    fc = avg_vec[3,1]
    return fa + fb + fc - 3*(fc^2 + (fb^2 - fc^2)*y + (fa^2 - fc^2)*x)/Fit(x,y,avg_vec)
end;

function Vlz(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    fa = avg_vec[1,1]
    fb = avg_vec[2,1]
    fc = avg_vec[3,1]
    Va = var_vec[1,1]
    Vb = var_vec[2,1]
    Vc = var_vec[3,1]
    return Va + Vb + Vc - 3*(fc*Vc + (fb*Vb - fc*Vc)*y + (fa*Va - fc*Vc)*x)/Fit(x,y,avg_vec)
end;

#########################

function Drx(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64},nu::Float64,N::Int64)
    return N*(ff(x,y,avg_vec) - Fit(x,y,avg_vec)^2)/Fit(x,y,avg_vec) + N*nu*flz(x,y,avg_vec) - (fV(x,y,avg_vec,var_vec) - Fit(x,y,avg_vec)*Var(x,y,var_vec))/Fit(x,y,avg_vec)^2
end;

function Dry(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64},nu::Float64,N::Int64)
    return N*(fV(x,y,avg_vec,var_vec) - Fit(x,y,avg_vec)*Var(x,y,var_vec))/Fit(x,y,avg_vec) + N*nu*Vlz(x,y,avg_vec,var_vec) - (VV(x,y,var_vec) - Var(x,y,var_vec)^2)/Fit(x,y,avg_vec)^2
end;

function Drx0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64},nu::Float64)
    return (ff(x,y,avg_vec) - Fit(x,y,avg_vec)^2)/Fit(x,y,avg_vec) + nu*flz(x,y,avg_vec)
end;

function Dry0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64},nu::Float64)
    return (fV(x,y,avg_vec,var_vec) - Fit(x,y,avg_vec)*Var(x,y,var_vec))/Fit(x,y,avg_vec) + nu*Vlz(x,y,avg_vec,var_vec) 
end;

function detA(avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    fa = avg_vec[1,1]
    fb = avg_vec[2,1]
    fc = avg_vec[3,1]
    Va = var_vec[1,1]
    Vb = var_vec[2,1]
    Vc = var_vec[3,1]
    return (fa - fc)*(Vb - Vc) - (fb - fc)*(Va - Vc)
end;

function xfunc(X::Float64, Y::Float64, avg_vec::Matrix{Float64}, var_vec::Matrix{Float64}) 
    fa = avg_vec[1,1]
    fb = avg_vec[2,1]
    fc = avg_vec[3,1]   
    Va = var_vec[1,1]
    Vb = var_vec[2,1]
    Vc = var_vec[3,1]
    return ((Vb - Vc)*X - (fb - fc)*Y)/detA(avg_vec,var_vec);
end;
    
function yfunc(X::Float64, Y::Float64, avg_vec::Matrix{Float64}, var_vec::Matrix{Float64})
    fa = avg_vec[1,1]
    fb = avg_vec[2,1]
    fc = avg_vec[3,1]   
    Va = var_vec[1,1]
    Vb = var_vec[2,1]
    Vc = var_vec[3,1]
    return (-(Va - Vc)*X + (fa - fc)*Y)/detA(avg_vec,var_vec);
end;

function drift_pt(pt::Point{2,Float32})
    X = Float64.(pt[1])
    Y = Float64.(pt[2])
    return Point2f(Drx(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec,nu,N),
                   Dry(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec,nu,N))
end;

###########################

function dff(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return (ff(x,y,avg_vec) - Fit(x,y,avg_vec)^2)*(1+Var(x,y,var_vec)/Fit(x,y,avg_vec)^2) + 
        (ffV(x,y,avg_vec,var_vec) - ff(x,y,avg_vec)*Var(x,y,var_vec) - 2*fV(x,y,avg_vec,var_vec)*Fit(x,y,avg_vec) + 
        2*Fit(x,y,avg_vec)^2*Var(x,y,var_vec))/Fit(x,y,avg_vec)^2
end;

function dvv(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return (VV(x,y,var_vec) - Var(x,y,var_vec)^2)*(1 + Var(x,y,var_vec)/Fit(x,y,avg_vec)^2) + 
        (VVV(x,y,var_vec) - 3*Var(x,y,var_vec)*VV(x,y,var_vec) + 2*Var(x,y,var_vec)^3)/Fit(x,y,avg_vec)^2
end;
    
function dfv(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return (fV(x,y,avg_vec,var_vec) - Fit(x,y,avg_vec)*Var(x,y,var_vec))*(1 - 
        Var(x,y,var_vec)/Fit(x,y,avg_vec)^2) + (fVV(x,y,avg_vec,var_vec) - VV(x,y,var_vec)*Fit(x,y,avg_vec))/(Fit(x,y,avg_vec)^2)
end;
function dff0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})  
    return ff(x,y,avg_vec) - Fit(x,y,avg_vec)^2
end;

function dvv0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})  
    return VV(x,y,var_vec) - Var(x,y,var_vec)^2
end;

function dfv0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})  
    return fV(x,y,avg_vec,var_vec) - Fit(x,y,avg_vec)*Var(x,y,var_vec)
end;

function Tr(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return dff(x,y,avg_vec,var_vec) + dvv(x,y,avg_vec,var_vec)
end;

function Det(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return dff(x,y,avg_vec,var_vec)*dvv(x,y,avg_vec,var_vec) - dfv(x,y,avg_vec,var_vec)^2
end;

function Tr0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return dff0(x,y,avg_vec,var_vec) + dvv0(x,y,avg_vec,var_vec)
end;

function Det0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return dff0(x,y,avg_vec,var_vec)*dvv0(x,y,avg_vec,var_vec) - dfv0(x,y,avg_vec,var_vec)^2
end;

function eigmax(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return 0.5*(Tr(x,y,avg_vec,var_vec) + sqrt(Tr(x,y,avg_vec,var_vec)^2 - 4*Det(x,y,avg_vec,var_vec)))
end;

function eigmin(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return 0.5*(Tr(x,y,avg_vec,var_vec) - sqrt(Tr(x,y,avg_vec,var_vec)^2 - 4*Det(x,y,avg_vec,var_vec)))
end;

function eigmax0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return 0.5*(Tr0(x,y,avg_vec,var_vec) + sqrt(Tr0(x,y,avg_vec,var_vec)^2 - 4*Det0(x,y,avg_vec,var_vec)))
end;

function eigmin0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return 0.5*(Tr0(x,y,avg_vec,var_vec) - sqrt(Tr0(x,y,avg_vec,var_vec)^2 - 4*Det0(x,y,avg_vec,var_vec)))
end;

function ufmax(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return abs(sqrt(eigmax(x,y,avg_vec,var_vec)^2/(1 + ((dff(x,y,avg_vec,var_vec) - eigmax(x,y,avg_vec,var_vec))/dfv(x,y,avg_vec,var_vec))^2)))
    #-sign((dff(x,y,avg_vec,var_vec) - eigmax(x,y,avg_vec,var_vec))/dfv(x,y,avg_vec,var_vec))*sqrt(eigmax(x,y,avg_vec,var_vec)^2/(1 + ((dff(x,y,avg_vec,var_vec) - eigmax(x,y,avg_vec,var_vec))/dfv(x,y,avg_vec,var_vec))^2));
end;

function ufmin(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return sign((dff(x,y,avg_vec,var_vec) - eigmin(x,y,avg_vec,var_vec))/dfv(x,y,avg_vec,var_vec))*sqrt(eigmin(x,y,avg_vec,var_vec)^2/(1 + ((dff(x,y,avg_vec,var_vec) - eigmin(x,y,avg_vec,var_vec))/dfv(x,y,avg_vec,var_vec))^2));
end;

function uvmax(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return ufmax(x,y,avg_vec,var_vec)*(dff(x,y,avg_vec,var_vec) - eigmax(x,y,avg_vec,var_vec))/dfv(x,y,avg_vec,var_vec) 
end;

function uvmin(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return ufmin(x,y,avg_vec,var_vec)*(dff(x,y,avg_vec,var_vec) - eigmin(x,y,avg_vec,var_vec))/dfv(x,y,avg_vec,var_vec)
end;

function ufmax0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return abs(sqrt(eigmax0(x,y,avg_vec,var_vec)^2/(1 + ((dff0(x,y,avg_vec,var_vec) - eigmax0(x,y,avg_vec,var_vec))/dfv0(x,y,avg_vec,var_vec))^2)))
end;

function ufmin0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return sign((dff0(x,y,avg_vec,var_vec) - eigmin0(x,y,avg_vec,var_vec))/dfv0(x,y,avg_vec,var_vec))*sqrt(eigmin0(x,y,avg_vec,var_vec)^2/(1 + ((dff0(x,y,avg_vec,var_vec) - eigmin0(x,y,avg_vec,var_vec))/dfv0(x,y,avg_vec,var_vec))^2));
end;

function uvmax0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return ufmax0(x,y,avg_vec,var_vec)*(dff0(x,y,avg_vec,var_vec) - eigmax0(x,y,avg_vec,var_vec))/dfv0(x,y,avg_vec,var_vec)
end;

function uvmin0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return ufmin0(x,y,avg_vec,var_vec)*(dff0(x,y,avg_vec,var_vec) - eigmin0(x,y,avg_vec,var_vec))/dfv0(x,y,avg_vec,var_vec)
end;

function diffmax_pt(pt::Point{2,Float32})
    X = Float64.(pt[1])
    Y = Float64.(pt[2])
    return Point2f(ufmax(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec),
            uvmax(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec))
end;

function diffmin_pt(pt::Point{2,Float32})
    X = Float64.(pt[1])
    Y = Float64.(pt[2])
    return Point2f(ufmin(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec),
            uvmin(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec))
end;

function diffmax0_pt(pt::Point{2,Float32})
    X = Float64.(pt[1])
    Y = Float64.(pt[2])
    return Point2f(ufmax0(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec),
                    uvmax0(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec))
end;

function diffmin0_pt(pt::Point{2,Float32})
    X = Float64.(pt[1])
    Y = Float64.(pt[2])
    return Point2f(ufmin0(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec),
        uvmin0(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec))
end;

function fixed(avg_vec::Matrix{Float64},var_vec::Matrix{Float64},nu::Float64,N::Int64)
    function drift_eq!(res,pt)
        X = Float64.(pt[1])
        Y = Float64.(pt[2])
        res[1] = Drx(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec,nu,N)
        res[2] = Dry(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec,nu,N)
    end
    return nlsolve(drift_eq!,[0.4,0.],ftol=1e-16,method=:newton)
end;

function fixed0(avg_vec::Matrix{Float64},var_vec::Matrix{Float64},nu::Float64)
    function drift_eq!(res,pt)
        X = Float64.(pt[1])
        Y = Float64.(pt[2])
        res[1] = Drx0(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec,nu)
        res[2] = Dry0(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec,nu)
    end
    return nlsolve(drift_eq!,[0.4,0.],ftol=1e-16,method=:newton)
end;

function norm_diff(avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    norm_drift = sqrt(Drx(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec,nu,N)^2 + 
                Dry(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec,nu,N)^2)
    norm_diff = sqrt(eigmax(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec)^2 + 
                    eigmin(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec)^2)
    return norm_drift - norm_diff
end;