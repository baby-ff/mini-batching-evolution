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
    
function Drx(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64},nu::Float64,N::Int64)
    return N*(ff(x,y,avg_vec) - Fit(x,y,avg_vec)^2)/Fit(x,y,avg_vec) + N*nu*flz(x,y,avg_vec) - (fV(x,y,avg_vec,var_vec) - Fit(x,y,avg_vec)*Var(x,y,var_vec))/Fit(x,y,avg_vec)^2
end;
    
function Dry(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64},nu::Float64,N::Int64)
    return N*(fV(x,y,avg_vec,var_vec) - Fit(x,y,avg_vec)*Var(x,y,var_vec))/Fit(x,y,avg_vec) + N*nu*Vlz(x,y,avg_vec,var_vec) - (VV(x,y,var_vec) - Var(x,y,var_vec)^2)/Fit(x,y,avg_vec)^2
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
    
function drift(X::Float64,Y::Float64)
    return Point2f(Drx(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec,nu,N),
                Dry(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec,nu,N))
end;
    
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
    return  0.5*(Tr0(x, y, avg_vec, var_vec) + sqrt(Tr0(x, y, avg_vec, var_vec)^2 - 4*Det0(x, y, avg_vec, var_vec)))
end;

function eigmin0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return  0.5*(Tr0(x, y, avg_vec, var_vec) - sqrt(Tr0(x, y, avg_vec, var_vec)^2 - 4*Det0(x, y, avg_vec, var_vec)))
end;

function ufmax(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return sign((dff(x,y,avg_vec,var_vec) - eigmax(x,y,avg_vec,var_vec))/dfv(x,y,avg_vec,var_vec))*sqrt(eigmax(x,y,avg_vec,var_vec)^2/(1 + ((dff(x,y,avg_vec,var_vec) - eigmax(x,y,avg_vec,var_vec))/dfv(x,y,avg_vec,var_vec))^2));
end;
    
function ufmin(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return sign((dff(x,y,avg_vec,var_vec) - eigmin(x,y,avg_vec,var_vec))/dfv(x,y,avg_vec,var_vec))*sqrt(eigmin(x,y,avg_vec,var_vec)^2/(1 + ((dff(x,y,avg_vec,var_vec) - eigmin(x,y,avg_vec,var_vec))/dfv(x,y,avg_vec,var_vec))^2));
end;
    
function uvmax(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    res = ufmax(x,y,avg_vec,var_vec)*(dff(x,y,avg_vec,var_vec) - eigmax(x,y,avg_vec,var_vec))/dfv(x,y,avg_vec,var_vec)
    return ufmax(x,y,avg_vec,var_vec)*(dff(x,y,avg_vec,var_vec) - eigmax(x,y,avg_vec,var_vec))/dfv(x,y,avg_vec,var_vec) 
end;
    
function uvmin(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return ufmin(x,y,avg_vec,var_vec)*(dff(x,y,avg_vec,var_vec) - eigmin(x,y,avg_vec,var_vec))/dfv(x,y,avg_vec,var_vec)
end;
    
function ufmax0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    sign((dff0(x, y, avg_vec, var_vec) - eigmax0(x, y, avg_vec, var_vec))/dfv0(x, y, avg_vec, var_vec))*sqrt(eigmax0(x,y, avg_vec, var_vec)^2/(1 + ((dff0(x, y, avg_vec, var_vec) - eigmax0(x, y, avg_vec, var_vec))/dfv0(x, y, avg_vec, var_vec))^2));
end;

function uvmax0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return ufmax0(x, y, avg_vec, var_vec)*(dff0(x, y, avg_vec, var_vec) - eigmax0(x, y, avg_vec, var_vec))/dfv0(x, y, avg_vec ,var_vec)
end

function ufmin0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return sign((dff0(x, y, avg_vec, var_vec) - eigmin0(x, y, avg_vec, var_vec))/dfv0(x, y, avg_vec, var_vec))*sqrt(eigmin0(x,y, avg_vec, var_vec)^2/(1 + ((dff0(x, y, avg_vec, var_vec) - eigmin0(x, y, avg_vec, var_vec))/dfv0(x, y, avg_vec, var_vec))^2));
end;

function uvmin0(x::Float64,y::Float64,avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    return ufmin0(x, y, avg_vec, var_vec)*(dff0(x, y, avg_vec, var_vec) - eigmin0(x, y, avg_vec, var_vec))/dfv0(x, y, avg_vec, var_vec)
end;

function diffmax_pt(pt::Point{2,Float32},avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    X = Float64.(pt[1])
    Y = Float64.(pt[2])
    return Point2f(ufmax(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec),
                uvmax(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec))
end;
    
function diffmin_pt(pt::Point{2,Float32},avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    X = Float64.(pt[1])
    Y = Float64.(pt[2])
    return Point2f(ufmin(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec),
                    uvmin(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec))
end;

function diffmax0_pt(pt::Point{2,Float32},avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    X = Float64.(pt[1])
    Y = Float64.(pt[2])
    return Point2f(ufmax0(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec),
                    uvmax0(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec))
end;

function diffmin0_pt(pt::Point{2,Float32},avg_vec::Matrix{Float64},var_vec::Matrix{Float64})
    X = Float64.(pt[1])
    Y = Float64.(pt[2])
    return Point2f(ufmin0(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec),
                    uvmin0(xfunc(X,Y,avg_vec,var_vec),yfunc(X,Y,avg_vec,var_vec),avg_vec,var_vec))
end;

function plot_land(land_vec::Matrix{Float64})
    ls = 0.07
    xstep = 0.05
    f = Figure(backgroundcolor = RGBf(1,1,1), size = (800, 600), fontsize = 20)
    gabc = f[1, 1] = GridLayout()
    gde = f[2, 1] = GridLayout()
    avg_vec = copy(transpose(mean(land_vec, dims=1)))
    var_vec = copy(transpose(var(land_vec, dims=1)))
    axh = [Axis(gabc[1, i], xlabel = "f", xgridvisible = false, ygridvisible = false) for i in 1:3]
    ad0 = Axis(gde[1, 1], xlabel = "f",ylabel = "V", xgridvisible = false, ygridvisible = false)
    adi = Axis(gde[1, 2], xlabel = "f",ylabel = "V", xgridvisible = false, ygridvisible = false)
    max = maximum(land_vec)
    hist!(axh[1],land_vec[:,1],normalization=:pdf)
    hist!(axh[2],land_vec[:,2],normalization=:pdf)
    hist!(axh[3],land_vec[:,3],normalization=:pdf)
    for k in 1:3
        xlims!(axh[k],0,max) 
        ylims!(axh[k],0,1.15) 
    end

    slope1 = (var_vec[2,1]-var_vec[3,1])/(avg_vec[2,1]-avg_vec[3,1])
    slope2 = (var_vec[1,1]-var_vec[3,1])/(avg_vec[1,1]-avg_vec[3,1])
    
    scatter!(adi,mean(land_vec,dims=1)[1,:].-avg_vec[3,1],var(land_vec,dims=1)[1,:].-var_vec[3,1],markersize=4,color=:black)
    lines!(adi,push!(avg_vec[:,1],avg_vec[1,1]).-avg_vec[3,1],push!(var_vec[:,1],var_vec[1,1]).-var_vec[3,1],color=:black,linewidth=0.5)
    
    scatter!(ad0,mean(land_vec,dims=1)[1,:].-avg_vec[3,1],var(land_vec,dims=1)[1,:].-var_vec[3,1],markersize=4,color=:black)
    lines!(ad0,push!(avg_vec[:,1],avg_vec[1,1]).-avg_vec[3,1],push!(var_vec[:,1],var_vec[1,1]).-var_vec[3,1],color=:black,linewidth=0.5)
    
    points_diff = reduce(vcat,[[Point2f(x, y) for y in range((slope1+0.02)*x,(slope2-0.02)*x,minimum([12,Int(floor((slope2-slope1)*x/0.001))]))] for x in 0:xstep:maximum(avg_vec)-avg_vec[3,1]])
    spmax = [diffmax_pt(x,avg_vec,var_vec) for x in points_diff]
    strength_max = [norm(spmax[i]) for i in 1:length(spmax)]
    spmin = [diffmin_pt(x,avg_vec,var_vec) for x in points_diff]
    strength_min = [norm(spmin[i]) for i in 1:length(spmin)]
    
    spmax0 = [diffmax0_pt(x,avg_vec,var_vec) for x in points_diff]
    strength_max0 = [norm(spmax0[i]) for i in 1:length(spmax)]
    spmin0 = [diffmin0_pt(x,avg_vec,var_vec) for x in points_diff]
    strength_min0 = [norm(spmin0[i]) for i in 1:length(spmin)]
    
    cmin = minimum([minimum(strength_min),minimum(strength_min0)])
    cmax = maximum([maximum(strength_max),maximum(strength_max0)])

    arrows!(adi,points_diff,spmax, arrowsize = 0, lengthscale = ls, arrowcolor = strength_max, linecolor = strength_max, 
        colormap = cgrad(:plasma),linewidth=2, align=:center, colorrange = (cmin, cmax))
    arrows!(adi,points_diff,spmin, arrowsize = 0, lengthscale = ls, arrowcolor = strength_min, linecolor = strength_min, 
        colormap = cgrad(:plasma),linewidth=1.5, align = :center, colorrange = (cmin, cmax))
    
    arrows!(ad0,points_diff,spmax0, arrowsize = 0, lengthscale = ls, arrowcolor = strength_max0, linecolor = strength_max0, 
        colormap = cgrad(:plasma),linewidth=2, align=:center, colorrange = (cmin, cmax))
    arrows!(ad0,points_diff,spmin0, arrowsize = 0, lengthscale = ls, arrowcolor = strength_min0, linecolor = strength_min0, 
            colormap = cgrad(:plasma),linewidth=1.5, align = :center, colorrange = (cmin, cmax))

    return f
end;