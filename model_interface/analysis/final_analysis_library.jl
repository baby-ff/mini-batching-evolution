function bin_to_int(str::Vector{Int})
#converts bin sequence to int
    str = reverse(join(string.(str)))
    x = parse(Int, str; base=2)
    return x
end;

function getsnap(filename::String)
    #reads snapshot data
        
    mat = readdlm(filename,header=true)[1]
    abund = mat[:,end]
    gen = mat[:,end-1]
    fit = mat[:,end-2]
    seq = mat[:,1]
    si = sortperm(fit)
    return (fit[si],cumsum(abund[si]),gen[si],abund[si],seq[si])
end;

function prf(folder::String,snaps::String,tini::Int,nsteps::Int,ncopies::Int,dt::Int,V::Vector{Any})
    popV = []
    popX = []
    popY = []
    
    inits = append!(["zeros/", "ones/"],["rand$(i)/" for i in 1:5])
    for init_cond in inits
        initfolder = folder*init_cond
        for n in 1:ncopies
            sfolder = initfolder*"/sample_$(n)/"
            for t in tini:nsteps
                x,_,y,ab,seq = getsnap(sfolder*snaps*"/snap_t$(t*dt).txt")
                v = [V[s+1] for s in seq]
                w = weights(Int.(ab))
                mpv = mean(v,w)
                mpx = mean(x,w)
                mpy = mean(y,w)
                push!(popV,mpv)
                push!(popX,mpx)
                push!(popY,mpy)
            end
        end
    end

    popX = round.(popX,digits=12)
    popV = round.(popV,digits=12)
    
    return popX,popV
end;


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

																									
function extract_potential(Fvec::Vector{Float64},Dvec::Vector{Float64},Dvec0::Vector{Float64},Tvec::Vector{Float64},Tvec0::Vector{Float64},adaptive_bin::Int64)
    perm = sortperm(Fvec)
    D = []
    T = []
    T0 = []
    Fx = []
    D0 =[ ]
    for k = 1:Int(size(Fvec)[1]/adaptive_bin)-1
        idxs = perm[k*adaptive_bin+1:(k+1)*adaptive_bin]
        push!(D,mean(Dvec[idxs]))
        push!(T,mean(Tvec[idxs]))
        push!(T0,mean(Tvec0[idxs]))
        push!(Fx,mean(Fvec[idxs]))
        push!(D0,mean(Dvec0[idxs]))
    end

    midFx = 0.5*(Fx[1:end-1]+Fx[2:end])
    dFx = append!([midFx[1]-Fx[1]],push!(midFx[2:end]-midFx[1:end-1],Fx[end]-midFx[end]))

    Pot = -(cumsum(dFx.*D)) 
    return Float64.(D),Float64.(D0),Float64.(T),Float64.(T0),Float64.(Fx),Float64.(Pot)
end;

function get_FVG(means::Matrix{Float64},cij::Matrix{Float64},L::Int64,K::Int64)
    test = reshape(append!(ones(K),zeros(L-K)),(1,L))
    F = []
    V = []
    G = []
    for spec in 0:2^L-1
        seq = digits(spec,base=2,pad=L)
        norm = sum(seq)
        if norm==0
            push!(G,K/L)
            push!(F,mean(means))
            push!(V,ones(L)'*cij*ones(L)./(L^2))
        else
            var = (seq'*cij*seq)./(norm^2)
            gen = (test*seq)[1]/norm
            fit = (means*seq)[1]/norm
            push!(V,var)
            push!(G,gen)
            push!(F,fit)
        end
    end
    
    return F,V,G
end;

function get_muller_tot(initfolder::String,snaps::String,ncopies::Int,L::Int,nsteps::Int,tini::Int,dt::Int)
    Ztr = mullerdata(initfolder*"sample_1/"*snaps,L,nsteps,tini,dt)
    for n in 2:ncopies
        Ztr = hcat(Ztr,mullerdata(initfolder*"sample_$(n)/"*snaps,L,nsteps,tini,dt))
    end
    return Ztr
end;

function mullerdata(f::String,L::Int,nsteps::Int,tini::Int,dt::Int)
#reads pop data from snaps
    
    pop = zeros(2^L,nsteps-tini+1)
    file_paths = [f * "snap_t$(t*dt).txt" for t in tini:nsteps]
    
    for (i, path) in enumerate(file_paths)
        _, _, _, ab, seq = getsnap(path)
        uni, ab_uni = unify(ab, seq)
        pop[uni .+ 1, i] = ab_uni
    end
    return pop
end;

function unify(ab::Vector{Any},seq::Vector{Any})
#unique sequence-ab occurrence
    
    uni = unique(seq)
    ab_uni  = []
    for u in uni
        idx = findall(x->x==u,seq)
        push!(ab_uni,sum(ab[idx]))
    end
    return (uni, ab_uni)
end;

function adjacency_hypercube(L::Int)
#adjacency matrix of hypercube 
    
    A = zeros(2^L,2^L)
    for spec in 0:2^L-1
        seq = digits(spec,base=2,pad=L)
        onehot = Matrix(I, L, L)
        neighs = [(seq .+ onehot[:,k]).%2 for k in 1:L]
        ngh = [bin_to_int(neighs[k]) for k in 1:L]
        A[spec+1,ngh.+1] = ones(L)
    end
    return A
end;

function lap_hypercube(L::Int)
#laplacian of hypercube
    
    lap = -L.*Matrix(I, 2^L, 2^L) + adjacency_hypercube(L)
    return lap
end;		

function get_obs_track(sfolder::String)
    bWF = readdlm(sfolder*"track_evo.txt",header=true)[1];
    mWF = readdlm(sfolder*"track_evo_ref.txt",header=true)[1];
return bWF[:,3],bWF[:,5],mWF[:,3],mWF[:,5]
end
    
function get_means_distr(initfolder::String,ncopies::Int,tini::Int)
    
    avg1,avg2,avg3,avg4 = get_obs_track(initfolder*"sample_1/")
    avg1_vec = copy(avg1[tini:end])
    avg2_vec = copy(avg2[tini:end])
    avg3_vec = copy(avg3[tini:end])
    avg4_vec = copy(avg4[tini:end])
    for n in 2:ncopies
        f,g,fref,gref = get_obs_track(initfolder*"sample_$(n)/")
        avg1 .+= f
        avg2 .+= g
        avg3 .+= fref
        avg4 .+= gref
        append!(avg1_vec,f[tini:end])
        append!(avg2_vec,g[tini:end])
        append!(avg3_vec,fref[tini:end])
        append!(avg4_vec,gref[tini:end])
    end
    return (mean(avg1[tini:end])/ncopies, mean(avg2[tini:end])/ncopies, mean(avg3[tini:end])/ncopies, mean(avg4[tini:end])/ncopies, 
        avg1_vec, avg2_vec, avg3_vec, avg4_vec) 
end;

function Nproj_driff(xtr::Matrix{Float64}, reF::Matrix{Float64}, reV::Matrix{Float64}, nu::Float64, 
        lap::Matrix{Float64}, Narr::Vector{Int}, v1::Matrix{Float64})
    tmax = size(xtr)[2]
    F_tr = reshape(reF*xtr,(tmax,1))
    V_tr = reshape(reV*xtr,(tmax,1))
    
    v1_tr = reshape(v1*xtr,(tmax,1))
    Fv1_tr = reshape((reF.*v1)*xtr,(tmax,1))
    Vv1_tr = reshape((reV.*v1)*xtr,(tmax,1))
    v12_tr = reshape((v1.^2)*xtr,(tmax,1))
    v12V_tr = reshape(((v1.^2).*reV)*xtr,(tmax,1))
    
    fvec = reduce(vcat,reF)
    Lvec = lap*(fvec.*xtr)
    
    dr_pc1 = Narr[1].*(Fv1_tr .- F_tr.*v1_tr .+ nu.* transpose(v1*Lvec))./(F_tr) .- (Vv1_tr .- V_tr.*v1_tr )./(F_tr.^2)

    for N in Narr[2:end]
        dr_pc1 = hcat(dr_pc1,N.*(Fv1_tr .- F_tr.*v1_tr .+ nu.* transpose(v1*Lvec))./(F_tr) .- (Vv1_tr .- V_tr.*v1_tr )./(F_tr.^2))
    end
    
    diff_pc1_0 = v12_tr .- (v1_tr.^2)
    diff_pc1 = diff_pc1_0.*(1 .+V_tr./(F_tr.^2)) .+  (v12V_tr .- v12_tr.*V_tr .- 2 .*(v1_tr.*Vv1_tr .- v12_tr.*V_tr))./(F_tr.^2)
    
    return dr_pc1, diff_pc1[:,1], diff_pc1_0[:,1]
end;

function proj_drift_diff(xtr::Matrix{Float64}, reF::Matrix{Float64}, reV::Matrix{Float64}, reG::Matrix{Float64}, nu::Float64, 
    lap::Matrix{Float64}, N::Int, v1::Matrix{Float64})
    tmax = size(xtr)[2]
    F_tr = reshape(reF*xtr,(tmax,1))
    V_tr = reshape(reV*xtr,(tmax,1))
    G_tr = reshape(reG*xtr,(tmax,1))

    V2_tr = reshape((reV.^2)*xtr,(tmax,1))
    F2_tr = reshape((reF.^2)*xtr,(tmax,1))
    FV_tr = reshape((reF.*reV)*xtr,(tmax,1))

    F2V_tr = reshape(((reF.^2).*reV)*xtr,(tmax,1))
    V3_tr = reshape((reV.^3)*xtr,(tmax,1))

    v1_tr = reshape(v1 * xtr, (tmax, 1))
    Fv1_tr = reshape((reF .* v1) * xtr, (tmax, 1))
    Vv1_tr = reshape((reV .* v1) * xtr, (tmax, 1))
    v12_tr = reshape((v1 .^ 2) * xtr, (tmax, 1))
    v12V_tr = reshape(((v1 .^ 2) .* reV) * xtr, (tmax, 1))

    Lvec = lap * (reduce(vcat, reF) .* xtr)

    dr_F0 = N.*(F2_tr .- (F_tr.^2) .+ nu.* transpose(reF*Lvec))./(F_tr)
    dr_F = dr_F0 - (FV_tr .- F_tr.*V_tr)./(F_tr.^2)

    dr_V0 = N.*(FV_tr .- F_tr.*V_tr .+ nu.* transpose(reV*Lvec))./(F_tr)
    dr_V = dr_V0 - (V2_tr .- V_tr.^2)./(F_tr.^2)

    dr_pc1_0 = N.*(Fv1_tr .- F_tr.*v1_tr .+ nu.* transpose(v1*Lvec))./(F_tr) 
    dr_pc1 = copy(dr_pc1_0) .- (Vv1_tr .- V_tr.*v1_tr )./(F_tr.^2)

    diff_F0 = F2_tr .- (F_tr.^2)
    diff_F = diff_F0.*(1 .+V_tr./(F_tr.^2)) .+ (F2V_tr .- F2_tr.*V_tr .- 2 .*F_tr.*(FV_tr .- F_tr.*V_tr))./(F_tr.^2)

    diff_V0 = V2_tr .- (V_tr.^2)
    diff_V = diff_V0.*(1 .+V_tr./(F_tr.^2)) .+ (V3_tr .- V_tr.^3 .- 3 .*V_tr.*(V2_tr .- V_tr.^2))./(F_tr.^2)

    diff_pc1_0 = v12_tr .- (v1_tr.^2)
    diff_pc1 = diff_pc1_0.*(1 .+V_tr./(F_tr.^2)) .+  (v12V_tr .- v12_tr.*V_tr .- 2 .*(v1_tr.*Vv1_tr .- v12_tr.*V_tr))./(F_tr.^2)

    return dr_F[:,1], dr_V[:,1], dr_pc1[:,1], dr_pc1_0[:,1], diff_F0[:,1], diff_F[:,1], diff_V0[:,1], diff_V[:,1], diff_pc1[:,1], diff_pc1_0[:,1], F_tr[:,1], V_tr[:,1], G_tr[:,1]
end

function vector_field(xtr::Matrix{Float64}, reF::Matrix{Float64}, reV::Matrix{Float64}, nu::Float64, 
        lap::Matrix{Float64}, N::Int, v1::Matrix{Float64}, v2::Matrix{Float64}, v3::Matrix{Float64})
    tmax = size(xtr)[2]
    F_tr = reshape(reF*xtr,(tmax,1))
    V_tr = reshape(reV*xtr,(tmax,1))
    
    v1_tr = reshape(v1*xtr,(tmax,1))
    Fv1_tr = reshape((reF.*v1)*xtr,(tmax,1))
    Vv1_tr = reshape((reV.*v1)*xtr,(tmax,1))
    
    v2_tr = reshape(v2*xtr,(tmax,1))
    Fv2_tr = reshape((reF.*v2)*xtr,(tmax,1))
    Vv2_tr = reshape((reV.*v2)*xtr,(tmax,1))
    
    v3_tr = reshape(v3*xtr,(tmax,1))
    Fv3_tr = reshape((reF.*v3)*xtr,(tmax,1))
    Vv3_tr = reshape((reV.*v3)*xtr,(tmax,1))
    
    fvec = reduce(vcat,reF)
    Lvec = lap*(fvec.*xtr)

    dr_pc1 = N.*(Fv1_tr .- F_tr.*v1_tr .+ nu.* transpose(v1*Lvec))./(F_tr) .- (Vv1_tr .- V_tr.*v1_tr )./(F_tr.^2)
    dr_pc2 = N.*(Fv2_tr .- F_tr.*v2_tr .+ nu.* transpose(v2*Lvec))./(F_tr) .- (Vv2_tr .- V_tr.*v2_tr )./(F_tr.^2)
    dr_pc3 = N.*(Fv3_tr .- F_tr.*v3_tr .+ nu.* transpose(v3*Lvec))./(F_tr) .- (Vv3_tr .- V_tr.*v3_tr )./(F_tr.^2)
    
    return dr_pc1, dr_pc2, dr_pc3
end

function get_ipr_traj(muller_traj::Matrix{Float64})
    return sum(muller_traj.^2,dims=1)[1,:]
end

function get_genbins(K::Int,L::Int,nbins::Int,consmatrix::Matrix{Float64},tol::Float64)
    G = []
    for spec in 0:2^L-1
        seq = digits(spec,base=2,pad=L)
        norm = sum(seq)
        if norm==0
            gen = K/L
        else
            gen = (consmatrix*seq)[1]/norm
        end
        push!(G,gen)
    end
    unigen = unique(G)
    
    unigen = Float64.(unigen)
    
    dist = LowerTriangular(pairwise(Euclidean(), unigen, unigen))
    to_merge = intersect(findall(x->x<tol,dist),findall(x->x>0,dist))
    for i in reverse(1:length(to_merge))
        deleteat!(unigen,to_merge[i][1])
    end
    
    if length(unigen)<nbins
        println("nbins greater than number of distinct gen values. Reset to $(length(unigen)).")
        nbins = length(unigen)
    end
    return sort(unigen,rev=true)[1:nbins] 
end;

function getbin(genvalue::Float64,genbins::Vector{Float64},tol::Float64)
    bool = abs.(genvalue.-genbins).<tol
    if sum(bool)==0
        println(genvalue)
        println(genbins)
        idx = length(genbins)+1
    else 
        idx = findfirst(x->x==1,bool)
    end
    return idx
end;

function mullerplot_gen(folder::String,snaps::String,tini::Int,ncopies::Int,nsteps::Int,N::Int,dt::Int,
    genbins::Vector{Float64},consmatrix::Matrix{Float64},nbins::Int64)

    tol = 10^-16

    inits = append!(["zeros/", "ones/"],["rand$(i)/" for i in 1:5])

    gentable = zeros(nsteps-tini+1,minimum([size(genbins)[1],nbins+1]))

    for  init_cond in inits
        initfolder = folder*init_cond
        for n in 1:ncopies
            sfolder = initfolder*"/sample_$(n)/"
            for t in tini:nsteps
                x,_,y,ab,seq = getsnap(sfolder*snaps*"/snap_t$(t*dt).txt")
                k = 1
                for s in seq
                    d = digits(s,base=2,pad=L)
                    norm = sum(d)
                    if norm==0
                        gen = K./L
                    else 
                        gen = (consmatrix*d)[1]/norm
                    end
                    gentable[t-tini+1,getbin(gen,genbins,tol)] += ab[k]
                    k+=1
                end
            end
        end
    end
    return gentable./(N*ncopies*length(inits))
end;

function get_dr2(reF::Matrix{Float64}, reV::Matrix{Float64}, xtr::Matrix{Float64}, lap::Matrix{Float64}, N::Int, nu::Float64)
    Ftr = (reF*xtr)[1,:]
    Vtr = (reV*xtr)[1,:]
    Lvec = transpose(lap * (reduce(vcat, reF) .* xtr))
    drift = ((N.*(reF.-Ftr) .-(reV.-Vtr))./Ftr.*transpose(xtr) .+ (N*nu).*Lvec)./Ftr
    drift2 = sum(drift.^2,dims=2)
    return drift2
end

function get_flatness(F::Vector{Float64}, L::Int)
    adj = lap_hypercube(L).+(L.*Matrix(I, 2^L, 2^L))
    invflat = sum(abs.(F.-transpose(adj.*F)).*adj,dims=2)
return invflat
end

function update_abs(ab_vec::Vector{Float64},ab_ref_vec::Vector{Float64},uni::Vector{Float64},startlist::Vector{Float64},
    xm::Matrix{Float64},zm::Matrix{Float64})
    i = 1
    for u in uni
        ab_vec[i] += sum(xm[findall(x -> x== u,startlist)])
        ab_ref_vec[i] += sum(zm[findall(x -> x== u,startlist)])
        i+=1
    end

    return ab_vec, ab_ref_vec
end

function plotcum(ax_tmp,absorted::Vector{Float64},ref_absorted::Vector{Float64},colsort::Vector{Float64})
    band!(ax_tmp[1],[0,1],zeros(2),ones(2).*absorted[1],color=colsort[1],colorrange=(colsort[1],colsort[end]),colormap=:magma)
    band!(ax_tmp[2],[0,1],zeros(2),ones(2).*ref_absorted[1],color=colsort[1],colorrange=(colsort[1],colsort[end]),colormap=:magma)
    for i in 2:length(absorted)
        band!(ax_tmp[1],[0,1],ones(2).*absorted[i-1],ones(2).*absorted[i].+0.01,color=colsort[i],colorrange=(colsort[1],colsort[end]),colormap=:magma,overdraw=true)
        band!(ax_tmp[2],[0,1],ones(2).*ref_absorted[i-1],ones(2).*ref_absorted[i].+0.01,color=colsort[i],colorrange=(colsort[1],colsort[end]),colormap=:magma,overdraw=true)
    end
    [ylims!(ax_tmp[k],0,1) for k in 1:2]
    [xlims!(ax_tmp[k],0,1) for k in 1:2]
    return 
end

function panelF(folder::String,F,G,N,ncopies,L,nsteps,tini,dt)
    init="zeros/"

    Xtr = get_muller_tot(folder*init,"saved_snaps/",ncopies,L,nsteps,tini,dt)./N
    Ztr = get_muller_tot(folder*init,"saved_snaps_ref/",ncopies,L,nsteps,tini,dt)./N

    invflat = get_flatness(Float64.(F),L)[:,1]

    xm = mean(Xtr,dims=2) 
    zm = mean(Ztr,dims=2) 

    uni_list = unique(invflat)

    gen_uni = unique(G)


    ab = zeros(length(uni_list))
    ab_ref = zeros(length(uni_list))

    genab = zeros(length(gen_uni))
    genab_ref = zeros(length(gen_uni))

    i = 1
    for u in uni_list
        ab[i] += sum(xm[findall(x -> x==u,invflat)])
        ab_ref[i] += sum(zm[findall(x -> x== u,invflat)])
        i+=1
    end

    i = 1
    for u in gen_uni
        genab[i] += sum(xm[findall(x -> x== u,G)])
        genab_ref[i] += sum(zm[findall(x -> x== u,G)])
        i+=1
    end


    for init in append!(["ones/"],["rand$(i)/" for i in 1:5])
        Xtr = get_muller_tot(folder*init,"saved_snaps/",ncopies,L,nsteps,tini,dt)./N
        Ztr = get_muller_tot(folder*init,"saved_snaps_ref/",ncopies,L,nsteps,tini,dt)./N
        xm = mean(Xtr,dims=2) 
        zm = mean(Ztr,dims=2) 
        i=1
        for u in uni_list
            ab[i] += sum(xm[findall(x -> x== u,invflat)])
            ab_ref[i] += sum(zm[findall(x -> x== u,invflat)])
            i+=1
        end
        i=1
        for u in gen_uni
            genab[i] += sum(xm[findall(x -> x== u,G)])
            genab_ref[i] += sum(zm[findall(x -> x== u,G)])
            i+=1
        end
    end

    return ab, ab_ref, genab, genab_ref, uni_list, gen_uni
end;

function get_data(w::Float64,Nmax::Int64,V0::Float64)
    inside_fp = Float64[]
    inside_ext = Float64[]
    Narr_fp = Float64[]
    Narr_ext = Float64[]
    for N in range(2,Nmax)
        #fixed points
        c0 = N*s - w 
        c1 = N*s^2 - 8*N*nu 
        c2 = -4*N*s + 4*w - 8*N*s*nu
        c3 = -4*N*s^2

        R = roots([c0,c1,c2,c3])
        r = real.(R[findall(x->isreal(x),R)])
        rapp = r[findall(x->x^2<=1/4,r)]
        append!(inside_fp,rapp)
        append!(Narr_fp,N*ones(length(rapp)))

        #extrema
        q0 = N*s - s*V0 - w
        q1 = 4 + 2*N*s^2 + 4*V0 - (-2 + s)*w - 8*N*nu
        q2 = -4*N*s + N*s^3 + 4*(s*(3 + 2*V0) + w) - 16*N*s*nu
        q3 = -8*N*s^2 + 4*(3*s^2 + (-4 + s)*w) - 8*N*s^2*nu
        q4 = -4*N*s^3 + 4*s*(s^2 - 2*w)

        Q = roots([q0,q1,q2,q3,q4])
        q = real.(Q[findall(x->isreal(x),Q)])
        qapp = q[findall(x->x^2<=1/4,q)]
        append!(inside_ext,qapp)
        append!(Narr_ext,N*ones(length(qapp)))
    end
    return inside_fp, inside_ext,Narr_fp, Narr_ext
end; 

function num(y::Float64,N::Int)
    return (1 + s*y)^(2*N*(1 - 8*nu/(4 - s^2)))*(1 - 2*y)^(4*N*nu/(2 + s) - 1)*(1 + 2*y)^(4*N*nu/(2 - s) - 1) 
end;

function P(y::Float64,N::Int)
    x, w = gausslegendre(3)
    x = x.*0.5
    den = dot(w,num.(x,N))
    return num(y,N)/den
end;

function get_avg_trajs(initfolder::String,ncopies::Int)
    F, _, Fref, _ = get_obs_track(initfolder*"sample_1/")
    for n in 2:ncopies
        f, _, fref, _ = get_obs_track(initfolder*"sample_$(n)/")
        if length(f)!=length(F)
            println(initfolder)
        end
        F.+=f
        Fref.+=fref
    end
    return F./ncopies, Fref./ncopies
end;