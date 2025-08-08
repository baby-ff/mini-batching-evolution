function replicate(pop::Population, N::Int, L::Int)
    f = (pop.fits).*(pop.abund)
    p = f./sum(f)
    copies = rand(Multinomial(N,p))  
    deads = findall(x->x==0,copies)
    if length(deads)>0
        deleteat!(pop.species,deads)
        deleteat!(pop.gens,deads)
        deleteat!(pop.fits,deads)
        deleteat!(copies,deads)
    end
    #DO NOT DELETE AVG_FITS - PROBLEMS DEPENDENCY ON FITS BUT ONLY SOMETIMES!! 
    pop.abund = copies
    return pop
end;


function trim(pop::Population)
    del_idx = findall(x->x==0,pop.abund)
    deleteat!(pop.abund,del_idx)
    deleteat!(pop.species,del_idx)
    deleteat!(pop.fits,del_idx)
    deleteat!(pop.gens,del_idx)
    return pop
end


function mutate(pop::Population,L::Int, nu::Float64)
    nspecies = length(pop.abund)
    muts = Vector{Guy}()
    for i in 1:nspecies
        mi = 0
        s = pop.species[i].seq
        abi = pop.abund[i]
        for k in 1:abi
            mut_site = [rand(Binomial(1,nu)) for j in 1:L]
            if sum(mut_site)!=0
                pseq = (s + mut_site).%2
                guy = Guy(pseq,pop.species[i].batch)
                push!(muts,guy)
                mi += 1
            end
        end
        pop.abund[i] -= mi
    end
    lm = length(muts)
    append!(pop.species,muts)
    append!(pop.abund,ones(lm))
    append!(pop.fits,zeros(lm))
    append!(pop.gens,zeros(lm))
    
    pop = trim(pop)
    
    return pop
end;


function guy_array(x::Sequence)
    return Guy(x,[0])
end


function recollect(pop::Population,L::Int)
    #groups by sequence, cancel env
    nspecies = length(pop.abund)
    vec = Vector{Int}() 
    for i in 1:nspecies
        x = bin_to_int(pop.species[i].seq)
        xvec = ones(pop.abund[i]).*x
        append!(vec,xvec)
    end
    cmap = countmap(vec)
    seqs = collect(keys(cmap))
    sp = digits.(seqs,base=2,pad=L)
    ab = collect(values(cmap))
    species = [Guy(s,[0]) for s in sp]
    return (species, ab)
end;


function recollect_by_env(s::Int,envs::Vector{Vector{Int}},set::Settings)
    #groups by batch
    fieldmatrix, consmatrix, L, K, B = set.fieldmatrix, set.consmatrix, set.L, set.K, set.B
    
    seq = digits.(s,base=2,pad=L)
    
    fits = Vector{Float64}()
    avg_fits = Vector{Float64}()
    gens = Vector{Float64}()
    sp = Vector{Guy}()

    cmap = countmap(envs)
    env = collect(keys(cmap))
    ab = collect(values(cmap))
    for u in env
        guy = Guy(seq,u)
        push!(sp,guy)
        push!(fits,compute_fitness(fieldmatrix,guy,L,B))
        push!(avg_fits,compute_avg_fitness(fieldmatrix,guy,L))
        push!(gens,compute_gen(L,K,guy))
    end
    return (sp, ab, fits, avg_fits, gens)
end;


function update(pop::Population, set::Settings)
    N, p_env, B, fieldmatrix, consmatrix, p_inherit = set.N, set.p_env, set.B, set.fieldmatrix, set.consmatrix, set.p_inherit 
    
    xvec = Vector{Int}()
    avec = Vector{Int}()
    bvec = Vector{Vector{Int}}()
    nspecies = length(pop.abund)
    for i in 1:nspecies
        x = bin_to_int(pop.species[i].seq)
        push!(xvec,x)
        push!(avec,pop.abund[i])
        push!(bvec,pop.species[i].batch)
    end

    #rebatching and grouping by batch
    tmp_sp = Vector{Guy}()
    tmp_fits = Vector{Float64}()
    tmp_avg_fits = Vector{Float64}()
    tmp_gens = Vector{Float64}()
    tmp_ab = Vector{Int}()
    for u in unique(xvec)
        idx = findall(x->x==u,xvec)
        envs = Vector{Vector{Int}}()
        for i in idx
            b = [batching(rand(Binomial(1,p_inherit),B),set,bvec[i]) for n in 1:avec[i]]
            append!(envs,b)
        end
        sp, ab, fits, avg_fits, gens = recollect_by_env(u,envs,set)
        append!(tmp_sp,sp)
        append!(tmp_fits,fits)
        append!(tmp_avg_fits,avg_fits)
        append!(tmp_gens,gens)
        append!(tmp_ab,ab)
    end

    pop.species = deepcopy(tmp_sp)
    pop.fits = copy(tmp_fits)
    pop.avg_fits = copy(tmp_avg_fits)
    pop.gens = copy(tmp_gens)
    pop.abund = copy(tmp_ab)
    return pop
end;


function ref_update(pop_ref::Population,fieldmatrix::Matrix{Float64},L::Int,K::Int)
    sp, ab = recollect(pop_ref,L)
    pop_ref.species = sp
    pop_ref.abund = ab
    fits = Vector{Float64}()
    gens = Vector{Float64}()
    for guy in pop_ref.species
        if guy.norm==0
            push!(fits, mean(fieldmatrix))
        else
            push!(fits, (mean(fieldmatrix,dims=1)*guy.seq)[1]/guy.norm)
        end
        push!(gens, compute_gen(L,K,guy))
    end
    pop_ref.fits = copy(fits)
    pop_ref.avg_fits = copy(fits)
    pop_ref.gens = copy(gens)
    
    return pop_ref
end;


function main_evolve(n::Int, set::Settings, pop_ini::Population, pop_ref_ini::Population, nsteps::Int, save_every::Int, save_after::Int, init_cond::String)    
    pop = deepcopy(pop_ini);
    pop_ref = deepcopy(pop_ref_ini);
    N, L, K, nu, maxc, fieldmatrix, consmatrix, p_inherit = set.N, set.L, set.K, set.nu, set.maxc, set.fieldmatrix, set.consmatrix, set.p_inherit
    
    dir, fout, fout_ref = setup_outdir(n,set,init_cond);
    
    fout, fout_ref = save_track(0, pop, pop_ref, maxc, N, fout, fout_ref);
    
    for t in 1:nsteps
        #mini-batching
        pop = replicate(pop,N,L);
        pop = mutate(pop,L,nu);  #fitness not computed yet
        pop = update(pop,set); #fitness updated
        
        #reference
        pop_ref = replicate(pop_ref,N,L); 
        pop_ref = mutate(pop_ref,L,nu);
        pop_ref = ref_update(pop_ref,fieldmatrix,L,K);
        
        if t%save_every==0 && t>save_after
            fout, fout_ref = save_track(t, pop, pop_ref, maxc, N, fout, fout_ref);
            savesnap(t,pop,dir*("/saved_snaps"));
            savesnap(t,pop_ref,dir*("/saved_snaps_ref"));
        end
    end
    
    close(fout);
    close(fout_ref);
    
    return 
end