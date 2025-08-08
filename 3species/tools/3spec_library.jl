#types
mutable struct Guy
    gtype::Int
    batch::Vector{Int}
    keep::Vector{Int}
end;

Guy(s,b) = Guy(s,b,zeros(length(b)));

mutable struct Population
    species::Vector{Guy}
    fits::Vector{Float64}
    avg_fits::Vector{Float64}
    vars::Vector{Float64}
    abund::Vector{Int}
end

mutable struct Settings
    N::Int
    B::Int
    M::Int
    land_vec::Matrix{Float64}
    avg_vec::Matrix{Float64}
    var_vec::Matrix{Float64}
end
Settings(N,B,M,land_vec) = Settings(N,B,M,land_vec,transpose(mean(land_vec,dims=2)),transpose(var(land_vec,dims=2)));

mutable struct EvoSettings
    nu::Float64
    nsteps::Int
    ncopies::Int
    save_every::Int
    save_after::Int
end

#landscapes
function get_m_s(f::Float64,V::Float64)
    sigma = sqrt(log(1+V/(f^2)))
    return log(f)-sigma^2/2, sigma
end

function generate_landscapes(M::Int64,m::Float64,s::Float64)
    return exp.(rand(Normal(m,s),M))
end

function get_landvec(M::Int,fa::Float64,fb::Float64,fc::Float64,Va::Float64,Vb::Float64,Vc::Float64)
    m_a, s_a = get_m_s(fa,Va)
    m_b, s_b = get_m_s(fb,Vb)
    m_c, s_c = get_m_s(fc,Vc)
    land_vec = hcat(generate_landscapes(M,m_a,s_a), generate_landscapes(M,m_b,s_b), generate_landscapes(M,m_c,s_c))
    return land_vec
    
end


#batching
function batching(keep::Vector{Int},B::Int,old_batch::Vector{Int})
    batch = copy(old_batch)
    new_batch = findall(x->x==0,keep)
    if length(new_batch)>0
        batch[new_batch] = rand(1:M,B)  #simpler than main sim - can modify main sim too
    end
    return batch
end;


#initialization
function initialize_pop(set::Settings,init_cond::String)
    N, B, land_vec, avg_vec, var_vec = set.N, set.B, set.land_vec, set.avg_vec, set.var_vec
    
    if init_cond=="balanced"
        s = rand(1:3, N) 
    elseif init_cond=="all_a"
        s = Int.(ones(N))
    elseif init_cond=="all_b"
        s = Int.(2 .* ones(N))
    elseif init_cond=="all_c"
        s = Int.(3 .* ones(N))
    end
    envs = [batching(Int.(zeros(B)),B,Int.(zeros(B))) for i in 1:N]
    
    if occursin("all",init_cond)
        #monoclonal initial condition
        species, cnt, fits, avg_fits, vars = recollect_by_env(s[1],envs,land_vec[s,:],avg_vec[s],var_vec[s])
        guy = Guy(s[1],[0])
        pop_ref = Population([guy],[],[],[],[set.N])
    else
        #polyclonal initial condition
        cm = countmap(s)
        species = Vector{Guy}()
        species_ref = Vector{Guy}()
        cnt = Vector{Int}()
        fits = Vector{Float64}()
        avg_fits = Vector{Float64}()
        vars = Vector{Float64}()
        ab_values = Int.(values(cm))
        for s_value in keys(cm)
            env_value = envs[findall(x->x==s_value,s)]
            sp, ct, f_sp, avg_sp, v_sp = recollect_by_env(s_value,env_value,land_vec[s_value,:],avg_vec[s_value],var_vec[s_value])
            append!(species,sp)
            append!(cnt,ct)
            append!(fits,f_sp)
            append!(avg_fits,avg_sp)
            append!(vars,v_sp)
            push!(species_ref,Guy(s_value,[0]))
        end
        pop_ref = Population(species_ref,[],[],[],ab_values)
    end
    pop = Population(species,fits,avg_fits,vars,cnt)
    pop_ref = ref_update(pop_ref,set)
    return (pop, pop_ref)
end;

function ref_update(pop_ref::Population,set::Settings)
    #computes fitnesses for reference population 
    avg_vec, var_vec = set.avg_vec, set.var_vec
    pop_ref.fits = copy([avg_vec[1,species.gtype] for species in pop_ref.species])
    pop_ref.avg_fits = copy(pop_ref.fits)
    pop_ref.vars = copy([var_vec[1,species.gtype] for species in pop_ref.species])
    return pop_ref
end


#grouping by batch
function recollect_by_env(s_value::Int,envs::Vector{Vector{Int}},land::Vector{Float64},avg_land::Float64,vars_land::Float64)
    
    fits = Vector{Float64}()
    vars = Vector{Float64}()
    avg_fits = Vector{Float64}()
    sp = Vector{Guy}()
    cmap = countmap(envs)
    env = collect(keys(cmap))
    ab = collect(values(cmap))
    for u in env
        guy = Guy(s_value,u)
        push!(sp,guy)
        push!(fits,land[u[1]])
        push!(avg_fits,avg_land)
        push!(vars,vars_land)
    end
    return (sp, ab, fits, avg_fits, vars)
end;


#evolve
function main_evolve(n::Int,set::Settings,evoset::EvoSettings,pop_ini::Population,pop_ref_ini::Population,int_cond::String,dir::String)
    pop = deepcopy(pop_ini)
    pop_ref = deepcopy(pop_ref_ini)
    
    N, nu = set.N, evoset.nu

    fout, fout_ref = setup_copydir(n,dir,init_cond);
    
    fout, fout_ref = save_track(0, pop, pop_ref,fout,fout_ref)
    
    for t in 1:nsteps
        #mini-batching
        pop = replicate(pop,N);
        pop = mutate(pop,nu);  #fitness not computed yet
        pop = update(pop,set); #fitness updated
        
        #reference
        pop_ref = replicate(pop_ref,N); 
        pop_ref = mutate(pop_ref,nu);
        pop_ref = ref_update(pop_ref,set);
        
        if t%save_every==0 && t>save_after
            fout, fout_ref = save_track(t, pop, pop_ref, fout, fout_ref);
        end
    end
    
    close(fout);
    close(fout_ref);
    
    return
end;

function replicate(pop::Population,N::Int)
    f = (pop.fits).*(pop.abund)
    p = f./sum(f)
    copies = rand(Multinomial(N,p))  
    deads = findall(x->x==0,copies)
    if length(deads)>0
        deleteat!(pop.species,deads)
        deleteat!(pop.vars,deads)
        deleteat!(pop.fits,deads)
        deleteat!(pop.avg_fits,deads)
        deleteat!(copies,deads)
    end
    pop.abund = copies
    return pop
end

function mutate(pop::Population,nu::Float64)
    nspecies = length(pop.abund)
    muts = zeros(3)
    for i in 1:nspecies
        abi = pop.abund[i]
        m1 = rand(Binomial(1,nu),abi)
        m2 = rand(Binomial(1,nu),abi)
        gt1 = (pop.species[i].gtype%3)+1
        gt2 = ((pop.species[i].gtype+1)%3)+1
        
        muts[gt1] += sum(m1)
        muts[gt2] += sum(m2)
        pop.abund[i] -= sum(m1)+sum(m2)    
    end
    
    for k in 1:3
        guy = Guy(k,[0])
        push!(pop.species,deepcopy(guy))
        push!(pop.fits,0)
        push!(pop.avg_fits,0)
        push!(pop.vars,0)
        push!(pop.abund,muts[k])
    end
    pop = trim(pop)
    return pop
end

function trim(pop::Population)
    del_idx = findall(x->x==0,pop.abund)
    deleteat!(pop.abund,del_idx)
    deleteat!(pop.species,del_idx)
    deleteat!(pop.fits,del_idx)
    deleteat!(pop.avg_fits,del_idx)
    deleteat!(pop.vars,del_idx)
    return pop
end

function update(pop::Population,set::Settings)
    B, land_vec, avg_vec, var_vec = set.B, set.land_vec, set.avg_vec, set.var_vec 
    #rebatching and grouping by batch
    tmp_sp = Vector{Guy}()
    tmp_fits = Vector{Float64}()
    tmp_avg_fits = Vector{Float64}()
    tmp_vars = Vector{Float64}()
    tmp_ab = Vector{Int}()

    for k in 1:3
        idx = findall(x->x.gtype==k,pop.species)
        envs = Vector{Vector{Int}}()
        for i in idx
            b = [batching(rand(Binomial(1,0.),B),set.B,[0]) for n in 1:pop.abund[i]]
            append!(envs,b)
        end
        sp, ab, fits, avg_fits, vars = recollect_by_env(k,envs,land_vec[k,:],avg_vec[k],var_vec[k])
        
        append!(tmp_sp,sp)
        append!(tmp_fits,fits)
        append!(tmp_avg_fits,avg_fits)
        append!(tmp_vars,vars)
        append!(tmp_ab,ab)
    end
    pop.species = deepcopy(tmp_sp)
    pop.fits = copy(tmp_fits)
    pop.avg_fits = copy(tmp_avg_fits)
    pop.vars = copy(tmp_vars)
    pop.abund = copy(tmp_ab)
    return pop
end

function ref_update(pop_ref::Population,set::Settings)
    avg_vec, var_vec = set.avg_vec, set.var_vec
    tmp_ab = zeros(3)
    tmp_sp = Vector{Guy}()
    for k in 1:3
        idx = findall(x->x.gtype==k,pop_ref.species)
        tmp_ab[k] = sum(pop_ref.abund[idx])
        push!(tmp_sp,Guy(k,[0]))
    end
    pop_ref.species = deepcopy(tmp_sp)
    pop_ref.abund = copy(tmp_ab)
    pop_ref.fits = copy(avg_vec[1,:])
    pop_ref.avg_fits = copy(avg_vec[1,:])
    pop_ref.vars = copy(var_vec[1,:])
    
    return pop_ref
end


#IO
function setup_outdir(set::Settings,evoset::EvoSettings)
    N,M,nu = set.N, set.M, evoset.nu
    land_vec, avg_vec, var_vec = set.land_vec, set.avg_vec, set.var_vec 
    fa, fb ,fc = avg_vec[1], avg_vec[2], avg_vec[3]
    Va, Vb, Vc = var_vec[1], var_vec[2], var_vec[3]
    #create directory
    dir = pwd()*"/data/N$(N)_M$(M)_nu$(nu)_fa$(round(fa,sigdigits=3))_fb$(round(fb,sigdigits=3))_fc$(round(fc,sigdigits=3))_Va$(round(Va,sigdigits=3))_Vb$(round(Vb,sigdigits=3))_Vc$(round(Vc,sigdigits=3))/"
    mkpath(dir)
    #write landscapes
    fout = open(dir*"landvec.txt","w")
    for i in 1:M
        print(fout, land_vec[1,i])
        print(fout, "\t")
        print(fout, land_vec[2,i])
        print(fout, "\t")
        print(fout, land_vec[3,i])
        print(fout, "\n")
    end
    close(fout)
    return dir
end

function setup_copydir(n::Int,dir::String,init_cond::String)
    #create directory
    copydir = dir*init_cond*"/copy_$(n)/"
    mkpath(copydir)
    #create files
    fout = open(copydir*"pop.txt","w")
    fout_ref = open(copydir*"pop_ref.txt","w")
    
    return fout,fout_ref
end

function save_track(t::Int,pop::Population,pop_ref::Population,fout::IOStream,fout_ref::IOStream)
    print(fout,t,"\t")
    print(fout_ref,t,"\t")
    nvec = zeros(3)
    nvec_ref = zeros(3)
    for k in 1:3
        idx = findall(x->x.gtype==k,pop.species)
        idx_ref = findall(x->x.gtype==k,pop_ref.species)
        nvec[k] = sum(pop.abund[idx])
        nvec_ref[k] = sum(pop_ref.abund[idx_ref])
    end
    print(fout,nvec[1],"\t",nvec[2],"\t",nvec[3],"\n")
    print(fout_ref,nvec_ref[1],"\t",nvec_ref[2],"\t",nvec_ref[3],"\n")
    return fout, fout_ref
end

