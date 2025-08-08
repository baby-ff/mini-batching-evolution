function initialize_pop(set::Settings,init_cond::String)
    #set initial condition population
    N, B, L, K, fieldmatrix, consmatrix = set.N, set.B, set.L, set.K, set.fieldmatrix, set.consmatrix
    
    if occursin("rand",init_cond)
        s = randseq(set.L)
    elseif init_cond=="zeros"
        s = Int.(zeros(set.L))
    elseif init_cond=="ones"
        s = Int.(ones(set.L))
    end

    envs = [batching(Int.(zeros(B)),set,Int.(zeros(B))) for i in 1:N]
    species, cnt, fits, avg_fits, gens = recollect_by_env(bin_to_int(s),envs,set)
    pop = Population(species,fits,avg_fits,gens,cnt)
    
    guy = Guy(s,[0])
    pop_ref = Population([guy],[],[],[],[set.N])
    pop_ref = ref_update(pop_ref,set.fieldmatrix,set.L,set.K)
    return (pop, pop_ref)
end;


