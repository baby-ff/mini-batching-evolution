function setup_outdir(s::Int,set::Settings,init_cond::String)
    Delta, F0, nu, N, M, B, L, K, pk, p_inherit = set.Delta, set.F0, set.nu, set.N, set.M, set.B, set.L, set.K, set.pk, set.p_inherit 
    consmatrix, fieldmatrix = set.consmatrix, set.fieldmatrix
    
    dir = "model_interface/data/persistent/"
    dir = dir*("/Delta$(Delta)_F0$(F0)_nu$(nu)_N$(N)_M$(M)_B$(B)_L$(L)_pers$(p_inherit)/half_fixed/long_K$(K)_pk$(pk)/")
    dir = dir*init_cond
    dir = dir*"/sample_$(s)"
    isdir(dir) || mkpath(dir)
    
    # print fields
    f = open(dir*("/fieldmatrix.txt"), "w")
    println(f,"#m \t fields-seq...")
    for m in 1:M
        print(f,m)
        print(f,"\t")
        for i in 1:L
            print(f,fieldmatrix[m,i])
            print(f,"\t")
        end
        print(f,"\n")
    end
    close(f)

    # print conserved part
    f = open(dir*"/consmatrix.txt", "w")
    println(f,"#fields-seq...")
    for i in 1:L
        print(f,consmatrix[1,i])
        print(f,"\t")
    end
    close(f)

    fout = open(dir*"/track_evo.txt", "w")
    println(fout,"#t \t F \t avg_F \t G")
    fout_ref = open(dir*"/track_evo_ref.txt", "w")
    println(fout_ref,"#t \t F \t (F) \t G")
    
    return (dir, fout, fout_ref)
end;


function track_evo(t::Int,F::Float64,avgF::Float64,G::Float64,nG::Float64,file::IOStream)
    print(file,t)
    print(file,"\t")
    print(file,F)
    print(file,"\t")
    print(file,avgF)
    print(file,"\t")
    print(file,G)
    print(file,"\t")
    println(file,nG)
    return file
end;


function save_track(t::Int, pop::Population, pop_ref::Population, maxc::Float64, N::Int, fout::IOStream, fout_ref::IOStream)
    F = transpose(pop.fits)*pop.abund
    avgF = transpose(pop.avg_fits)*pop.abund
    G = transpose(pop.gens)*pop.abund
    norm_G = G/maxc
    F_ref = transpose(pop_ref.fits)*pop_ref.abund
    G_ref = transpose(pop_ref.gens)*pop_ref.abund
    norm_G_ref = G_ref/maxc 
    fout = track_evo(t,F,avgF,G,norm_G,fout)
    fout_ref = track_evo(t,F_ref,F_ref,G_ref,norm_G_ref,fout_ref) 
    
    return (fout, fout_ref)
end;


function bin_to_int(str::Vector{Int})
    str = reverse(join(string.(str)))
    x = parse(Int, str; base=2)
    return x
end


function savesnap(t::Int,pop::Population,dir::String)
    isdir(dir) || mkpath(dir)
    f = open(dir*"/snap_t$(t).txt", "w")
    println(f,"#seq\t batch \t fit \t gen \t abund")
    for i in 1:length(pop.abund)
        print(f,bin_to_int(pop.species[i].seq))
        print(f,"\t")
        print(f,pop.species[i].batch)
        print(f,"\t")
        print(f,pop.fits[i])
        print(f,"\t")
        print(f,pop.avg_fits[i])
        print(f,"\t")
        print(f,pop.gens[i])
        print(f,"\t")
        print(f,pop.abund[i])
        println(f,);
    end
    close(f)
    return
end;