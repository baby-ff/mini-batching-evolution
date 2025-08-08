module SWF

using Random, Distributions, Statistics, DelimitedFiles, StatsBase;

push!(LOAD_PATH,"./")
include("SWF/types.jl");
include("SWF/random.jl");
include("SWF/batch.jl");
include("SWF/fitness.jl");
include("SWF/initialize.jl");
include("SWF/evolve.jl");
include("SWF/IO.jl");

export
    #from dataTypes
    Guy, Sequence, Population, Settings, EvoSettings,

    #from random
    randseq, generate_fields, copymatrix,

    #from batch
    sample_env, batching,

    #from fitness 
    compute_fitness, compute_gen, compute_avg_fitness,

    #from initialize
    initialize_pop,

    #from evolve
    replicate, mutate, update, ref_update, main_evolve,

    #from IO
    setup_outdir, track_evo , save_track, savesnap, bin_to_int
end