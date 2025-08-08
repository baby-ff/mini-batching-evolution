mutable struct Guy
    seq::Vector{Int}
    norm::Int
    batch::Vector{Int}
    keep::Vector{Int}
end;

Guy(s,b) = Guy(s,sum(s),b,zeros(length(b)));


mutable struct Population
    species::Vector{Guy}
    fits::Vector{Float64}
    avg_fits::Vector{Float64}
    gens::Vector{Float64}
    abund::Vector{Int}
end


mutable struct Settings
    Delta::Float64
    F0::Float64
    nu::Float64
    N::Int
    M::Int
    B::Int
    p_inherit::Float64
    L::Int
    K::Int
    pk::Float64
    p_env::Vector{Float64}
    consmatrix::Matrix{Float64}
    fieldmatrix::Matrix{Float64}
    maxc::Float64
end;


mutable struct EvoSettings
    nsteps::Int
    save_every::Int
    save_after::Int
    ncopies::Int
end;

# Define alias
const Sequence = Vector{Int}