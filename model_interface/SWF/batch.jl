function sample_env(p_env::Array{Float64},B::Int)
    multi = rand(Multinomial(B,p_env))
    env_list = Vector{Int}()
    for j in 1:maximum(multi)
        append!(env_list,findall(multi.>=j))
    end
    return env_list
end;


function batching(keep::Vector{Int},set::Settings,old_batch::Vector{Int})
    p_env, B = set.p_env, set.B

    batch = copy(old_batch)
    new_batch = findall(x->x==0,keep)
    if length(new_batch)>0
        batch[new_batch] = sample_env(p_env,length(new_batch))
    end
    return batch
end;