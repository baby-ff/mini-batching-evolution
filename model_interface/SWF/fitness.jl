function compute_fitness(fieldmatrix::Matrix{Float64}, guy::Guy, L::Int, B::Int)
    f = 0
    if guy.norm==0
        for idx in guy.batch
            f += mean(fieldmatrix[idx,:])
        end
    else
        for idx in guy.batch
            f += (fieldmatrix[idx,:]'*guy.seq)/guy.norm
        end
    end
    return f/B
end;


function compute_avg_fitness(fieldmatrix::Matrix{Float64}, guy::Guy, L::Int)
    if guy.norm==0
        f = mean(fieldmatrix)
    else
        f = (mean(fieldmatrix,dims=1)*guy.seq)[1]/guy.norm
    end
    return f
end;


function compute_gen(L::Int, K::Int, guy::Guy)
    if guy.norm==0
        g = K/L
    else
        g = sum(guy.seq[1:K])/guy.norm
    end
    return g
end;

