function generate_fields(M::Int, L::Int, K::Int, pk::Float64, Delta::Float64, F0::Float64)
    c = reshape(append!(ones(K),zeros(L-K)),(1,L))
    c = c.+ (0.001.*ones(1,L))
    common = ones((M,1))*c
    variable = hcat(zeros(M,K),Delta*rand(Binomial(1,pk),(M,L-K)))
    return (F0.*c, F0.*(common + variable))
end;


function randseq(L::Int)
    #seq = 2 .* ([rand(Binomial(1)) for i in 1:L] .- 0.5)  #needs adjustment of fitness definition
    #seq = convert.(Int,seq)
    seq = rand(Binomial(1),L)
    return seq
end;


function copymatrix(folder::String)
    sfolder = folder*"ones/sample_1/"
    copymatrix = readdlm(sfolder*"consmatrix.txt",header=true)[1][1:end-1]
    fieldmatrix = readdlm(sfolder*"fieldmatrix.txt",header=true)[1][:,2:end]
    return copymatrix, fieldmatrix
end
