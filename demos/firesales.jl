using FinNetValu

N = 3  # banks
M = 5  # assets
A = [0 1 1 0 0;
     0 1 0 1 0;
     1 0 0 0 0]
aᵉ = [0.8, 0.8, 0.8]
aⁱ = [0.2, 0.2, 0.2]
d = [0.9, 0.9, 0.9]

loss = 0.3
α = 1.0536

aⁱ[vec(sum(A, dims=2) .== 0)] .= 1

fsmodel = FSModel(A, aᵉ, aⁱ, d, α=α)

for idx in 1:M
    a = ones(M)
    a[idx] *= loss

    fp = fixvalue(fsmodel, a, m=0)

    println("shocked asset index: ", idx,", failed banks: " , N - sum(solvent(fsmodel, fp)))
end
