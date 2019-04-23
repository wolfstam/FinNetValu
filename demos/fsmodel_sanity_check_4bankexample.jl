using FinNetValu

# N = 4  # banks
# M = 7  # assets
# A = [0 0 1 1 0 0 0;
#      0 1 0 0 1 0 1;
#      1 0 1 0 1 1 0;
#      1 0 0 1 0 0 1]
#
# aᵉ = [0.8, 0.8, 0.8, 0.8]
# aⁱ = [0.2, 0.2, 0.2, 0.2]
# d = [0.9, 0.9, 0.9, 0.9]

# N = 3  # banks
# M = 5  # assets
# A = [0 1 0 0 0;
#     0 0 0 1 0;
#     1 0 0 0 0]
# aᵉ = [0.8, 0.8, 0.8]
# aⁱ = [0.2, 0.2, 0.2]
# d = [0.9, 0.9, 0.9]

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
    println("Idx: ", idx)
    a = ones(M)
    a[idx] *= loss
    println("shocked asset prices: ", a)

    fp = fixvalue(fsmodel, a, m=0)
    # x = FinNetValu.init(fsmodel, a)
    # fp = ones(size(x))
    # valuation!(fp, fsmodel, x, a)
    # valuation!(fp, fsmodel, fp, a)
    # valuation!(fp, fsmodel, fp, a)
    # # valuation!(fp, fsmodel, fp, a)
    # # valuation!(fp, fsmodel, fp, a)

    println("Failed banks: " , N - sum(solvent(fsmodel, fp)))
    println("------------------------------------")
end
