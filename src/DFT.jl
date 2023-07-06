using DelimitedFiles, LinearAlgebra, Tullio

#=
QCPractice.DFT

Tiny snippet related to Density Functional Theory
=#

# Input file name (194.Density.dat, 302.Density.dat, 590.Density.dat)
grid_path = "src/Tests/194.Density.dat"

# Each .dat file has the following format, delimited by space
# Idx 1 2 3 4      5   6   7      8      9      10     11     12
#     x y z weight ρ^A ρ^B ∇ρ^A_x ∇ρ^A_y ∇ρ^A_z ∇ρ^B_x ∇ρ^B_y ∇ρ^B_z
input = readdlm(grid_path)

# Sanity check: make sure that ∫(ρ^A+ρ^B)dr = n_elec
@tullio numint := input[i,4] * (input[i,5] + input[i,6])
@assert isapprox(numint, 14.0, atol = 1e-5) "Numerically integrating ρ over all grid points did not sum up to 14."

# Applying LDA Exchange Functional from Parr & Yang Eqn. 7.4.5,
# ε_x(ρ) = -(3/4 * (3/π)^(1/3)) * ∫ρ(r)^(4/3) d³r
@tullio ldax := -(0.75 * (3/3.141592)^(1/3)) * input[i,4]*(input[i,5] + input[i,6])^(4/3)

# Using B88 Exchange Functional from Becke, Phys. Rev. A, 38 (6), 1988 Eqn. 8,
# ε'_x(ρ) = ε_{LDAx}(ρ) - β ∫ρ^(4/3) * x^2 / (1 + 6βx sinh^(-1)x) d³r
# x = |∇ρ|/(ρ^(4/3)), β = 0.0042
# Sieve out grids that will cause divide-by-zero
validgrids = .~isapprox.(0.,input[:,5].^(4/3))
input = input[validgrids,:]
@tullio b88x_correction := -(0.0042) * input[i,4]*(input[i,5])^(4/3) * (sqrt(input[i,7]^2 + input[i,8]^2 + input[i,9]^2)/input[i,5]^(4/3))^2 / (1 + 6*0.0042*(sqrt(input[i,7]^2 + input[i,8]^2 + input[i,9]^2)/input[i,5]^(4/3))*asinh((sqrt(input[i,7]^2 + input[i,8]^2 + input[i,9]^2)/input[i,5]^(4/3)))) - (0.0042) * input[i,4]*(input[i,6])^(4/3) * (sqrt(input[i,10]^2 + input[i,11]^2 + input[i,12]^2)/input[i,6]^(4/3))^2 / (1 + 6*0.0042*(sqrt(input[i,10]^2 + input[i,11]^2 + input[i,12]^2)/input[i,6]^(4/3))*asinh((sqrt(input[i,10]^2 + input[i,11]^2 + input[i,12]^2)/input[i,6]^(4/3))))
b88x = ldax + b88x_correction
