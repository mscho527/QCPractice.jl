using DelimitedFiles, LinearAlgebra, Tullio
export RHF

#=
QCPractice.RHF

Restricted Hartree Fock Calculation
=#

mutable struct RHF
    """
    Main object for Restricted Hartree Fock Calculation
    
    * Inputs *
        mol     QCPractice.Molecule Object that includes atomic symbols, coordinates, and basis set information
    * Outputs *
        
    """
    mol::Molecule # Input Molecule Object
    S::Hermitian # Overlap Matrix
    X::Matrix{Float64} # Transformation Matrix
    hcore::Matrix{Float64} # Core Hamiltonian
    F::Matrix{Float64} # Fock Matrix
    C::Matrix{Float64} # Coefficient Matrix
    P::Matrix{Float64} # Density Matrix
    ijkl::Array{Float64} # Two-electron Integral Matrix
    E::Float64 # RHF Energy
    ϵ::Vector{Float64} # Orbital Energies
    Nuc::Float64 # Nuclear-nuclear Repulsion Energy (Hartrees)

    function RHF(mol::Molecule)
        new(mol)
    end
end

function run_rhf!(rhfobj::RHF, istroyint::Bool=false)
    """
    Performs Hartree-Fock calculation, based on the procedure introduced in Szabo & Ostlund, pg. 146
    Symmetric Orthogonalization is used for the calculation.
    """
    # Check n_{occupied} is even
    try Int(rhfobj.mol.nelec/2)
    catch InexactError
        throw(ArgumentError("Restricted Hartree Fock calculation does not accept odd number of electrons."))
    end

    if istroyint ≠ true
        # For libcint binding (Default)
        # Obtain transformation matrix X for symmetric orthogonalization
        # X = Us^(-1/2)U
        # One could also eigendecompose, and do X = S^(-1/2) = U*s^(-1/2)*Uᵀ
        # However, @btime runs show that S^(-1/2) is faster
        # (by using Symmetric type-specific routines)
        # S_eigen = eigen(rhfobj.S)
        # rhfobj.X = S_eigen.vectors * diagm(S_eigen.values.^(-1/2)) * transpose(S_eigen.vectors)
        rhfobj.S = Hermitian(libcint_S(rhfobj.mol)) # overlap matrix from libcint
                                                    # Hermitian tag improves efficiency
        rhfobj.X = rhfobj.S^(-1/2)

        # TODO: Add sparse NTuple version
        # Obtain two-electron integrals and store for future use
        rhfobj.ijkl = libcint_ijkl(rhfobj.mol)

        rhfobj.hcore = libcint_T(rhfobj.mol) + libcint_V(rhfobj.mol)

        # TODO: Implement Nuc-Nuc energy
        rhfobj.Nuc = 0.0
    else
        # For troyint, mol.basis contains path to the input file
        # Line 1 // Nuc-Nuc Repulsion Energy
        # Line 2 // Int::nbasis, Int::nocc/2 (Molecule.jl reads this line)
        # Line 3 // Vector{Float64}::S
        # Line 4 // Vector{Float64}::h
        # Line 5 // Vector{Float64}::ijkl
        input = readdlm(rhfobj.mol.basis)
        nbasis = input[2,1]
        rhfobj.S = Hermitian(reshape(convert(Vector{Float64}, input[3,1:nbasis^2]), nbasis, nbasis))
        rhfobj.X = rhfobj.S^(-1/2)

        rhfobj.hcore = reshape(convert(Vector{Float64}, input[4,1:nbasis^2]), nbasis, nbasis)
        rhfobj.Nuc = input[1]

        # Unpack packed input matrix
        # For a given ijkl, set i = max(i,j), j = min(i,j)
        #                       k = max(k,l), l = min(k,l)
        # then assume column-major ordering of 4-dim tensor
        # TODO: With Sparse N-Tuple above, add implementation
        # that does not unpack the matrix (reduce declarations)
        rhfobj.ijkl = Array{Float64}(undef, nbasis, nbasis, nbasis, nbasis) # preallocate
        for i in 1:nbasis
            for j in 1:nbasis
                for k in 1:nbasis
                    @simd for l in 1:nbasis # use simd directive
                        ij = (max(i, j) - 1)*(max(i, j)) ÷ 2 + min(i, j)
                        kl = (max(k, l) - 1)*(max(k, l)) ÷ 2 + min(k, l)
                        rhfobj.ijkl[i,j,k,l] = input[5, (max(ij, kl) - 1)*(max(ij, kl)) ÷ 2 + min(ij, kl)]
                    end
                end
            end
        end
    end

    # Using zero charge density P_{initial}, calculate the Fock matrix
    rhfobj.F = rhfobj.hcore # initial guess for the Fock matrix := hcore

    # Extended Huckel for Initial Fock Matrix
    #rhfobj.F = Array{Float64}(undef, size(rhfobj.hcore)...) # preallocate
    #for i = 1:size(rhfobj.hcore)[1]
    #    @simd for j = 1:size(rhfobj.hcore)[1]
    #        rhfobj.F[i,j] = 1.75 * rhfobj.S[i,j] * (rhfobj.hcore[i,i] * rhfobj.hcore[j,j]) / 2.0
    #    end
    #end

    # Using the Fock matrix, find initial ϵ and C
    F_trans = transpose(rhfobj.X)*rhfobj.F*rhfobj.X
    F_trans_eigen = eigen(F_trans)
    rhfobj.C = rhfobj.X*F_trans_eigen.vectors
    rhfobj.ϵ = F_trans_eigen.values

    # Form initial guess for the density matrix
    rhfobj.P = build_density_matrix(rhfobj.C, rhfobj.mol.nelec)
    rhfobj.E = evaluate_energy(rhfobj.P, rhfobj.hcore, rhfobj.F) + rhfobj.Nuc
    print("iter = 0, E = "*string(rhfobj.E)*"\n")

    iter = 1
    # SCF Loop
    while true
        # Build Fock Matrix using Szabo Eqn. 3.154
        # F[μ,ν] = Hcore[μ,ν]+∑_{λ,σ}(P[λ,σ]*((μν|σλ) - 0.5*(μλ|σν)))
        @tullio G[μ,ν] := rhfobj.P[λ,σ]*(rhfobj.ijkl[μ,ν,σ,λ] - 0.5*rhfobj.ijkl[μ,λ,σ,ν])
        rnd = 0.4 + 0.3 * rand(Float64)
        rhfobj.F = rnd*rhfobj.F + (1-rnd) * (rhfobj.hcore + G) # prevent oscillation
        F_trans = transpose(rhfobj.X)*rhfobj.F*rhfobj.X

        F_trans_eigen = eigen(F_trans)
        rhfobj.C = real.(rhfobj.X*F_trans_eigen.vectors)
        rhfobj.ϵ = real.(F_trans_eigen.values)

        oldp = copy(rhfobj.P)
        rhfobj.P = build_density_matrix(rhfobj.C, rhfobj.mol.nelec)
        prevE = rhfobj.E
        rhfobj.E = evaluate_energy(rhfobj.P, rhfobj.hcore, rhfobj.F) + rhfobj.Nuc
        print("iter = "*string(iter)*", E = "*string(rhfobj.E)*"\n")
        if prevE ≈ rhfobj.E && oldp ≈ rhfobj.P # Convergence
            break
        end
        iter = iter + 1
    end
end

function build_density_matrix(C, nelec)
    """
    Builds new density matrix P using given coefficient matrix C
    From Szabo Eqn. 3.145,
    P_{μ,ν} = 2 ∑_{a ∈ 1:N/2} (C_{μ,a} C_{ν,a})
    """
    Coef = copy(C)
    Coef[:, nelec÷2+1:end] .= zero(C[1,1])
    P = Matrix{Float64}(undef,size(C)...)
    @tullio P[μ,ν] := 2*Coef[μ,a]*Coef[ν,a]'

    return P
end

function evaluate_energy(P, hcore, F)
    """
    Evaluate Hartree Fock Energy from the given inputs
    From Szabo Eqn. 3.184,
    E_0 = 0.5*∑_{μ,ν}(P[ν,μ]*(hcore[μ,ν]+F[μ,ν]))
    """
    @tullio E := 0.5 * P[ν,μ]*(hcore[μ,ν]+F[μ,ν])
    return E
end
