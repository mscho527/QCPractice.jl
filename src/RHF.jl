using LinearAlgebra, Tullio
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

    function RHF(mol::Molecule)
        new(mol)
    end
end

function run_rhf!(rhfobj::RHF)
    """
    Performs Hartree-Fock calculation, based on the procedure introduced in Szabo & Ostlund, pg. 146
    Symmetric Orthogonalization is used for the calculation.
    """
    # Check n_{occupied} is even
    try Int(rhfobj.mol.nelec/2)
    catch InexactError
        throw(ArgumentError("Restricted Hartree Fock calculation does not accept odd number of electrons."))
    end

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

    # TODO: Add other initial SCF guesses later
    # Using zero charge density P_{initial}, calculate the Fock matrix
    rhfobj.hcore = libcint_T(rhfobj.mol) + libcint_V(rhfobj.mol)
    rhfobj.F = rhfobj.hcore # initial guess for the Fock matrix := hcore
    F_trans = transpose(rhfobj.X)*rhfobj.F*rhfobj.X

    # Using the Fock matrix, find initial ϵ and C
    F_trans_eigen = eigen(F_trans)
    rhfobj.C = rhfobj.X*F_trans_eigen.vectors
    rhfobj.ϵ = F_trans_eigen.values

    # Form initial guess for the density matrix
    rhfobj.P = build_density_matrix(rhfobj.C, rhfobj.mol.nelec)

    # SCF Loop
    while true
        # Build Fock Matrix using Szabo Eqn. 3.154
        # F[μ,ν] = Hcore[μ,ν]+∑_{λ,σ}(P[λ,σ]*((μν|σλ) - 0.5*(μλ|σν)))
        @tullio G[μ,ν] := rhfobj.P[λ,σ]*(rhfobj.ijkl[μ,ν,σ,λ] - 0.5*rhfobj.ijkl[μ,λ,σ,ν])
        rhfobj.F = rhfobj.hcore + G
        F_trans = transpose(rhfobj.X)*rhfobj.F*rhfobj.X

        F_trans_eigen = eigen(F_trans)
        rhfobj.C = rhfobj.X*F_trans_eigen.vectors
        rhfobj.ϵ = F_trans_eigen.values

        Pnew = build_density_matrix(rhfobj.C, rhfobj.mol.nelec)
        if Pnew ≈ rhfobj.P # Convergence
            break
        end
        rhfobj.P = Pnew
    end
    rhfobj.E = evaluate_energy(rhfobj.P, rhfobj.hcore, rhfobj.F)
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
