using DelimitedFiles, LinearAlgebra, Tullio
export MP

#=
QCPractice.MP

Moller-Plesset Perturbation-type Calculation
=#

mutable struct MP
    """
    Main object for Moller-Plesset Perturbation-type Calculation
    
    * Inputs *
        rhf     QCPractice.RHF Object that includes restricted Hartree-Fock calculation
    * Outputs *
        
    """
    rhf::RHF # Input RHF Object
    Ecorr_mp2::Float64 # MP2 Correlation Energy

    function MP(rhf::RHF)
        new(rhf)
    end
end

function run_mp2!(mp::MP)
    """
    Performs MP2 calculation based on the procedure introduced in Szabo & Ostlund, pg. 352
    """

    # Check n_{occupied} is even
    try Int(mp.rhf.mol.nelec/2)
    catch InexactError
        throw(ArgumentError("Restricted Hartree Fock calculation does not accept odd number of electrons."))
    end

    # From Equation 6.73,
    # E_MP2 = (1/2) \sum_{abrs} <ab|rs><rs|ab> / (e_a + e_b - e_r - e_s)
    #        -(1/2) \sum_{abrs} <ab|rs><rs|ba> / (e_a + e_b - e_r - e_s)
    # a, b loop over occupied and r, s loop over virtual
    occidx = 1:mp.rhf.mol.nelec # two spin orbitals per spatial orbitals,
                                # thus spin orbital per electron
    virtidx = mp.rhf.mol.nelec+1:2*size(mp.rhf.S)[1]

    # Integral Transformation from AO basis to MO basis
    # (pr|qs) = \sum_{μ,ν,λ,σ} (μν|λσ) C_{μp} C_{νr} C_{λq} C_{σs} = <pq|rs>
    @tullio temp1[p,ν,λ,σ] := mp.rhf.ijkl[μ,ν,λ,σ] * mp.rhf.C[μ,p]
    @tullio temp2[p,r,λ,σ] := temp1[p,ν,λ,σ] * mp.rhf.C[ν,r]
    @tullio temp3[p,r,q,σ] := temp2[p,r,λ,σ] * mp.rhf.C[λ,q]
    @tullio moint[p,r,q,s] := temp3[p,r,q,σ] * mp.rhf.C[σ,s]

    mp.Ecorr_mp2 = 0
    for a in occidx
        for b in occidx
            for r in virtidx
                for s in virtidx
                    abrs = moint[(a+1)÷2,(r+1)÷2,(b+1)÷2,(s+1)÷2] * (a%2 == r%2) * (b%2 == s%2) - moint[(a+1)÷2,(s+1)÷2,(b+1)÷2,(r+1)÷2] * (a%2 == s%2) * (b%2 == r%2)
                    mp.Ecorr_mp2 += 0.25 * abrs * abrs / (mp.rhf.ϵ[(a+1)÷2]+mp.rhf.ϵ[(b+1)÷2]-mp.rhf.ϵ[(r+1)÷2]-mp.rhf.ϵ[(s+1)÷2])
                end
            end
        end
    end
    print("MP2 Correlation Energy = "*string(mp.Ecorr_mp2)*"\n")
end

