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
    moint::Array{Float64} # Integrals in MO Basis
    Ecorr_mp2::Float64 # MP2 Correlation Energy
    Ecorr_mp3::Float64 # MP3 Correlation Energy

    function MP(rhf::RHF)
        # Integral Transformation from AO basis to MO basis
        # (pr|qs) = \sum_{μ,ν,λ,σ} (μν|λσ) C_{μp} C_{νr} C_{λq} C_{σs} = <pq|rs>
        @tullio temp1[p,ν,λ,σ] := rhf.ijkl[μ,ν,λ,σ] * rhf.C[μ,p]
        @tullio temp2[p,r,λ,σ] := temp1[p,ν,λ,σ] * rhf.C[ν,r]
        @tullio temp3[p,r,q,σ] := temp2[p,r,λ,σ] * rhf.C[λ,q]
        @tullio moint[p,r,q,s] := temp3[p,r,q,σ] * rhf.C[σ,s]
        new(rhf, moint)
    end
end

function run_mp2!(mp::MP)
    """
    Performs MP2 calculation based on the procedure introduced in Szabo & Ostlund, pg. 352
    """
    # From Equation 6.73,
    # E_MP2 = (1/4) \sum_{abrs} |<ab||rs>|^2 / (e_a + e_b - e_r - e_s)
    # a, b loop over occupied and r, s loop over virtual
    occidx = 1:mp.rhf.mol.nelec # two spin orbitals per spatial orbitals,
                                # thus spin orbital per electron
    virtidx = mp.rhf.mol.nelec+1:2*size(mp.rhf.S)[1]

    mp.Ecorr_mp2 = 0.
    for a in occidx
        for b in occidx
            for r in virtidx
                @simd for s in virtidx
                    abrs = mp.moint[(a+1)÷2,(r+1)÷2,(b+1)÷2,(s+1)÷2] * (a%2 == r%2) * (b%2 == s%2) - mp.moint[(a+1)÷2,(s+1)÷2,(b+1)÷2,(r+1)÷2] * (a%2 == s%2) * (b%2 == r%2)
                    mp.Ecorr_mp2 += 0.25 * abrs * abrs / (mp.rhf.ϵ[(a+1)÷2]+mp.rhf.ϵ[(b+1)÷2]-mp.rhf.ϵ[(r+1)÷2]-mp.rhf.ϵ[(s+1)÷2])
                end
            end
        end
    end
    print("MP2 Correlation Energy = "*string(mp.Ecorr_mp2)*"\n")
end

function run_mp3!(mp::MP)
    """
    Performs MP3 calculation based on the procedure introduced in Szabo & Ostlund, pg. 352
    """
    # From Equation 6.75,
    # E_MP3 = (1/8) \sum_{abcdrs} <ab||rs><cd||ab><rs||cd> / ((e_a + e_b - e_r - e_s)(e_c + e_d - e_r - e_s))
    #        +(1/8) \sum_{abrstu} <ab||rs><rs||tu><tu||ab> / ((e_a + e_b - e_r - e_s)(e_a + e_b - e_t - e_u))
    #        + \sum_{abcrst} <ab||rs><cs||tb><rt||ac> / ((e_a + e_b - e_r - e_s)(e_a + e_c - e_r - e_t))
    # a, b, c, d loop over occupied and r, s, t, u loop over virtual
    occidx = 1:mp.rhf.mol.nelec # two spin orbitals per spatial orbitals,
                                # thus spin orbital per electron
    virtidx = mp.rhf.mol.nelec+1:2*size(mp.rhf.S)[1]

    mp.Ecorr_mp3 = 0.
    for a in occidx
        for b in occidx
            for c in occidx
                for d in occidx
                    for r in virtidx
                        @simd for s in virtidx
                            mp.Ecorr_mp3 += 0.125 * antisym(mp.moint, a, b, r, s)*antisym(mp.moint, c, d, a, b)*antisym(mp.moint, r, s, c, d) / ((mp.rhf.ϵ[(a+1)÷2]+mp.rhf.ϵ[(b+1)÷2]-mp.rhf.ϵ[(r+1)÷2]-mp.rhf.ϵ[(s+1)÷2])*(mp.rhf.ϵ[(c+1)÷2]+mp.rhf.ϵ[(d+1)÷2]-mp.rhf.ϵ[(r+1)÷2]-mp.rhf.ϵ[(s+1)÷2]))
                        end
                    end
                end
            end
        end
    end

    for a in occidx
        for b in occidx
            for r in virtidx
                for s in virtidx
                    for t in virtidx
                        @simd for u in virtidx
                            mp.Ecorr_mp3 += 0.125 * antisym(mp.moint, a, b, r, s)*antisym(mp.moint, r, s, t, u)*antisym(mp.moint, t, u, a, b) / ((mp.rhf.ϵ[(a+1)÷2]+mp.rhf.ϵ[(b+1)÷2]-mp.rhf.ϵ[(r+1)÷2]-mp.rhf.ϵ[(s+1)÷2])*(mp.rhf.ϵ[(a+1)÷2]+mp.rhf.ϵ[(b+1)÷2]-mp.rhf.ϵ[(t+1)÷2]-mp.rhf.ϵ[(u+1)÷2]))
                        end
                    end
                end
            end
        end
    end

    for a in occidx
        for b in occidx
            for c in occidx
                for r in virtidx
                    for s in virtidx
                        @simd for t in virtidx
                            mp.Ecorr_mp3 += antisym(mp.moint, a, b, r, s)*antisym(mp.moint, c, s, t, b)*antisym(mp.moint, r, t, a, c) / ((mp.rhf.ϵ[(a+1)÷2]+mp.rhf.ϵ[(b+1)÷2]-mp.rhf.ϵ[(r+1)÷2]-mp.rhf.ϵ[(s+1)÷2])*(mp.rhf.ϵ[(a+1)÷2]+mp.rhf.ϵ[(c+1)÷2]-mp.rhf.ϵ[(r+1)÷2]-mp.rhf.ϵ[(t+1)÷2]))
                        end
                    end
                end
            end
        end
    end

    print("MP3 Correlation Energy = "*string(mp.Ecorr_mp3)*"\n")
end

function antisym(moint, p, q, r, s)
    """
    Helper function that returns antisymmetrized two electron integral <pq||rs>
    """
    return moint[(p+1)÷2,(r+1)÷2,(q+1)÷2,(s+1)÷2] * (p%2 == r%2) * (q%2 == s%2) - moint[(p+1)÷2,(s+1)÷2,(q+1)÷2,(r+1)÷2] * (p%2 == s%2) * (q%2 == r%2)
end