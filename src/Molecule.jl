#= QCPractice.Molecule
   Submodule that contains Julia struct Molecule =#

using DelimitedFiles
struct Molecule
    xyzstring::String # xyz string in the format below:
    # "H 0.00 0.00 0.00\nH 0.70 0.00 0.00" (unit = angstroms = 1.889 bohr)
    basis::String # basis set (if troyint, path to data file)
    nelec::Int # number of electrons
end

function Molecule(element::Array{String}, atomxyz::AbstractArray, basis::String, nelec::Int)
    """
    TODO: Make this more generic using AbstractTypes
    Constructor function
        Takes in necessary inputs and creates xyzstring to prepare for libcint calls
        element     atomic symbols of constituent atoms
        atomxyz     atomic coordinates in the Cartesian system
    """
    geometry = hcat(element, atomxyz)
    xyzstring = join([join(view(geometry, i, :), " ") for i in 1:size(geometry)[1]], "\n")
    return Molecule(xyzstring, basis, nelec)
end
    
function Molecule(xyzstring::String, basis::String, nelec::Int)
    """
    Constructor function
        Takes in formatted xyzstring to initiate Molecule struct
    """
    return Molecule(xyzstring, basis, nelec)
end

function Molecule(xyzstring::String, troyint::String)
    """
    Constructor function
        Reads pre-calculated integrals from the specified file
        and initiate Molecule struct
    """
    input = readdlm(troyint)
    # Line 1 // Nuc-Nuc Repulsion Energy (Hartrees; RHF.jl will read this)
    # Line 2 // Int::nbasis, Int::nocc/2
    # Line 3 // Vector{Float64}::S
    # Line 4 // Vector{Float64}::h
    # Line 5 // Vector{Float64}::ijkl
    return Molecule(xyzstring, troyint, input[2,2]*2)
end
