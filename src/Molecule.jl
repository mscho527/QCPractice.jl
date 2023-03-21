#= QCPractice.Molecule
   Submodule to that contains Julia struct Molecule =#

struct Molecule
    xyzstring::String # xyz string in the format below:
    # "H 0.00 0.00 0.00\nH 0.70 0.00 0.00" (unit = angstroms = 1.889 bohr)
    basis::String # basis set
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

