"""
QCPractice.Molecule

Submodule to that contains Julia struct Molecule
"""

struct Molecule
    atomnos::Array{String} # atomic symbols of constituent atoms
    atomxyz::Array{Vector} # atomic coordinates in the Cartesian system
    basis::String # basis set
    nelec::Int # number of electrons
    xyzstring::String # xyz string in the format below:
                      # "H 0.00 0.00 0.00\nH 0.70 0.00 0.00"

    function Molecule(atomnos::Array{String}, atomxyz::Array{Vector}, basis::String, nelec::Int)
        """
        Constructor function
            Takes in necessary inputs and creates xyzstring to prepare for libcint calls
        """
        geometry = hcat(atomnos, atomxyz)
        xyzstring = join([join(view(geometry, i, :), " ") for i in 1:size(geometry)[1]], "\n")
        new(atomnos, atomxyz, basis, nelec, xyzstring)
    end

end