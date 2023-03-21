using GaussianBasis # Julia bindings to libcint
include("Molecule.jl")

#using PyCall # Allows calls to Python (for pyscf bindings)
#@pyimport pyscf

#= QCPractice.Inputs

Submodule to manage IO operations and bindings to other packages

# Note
PyCall should be rebuilt with Conda activated if `pyscf` was installed within a Conda environment.
Ideally, submodules for each package should be separate (to simplify the dependencies)
=#

### GaussianBasis.jl Wrappers
### Functions below simply provide interface to GaussianBasis.jl
### (Julia libcint bindings)

function libcint_S(mol::Molecule) # Overlap Integrals
    bset = BasisSet(mol.basis, mol.xyzstring)
    return overlap(bset)
end

function libcint_T(mol::Molecule) # Kinetic Integrals
    bset = BasisSet(mol.basis, mol.xyzstring)
    return kinetic(bset)
end

function libcint_V(mol::Molecule) # Nuclear Attraction Integrals
    bset = BasisSet(mol.basis, mol.xyzstring)
    return nuclear(bset)
end

function libcint_ijkl(mol::Molecule, sparse::Bool=false) # (ij|kl) 4-Centered ERI
    bset = BasisSet(mol.basis, mol.xyzstring)
    sparse ? (return sparseERI_2e4c(bset)) : (return ERI_2e4c(bset))
end

function libcint_df(mol::Molecule, auxbasis::String) # (ij|K) 3-Centered ERI from Density Fitting
    bset = BasisSet(mol.basis, mol.xyzstring)
    return ERI_2e3c(bset, auxbasis)
end

function libcint_ij(mol::Molecule) # (i|j) 2-Centered ERI
    bset = BasisSet(mol.basis, mol.xyzstring)
    return ERI_2e2c(bset)
end
