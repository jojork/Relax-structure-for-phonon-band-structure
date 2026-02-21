from ase.io import read, write
from ase.optimize import BFGS
from ase.filters import UnitCellFilter
from gpaw import GPAW,PW, FermiDirac

# Read structure
atoms = read(" put your vasp or xyz structure file to call")

print("Structure:", atoms.get_chemical_formula())

# ---------- Set magnetic moments ----------
magmoms = []
for atom in atoms:
    if atom.symbol == "Fe":
        magmoms.append(2.5)   # typical Fe moment
    else:
        magmoms.append(0.0)

atoms.set_initial_magnetic_moments(magmoms)

# ---------- GPAW Calculator ----------
calc = GPAW(
    mode=PW(700),
    xc="PBE",
    kpts=(12, 12, 1),   # Use (12,12,1) if 2D
    occupations=FermiDirac(0.1),
    spinpol=True,
    symmetry="off",
    txt="relax.txt",
    convergence={'energy': 1e-6}
)

atoms.calc = calc

# ---------- Relax BOTH atoms + cell ----------
ucf = UnitCellFilter(atoms)

print("\nStarting full relaxation (cell + atoms)...")

opt = BFGS(ucf, trajectory="relax.traj", logfile="relax.log")
opt.run(fmax=0.03)

print("\nRelaxation complete!")

# ---------- Save relaxed structure ----------
write("relaxed.vasp", atoms, vasp5=True)
print("Relaxed structure saved as relaxed.vasp")
