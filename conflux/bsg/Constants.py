from scipy import constants

######################
# Declaring constants

ELECTRON_MASS_MEV = constants.m_e/constants.eV*constants.speed_of_light**2/1E6
ELECTRON_MASS_KEV = ELECTRON_MASS_MEV*1E3
ELECTRON_MASS_EV = ELECTRON_MASS_MEV*1E6
PROTON_MASS_EV = constants.m_p/constants.eV*constants.speed_of_light**2
PROTON_MASS_KEV = PROTON_MASS_EV/1E3
NUCLEON_MASS_EV = (constants.m_p+constants.m_n)/2.0/constants.eV*constants.speed_of_light**2
NUCLEON_MASS_KEV = NUCLEON_MASS_EV/1E3
NUCLEON_MASS_W = NUCLEON_MASS_EV/(1.0*ELECTRON_MASS_EV)
PROTON_MASS_W = PROTON_MASS_EV/(1.0*ELECTRON_MASS_EV)
AMU_MASS_KEV = constants.physical_constants['atomic mass constant energy equivalent in MeV'][0]*1E3
NATURAL_LENGTH = (constants.hbar*constants.speed_of_light/constants.e)/ELECTRON_MASS_EV
ALPHA = constants.alpha

####################
# Nucleon coupling constants

GA = 1.27
GM = 4.706

####################
# Element names

atoms = ["N", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt"]
