using DEBmicroTrait, Test

el = "C7H8"
chemFormBiom = [1, 1.8, 0.2, 0.5, 0, 0, 0]

chemical_indices = DEBmicroTrait.extract_composition(el)
a = chemical_indices[1] #C
b = chemical_indices[2] #H
c = chemical_indices[3] #N
d = chemical_indices[4] #O
e = chemical_indices[5] #S
f = chemical_indices[6] #P
z = 0

ne = -z+4*a+b-3*c-2*d+5*e-2*f  # number of electrons transferred in D
nosc = -(ne/a) + 4
γ_c = 4 - nosc


chemical_composition = DEBmicroTrait.extract_composition(el) # CHNOSP

stoich_electron_donor = DEBmicroTrait.get_stoich_electron_donor(el)

stoich_electron_acceptor = DEBmicroTrait.get_stoich_electron_acceptor()

stoich_cat_rxns = DEBmicroTrait.get_stoich_catabolic_reaction(el)

stoich_anabolic_O2, stoich_anabolic_HCO3 = DEBmicroTrait.get_stoich_anabolic_reaction(el, chemFormBiom)

out = DEBmicroTrait.get_lambda(el, chemFormBiom)

bio = out[4][10]

sub = out[4][1]

(bio)/(sub*a)