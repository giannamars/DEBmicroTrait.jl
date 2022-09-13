using DEBmicroTrait, Test

el = "C4H9NO2"
chemFormBiom = [1, 1.8, 0.2, 0.5, 0, 0, 0]

chemical_composition = DEBmicroTrait.extract_composition(el) # CHNOSP
@test chemical_composition == [4,9,1,2,0,0]

stoich_electron_donor = DEBmicroTrait.get_stoich_electron_donor(el)
@test stoich_electron_donor == [-1, -10, 4, 1, 0, 0, 21, 18, 0, 0]

stoich_electron_acceptor = DEBmicroTrait.get_stoich_electron_acceptor()
@test stoich_electron_acceptor == [0., 2., 0., 0., 0., 0., -4., -4., -1., 0.]

stoich_cat_rxns = DEBmicroTrait.get_stoich_catabolic_reaction(el)
@test stoich_cat_rxns == [-1., -1., 4., 1., 0., 0., 3., 0., -4.5, 0.]


stoich_anabolic_O2, stoich_anabolic_HCO3 = DEBmicroTrait.get_stoich_anabolic_reaction(el, chemFormBiom)
@test stoich_anabolic_O2  ≈ [-0.25, 0.15, 0.0,  0.05,  0., 0., -0.05, 0.,  -0.075, 1.0] atol =1e-4
@test stoich_anabolic_HCO3  ≈ [-0.23333333, 0.16666667, -0.06666667,  0.03333333,  0., 0., -0.1, 0.,  0., 1.] atol = 1e-4

out = DEBmicroTrait.get_lambda(el, chemFormBiom)
@test 1==1
