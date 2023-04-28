using DEBmicroTrait

u_c = [171.0]
l_c = [6.0]
ν_kin = 1e6

Re = DEBmicroTrait.reynolds_number(u_c, l_c, ν_kin)
