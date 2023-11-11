abstract type AbstractForcing end

struct Forcing{RootC,gRoot,ExC} <: AbstractForcing
    root_C::RootC      # Root carbon [C-mol/m3]
    γ_root::gRoot      # First-order root turnover rate [1/hr]
    exud_C::ExC        # Root exudation rate [mol/m3]
end

for fn in fieldnames(Forcing)
    @eval $fn(p::Forcing) = p.$fn
end

root_exudation(p::AbstractForcing) = begin
    x = collect(1:size(exud_C(p),1))
    y = exud_C(p)
    spl_exudC = Spline1D(x, y)
    return spl_exudC
end

root_decay(p::AbstractForcing) = begin
    x = collect(1:size(root_C(p),1))
    y = root_C(p).*γ_root(p)
    spl_rootC = Spline1D(x, y)
    return spl_rootC
end