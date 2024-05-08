using DEBmicroTrait
using Dates
using CSV, DataFrames, Statistics
using JLD2
using OrdinaryDiffEq
using Plots; plotlyjs()

# I/O
dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd();
avena_ex                = load(joinpath(dir, "files/avena/avena_exudation_2006_2022.jld2"));
avena_root              = load(joinpath(dir, "files/avena/avena_root_2006_2022.jld2"));
γ_root                  = 5e-5
ex_rate                 = avena_ex["exudation_rate"]

# ODE
p_force = Forcing(avena_root["root_C"], γ_root, avena_ex["exudation_rate"])
function test!(du,u,p_force,t)
    spl_ex   = DEBmicroTrait.root_exudation(p_force)
    spl_root = DEBmicroTrait.root_decay(p_force)
    #
    du[1] = spl_ex(t)
    du[2] = spl_root(t)
end
u0 = [0.0, 0.0]

# 2018
week_hours = 24*7
month_hours = 365.25/12*24
months = 12
nyears  = 2019-2006
tstop = month_hours*months*nyears
MONTHS = [Dates.LOCALES["english"].months... Dates.LOCALES["english"].months...]
xticks_month=(repeat((0:month_hours:tstop),nyears), repeat(MONTHS[1:months], nyears))
week3start = 8760*(nyears-1)+1+345
start2018 = week3start
# total C 
tspan = [8760*(nyears-1)+1,Int(tstop)]
prob = ODEProblem(test!,u0,tspan,p_force)
sol = solve(prob, Tsit5());
plot(sol, idxs=2, xticks=xticks_month, xrotation=90)
week3C_2018 = sol(week3start+week_hours*3, idxs=1) - sol(week3start, idxs=1);
week6C_2018 = sol(week3start+week_hours*6, idxs=1) - sol(week3start+week_hours*3+1, idxs=1);
week9C_2018 = sol(week3start+week_hours*9, idxs=1) - sol(week3start+week_hours*6+1, idxs=1);
week12C_2018 = sol(week3start+week_hours*12, idxs=1) - sol(week3start+week_hours*9+1, idxs=1);
totalC_2018 = sol(tstop, idxs=1);

week3C_2018_root = sol(week3start+week_hours*3, idxs=2);
week6C_2018_root = sol(week3start+week_hours*6, idxs=2);
week9C_2018_root = sol(week3start+week_hours*9, idxs=2);
week12C_2018_root = sol(week3start+week_hours*12, idxs=2);
totalC_2018_root = sol(tstop, idxs=2);

# rate
#plot(ex_rate[8760*(nyears-1)+1:Int(tstop)], xticks=xticks_month, xrotation=90, ylabel="molC/m3/hr", label="2022")
week3rate_2018 = median(ex_rate[week3start:week3start+week_hours*3])
week6rate_2018 = median(ex_rate[week3start+week_hours*3+1:week3start+week_hours*6])
week9rate_2018 = median(ex_rate[week3start+week_hours*6+1:week3start+week_hours*9])
week12rate_2018 = median(ex_rate[week3start+week_hours*9+1:week3start+week_hours*12])

rateseries_2018 = ex_rate[8760*(nyears-1)+1:Int(tstop)]
totalCseries_2018 = collect(sol(8760*(nyears-1)+1:Int(tstop), idxs=1))
jldsave("files/output/timeseries_avena2018.jld2"; rateseries_2018, totalCseries_2018)

# 2012
week_hours = 24*7
month_hours = 365.25/12*24
months = 12
nyears  = 2013-2006
tstop = month_hours*months*nyears
MONTHS = [Dates.LOCALES["english"].months... Dates.LOCALES["english"].months...]
xticks_month=(repeat((0:month_hours:tstop),nyears), repeat(MONTHS[1:months], nyears))
week3start = 8760*(nyears-1)+1+1007
start2012 = week3start
# total C 
tspan = [8760*(nyears-1)+1,Int(tstop)]
prob = ODEProblem(test!,u0,tspan,p_force)
sol = solve(prob, Tsit5());
#plot(sol, vars=1, xticks=xticks_month, xrotation=90)
week3C_2012 = sol(week3start+week_hours*3, idxs=1) - sol(week3start, idxs=1)
week6C_2012 = sol(week3start+week_hours*6, idxs=1) - sol(week3start+week_hours*3+1, idxs=1)
week9C_2012 = sol(week3start+week_hours*9, idxs=1) - sol(week3start+week_hours*6+1, idxs=1)
week12C_2012 = sol(week3start+week_hours*12, idxs=1) - sol(week3start+week_hours*9+1, idxs=1)
totalC_2012 = sol(tstop, idxs=1);

week3C_2012_root = sol(week3start+week_hours*3, idxs=2);
week6C_2012_root = sol(week3start+week_hours*6, idxs=2);
week9C_2012_root = sol(week3start+week_hours*9, idxs=2);
week12C_2012_root = sol(week3start+week_hours*12, idxs=2);
totalC_2012_root = sol(tstop, idxs=2);

# rate
#plot(ex_rate[8760*(nyears-1)+1:Int(tstop)], xticks=xticks_month, xrotation=90, ylabel="molC/m3/hr", label="2022")
week3rate_2012 = median(ex_rate[week3start:week3start+week_hours*3])
week6rate_2012 = median(ex_rate[week3start+week_hours*3+1:week3start+week_hours*6])
week9rate_2012 = median(ex_rate[week3start+week_hours*6+1:week3start+week_hours*9])
week12rate_2012 = median(ex_rate[week3start+week_hours*9+1:week3start+week_hours*12])

rateseries_2012 = ex_rate[8760*(nyears-1)+1:Int(tstop)]
totalCseries_2012 = collect(sol(8760*(nyears-1)+1:Int(tstop), idxs=1))
jldsave("files/output/timeseries_avena2012.jld2"; rateseries_2012, totalCseries_2012)


# 2013
week_hours = 24*7
month_hours = 365.25/12*24
months = 12
nyears  = 2014-2006
tstop = month_hours*months*nyears
MONTHS = [Dates.LOCALES["english"].months... Dates.LOCALES["english"].months...]
xticks_month=(repeat((0:month_hours:tstop),nyears), repeat(MONTHS[1:months], nyears))
week3start = 8760*(nyears-1)+1+1007
start2013 = week3start
# total C 
tspan = [8760*(nyears-1)+1,Int(tstop)]
prob = ODEProblem(test!,u0,tspan,p_force)
sol = solve(prob, Tsit5());
#plot(sol, vars=1, xticks=xticks_month, xrotation=90)
week3C_2013 = sol(week3start+week_hours*3, idxs=1) - sol(week3start, idxs=1)
week6C_2013 = sol(week3start+week_hours*6, idxs=1) - sol(week3start+week_hours*3+1, idxs=1)
week9C_2013 = sol(week3start+week_hours*9, idxs=1) - sol(week3start+week_hours*6+1, idxs=1)
week12C_2013 = sol(week3start+week_hours*12, idxs=1) - sol(week3start+week_hours*9+1, idxs=1)
totalC_2013 = sol(tstop, idxs=1);

week3C_2013_root = sol(week3start+week_hours*3, idxs=2);
week6C_2013_root = sol(week3start+week_hours*6, idxs=2);
week9C_2013_root = sol(week3start+week_hours*9, idxs=2);
week12C_2013_root = sol(week3start+week_hours*12, idxs=2);
totalC_2013_root = sol(tstop, idxs=2);


# rate
#plot(ex_rate[8760*(nyears-1)+1:Int(tstop)], xticks=xticks_month, xrotation=90, ylabel="molC/m3/hr", label="2022")
week3rate_2013 = median(ex_rate[week3start:week3start+week_hours*3])
week6rate_2013 = median(ex_rate[week3start+week_hours*3+1:week3start+week_hours*6])
week9rate_2013 = median(ex_rate[week3start+week_hours*6+1:week3start+week_hours*9])
week12rate_2013 = median(ex_rate[week3start+week_hours*9+1:week3start+week_hours*12])

rateseries_2013 = ex_rate[8760*(nyears-1)+1:Int(tstop)]
totalCseries_2013 = collect(sol(8760*(nyears-1)+1:Int(tstop), idxs=1))
jldsave("files/output/timeseries_avena2013.jld2"; rateseries_2013, totalCseries_2013)

######## scale exudation rate
p_force = Forcing(avena_root["root_C"], γ_root, avena_ex["exudation_rate"]./10)
# 2006
week_hours = 24*7
month_hours = 365.25/12*24
months = 12
nyears  = 2007-2006
tstop = month_hours*months*nyears
MONTHS = [Dates.LOCALES["english"].months... Dates.LOCALES["english"].months...]
xticks_month=(repeat((0:month_hours:tstop),nyears), repeat(MONTHS[1:months], nyears))
week3start = 8760*(nyears-1)+1+804
start2006 = week3start;
# total C 
tspan = [8760*(nyears-1)+1,Int(tstop)]
prob = ODEProblem(test!,u0,tspan,p_force)
sol = solve(prob, Tsit5());
#plot(sol, vars=1, xticks=xticks_month, xrotation=90)
week3C_2006 = sol(week3start+week_hours*3, idxs=1) - sol(week3start, idxs=1)
week6C_2006 = sol(week3start+week_hours*6, idxs=1) - sol(week3start+week_hours*3+1, idxs=1)
week9C_2006 = sol(week3start+week_hours*9, idxs=1) - sol(week3start+week_hours*6+1, idxs=1)
week12C_2006 = sol(week3start+week_hours*12, idxs=1) - sol(week3start+week_hours*9+1, idxs=1)
totalC_2006 = sol(tstop, idxs=1);

week3C_2006_root = sol(week3start+week_hours*3, idxs=2);
week6C_2006_root = sol(week3start+week_hours*6, idxs=2);
week9C_2006_root = sol(week3start+week_hours*9, idxs=2);
week12C_2006_root = sol(week3start+week_hours*12, idxs=2);
totalC_2006_root = sol(tstop, idxs=2);

# rate
#plot(ex_rate[8760*(nyears-1)+1:Int(tstop)], xticks=xticks_month, xrotation=90, ylabel="molC/m3/hr", label="2022")
week3rate_2006 = median(ex_rate[week3start:week3start+week_hours*3])
week6rate_2006 = median(ex_rate[week3start+week_hours*3+1:week3start+week_hours*6])
week9rate_2006 = median(ex_rate[week3start+week_hours*6+1:week3start+week_hours*9])
week12rate_2006 = median(ex_rate[week3start+week_hours*9+1:week3start+week_hours*12])

rateseries_2006 = ex_rate[8760*(nyears-1)+1:Int(tstop)]
totalCseries_2006 = collect(sol(8760*(nyears-1)+1:Int(tstop), idxs=1))
jldsave("files/output/timeseries_avena2006.jld2"; rateseries_2006, totalCseries_2006)

# I/O

df_out = DataFrame()
df_out.weeks = vcat(repeat(["week3"], 4), repeat(["week6"], 4), repeat(["week9"], 4), repeat(["week12"], 4))
df_out.weekC = vcat(week3C_2018, week3C_2013, week3C_2012, week3C_2006, week6C_2018, week6C_2013, week6C_2012, week6C_2006, week9C_2018, week9C_2013, week9C_2012, week9C_2006, week12C_2018, week12C_2013, week12C_2012, week12C_2006)
df_out.weekCroot = vcat(week3C_2018_root, week3C_2013_root, week3C_2012_root, week3C_2006_root, week6C_2018_root, week6C_2013_root, week6C_2012_root, week6C_2006_root, week9C_2018_root, week9C_2013_root, week9C_2012_root, week9C_2006_root, week12C_2018_root, week12C_2013_root, week12C_2012_root, week12C_2006_root)
df_out.rateC = vcat(week3rate_2018, week3rate_2013, week3rate_2012, week3rate_2006, week6rate_2018, week6rate_2013, week6rate_2012, week6rate_2006, week9rate_2018, week9rate_2013, week9rate_2012, week9rate_2006, week12rate_2018, week12rate_2013, week12rate_2012, week12rate_2006)
df_out.years = vcat(repeat(["2018"], 4), repeat(["2013"], 4), repeat(["2012"], 4), repeat(["2006"], 4))
df_out.totalC = vcat(repeat([totalC_2018], 4), repeat([totalC_2013], 4), repeat([totalC_2012], 4), repeat([totalC_2006], 4))
df_out.starthr =  vcat(repeat([start2018], 4), repeat([start2013], 4), repeat([start2012], 4), repeat([start2006], 4))
df_out.deadhr = vcat(repeat([4384], 4), repeat([3740], 4), repeat([4578], 4), repeat([4214], 4))
CSV.write(joinpath(dir, "files/output/stats_avena_2018_2013_2012_2006.csv"), df_out)

# Extend times



#= week3C = zeros(17);
week6C = zeros(17);
week9C = zeros(17);
week12C = zeros(17);
totalC = zeros(17);
week3rate = zeros(17);
week3rateday = zeros(17);
week3ratenight = zeros(17);
week6rate = zeros(17);
week6rateday = zeros(17);
week6ratenight = zeros(17);
week9rate = zeros(17);
week9rateday = zeros(17);
week9ratenight = zeros(17);
week12rate = zeros(17);
week12rateday = zeros(17);
week12ratenight = zeros(17);

for i in 1:17
    nyears = 2006+i-2006
    tstop = month_hours*months*nyears
    tspan = [8760*(nyears-1)+1,Int(tstop)]
    prob = ODEProblem(test!,u0,tspan,p_force)
    sol = solve(prob, Tsit5());
    #
    week3start = 8760*(nyears-1)+1+345
    #week3start = 8760*(nyears-1)+1+week_hours*3
    week3C[i] = sol(week3start+week_hours*3, idxs=1) - sol(week3start, idxs=1);
    week6C[i] = sol(week3start+week_hours*6, idxs=1) - sol(week3start+week_hours*3+1, idxs=1);
    week9C[i] = sol(week3start+week_hours*9, idxs=1) - sol(week3start+week_hours*6+1, idxs=1);
    week12C[i] = sol(week3start+week_hours*12, idxs=1) - sol(week3start+week_hours*9+1, idxs=1);
    totalC[i] = sol(tstop, idxs=1);
    #
    week3rate[i] = median(ex_rate[week3start:week3start+week_hours*3])
    week3rateday[i] = median(ex_rate[week3start:24:week3start+week_hours*3])
    week3ratenight[i] = median(ex_rate[week3start+12:24:week3start+week_hours*3])
    week6rate[i] = median(ex_rate[week3start+week_hours*3+1:week3start+week_hours*6])
    week6rateday[i] = median(ex_rate[week3start+week_hours*3+1:24:week3start+week_hours*6])
    week6rateday[i] = median(ex_rate[week3start+week_hours*3+1+12:24:week3start+week_hours*6])
    week9rate[i] = median(ex_rate[week3start+week_hours*6+1:week3start+week_hours*9])
    week9rateday[i] = median(ex_rate[week3start+week_hours*6+1:24:week3start+week_hours*9])
    week9ratenight[i] = median(ex_rate[week3start+week_hours*6+1+12:24:week3start+week_hours*9])
    week12rate[i] = median(ex_rate[week3start+week_hours*9+1:week3start+week_hours*12])
    week12rateday[i] = median(ex_rate[week3start+week_hours*9+1:24:week3start+week_hours*12])
    week12ratenight[i] = median(ex_rate[week3start+week_hours*9+1+12:24:week3start+week_hours*12])
end
 =#






#= plot(sol, vars=1, xticks=xticks_month, xrotation=90)


week3start = 8760*(nyears-1)+1+345
week3C = sol(week3start+week_hours*3, idxs=1) - sol(week3start, idxs=1);
week6C = sol(week3start+week_hours*6, idxs=1) - sol(week3start+week_hours*3+1, idxs=1);
week9C = sol(week3start+week_hours*9, idxs=1) - sol(week3start+week_hours*6+1, idxs=1);
week12C = sol(week3start+week_hours*12, idxs=1) - sol(week3start+week_hours*9+1, idxs=1);
totalC = sol(tstop, idxs=1);

week3rate = median(ex_rate[week3start:week3start+week_hours*3])
week3rateday = median(ex_rate[week3start:24:week3start+week_hours*3])
week3ratenight = median(ex_rate[week3start+12:24:week3start+week_hours*3])
week6rate = median(ex_rate[week3start+week_hours*3+1:week3start+week_hours*6])
week6rateday = median(ex_rate[week3start+week_hours*3+1:24:week3start+week_hours*6])
week6rateday = median(ex_rate[week3start+week_hours*3+1+12:24:week3start+week_hours*6])
week9rate = median(ex_rate[week3start+week_hours*6+1:week3start+week_hours*9])
week9rateday = median(ex_rate[week3start+week_hours*6+1:24:week3start+week_hours*9])
week9ratenight = median(ex_rate[week3start+week_hours*6+1+12:24:week3start+week_hours*9])
week12rate = median(ex_rate[week3start+week_hours*9+1:week3start+week_hours*12])
week12rateday = median(ex_rate[week3start+week_hours*9+1:24:week3start+week_hours*12])
week12ratenight = median(ex_rate[week3start+week_hours*9+1+12:24:week3start+week_hours*12])


# Benchmark 
avena_TOC               = CSV.read(joinpath(dir, "files/avena/TOC_data_Hopland_Avena.csv"), DataFrame, missingstring="N/A")
avena_TOC_mM            = avena_TOC.TOC./12.011

median(avena_TOC_mM[13:end])
median(avena_TOC_mM[1:4])
median(avena_TOC_mM[5:8])
median(avena_TOC_mM[9:12])

avena_TOC_mM[9]

1056429+168*3

8760*(nyears-1)+1+month_hours =#