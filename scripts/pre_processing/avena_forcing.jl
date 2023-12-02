using JLD2

# I/O
dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd();
root_2006               = load(joinpath(dir, "files/avena/root_2006.jld2"));
root_2007               = load(joinpath(dir, "files/avena/root_2007.jld2"));
root_2008               = load(joinpath(dir, "files/avena/root_2008.jld2"));
root_2009               = load(joinpath(dir, "files/avena/root_2009.jld2"));
root_2010               = load(joinpath(dir, "files/avena/root_2010.jld2"));
root_2011               = load(joinpath(dir, "files/avena/root_2011.jld2"));
root_2012               = load(joinpath(dir, "files/avena/root_2012.jld2"));
root_2013               = load(joinpath(dir, "files/avena/root_2013.jld2"));
root_2014               = load(joinpath(dir, "files/avena/root_2014.jld2"));
root_2015               = load(joinpath(dir, "files/avena/root_2015.jld2"));
root_2016               = load(joinpath(dir, "files/avena/root_2016.jld2"));
root_2017               = load(joinpath(dir, "files/avena/root_2017.jld2"));
root_2018               = load(joinpath(dir, "files/avena/root_2018.jld2"));
root_2019               = load(joinpath(dir, "files/avena/root_2019.jld2"));
root_2020               = load(joinpath(dir, "files/avena/root_2020.jld2"));
root_2021               = load(joinpath(dir, "files/avena/root_2021.jld2"));
root_2022               = load(joinpath(dir, "files/avena/root_2022.jld2"));


# root exudation rate 
r_soil                  = 2.3e-3
#r_soil                  = 5e-3
V_soil                  = pi*(r_soil)^3
exudC_2006              = abs.(root_2006["exudC"].*root_2006["rootEC"]/V_soil)
exudC_2007              = abs.(root_2007["exudC"].*root_2007["rootEC"]/V_soil)
exudC_2008              = abs.(root_2008["exudC"].*root_2008["rootEC"]/V_soil)
exudC_2009              = abs.(root_2009["exudC"].*root_2009["rootEC"]/V_soil)
exudC_2010              = abs.(root_2010["exudC"].*root_2010["rootEC"]/V_soil)
exudC_2011              = abs.(root_2011["exudC"].*root_2011["rootEC"]/V_soil)
exudC_2012              = abs.(root_2012["exudC"].*root_2012["rootEC"]/V_soil)
exudC_2013              = abs.(root_2013["exudC"].*root_2013["rootEC"]/V_soil)
exudC_2014              = abs.(root_2014["exudC"].*root_2014["rootEC"]/V_soil)
exudC_2015              = abs.(root_2015["exudC"].*root_2015["rootEC"]/V_soil)
exudC_2016              = abs.(root_2016["exudC"].*root_2016["rootEC"]/V_soil)
exudC_2017              = abs.(root_2017["exudC"].*root_2017["rootEC"]/V_soil)
exudC_2018              = abs.(root_2018["exudC"].*root_2018["rootEC"]/V_soil)
exudC_2019              = abs.(root_2019["exudC"].*root_2019["rootEC"]/V_soil)
exudC_2020              = abs.(root_2020["exudC"].*root_2020["rootEC"]/V_soil)
exudC_2021              = abs.(root_2021["exudC"].*root_2021["rootEC"]/V_soil)
exudC_2022              = abs.(root_2022["exudC"].*root_2022["rootEC"]/V_soil)

#plot(abs.(root_2018["exudC"].*root_2018["rootEC"]*12*24/8e-3))

exudation_rate          = vcat(exudC_2006, exudC_2007, exudC_2008, exudC_2009, exudC_2010, exudC_2011,
                               exudC_2012, exudC_2013, exudC_2014, exudC_2015, exudC_2016, exudC_2017,
                               exudC_2018, exudC_2019, exudC_2020, exudC_2021, exudC_2022);
@save joinpath(dir,"files/avena/avena_exudation_2006_2022.jld2") exudation_rate

# root decay rate
r_soil                  = 5e-3
A_soil                  = pi*(r_soil)^2
L_root                  = 40*1e-2
V_soil                  = L_root*A_soil
rootC_2006              = root_2006["rootVC"]./V_soil
rootC_2007              = root_2007["rootVC"]./V_soil
rootC_2008              = root_2008["rootVC"]./V_soil
rootC_2009              = root_2009["rootVC"]./V_soil
rootC_2010              = root_2010["rootVC"]./V_soil
rootC_2011              = root_2011["rootVC"]./V_soil
rootC_2012              = root_2012["rootVC"]./V_soil
rootC_2013              = root_2013["rootVC"]./V_soil
rootC_2014              = root_2014["rootVC"]./V_soil
rootC_2015              = root_2015["rootVC"]./V_soil
rootC_2016              = root_2016["rootVC"]./V_soil
rootC_2017              = root_2017["rootVC"]./V_soil
rootC_2018              = root_2018["rootVC"]./V_soil
rootC_2019              = root_2019["rootVC"]./V_soil
rootC_2020              = root_2020["rootVC"]./V_soil
rootC_2021              = root_2021["rootVC"]./V_soil
rootC_2022              = root_2022["rootVC"]./V_soil

root_C          = vcat(rootC_2006, rootC_2007, rootC_2008, rootC_2009, rootC_2010, rootC_2011,
                               rootC_2012, rootC_2013, rootC_2014, rootC_2015, rootC_2016, rootC_2017,
                               rootC_2018, rootC_2019, rootC_2020, rootC_2021, rootC_2022);
@save joinpath(dir,"files/avena/avena_root_2006_2022.jld2") root_C
