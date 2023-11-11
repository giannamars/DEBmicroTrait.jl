using CSV, DataFrames, Statistics
using JLD

dir                     = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
df_isolates             = CSV.read(joinpath(dir, "files/input/isolates2traits.csv"), DataFrame, missingstring="N/A");
df_metabolites          = CSV.read(joinpath(dir, "files/input/root_exudates.csv"), DataFrame, missingstring="N/A");

assimilation            = load(joinpath(dir, "files/output/isolates_assimilation.jld"));

id_sugars               = findall(df_metabolites.Ontology .== "Sugars")
id_aminos               = findall(df_metabolites.Ontology .== "Amino acids")
id_organics             = findall(df_metabolites.Ontology .== "Organic acids")
id_fattys               = findall(df_metabolites.Ontology .== "Fatty acids")
id_auxins               = findall(df_metabolites.Ontology .== "Auxins")
id_nucleos              = findall(df_metabolites.Ontology .== "Nucleotides")

N_SB_plant = vcat(median(assimilation["NSB"][id_sugars,:], dims=1), median(assimilation["NSB"][id_aminos,:], dims=1),
                    median(assimilation["NSB"][id_organics,:], dims=1), median(assimilation["NSB"][id_fattys,:], dims=1),
                    median(assimilation["NSB"][id_auxins,:], dims=1), median(assimilation["NSB"][id_nucleos,:], dims=1))

K_D_plant = vcat(median(assimilation["KD"][id_sugars,:], dims=1), median(assimilation["KD"][id_aminos,:], dims=1),
                    median(assimilation["KD"][id_organics,:], dims=1), median(assimilation["KD"][id_fattys,:], dims=1),
                    median(assimilation["KD"][id_auxins,:], dims=1), median(assimilation["KD"][id_nucleos,:], dims=1))

y_DE_plant = vcat(median(assimilation["yDE"][id_sugars,:], dims=1), median(assimilation["yDE"][id_aminos,:], dims=1),
                    median(assimilation["yDE"][id_organics,:], dims=1), median(assimilation["yDE"][id_fattys,:], dims=1),
                    median(assimilation["yDE"][id_auxins,:], dims=1), median(assimilation["yDE"][id_nucleos,:], dims=1))

N_C_plant   = vcat([median(assimilation["NC"][id_sugars])], [median(assimilation["NC"][id_aminos])], 
                    [median(assimilation["NC"][id_organics])], [median(assimilation["NC"][id_fattys])],
                    [median(assimilation["NC"][id_auxins])], [median(assimilation["NC"][id_nucleos])])

id_sugars               = [21, 31]  # only sugars detected in week 3
LC_MS_3     = vcat([median(df_metabolites.week3[id_sugars])], [median(df_metabolites.week3[id_aminos])],
                   [median(df_metabolites.week3[id_organics])], [median(df_metabolites.week3[id_fattys])],
                   [median(df_metabolites.week3[id_auxins])], [median(df_metabolites.week3[id_nucleos])]) 

LC_MS_6     = vcat([median(df_metabolites.week6[id_sugars])], [median(df_metabolites.week6[id_aminos])],
                   [median(df_metabolites.week6[id_organics])], [median(df_metabolites.week6[id_fattys])],
                   [median(df_metabolites.week6[id_auxins])], [median(df_metabolites.week6[id_nucleos])])

LC_MS_9     = vcat([median(df_metabolites.week9[id_sugars])], [median(df_metabolites.week9[id_aminos])],
                   [median(df_metabolites.week9[id_organics])], [median(df_metabolites.week9[id_fattys])],
                   [median(df_metabolites.week9[id_auxins])], [median(df_metabolites.week9[id_nucleos])])

LC_MS_12     = vcat([median(df_metabolites.week12[id_sugars])], [median(df_metabolites.week12[id_aminos])],
                   [median(df_metabolites.week12[id_organics])], [median(df_metabolites.week12[id_fattys])],
                   [median(df_metabolites.week12[id_auxins])], [median(df_metabolites.week12[id_nucleos])])

LC_MS        = hcat(LC_MS_3, LC_MS_6, LC_MS_9, LC_MS_12)


y_EM                    = ones(39)
# I/O
save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/output/isolates_assimilation_plant.jld", "NSB", N_SB_plant, "KD", K_D_plant, "yEM", y_EM, "yDE", y_DE_plant, "NC", N_C_plant, "LCMS", LC_MS)





