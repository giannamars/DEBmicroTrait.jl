# DEBmicroTrait

**DEBmicroTrait** is a genome-informed, trait-based dynamic energy budget (DEB) modeling package for microbial systems, implemented in Julia. The framework links genome-derived traits to microbial energy and resource acquisition, allocation, and growth, enabling metabolic trade-offs and ecological strategies to emerge from first principles under varying environmental conditions.

DEBmicroTrait is designed to support **machine learning–enabled analysis of emergent microbial strategies**. Model outputs and genome-derived trait spaces can be analyzed using **unsupervised learning** to identify dominant metabolic and life-history strategies, cluster taxa based on variance in genome-informed traits, and map emergent behaviors across environmental gradients.

# Key Features

**Genome-informed trait integration**
Uses genome-derived traits (via [microTrait](https://github.com/ukaraoz/microtrait), qSIP, iRep) to constrain growth, uptake kinetics, maintenance, enzyme production, and mortality.

**Unsupervised learning of emergent strategies**
Supports clustering, dimensionality reduction, and manifold learning to identify dominant metabolic and life-history strategies emerging from simulations.

Genome-scale trait variance analysis
Enables clustering of taxa based on variance and covariance in genome-derived trait distributions.

Trait-based DEB framework
Explicit representation of energy and resource allocation allows trade-offs to emerge rather than be imposed.

Scalable reaction-network formulation
SUPECA-based kinetics ensure scale invariance across complex substrate–consumer and enzyme-mediated reaction networks.

Flexible experimental configurations
Batch, chemostat, and time-dependent forcing scenarios for arbitrary numbers of substrates, taxa, and enzymes.

Environmental coupling
Optional integration of soil microscale processes and environmental controls, such as temperature, moisture, texture, and mineral matrix composition.



**DEBmicroTrait** is a genome-informed trait-based dynamic energy budget modeling package developed for trait-based microbial modeling in Julia. The dynamic energy budget approach captures interacting strategies for energy and resource acquisition and allocation, letting the shape of trade-offs and trait variation at population or community level emerge in response to environmental conditions.


 In the current version, **DEBmicroTrait** can be built in batch mode, chemostat mode, or with time-dependent substrate forcings for any number of substrates (polymers, monomers), microbes, and enzymes. The model structure is based on synthesizing unit equilibrium chemistry approximation kinetics ([SUPECA](https://gmd.copernicus.org/articles/10/3277/2017/)) which guarantees scale-invariance across an arbitrary number of complex substrate-consumer reactions. Representations of the soil mineral matrix and microscale-environment feedback on process rates (soil temperature, soil saturation and texture) can be added based on sample metadata availability.

**DEBmicroTrait** enables integration of [microTrait](https://github.com/ukaraoz/microtrait) information. **microTrait** is an R package that provides a workflow to extract fitness traits from microbial genome sequences, and serves as a basis to constrain model parameters, e.g., for maximum specific growth rate, substrate uptake kinetics, ribosome biosynthesis potential, basal maintenance and extracellular enzyme synthesis. The choice of specific traits and their representational granularity depend on the research question of interest and can be toggled across a hierarchical trait integration framework spanning genome-derivable resource acquisition, resource use, biophysical/life history, and stress tolerance traits underlying ecological strategies.

**DEBmicroTrait** can directly query trait ontologies for taxa and implement traits from different sources (microTrait, [qSIP](https://github.com/bramstone/qsip), iRep), enabling model benchmarking of emergent growth rates and growth efficiencies with qSIP/iRep data, as well as integration of qSIP-derived mortality rates.

An iterative workflow for trait-based model development using the **DEBmicroTrait** pipeline is outlined in our Review on *Life and Death in the Soil Microbiome* ([Box 2](https://www.nature.com/articles/s41579-022-00695-z)). The scripts on the main branch accompany a case study in the rhizosphere revealing how *Life History Strategies and Niches of Soil Bacteria Emerge from Interacting Thermodynamic, Biophysical, and Metabolic traits*  ([Preprint](https://www.biorxiv.org/content/10.1101/2022.06.29.498137v1.abstract)).

 **microTrait** and **DEBmicroTrait** app [development](https://www.kbase.us/research/pett-ridge-sfa/) and documentation continues on [KBase](https://www.kbase.us) under U.S. Department of Energy contract number DE-577AC02-05CH1123.   
