# DEBmicroTrait

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://giannamars.github.io/DEBmicroTrait.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://giannamars.github.io/DEBmicroTrait.jl/dev/)


**DEBmicroTrait** is a genome-informed trait-based dynamic energy budget modeling package developed for trait-based microbial modeling in Julia. The dynamic energy budget approach captures interacting strategies for energy and resource acquisition and allocation, letting the shape of trade-offs and trait variation at population or community level emerge in response to environmental conditions.

<img src="/files/EESA21-042.png" width="150">

 In the current version, **DEBmicroTrait** can be built in batch mode, chemostat mode, or with time-dependent substrate forcings for any number of substrates (polymers, monomers), microbes, and enzymes. Representations of the soil mineral matrix and microscale-environment feedback on process rates (soil temperature, soil saturation and texture) can be added based on sample metadata availability.

**DEBmicroTrait** enables integration of [microTrait](https://github.com/ukaraoz/microtrait) information. **microTrait** is an R package that provides a workflow to extract fitness traits from microbial genome sequences, and serves as a basis to constrain model parameters for maximum specific growth rate, substrate uptake kinetics, ribosome biosynthesis potential and extracellular enzyme synthesis. The model structure is based on synthesizing unit equilibrium chemistry approximation kinetics which guarantees scale-invariance across an arbitrary number of complex substrate-consumer reactions.

**DEBmicroTrait** can directly query trait ontologies for taxa and implement traits from different sources (microTrait, qSIP, iRep), enabling model benchmarking of emergent growth rates and growth efficiencies with qSIP/iRep data, as well as integration of qSIP-derived mortality rates.

An iterative workflow for trait-based model development using the **DEBmicroTrait** pipeline is outlined in our Review on *Life and Death in the Soil Microbiome* ([Box 2](https://www.nature.com/articles/s41579-022-00695-z)). The scripts on the main branch accompany a case study in the rhizosphere of a Mediterranean grassland showing how *Life History Strategies and Niches of Soil Bacteria Emerge from Interacting Thermodynamic, Biophysical, and Metabolic traits*  ([Preprint](https://www.biorxiv.org/content/10.1101/2022.06.29.498137v1.abstract)).

 **microTrait** and **DEBmicroTrait** app development and documentation continues on [KBase](https://www.kbase.us/) under U.S. Department of Energy contract number DE-577AC02-05CH1123.   
