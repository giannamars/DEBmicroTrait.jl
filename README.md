# DEBmicroTrait
## Genome-informed, trait-based dynamic energy budget (DEB) modeling of microbes, with support for unsupervised learning and ML workflows to analyze emergent metabolic and ecological strategies.


**DEBmicroTrait** is a genome-informed, trait-based dynamic energy budget (DEB) modeling package for microbial systems, implemented in Julia. The framework links genome-derived traits to microbial energy and resource acquisition, allocation, and growth, enabling metabolic trade-offs and ecological strategies to emerge from first principles under varying environmental conditions.

DEBmicroTrait is designed to support **machine learning–enabled analysis of emergent microbial strategies**. Model outputs and genome-derived trait spaces can be analyzed using **unsupervised learning** to identify dominant metabolic and life-history strategies, cluster taxa based on variance in genome-informed traits, and map emergent behaviors across environmental gradients.

## Key Features

**Genome-informed trait integration**
- Uses genome-derived traits (via [microTrait](https://github.com/ukaraoz/microtrait), qSIP, iRep) to constrain growth, uptake kinetics, maintenance, enzyme production, and mortality.

**Unsupervised learning of emergent strategies**
- Supports clustering, dimensionality reduction, and manifold learning to identify dominant metabolic and life-history strategies emerging from simulations.

**Genome-scale trait variance analysis**
- Enables clustering of taxa based on variance and covariance in genome-derived trait distributions.

**Trait-based DEB framework**
- Explicit representation of energy and resource allocation allows trade-offs to emerge rather than be imposed.

**Scalable reaction-network formulation**
- SUPECA-based kinetics ensure scale invariance across complex substrate–consumer and enzyme-mediated reaction networks.

**Flexible experimental configurations**
- Batch, chemostat, and time-dependent forcing scenarios for arbitrary numbers of substrates, taxa, and enzymes.

**Environmental coupling**
- Optional integration of soil microscale processes and environmental controls, such as temperature, moisture, texture, and mineral matrix composition.

 **microTrait** and **DEBmicroTrait** app [development](https://www.kbase.us/research/pett-ridge-sfa/) and documentation continues on [KBase](https://www.kbase.us) under U.S. Department of Energy contract number DE-577AC02-05CH1123.   
