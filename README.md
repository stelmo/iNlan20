# iNlan20
 Reconstruction and analysis of a genome-scale metabolic model of *Neocallimastix lanati*. 
 Memote snapshot of the model can be found here: [iNlan20](https://gitcdn.link/repo/stelmo/iNlan20/main/iNlan20_memote_report.html).
 
 ### Introduction
Anaerobic gut fungi, from the phylum Neocallimastigomycota, are a clade of early diverging, obligately anaerobic fungi that specialize in lignocellulosic plant biomass degradation. These fungi are primarily found in the stomachs of large herbivores where they play a crucial role in the digestive system of their hosts. Their primary carbon source is lignocellulose, which is an energy rich, abundant and renewable biopolymer that is under-utilized due to the recalcitrant properties of lignin that slow down decomposition. In contrast to most industrially used microbes, anaerobic gut fungi excel at decomposing lignocellulose. By leveraging the high-quality genome of *N. lanati*, the first genome-scale metabolic model of an anaerobic gut fungus is introduced here. The model captures key aspects of its primary metabolism and reveals new metabolic capabilities present in the energy generating organelles of anaerobic gut fungi. Further, coupling the model with flux balance analysis provides a route to predict growth rates and intra-cellular metabolic fluxes in silico, which are experimentally validated using 13C tracer experiments. These simulations can be used to systematically find metabolic connection points for model-based consortia design, as well as highlight understudied aspects of gut fungal metabolism that require further investigation.
 
![N. lanati](/Miscellaneous/nlanimg.jpg)

For more information see: 
1) [Harnessing Nature's Anaerobes for Biotechnology and Bioprocessing](https://doi.org/10.1146/annurev-chembioeng-060718-030340)
2) [The importance of sourcing enzymes from non-conventional fungi for metabolic engineering and biomass breakdown](https://doi.org/10.1016/j.ymben.2017.09.008)
3) [Early-branching gut fungi possess a large, comprehensive array of biomass-degrading enzymes](https://doi.org/10.1126/science.aad1431)

### Outline of the repo
1) The model is supplied in two formats, SBML (.xml) and cobrapy compatible JSON (.json), `iNlan20.xml` and `iNlan20.json` respectively.
2) The Memote report of the model is supplied as `iNlan20_memote_report.html`.
3) The entire reconstruction pipeline is supplied in `/ReconstructionAndAnalysis/`, which is further explained below.
4) The omics data (genomics, transcriptomics and reference databases) are supplied in `/OmicsData/`.
5) The bidirectional blast annotation results for *N. lanati* are supplied in `/MetabolicTables/`.
6) All the other folders have self explanatory names.

### Reconstruction and Analysis
The model can be reconstructed from scratch using a combination of manually curated files and scripts that are collected in `/ReconstructionAndAnalysis/`. To create the basic model run the notebook `/ReconstructionAndAnalysis/Reconstruction of N. lanati GEM.ipynb`. This requires at least version 1.5 of the [Julia language](https://julialang.org/), as well as the packages: JSON, BioSequences, ExcelReaders, DataValues, Statistics, FASTX, StatsBase and DataFrames. Once the basic model has been built, it is necessary to open the resulting `iNlan20.json` file in cobrapy and save it as both a JSON and SBML model, please use the notebook `/ReconstructionAndAnalysis/Write model.ipynb` for this. 

After the preceding steps the model can be run. We supply example scripts: `/ReconstructionAndAnalysis/iNlan20 Example Simulations.ipynb` to make this easy. To use this script Python 3.8.3 is required, the cobrapy package, pandas, json, as well as Gurobi's linear solver.

### Notes
Please note that each fungus is referred to by a shorthand code in the `OmicsData` section. This is described here:
 
| Fungus name | Short name |
| --- | --- |
|Neocallimastix lanati | Neosp3 |
|Neocallimastix californiae | Neosp1 |
|Piromyces finnis | Pirfi3 |
|Anaeromyces robustus | Anasp1 |
|Piromyces sp. E2 | PirE2 |
