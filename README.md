# :wave: About
This repository contains the code and data for the paper: "Contextual evidence for categorising interactions and guiding their discovery".


# :page_facing_up: Paper and citing
Abramov, K., Nehoray, S. M., Puzis, R., & Pilosof, S. (2026). **Contextual evidence for categorising interactions and guiding their discovery**. 


# Abstract:
A surge of recent methods tackles incomplete ecological interaction data through link prediction. However, validating predictions means field-sampling an immense number of candidate interactions, possibly with a poorly suited method. This costly endeavour lacks operative guides grounded in strong ecological and statistical theory. We present a guided-sampling framework combining link prediction with within-system variability across replicated networks as contextual evidence, generating a fine-resolution, ecologically-aware link taxonomy. It resolves categories ecologists have long struggled to separate, such as forbidden interactions versus those merely missing from the local sample, pinpoints the most probable unknowns, and directs targeted, cost-efficient discovery. We formalise how confidence in each assignment is quantified and grows as evidence accumulates. Treating validation as evidence accumulation makes it a tractable sampling problem rather than a permanent limit on knowledge.

<p align="center"><img src="results/figures/framework_flow.png" alt="Framework overview" width="600"/></p>

# :bar_chart: Section Bayesian Inference & SI, Section Empirical Demonstration — interactive explorers

Both ship with a standalone interactive explorer, live at <https://lpguide.ecomplab.com/>.

**Bayesian reading of the link taxonomy.** Accompanies the development of Section Bayesian Inference & Supplementary Information. Figure script, rendered figure, the interactive explorer, and a math/assumptions context doc.

**The evidence accumulation framework in action.** Accompanies Section Empirical Demonstration: the framework applied to plant–pollinator networks at six sites: catrgorisation and prediction scripts, rendered alluvial figure, and the interactive alluvial explorer.

# :computer: Code:
Instructions for running the code and reproducing the results are in the repository Wiki under ["Code"](https://github.com/Ecological-Complexity-Lab/LP_validation/wiki/Code).

# :file_cabinet: Data:
Data used for empirical demonstration is detailed in the repository Wiki under ["Data"](https://github.com/Ecological-Complexity-Lab/LP_validation/wiki/Data).
