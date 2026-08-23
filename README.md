# :wave: About
This repository contains the code and data for the paper: "Contextual evidence for categorising interactions and guiding their discovery".


# :page_facing_up: Paper and citing
Abramov, K., Nehoray, S. M., Puzis, R., Mucha, P. J., & Pilosof, S. (2026). **Contextual evidence for categorising interactions and guiding their discovery**. EcoEvoRxiv. [https://doi.org/10.32942/X2JD60](https://doi.org/10.32942/X2JD60)


# Abstract:
Ecological interactions shape community dynamics and stability, yet many go unrecorded. Link prediction methods attempt to tackle incomplete knowledge, but predictions are only conjectures, and validating them all in the field is neither feasible nor desirable. Furthermore, interactions differ in the action needed for their detection. Validation therefore must be guided by ecological and statistical theory. We present a guided-sampling framework that combines predictions with local and regional observations, using within-system variability across replicated networks as contextual evidence. It assigns each link to a fine-resolution, ecologically-aware taxonomy that separates conflated categories (e.g., forbidden versus undersampled interactions). Each category then implies a concrete action (e.g., sample more, change method). A Bayesian formulation quantifies confidence in each assignment and incorporates researchers' knowledge. We include a worked empirical example, and an interactive Bayesian version online ([http://lpguide.ecomplab.com](http://lpguide.ecomplab.com)). Contextual evidence turns link prediction from a source of untested hypotheses into a plan for fieldwork.

<p align="center"><img src="results/figures/framework_flow.png" alt="Framework overview" width="600"/></p>

# :bar_chart: Section Bayesian Inference & SI, Section Empirical Demonstration — interactive explorers

Both ship with a standalone interactive explorer, live at <https://lpguide.ecomplab.com/>.

**Bayesian reading of the link taxonomy.** Accompanies the development of Section Bayesian Inference & Supplementary Information. Figure script, rendered figure, the interactive explorer, and a math/assumptions context doc.

**The evidence accumulation framework in action.** Accompanies Section Empirical Demonstration: the framework applied to plant–pollinator networks at six sites: catrgorisation and prediction scripts, rendered alluvial figure, and the interactive alluvial explorer.

# :computer: Code:
Instructions for running the code and reproducing the results are in the repository Wiki under ["Code"](https://github.com/Ecological-Complexity-Lab/LP_validation/wiki/Code).

# :file_cabinet: Data:
Data used for empirical demonstration is detailed in the repository Wiki under ["Data"](https://github.com/Ecological-Complexity-Lab/LP_validation/wiki/Data).
