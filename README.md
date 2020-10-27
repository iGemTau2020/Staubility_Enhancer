# Staubility_Enhancer
Welcome to the Staubility Enhancer! We hope that using our product, you will be able to easily design synthetic sequences with higher stability and expression level.  A key challenge in the field of synthetic biology is genomic instability of introduced genes. Once a gene is inserted into a host organism, it can cause an additional metabolic load, significantly reducing host's fitness. Mutations that damage the introduced gene are therefore likely to be selected for, diminishing its expression. These mutations could render synthetic-biology products obsolete and require constant maintenance. 

We propose interlocking a target gene to the N-terminus of an essential gene in the host’s genome, under the same promoter. This way, mutations on the target gene are likely to affect the expression of the essential gene, leading to the mutated host's mortality. 

We developed a software called the sTAUbility Enhancer, presented here. Given a target gene, it will provide you with the best-fitting conjugated gene candidates. After you select the gene and supply your construct’s purpose, the software will provide you with the best linker candidates. After selection of a linker, the software will scan the construct for mutational and epigenetic hotspots. Using these sites, it will optimize the combined construct for efficient gene expression and increased stability.

The Enhancer directory contains the compiled product for our users, while the other directories contain the models, code and data used in its construction.


Unfortunately, our code combines many frameworks, some of whom only work on certain operating systems, and some of whom require highly nontrivial installations and dependencies. Thus, running our code as a script is extremely tricky. This is part of the reason we built the GUI framework, to avoid having our users go through this process.

The code is still presented, and many internal functions from it can be scavenged for future projects. For any help, feel free to contact us at igem.tau.2019@gmail.com.
