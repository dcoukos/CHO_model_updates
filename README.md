# CHO_network
This project was done when I was working at the Center for Molecular Immunology in Havana, Cuba, in 2018-2019. Jorge Fernandez-de-Cossio-Diaz, Prof. Kalet Léon, and Robert Mulet had just published a paper which described the production dynamics of CHO cells in bioreactors, based on the flow in the bioreactor. This model was an important advance, since the pharmaceuticals produced by the Center for Molecular Immunology were produced in bioreactors, and this work promised to make production in continous-flow bioreactors much more reliable. You can find their work [here](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005835).

This work relies on the principle of mass-balance analysis, but uses average parameters for CHO cell proteins. The idea behind my project was that by using precise biochemical data for the proteins in CHO cells, the model could be optimized to make manufacture even more reliable.

Through the use of custom, appropriate scripts, this work communicated with BRENDA16, BiGG17, the Chemical Translation Service18, and the new addendum to the genes category in KEGG19 to collect such data.

This repo is organized in a flat hierarchy with files representing individual runfiles which collect data or process it. This project taught me a lot about I/O in python, making scripts that run on multi CPU cores, web scraping, and interacting with databases using REST protocols from python.
