# SAINTS-CF 

"Parallel metagenomic- and culture-based approaches show nasal swabs are a good proxy for broncho-alveolar lavage in children with cystic fibrosis"

Code for the research carried out with Teagasc, CHI-Crumlin, and RCSI.

Nominally, code should be run in the following order, but hopefully it would only be becessary to run the fig_manus script to re-create the material in the paper (once you have already run everything once...)


 - `mtu__paedcf__assembly_0.2` - clean, filter, decontaminate, assign taxonomy, and compute k-mer diversity, all the real work

 - `mtu__paedcf__spreadsheeting__shortcuts.R` - sets up the data (created by Kraken2 and NonPareil) for R
 - `mtu__paedcf__chunk__diversity__alpha.R` - alpha diversity
 - `mtu__paedcf__chunk__diversity__beta.R` - beta diversity
 - `mtu__paedcf__chunk__microbiology.R` - checking the pathogen load
 - `mtu__paedcf__chunk__DA-LM-EMM__plots.R` - differential abundance testing
 - `mtu__paedcf__figs__manus_0.3.5__noHsap.R` - presenting findings and tweaking the edges