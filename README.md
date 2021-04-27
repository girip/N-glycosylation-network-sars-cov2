# N-glycan biosynthesis network for SARS-CoV-2 spike S protein



# Instructions for using the scripts to generate N-Glycosylation pathway
1. Download GNAT toolbox from http://gnatmatlab.sourceforge.net/
2. Replace the db folder in GNAT/examplefiles/ subfolder. 
3. Install GNAT in matlab. Follow instructions for installing GNAT from the instruction manual.
4. Add path to glycan folder provided in this git.
5. Run infer_from_mspec_covid file from scripts folder. Change the path for using different initial glycans.

## Plot glycosylation network
1. Run plot_network.m after infer_from_mspec_covid script to generate plot of full network.
2. To block specific enzyme and replot network, use block_enzyme_and_plot.m script.

## Contructed glycosylation network
Reconstructed network along with enzyme names are in contructed_network folder.
Inferred glycans can be downloaded from inferred_glycans in the glycan folder.
