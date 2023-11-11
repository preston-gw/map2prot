# map2prot
Tool for mapping peptides onto a protein sequence

## About
The map2prot tool is an R script that maps a set of peptides onto a protein sequence. It generates a chart in which the peptides appear as coloured blocks underneath the protein sequence. Sequence coverage (%) and various metadata are printed on the chart.

The inputs are:
* a text file containing the sequences of the peptides;
* a *.fasta file containing the sequence of the protein.

If the file containing the peptide sequences is called 'peptides.txt', it will be interpreted as the [MaxQuant](https://www.maxquant.org/) output file of the same name. If the file is called something else, it will be interpreted as a plain list (as would be uploaded for a [UniProt peptide search](https://www.uniprot.org/peptide-search), for example).

In the case of the *.fasta file, you just need to know the index number of the protein sequence you want to use (e.g., 1 for the first or only sequence, 2 for the second sequence, and so on).

map2prot was written in [R for Windows](https://www.R-project.org/) and requires R package [seqinR](https://cran.r-project.org/web/packages/seqinr/index.html). It was last tested on files generated by MaxQuant version 1.6.0.1, using R version 4.2.1 and seqinR version 4.2-16.

map2prot was adapted from my earlier work on [dependent peptides](https://github.com/preston-gw/dependent-peptides). Please see the [citation metadata](https://github.com/preston-gw/map2prot/blob/main/CITATION.cff) for more information.

## Files
#### map2prot_v1.0.0.R
This is the tool itself. As well as the R code, the file contains usage notes and instructions.  
#### 4f5s_A.fasta
This is an example of a *.fasta file. It contains the sequence of mature bovine serum albumin (A190T variant). The file itself originates from PRIDE project [PXD013040](https://www.ebi.ac.uk/pride/archive/projects/PXD013040), and the sequence within it originates from Protein Data Bank accession [4F5S](https://www.rcsb.org/structure/4f5s). 
#### peptides.txt
This is an example of a MaxQuant 'peptides.txt' file. It can be used together with 4f5s_A.fasta to test that the script is working properly. The file originates from PRIDE project [PXD013040](https://www.ebi.ac.uk/pride/archive/projects/PXD013040).

## Getting started
1. Download the script (**map2prot_v1.0.0.R**) and the example input files (**4f5s_A.fasta** and **peptides.txt**).
2. Follow the instructions in the script header.
3. Review the text output in the R console, checking for any warnings.
4. Review the graphical output. It should look like **example_output.jpg**.

## Acknowledgements
I thank Dr Stuart Warriner (University of Leeds) and Professor Andrew Wilson (University of Birmingham) for supporting this project. My work at the University of Leeds has been funded by the [Engineering and Physical Sciences Research Council](https://www.ukri.org/councils/epsrc/).
