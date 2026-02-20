Hello, this pipeline was created for the identification of highly confident horizontal transfer genes (hcHTGs) in a genome by comparing the Codon Usage Bias (CUB) and GC content of each gene with highly conserved and expressed genes, such as ribosomal genes.
The pipeline was developed during the analyzes of Selleri et al 2026, mSystems "Assessment of genome evolution in Bifidobacterium adolescents indicates genetic adaptation to the human gut" (doi.org/10.1128/msystems.01173).
The pipeline is executed automatically by the HGT_analysis.sh script.

REQUIREMENTS: 
A) A folder in the working directory containing multifasta files for every genome (.ffn files, PROKKA output-like format).
B) A conda environment with all dependencies. Required packages can be found in the cub.yml file.
C) The three Python scripts located in the working folder.

PIPELINE OVERVIEW
    1) cub_calculator.py 	-> This script carries out the CUB calculation starting from the .ffn files. It analyzes GC content and CODON USAGE BIAS (CUB) using the following indices:
	- ENC (Effective Number of Codons): A measure of codon bias. Values range from 20 (maximum bias) to 61 (minimum bias/random usage).
	- CAI (Codon Adaptation Index): Measures codon usage compared to reference genes (ribosomal). A higher CAI indicates greater similarity to highly expressed host genes and a lower probability of HGT.
	- RSCU (Relative Synonymous Codon Usage): Codon usage compared to the theoretical maximum.
    The script also evaluates gene quality in the 'Quality_note' column, flagging short genes or non-coding RNAs.

    2) HGT_identification.py 	-> This script identifies potential HGT candidates through statistical analysis of the host genome. It uses indirect criteria such as low CAI, high ENC, high GC Z-score, and high Codon_Avg_Dist_self.

    3) HGT_extractor.py 	-> This script filters the previous output to retain only the predicted HGT genes. A gene is classified as 'HGT_predicted = True' if it meets at least 2 out of the 4 statistical factors (based on the extreme percentiles of the genome). Finally, it creates a single combined file with all the results for further processing.
