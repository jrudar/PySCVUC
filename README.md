# PySCVUC v1.1
A Python3 implementation of the SCVUC amplicon sequencing pipeline which is maintained by Dr. Teresita M. Porter.

Based upon the SCVUC repository by Dr. Teresita M. Porter (https://github.com/EcoBiomics-Zoobiome/SCVUC_COI_metabarcode_pipeline and https://github.com/terrimporter/SCVUS_18S_metabarcode_pipeline), PySCVUC uses the same programs, algorithms, and reference dataset outlined in the aformentioned work. The goal of this project is to unify the aformentioned pipeline under a common but modular code base. For further information, please refer to the link.

This pipeline denoises and classifies COI and 18S metabarcode data produced using the Illumina MiSeq platform. It currently supports only paired end reads, however support for single end reads will be included in the future. Ultimatly, this pipeline produces a set of Exact Sequence Variants (ESVs) which are then identified down to the species rank using the RDP Classifier and the COI and 18S V3.2 Training Set (Porter & Hajibabaei, 2018 Sci Rep). This pipeline can also automatically demultiplex (the "primer" in this case needs to include the unique tag). In order to use this pipeline SeqPrep, CutAdapt 2.1, USearch 10.0.240, VSearch 2.8.2, and the RDP Classifier 2.12 must be installed and the location of each of the programs must be in your path. In addition, this pipeline has been tested with with Python 3.6+ and the following Python dependencies are needed: NumPy, Pandas, BioPython, and Scipy. 

## Install:
1) Clone this repository to your home directory.
2) Install the above programs.
3) Install the above Python libraries.
4) Edit your .bashrc file and export the paths to external dependencies (USearch, VSearch).
5) Ensure that you export the path to the location of the PySCVUC directory and make the script executable.
6) Install the following Python libraries: NumPy, SciPy, Pandas, and BioPython
7) Enter the PySCVUC directory and run ./configure.sh.
8) The latest version of the COI Training Files should be downloaded. Extract these files to the 'PySCVUC/training_files/COI/' and 'PySCVUC/training_files/18S/' directories.

## Example Usage: 
python3 /home/user/PySCVUC/PySCVUC.py --primers GTCAGAGGTTCGAAGGCG-GATCCTTCCGCAGGTTCACC-AGGGCAAKYCTGGTGCCAGC-GRCGGTATCTRATCGYCTT --amplicon_names 18SEukF-18SEukR-Uni18SF-Uni18SR  --indices 1-2 --input_dir /home/user/Analysis/18S/Data --results_dir /home/user/Analysis/18S --pyscvuc_dir /home/user --rdpclf_dir /home/user/bin/RDPTools/classifier.jar --classifier 18S 

## Command Line Options: 
##### --primers PRIMERS

````Required. This option specifies the primers which will be used. Each primer should be separated by a dash. This pipeline automatically reverse complements the reverse primer.````

````Example Usage: --primers ATCGATCG-ATCGATCG````

##### --amplicon_names AMPLICON_NAMES

````Required. This option specifies the name for each amplicon. Each name should be separated by a dash.````

````Example Usage: --amplicon_names ForwardPrimer-ReversePrimer````

##### --n N

````Optional. Default is 3. This option controls the number of mismatches that CutAdapt will allow.````

````Example Usage: --n 3````

##### --classifier CLASSIFIER

````Optional. Default is COI. The name of the training set which will be used. Currently, only the V3 COI Training Set is hard-coded. However, adding additonal training sets is very easy.````

````Example Usage: --classifier COI````

##### --threads THREADS

````Optional. Default is 10. This option controls the number of threads which will be used for SeqPrep, CutAdapt, and VSearch.````

````Example Usage: --threads 10````

##### --pe_names PE_NAMES

````Optional. Default is L001_R1_001-L001_R2_001. This option specifies the string which will identify each pair of reads.````

````Example Usage: ----pe_names L001_R1_001-L001_R2_001````

##### --indices INDICES

````Required. This option identifies the sample name from the file name of the fastq.gz files. For example, if the fastq.gz files are: COI-Sample-ID_S1_L001_R1_001.fastq.gz and COI-Sample-ID_S1_L001_R2_001.fastq.gz, the sample index would be 1-2. These correspond to the words Sample and ID. Currently, file naming conventions are very strict and all input files must end with L001_R1_001.fastq.gz or L001_R2_001.fastq.gz.````

````Example Usage: --indices 1-2````

##### --input_dir INPUT_DIR

````Required. This option specifies the path to the directory which contains the fastq.gz files.````

````Example Usage: /home/data_set/fastq````

##### --results_dir RESULTS_DIR

````Required. This option specifies the path to the directory which will store the results of the taxonomic assignment.````

````Example Usage: /home/data_set/results````

##### --PySCVUC_dir PYSCVUC_DIR

````Reqiored. The path to the directory containing the PySCVUC diretory.````
                        
````Example Usage: --PySCVUC_Path /home/scripts/PySCVUC````

##### --RDPClf_dir RDPCLF_DIR

````Required. The path to the directory containing the RDP Classifier.````
                        
````Example Usage: --RDPClf_Path /home/RDPTools/classifier.jar````

## Updated in v1.1a:

- Support For the 18S Training Set (https://github.com/terrimporter/18SClassifier)
- Numeric Columns Output Correctly
- Multiple Processes for SeqPrep and CutAdapt Now Created using the Python 3.6 Multiprocessing Module

## Future Work:

- Support Single End Reads
- Support Multiple Training Sets (16S, ITS, etc)
- Clean Up and Improve Readability of the Code

## References:

Porter, T. M., & Hajibabaei, M. (2018). Automated high throughput animal CO1 metabarcode classification. Scientific Reports, 8, 4226. Available from https://github.com/terrimporter/CO1Classifier

Edgar, R. C. (2016). UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing. BioRxiv. doi:10.1101/081257

St. John, J. (2016, Downloaded). SeqPrep. Retrieved from https://github.com/jstjohn/SeqPrep/releases

Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Applied and Environmental Microbiology, 73(16), 5261–5267. doi:10.1128/AEM.00062-07

Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584. doi:10.7717/peerj.2584

Tange, O. (2011). GNU Parallel - The Command-Line Power Tool. ;;Login: The USENIX Magazine, February, 42–47.

