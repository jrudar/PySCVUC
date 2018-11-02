# PySCVUC
A Python3 implementation of the SCVUC amplicon sequencing pipeline which is maintained by Dr. Teresita M. Porter.

Based upon the repository by Dr. Teresita M. Porter (https://github.com/EcoBiomics-Zoobiome/SCVUC_COI_metabarcode_pipeline), PySCVUC uses the same programs, algorithms, and reference dataset outlined in the aformentioned work. The goal of this project is to unify the aformentioned pipeline under a common but modular code base. For further information, please refer to the link. At the moment, this pipeline only works for paired end reads. An option will be available in the future to make use of single end reads.

Denoising and classification pipeline for paired amplicon sequencing data. In order to use this pipeline SeqPrep, CutAdapt, USearch 10.0.240, VSearch 2.8.2, and the RDP Tools must be installed and the location of each of the programsmust be in your path. In addition, this pipeline is compatable with Python 3.4+ and the following Python dependencies are needed: NumPy, Pandas, BioPython, and Scipy. 

# Install:
1) Clone this repository to your home directory.
2) Download and install the above programs.
3) Download and install the above Python libraries.
4) Edit your .bashrc file and export the paths to external dependencies (USearch, VSearch). The path to the RDP Tools .jar files should be as follows: export CLASSPATH=<PATH TO RDP TOOLS>/:$CLASSPATH
5) Ensure that you export the path to the location of the PySCVUC directory and make the script executable.
6) Install the following Python libraries: NumPy, SciPy, Pandas, and BioPython

# Command Line Options: 
--primers PRIMERS

*Required. This option specifies the primers which will be used. Each primer should be separated by a dash. Currently the reverse primer is reverse complimented.

*Example Usage: --primers ATCGATCG-ATCGATCG

--amplicon_names AMPLICON_NAMES

*Required. This option specifies the name for each amplicon. Each name should be separated by a dash.

*Example Usage: --amplicon_names ForwardPrimer-ReversePrimer

--n N

*Optional. Default is 3. This option controls the number of mismatches that CutAdapt will allow.

*Example Usage: --n 3

--classifier CLASSIFIER

*Optional. Default is V3. The name of the training set which will be used. Currently, only the V3 COI Training Set is hard-coded. However, adding additonal training sets is very easy.

*Example Usage: --classifier V3

--threads THREADS

*Optional. Default is 10. This option controls the number of threads which will be used for SeqPrep, CutAdapt, and VSearch.

*Example Usage: --threads 10

--indices INDICES

*Required. This option identifies the sample name from the file name of the fastq.gz files. For example, if the fastq.gz files are: COI-Sample-ID_S1_L001_R1_001.fastq.gz and COI-Sample-ID_S1_L001_R2_001.fastq.gz, the sample index would be 1-2. These correspond to the words Sample and ID. Currently, file naming conventions are very strict and all input files must end with L001_R1_001.fastq.gz or L001_R2_001.fastq.gz.

*Example Usage: --indices 1-2

--input_dir INPUT_DIR

*Required. This option specifies the path to the directory which contains the fastq.gz files.

*Example Usage: /home/data_set/fastq

--results_dir RESULTS_DIR

*Required. This option specifies the path to the directory which will store the results of the taxonomic assignment.

*Example Usage: /home/data_set/results

--PySCVUC_Path PYSCVUC_PATH

*Optional. Default is $HOME. The path to the directory containing the PySCVUC diretory.
                        
*Example Usage: --PySCVUC /home/scripts/
