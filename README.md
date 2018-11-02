# PySCVUC
A Python3 implementation of the SCVUC amplicon sequencing pipeline which is maintained by Dr. Teresita M. Porter.

Based upon the repository by Dr. Teresita M. Porter (https://github.com/EcoBiomics-Zoobiome/SCVUC_COI_metabarcode_pipeline), PySCVUC uses the same programs, algorithms, and reference dataset outlined in the aformentioned work. The goal of this project is to unify the aformentioned pipeline under a common but modular code base. For further information, please refer to the link. At the moment, this pipeline only works for paired end reads. An option will be available in the future to make use of single end reads.

Denoising and classification pipeline for paired amplicon sequencing data. In order to use this pipeline SeqPrep, CutAdapt, USearch 10.0.240, VSearch 2.8.2, and the RDP Tools must be installed and the location of each of the programsmust be in your path. In addition, this pipeline is compatable with Python 3.4+ and the following Python dependencies are needed: NumPy, Pandas, BioPython, and Scipy. 

Install:
1) Clone this repository to your home directory
2) Download and install the above programs
3) Download and install the above Python libraries
4) Edit your .bashrc file found in your home directory so the following paths are present:

export PATH=<PATH TO USEARCH10>/:$PATH
export PATH=<PATH TO VSEARCH 2.8.2>/bin/:$PATH
export PATH=<PATH TO RDP TOOLS>/:$PATH
export PATH=<PATH TO PySCVUC>/:$PATH
export CLASSPATH=/home/joe/bin/RDPTools/:$CLASSPATH
  
5) Install the specified Python libraries: NumPy, SciPy, Pandas, and BioPython

Command Line Options:
--primers PRIMERS       The primer sequences which will be trimmed by
                        CutAdapt. Each sequences should be separated by a
                        dash.
--amplicon_names AMPLICON_NAMES
                        The name of each primer sequence. Each name should be
                        separated by a dash.
--n N                 The stringency setting for CutAdapt. (Default: 3)
--classifier CLASSIFIER     
                        The name of the classifier to use (Default: V3).
--threads THREADS     The number of threads to spawn. (Default: 10)
--indices INDICES     The index values in the filename which correspond to
                      the sample name. Each index should be separated by a
                      dash.
--input_dir INPUT_DIR
                      The directory which contains the fastq.gz files.
--results_dir RESULTS_DIR
                      The output directory.
--PySCVUC_Path PYSCVUC_PATH
                        The path to the directory containing the PySCVUC
                        Directory (Default: Home directory).
                        
Example Usage: 

If the fastq.gz files are: COI-Sample-ID_S1_L001_R1_001.fastq.gz and COI-Sample-ID_S1_L001_R2_001.fastq.gz, the sample index (--indices option) would be 1-2. These correspond to the words Sample and ID.

PySCVUC.py --primers ATCGATCG-CGATCGAT --amplicon_names ForwardPrimer-ReversePrimer --classifier V3 --indices 1-2 --input_dir /home/seqdata --results_dir /home/seqresults
