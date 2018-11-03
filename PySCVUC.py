#!/usr/bin/python3

#########################################################################
"Project Information"
#########################################################################
"""
This script was designed to automate a bioinformatic pipeline
and consolidate a number of scripts under a single code base.
The general steps of the analysis were adapted by Josip Rudar 
from various ZSH and Perl scripts written by Terri Porter.

Citation: https://github.com/EcoBiomics-Zoobiome/SCVUC_COI_metabarcode_pipeline 

"""

#########################################################################
"Global Imports"
#########################################################################
import subprocess

import pandas as pd

import numpy as np

from scipy.stats import mode

from sys import argv, exit

import argparse

from os.path import isfile, join, exists
from os import getcwd, makedirs, listdir, chdir

import gzip

from re import compile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from multiprocessing import Pool

#########################################################################
"Command Line Arguments"
#########################################################################
parser = argparse.ArgumentParser(description = "Denoising and classification pipeline for paired amplicon sequencing data. \
                                                In order to use this pipeline SeqPrep, CutAdapt 1.14, USearch 10.0.240, \
                                                VSearch 2.8.2, and the RDP Tools must be installed and the location \
                                                of each of the programs must be in your path. \
                                                In addition, this pipeline is compatable with Python 3.4+ and the following \
                                                Python dependencies are needed: NumPy, Pandas, BioPython, and Scipy.  \
                                                Example Usage: PySCVUC.py --primers ATCGATCG-CGATCGAT --amplicon_names A1-A2 --classifier V3 --indices 2-3 --input_dir /home/seqdata --results_dir /home/seqresults")

parser.add_argument("--primers", help = "The primer sequences which will be trimmed by CutAdapt. Each sequences should be separated by a dash.", required = True)
parser.add_argument("--amplicon_names", help = "The name of each amplicon. Each name should be separated by a dash.", required = True)
parser.add_argument("--n", help = "The stringency setting for CutAdapt. (Default: 3)", default = "3")
parser.add_argument("--classifier", help = "The name of the training set to use. (Default: V3)", default = "V3")
parser.add_argument("--threads", help = "The number of threads to spawn. (Default: 10)", default = "10")
parser.add_argument("--pe_names", help = "The file ending for paired-end reads. Should be separated by a dash. Example: L001_R1_001 and L001_R2_001. (Default: L001_R1_001-L001_R2_001)", default = "L001_R1_001-L001_R2_001")
parser.add_argument("--indices", help = "The index values in the filename which correspond to the sample name. Each index should be separated by a dash.", required = True)
parser.add_argument("--input_dir", help = "The directory which contains the fastq.gz files.", required = True)
parser.add_argument("--results_dir", help = "The output directory.", required = True)
parser.add_argument("--PySCVUC_Path", help = "The path to the directory containing the PySCVUC Directory (Default: Home directory).", default = "$HOME")

#########################################################################
"Utility Functions"
#########################################################################
def change_inosine(primer):
    
    new_primer = "".join([x if x != "I" else "D" for x in primer])

    return new_primer

def reverse_complement(primer):

    new_primer = Seq("".join([x if x != "I" else "D" for x in primer]),
                     generic_dna)

    new_primer = str(new_primer.reverse_complement())
    
    return new_primer

def fix_headers(fa_file):
    "Rename each file and correct the headers of each FASTA file"

    seq_origin = "_".join(fa_file.split(".")[0].split("-"))

    header_components = [">", seq_origin, ";", ""]

    file_data = []
    fasta_file = []

    with open(fa_file, "r") as file:
        file_data = file.readlines()

        headers = file_data[0::2]
        seqs = file_data[1::2]

    for sequence_set in zip(headers, seqs):
        header_components[3] = sequence_set[0][1:]
        header = "".join(header_components)
        sequence = sequence_set[1]
        fasta_file.append(header)
        fasta_file.append(sequence)

    return fasta_file

def get_fnames(ftype):
    "Get a list of files of type 'ftype' in the current working directory."
    cwdpath = getcwd()

    fnames = [file 
              for file in listdir(cwdpath) 
              if isfile(join(cwdpath, file)) and ftype in file]

    return fnames

def subprocess_command(command):
    """
    Issue multiple commands to the command line.
    From: https://stackoverflow.com/questions/17742789/running-multiple-bash-commands-with-subprocess
    """
    process = subprocess.Popen(command,
                               stdout=subprocess.PIPE, 
                               shell=True)

    proc_stdout = process.communicate()[0].strip()

def fix_zotu(fa_file):

    counter = 1

    with open("cat.mod.denoised", "w") as file:
        with open("cat.original.denoised", "r") as infile:
            data = infile.readlines()

            for line in data:
                new_line = None

                if ">" in line:
                    new_line = ">Otu%s" %str(counter)
                    counter += 1

                else:
                    new_line = line.strip("\n")

                print (new_line, file = file)
    
#########################################################################
"Pipeline statistics object and related functions"
#########################################################################
class NGSPipelineStats:
    """
    This class calculates the statistics for each sample or the 
    concatenated and/or denoised data.
    """
    def prep_detailed_stats(self):
        "Writes out the detailed per-file statistics."

        writer = pd.ExcelWriter("detailed_statistics.xlsx", 
                                engine = "xlsxwriter")

        for detailed_stat in self.tables.keys():
            if detailed_stat != "taxonomic_assignments_raw":
                df = pd.DataFrame.from_records(self.tables[detailed_stat][1:], 
                                               columns = self.tables[detailed_stat][0])

                df.to_excel(writer, 
                            sheet_name = detailed_stat)
    
        writer.save()

    def prep_otu_summary(self):
        """
        This function will create a summary of taxonomic assignments.
        """
        self.tables["taxonomic_assignments_raw"][0].append("Strand")

        final_table = []
        for entry in self.tables["taxonomic_assignments_raw"][1:]:
            final_entry = entry[:]
            final_entry.append(" ")
            final_table.append(final_entry)

        tax_df = pd.DataFrame.from_records(final_table, 
                                           columns = self.tables["taxonomic_assignments_raw"][0])

        if self.trained_clf == "V3":
            cutoff_table = [
                   ["Rank", "500+bp", "400bp", "200bp", "100bp", "50bp"],
                   ["Superkingdom", "0", "0", "0", "0", "0"],
                   ["Kingdom", "0", "0", "0", "0", "0"],
                   ["Phylum", "0", "0", "0", "0", "0"],
                   ["Class", "0", "0", "0", "0", "30"],
                   ["Order", "0", "0", "0", "20", "70"],
                   ["Family", "0", "10", "20", "20", "70"],
                   ["Genus", "40", "30", "30", "40", "90"],
                   ["Species", "NA", "NA", "NA", "NA", "NA"],
                   [" "],
                   ["NA = No cutoff available will result in 99% correct assignments"],
                   ["Minimum Bootstrap Support Cutoff (Species Rank Classifier) (%)"],
                   ["Prescribed Cutoffs for the RDP Classifier 2.12 with the COIv3.2 Training Set and BR5 Amplicons."]
                   ]

            tax_df[["ESVsize", 
                    "skBP", 
                    "pBP", 
                    "cBP", 
                    "oBP", 
                    "fBP", 
                    "sBP"]] = tax_df[["ESVsize", 
                                      "skBP", 
                                      "pBP", 
                                      "cBP", 
                                      "oBP", 
                                      "fBP", 
                                      "sBP"]].apply(pd.to_numeric)

        cutoff_table_df = pd.DataFrame.from_records(cutoff_table[1:], 
                                                    columns = cutoff_table[0])  

        writer = pd.ExcelWriter("taxonomic_assignments_raw.xlsx", 
                                engine = "xlsxwriter")

        tax_df["COI_GlobalESV"] = tax_df[["Amplicon", "GlobalESV"]].apply(lambda entry: "_".join(entry), axis = 1)

        tax_df = tax_df[self.final_cols]

        tax_df.to_excel(writer, 
                        sheet_name = "RDP 2.12 Taxonomic Results")

        cutoff_table_df.to_excel(writer, 
                                 sheet_name = "RDP Classifier 2.12 Cutoffs")

        writer.save()

    def add_abundance(self, amplicon, index_set):
        "Add abundance data to each OTU from each amplicon."

        #Variables to store final data, denoised data, abundance 
        #data, etc
        final_data = []
        table_samples = [] 
        denoised_data = dict()
        table_data = dict()

        #Open the denoised file, use it to populate a hash
        with open("results.out", "r") as file:
            for line in file:
                temp = line.strip("\n")
                temp = temp.split("\t")
        
                denoised_data[temp[0]] = [temp[0], 
                                          "sample", 
                                          amplicon, 
                                          "abundance"]

                denoised_data[temp[0]].extend(temp[2:])

        #Open the abundance data, read it, and extract and organize 
        #the information
        with open("cat.fasta.table", "r") as file:
            full_data = file.readlines()
            otu_abd_data = full_data[1:]

            #Get the sample line
            samples = full_data[0][1:].strip("\n").split("\t")[1:] 
            ind_samples = []

            #Rename each sample (for building raw assignment table)
            for full_sample_name in samples:
                if len(index_set) > 0:
                    temp = full_sample_name.split("_")
        
                    sample_components = []
                    for i, component in enumerate(temp):
                        if i in index_set: sample_components.append(component)
                    table_samples.append("_".join(sample_components))

                else:
                    table_samples.append(full_sample_name)

            #Report raw reads
            final_data = []
            for otu_line in otu_abd_data: 
                    line = otu_line.strip("\n").split("\t")
                    if line[0] in denoised_data: 
                        for i, entry in enumerate(line[1:]):
                            if int(entry) > 2:
                                file_entry = denoised_data[line[0]][:]
                                file_entry[1] = table_samples[i]
                                file_entry[3] = entry
                                final_data.append(file_entry)

            self.tables["taxonomic_assignments_raw"].extend(final_data)

    def uni_stats(self, fa_file, denoised = False):
        """
        This function calculates the statistics on unique and denoised file.
        This function also creates a file in which multi-line sequences are stored as one line.
        """
        
        #Set up variables for holding data
        fa_data = []
        counts = 0

        #Get the key name
        sample_key = fa_file

        #Create a modified cat.uniques fasta file with each sequence on a single line
        if denoised == False:
            for record in SeqIO.parse(fa_file, "fasta"):
                key = record.description
                seq = str(record.seq)

                #Get the length of the sequence
                fa_data.append(len(seq.strip("\n")))

                counts += int(key.strip("\n").split("=")[1].strip(";"))

        else:
            for record in SeqIO.parse(fa_file, "fasta"):
                key = record.description
                seq = str(record.seq)

                #Get the length of the sequence
                fa_data.append(len(seq.strip("\n")))

                counts += 1

        #Determine the number of OTUs and stats
        num_seqs = len(fa_data)
        stats = self.stat_math(np.asarray(fa_data, dtype = np.int))
        
        return [sample_key, 
                num_seqs, 
                stats[0], 
                stats[1], 
                stats[2], 
                stats[3], 
                stats[4], 
                counts]

    def trimmed_stats(self, index, amp_index, stage):
        "This function calculates statistics on trimmed and denoised data."

        #Get the names of all FASTA files in that directory
        if stage == "trim1": dir_files = get_fnames("Ftrimmed.fastq.gz")
        else: dir_files = get_fnames("Rtrimmed.fasta")

        trimmed_stats_table = [["Sample", 
                                "NumSeqs", 
                                "MinLen", 
                                "MaxLen", 
                                "MeanLen", 
                                "MedianLen", 
                                "ModeLen"]]

        num_reads = 0

        if stage == 'trim1':
            trimmed_results = [self.raw_stats(file) for file in dir_files]
        else:
            trimmed_results = [self.raw_stats(file, True) for file in dir_files]

        for result in trimmed_results:
            modified_result = [result[0]]
            modified_result.extend(result[2:-1])
            num_reads += int(result[2])
            trimmed_stats_table.append(modified_result)

        table_key = str(index) + "_" + stage + "_" + amp_index
        self.tables[table_key] = trimmed_stats_table

        summary_key = "Number of Trimmed Reads (" + amp_index + ")"
        self.summary_dict[summary_key] = num_reads
        
        summary_key = "% Reads Trimmed (" + amp_index + ")" 

        self.summary_dict[summary_key] = num_reads / self.summary_dict["Total Paired Reads"]

    def raw_stats(self, fa_file, trim2 = False):
        "This function calculates the raw and paired statistics."

        #Get the name of the keys
        key_names = fa_file
        
        #Replace all dashes with underscores
        if "-" in key_names: 
            expr_sub = compile("-")
            key_names = expr_sub.sub("_", key_names)
    
        #Get the sample name and the read
        key_names = key_names.split("_")
    
        #Create the sample key given the indices        
        sample_key = [component 
                      for i, component in enumerate(key_names) 
                      if i in self.indices]

        sample_key = "_".join(sample_key)

        #Get the end
        if "R1" in set(key_names): end = "R1"
        elif "R2" in set(key_names): end = "R2"
        else: end = "Paired"

        #Open and read the file
        data = None
        try:
            with gzip.open(fa_file, "rt") as file: 
                data = file.readlines()
        except:
            with open(fa_file, "rt") as file:
                data = file.readlines()

        #Get every sequence
        if trim2 == False:
            data = data[1::4]
        else:
            data = data[1::2]

        #Count the number of entries in the file
        num_seqs = len(data) 
                                   
        #Create an array the size of num_seqs
        seq_lens = np.zeros((num_seqs,), dtype = np.int) 
      
        #Determine the length of the sequence
        for i, line in enumerate(data):
            seq_lens[i] = len(line.strip("\n")) 
        
        #Calculate the stats
        stats = self.stat_math(seq_lens)

        return [sample_key, 
                end, 
                num_seqs, 
                stats[0], 
                stats[1], 
                stats[2], 
                stats[3], 
                stats[4], 
                fa_file]

    @staticmethod
    def stat_math(np_array):
        """
        This function runs all the stats - min, max, mean, median, 
        and mode on a numpy array
        """
        try:
            min_len = np_array.min() #Find the smallest value
            max_len = np_array.max() #Find the largest value
            mean_len = np_array.mean() #Find the mean
            mode_len = int(mode(np_array)[0]) #Find the mode value
            median_len = np.median(np_array) #Find the median value

        except:
            min_len = 0
            max_len = 0
            mean_len = 0
            mode_len = 0
            median_len = 0

        return min_len, max_len, mean_len, median_len, mode_len 

    def __init__(self, samp_ind, clf_type):
        #Set the indices which will determine the sample names
        self.indices = samp_ind
        self.trained_clf = clf_type

        #Setup hash to hold results for calculating summary stats
        self.summary_dict = {}

        if clf_type == "V3":
            self.column_headers = ["GlobalESV", 
                                   "SampleName", 
                                   "Amplicon", 
                                   "ESVsize",
                                   "Root",
                                   "RootRank",
                                   "rBP",
                                   "SuperKingdom", 
                                   "SuperKingdomRank", 
                                   "skBP", 
                                   "Kingdom", 
                                   "KingdomRank", 
                                   "kBP", 
                                   "Phylum", 
                                   "PhylumRank", 
                                   "pBP", 
                                   "Class", 
                                   "ClassRank", 
                                   "cBP",
                                   "Order", 
                                   "OrderRank", 
                                   "oBP", 
                                   "Family", 
                                   "FamilyRank", 
                                   "fBP", 
                                   "Genus", 
                                   "GenusRank", 
                                   "gBP", 
                                   "Species", 
                                   "SpeciesRank", 
                                   "sBP"]

        self.final_cols = ["COI_GlobalESV", "SampleName", "ESVsize", "Strand"]
        self.final_cols.extend(self.column_headers[4:])

        #Hash to hold the tables produced by the analysis
        self.tables = {
                       "raw_stats_table": [["Sample", 
                                            "Read", 
                                            "NumSeqs", 
                                            "MinLen", 
                                            "MaxLen", 
                                            "MeanLen", 
                                            "MedianLen", 
                                            "ModeLen", 
                                            "Sample", 
                                            "Read", 
                                            "NumSeqs", 
                                            "MinLen", 
                                            "MaxLen", 
                                            "MeanLen", 
                                            "MedianLen", 
                                            "ModeLen"]],

                       "paired_stats_table": [["Sample", 
                                               "NumSeqs", 
                                               "MinLen", 
                                               "MaxLen", 
                                               "MeanLen", 
                                               "MedianLen", 
                                               "ModeLen"]],

                       "unique_stats_table": [["Sample", 
                                               "NumOTUs", 
                                               "MinLen", 
                                               "MaxLen",
                                               "MeanLen", 
                                               "MedianLen", 
                                               "ModeLen", 
                                               "Counts"]],

                       "denoised_stats_table": [["Sample", 
                                                 "NumOTUs", 
                                                 "MinLen", 
                                                 "MaxLen",
                                                 "MeanLen", 
                                                 "MedianLen", 
                                                 "ModeLen", 
                                                 "Counts"]],

                       "taxonomic_assignments_raw": [self.column_headers]
                      }

#########################################################################
"This section contains all the main steps of the pipeline."
#########################################################################
print (parser.parse_args())

args = parser.parse_args()

install_path = vars(args)["PySCVUC_Path"]

#Dictionary containing command to use each training set
classifier_dict = {"V3": "_JAVA_OPTIONS='Xmx8g';classifier.jar classify -t %s/PySCVUC/training_files/V3/rRNAClassifier.properties -o results.out cat.denoised" %install_path,}

primer_pairs = vars(args)["primers"].split("-")
amplicons = vars(args)["amplicon_names"].split("-")
classifier = vars(args)["classifier"]
n_val = vars(args)["n"]

n_threads = vars(args)["threads"]

pe_names = vars(args)["pe_names"].split("-")

sample_ind = {int(index) 
              for index in set(vars(args)["indices"].split("-"))}
    
input_dir = vars(args)["input_dir"]
results_dir = vars(args)["results_dir"]

chdir(results_dir)

#Create an instance of the statistics class
stats_class = NGSPipelineStats(sample_ind, 
                               classifier)

#Make the directories and copy files
if not exists("final"): makedirs("final")
chdir("final")
final_path = getcwd()
chdir("..")

if not exists("raw"): makedirs("raw")

command = "cp %s/*.* raw/" %input_dir
subprocess_command(command)

try:
    command = "unzip raw/*.zip -d raw/"
    subprocess_command(command)

    command = "rm -f raw/*.zip"
    subprocess_command(command)

except:
    pass

chdir("raw")

#########################################################################
"Calculate Raw Statistics"
#########################################################################
print (' ')

dir_files = get_fnames("fastq.gz")
workers = Pool(int(n_threads))
raw_results = workers.map(stats_class.raw_stats, dir_files) 
workers.close()
workers.join()

#Populate the hash table with the results for each sample
temp_table = {}
for result in raw_results:
    if result[0] not in temp_table:
        temp_table[result[0]] = list()

    temp_table[result[0]].append(result[0:-1])

r1_reads = 0
r2_reads = 0
for entry, result in temp_table.items():
    temp_list = None
    print (entry, result)
    if result[0][1] == "R1":
        temp_list = result[0][:]
        temp_list.extend(result[1][:])
        r1_reads += int(result[0][2])
        r2_reads += int(result[1][2])
    else:
        temp_list = result[1][:]
        temp_list.extend(result[0][:])
        r1_reads += int(result[1][2])
        r2_reads += int(result[0][2])

    stats_class.tables["raw_stats_table"].append(temp_list)

print (' ')
        
number_of_sites = len(stats_class.tables["raw_stats_table"]) - 1
stats_class.summary_dict["Number of Samples"] = number_of_sites
stats_class.summary_dict["Total Number of Reads (R1)"] = r1_reads
stats_class.summary_dict["Total Number of Reads (R2)"] = r2_reads
stats_class.summary_dict["Read Coverage"] = int(r1_reads / number_of_sites)

temp_table.clear()

#########################################################################
"Pair Reads and Calculate Paired Statistics"
#########################################################################
#Pair the reads and calculate the paired statistics
dir_files =  get_fnames("fastq.gz")

raw_file_hash = {}
for entry in dir_files:
    file_name = entry.split(".")[0]

    new_fname = file_name.replace(pe_names[0], '').replace(pe_names[1], '')

    if new_fname not in raw_file_hash:
        raw_file_hash[new_fname] = list()

    raw_file_hash[new_fname].append(entry)

print (" ")
print ("Getting ready to pair reads...")
for key, value in raw_file_hash.items():
    value.sort()

    print (key, value)

print (" ")

seq_commands = []
for key, value in raw_file_hash.items():
    file_1 = value[0]
    file_2 = value[1]

    seq_command = 'seqprep -f %s -r %s -1 %s.out -2 %s.out -q 20 -s %s.paired.fastq.gz -o 25' %(file_1, file_2, file_1, file_2, key)
    seq_commands.append(seq_command)

workers = Pool(int(n_threads))
workers.map(subprocess_command, seq_commands)
workers.close()
workers.join()

if not exists("paired"): makedirs("paired")
command = "rm -f *.out;mv *.paired.fastq.gz paired;rm -f *.fastq.gz"
subprocess_command(command)
chdir("paired")

dir_files =  get_fnames("fastq.gz")
workers = Pool(int(n_threads))
paired_results = workers.map(stats_class.raw_stats, dir_files)
workers.close()
workers.join()

reads_paired = 0
for entry in paired_results:
    reads_paired += int(entry[2])
    final_entry = [entry[0]]
    final_entry.extend(entry[2:-1])
    stats_class.tables["paired_stats_table"].append(final_entry)

stats_class.summary_dict["Total Paired Reads"] = reads_paired

paired_percentage = (reads_paired / stats_class.summary_dict["Total Number of Reads (R1)"]) * 100
paired_percentage = str(paired_percentage)

len_perc = len(paired_percentage)

if len_perc > 3:
    stats_class.summary_dict["Percentage Paired"] = paired_percentage[0:4] + '%'
else:
    stats_class.summary_dict["Percentage Paired"] = paired_percentage + '%'

#########################################################################
"""
Loop through each amplicon and do the following:
    1) Run Cutadapt
    2) Dereplicate and Denoise (Vsearch/Usearch)
    3) Create OTU Table (Vsearch)
    4) Classify (RDP Classifier) and Add Abundance Information
"""
#########################################################################
num_primers = len(primer_pairs)
primer_index = 0
while True:
    #########################################################################
    "Trim with Cutadapt"
    #########################################################################
    dir_name = amplicons[primer_index] + "_Trimmed"
    if not exists(dir_name): makedirs(dir_name)

    print (' ')
    print ("Trimming %s from sequences with Cutadapt... " %primer_pairs[primer_index])

    new_primer = change_inosine(primer_pairs[primer_index])

    dir_files =  get_fnames("paired.fastq.gz")
    cut_commands = []
    for file_name in dir_files:
        cut_commands.append("cutadapt -g %s -m 150 -q 20,20 --max-n=%s --discard-untrimmed %s -o %s.Ftrimmed.fastq.gz" %(new_primer, str(n_val), file_name, file_name))

    workers = Pool(int(n_threads))
    workers.map(subprocess_command, cut_commands)
    workers.close()
    workers.join()

    stats_class.trimmed_stats(primer_index, 
                              amplicons[primer_index], 
                              "trim1")
    
    command = "mv *.Ftrimmed.fastq.gz " + dir_name
    subprocess_command(command)

    chdir(dir_name)
    
    #Trim reverse primers if the sequences are paired
    primer_index += 1
    print (' ')
    print ("Trimming %s from sequences with Cutadapt... " %primer_pairs[primer_index])
    
    new_primer = reverse_complement(primer_pairs[primer_index])

    dir_name = amplicons[primer_index] + "_Trimmed"
    if not exists(dir_name): makedirs(dir_name)
     
    dir_files =  get_fnames("Ftrimmed.fastq.gz")
    cut_commands = []
    for file_name in dir_files:
        cut_commands.append("cutadapt -a %s -m 150 -q 20,20 --max-n=%s  --discard-untrimmed %s -o %s.Rtrimmed.fasta.gz" %(new_primer, str(n_val), file_name, file_name))

    workers = Pool(int(n_threads))
    workers.map(subprocess_command, cut_commands)
    workers.close()
    workers.join()

    command = "mv *.Rtrimmed.fasta.gz " + dir_name
    subprocess_command(command)

    chdir(dir_name)

    stats_class.trimmed_stats(primer_index, 
                              amplicons[primer_index], 
                              "trim2")

    print("\nCurrent Statistical Summary After Trimming Reads: ")
    for key, value in stats_class.summary_dict.items():
        print (key, value)
    print("\n")

    current_dir = getcwd()
    chdir("..")
    chdir("..")
    chdir("..")

    linked_dir_name = "QCdFastas-" + amplicons[primer_index]
    command = "ln -s " + current_dir + " " + linked_dir_name
    subprocess_command(command)

    chdir(linked_dir_name)

    #########################################################################
    "Dereplicate and Denoise"
    #########################################################################
    command = "ls | grep .gz | parallel -j %s gunzip {}" %n_threads
    subprocess_command(command)

    dir_files =  get_fnames("fasta")
    fixed_files = [fix_headers(file) for file in dir_files]

    final_file = [line for file in fixed_files for line in file]

    with open ("cat.fasta", "w") as fa_file:
        [print (line, end = "", file = fa_file) for line in final_file]

    command = "vsearch --threads %s --derep_fulllength cat.fasta --output cat.uniques --sizein --sizeout " %n_threads
    subprocess_command(command)

    if not exists("global_uniques"): makedirs("global_uniques")
    command = "mv cat.uniques global_uniques"
    subprocess_command(command)
    chdir("global_uniques")

    unique_results = stats_class.uni_stats("cat.uniques")

    stats_class.tables["unique_stats_table"].append(unique_results)
    stats_class.summary_dict["Unique Sequences" + amplicons[primer_index - 1]] = unique_results[1]

    command = "usearch10 -unoise3 cat.uniques -zotus cat.original.denoised -minsize 3"
    subprocess_command(command)

    if not exists("global_denoised"): makedirs("global_denoised")
    command = "mv cat.original.denoised global_denoised"
    subprocess_command(command)
    chdir ("global_denoised")

    denoised_results = stats_class.uni_stats("cat.original.denoised", True)
    stats_class.tables["denoised_stats_table"].append(unique_results)
    stats_class.summary_dict["Denoised Sequences" + amplicons[primer_index - 1]] = denoised_results[1]

    print("\nCurrent Statistical Summary After Denoising: ")
    for key, value in stats_class.summary_dict.items():
        print (key, value)
    print("\n")

    fix_zotu("cat.original.denoised")

    chdir("..")
    chdir("..")

    #########################################################################
    "Create OTU Table and Classify"
    #########################################################################
    command = "vsearch --usearch_global cat.fasta --db global_uniques/global_denoised/cat.mod.denoised --id 1.0 --otutabout cat.fasta.table --threads %s" %n_threads
    subprocess_command(command)

    if not exists("global_OTU_table"): makedirs("global_OTU_table")
    command = "mv *.table global_OTU_table/"
    subprocess_command(command)

    chdir("global_uniques")
    chdir("global_denoised")

    command = "mv cat.original.denoised cat.original.denoised.bak;mv cat.mod.denoised cat.denoised"
    subprocess_command(command)

    command = classifier_dict[classifier]
    subprocess_command(command)

    if not exists("clf_results"): makedirs("clf_results")
    command = "mv results.out clf_results/"
    subprocess_command(command)
    chdir("clf_results")

    command = "cp ../../../global_OTU_table/cat.fasta.table ."
    subprocess_command(command)

    stats_class.add_abundance(amplicons[primer_index-1], 
                              sample_ind)

    chdir("..")
    chdir("..")
    chdir("..")
    chdir("..")
    chdir("..")
    chdir("..")

    primer_index += 1
    if primer_index >= num_primers: break

    chdir("paired")

chdir("..")
chdir("final")


#########################################################################
"""
1) Prepare Taxonomic Summary Table (XLSX)
2) Prepare full statistics (XLSX)
"""
#########################################################################
stats_class.prep_otu_summary()
stats_class.prep_detailed_stats()

if not exists("Assignments"): makedirs("Assignments")
command = "mv *_assignments_raw.xlsx Assignments"
subprocess_command(command)

if not exists("Statistics"): makedirs("Statistics")
command = "mv *_statistics.xlsx Statistics"
subprocess_command(command)
