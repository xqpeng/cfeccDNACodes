# Introduction
Cell free eccDNA codes (cfeccDNACodes) contains scripts for linearized cell-free eccDNA (cf-eccDNA) fragments analysis, cf-eccDNA feature extraction and cf-eccDNA feature-based cancer detection models.
The functions provided in cf-eccDNACodes are illustrated as follows.
![codes](https://github.com/user-attachments/assets/5770489b-7ad9-4814-b872-ed6bbe63fa00)



## Getting Started
## Hardware requirements
`cf-eccDNACodes` requires enough RAM to support the in-memory operations.

### OS Requirements
This package is supported for *Linux* and *macOS*. The package has been tested on the following systems:
+ Linux: Ubuntu 16.04.7
+ macOS: ventura 13.2.1

### Prerequisites
Please ensure the following dependencies are installed:

**Python Packages:**
- Python (>=3.7)
- numpy
- pandas
- scikit-learn

**R Packages:**
- R (>=4.0)
- Bioconductor: BSgenome, BSgenome.Hsapiens.UCSC.hg38, GenomeInfoDb, GenomicAlignments, GenomicRanges, IRanges, Rsamtools
- CRAN: dplyr, mgsub, readxl, stringr, tidyr

**Other:**
- Perl
	
 
### Installation
1. **Necessary Step:** Download from Github:
```
git clone https://github.com/xqpeng/cfeccDNACodes.git
cd cfeccDNACodes
```
2. **Recommended Step:** (Conda users, Conda version: 4.12.0) Create your environment and activate it:
```
conda env create -f environment.yml
conda activate cfeccDNA_analysis
Rscript setup.R
```

## Usage
1. **sh map.sh *<input_dir>* *<output_dir>* *<genome_index_dir>***
   
   Function: Align the paired fastq files from the file directory *input_dir*, and output the sorted bam files to the file directory *output_dir*.
   
   Make sure you have built indexed reference genome and provided it in *genome_index_dir*.
   
   Input parameters:
   
          (1) <input_dir>: the file directory containing the paired fastq files
   
          (2) <output_dir>: the file directory to store the sorted bam files
   
          (3) <genome_index_dir>: the directory containing the indexed reference genome
   
    Output: Sorted bam files

2. **sh runSplit_read.sh *<input_dir>* *<output_dir>***
   
   Function: Extract split-aligned fragments from sorted BAM files.
   
   Input parameters:
   
         (1) <input_dir>: the file directory containing the sorted bam files
   
         (2) <output_dir>: the file directory to store the bam files of split-aligned fragments
   
   Output: The bam files of split-aligned fragments
   
3. **sh runDiscordant_read.sh *<input_dir>* *<output_dir>***

   Function: Extract discordantly-aligned fragments from sorted BAM files.
   
   Input parameters:
   
         (1) <input_dir>: the file directory containing the sorted bam files
   
         (2) <output_dir>: the file directory to store the bam files of split-aligned fragments
   
   Output: The bam files of split-aligned fragments
   
4. **perl SplitAnalysis.pl *<list_file>***

   Function:
   
          (1) Extract two ends of each fragment, following the sequencing direction 5'-3'.
   
          (2) Perform statistic analysis on the split-aligned fragments from each sample.
   
          (3) Perform length analysis on the split-aligned fragments from a group of samples.

          (4) Perform lengths, endpoint distances, and identical 6-mer end motifs analysis on the split-aligned fragments from each group.

   Input parameter:

           (1)  <list_file> is a file contain the split-aligned SAM filenames of a group of samples.
   
              For example: list_HCC
   
		      which contains three lines and there is a filename on each line, such as following:
   
				         HCC1.split_reads.sam
   
				         HCC2.split_reads.sam
   
				         HCC3.split_reads.sam
   
                         ...
   Output:

           (1) statics_proportion_length_list_HCC; # corresponding to  function 2
   
           (Belowing are optional)
   
		   (2) File of Fragment length distribution
   
		   (3) File of eccDNA molecules size distribution
   
		   (4) File of endpoint distance distribution of Same-Strand 5' Ends fragments
   
		   (5) File of endpoint distance distribution of Opposite-Strand 5' Ends fragments
   
		   (6) File of Identical 6-mer end motifs profile of Same-Strand 5' Ends fragments
   
		   (7) File of Identical 6-mer end motifs profile of split-aligned fragments fragments
   
	       (8) Files of the 5ends of fragments corresponding to the files in the input file. 
   
5. **perl DiscordantAnalysis.pl  *<list_file>***
   
   Function:
   
         (1) Perform statistic analysis on the discordant-aligned fragments from each sample.
   
         (2) Perform  identical 6-mer end motifs analysis on the discordant-aligned fragments from each group

   Input parameter:
   
         (1) list_file is a file contain the discordantly-aligned SAM filenames of a group of samples.
   
             For example: list_HCC_discordant_pairs
   
		      which contains three lines and there is a filename on each line, such as following:
   
				         HCC1_discordant_pairs.sam
   
				         HCC2_discordant_pairs.sam
   
				         HCC3_discordant_pairs.sam
                         ...
	 Output:

             (1) statics_discordant_fragments_list_HCC_discordant_pairs;

             (2) File of Identical 6-mer end motifs profile of discordantly-aligned fragments
   
7. **Rscript EccDNAFE.R *<feature1,feature2,...>*  *\<dir1\>*  *\<dir2\>* ... *\<dirn\>***

   Function: Extract eight types of feature profile, including BPM, EDM, JNM, SBM, OJM, OLR, CNV_onco, and CNV_im for each sample in the input directories, such as *dir1* *dir2* ...
   
   Make sure the dependencies have properly installed. 
   
   Input parameter:
   
   	     (1)<feature1,feature2,...>: available features include BPM, EDM, JNM, SBM, OJM, OLR, CNV_onco, and CNV_im.
		 
		 (2)<dir1>  <dir2> ... <dirn>: the file directories containing the bam files of split_fragments/discordant_fragments. For bam files with discordant_fragments, only BPM, EDM, CNV_onco, and CNV_im are extracted. For bam files with split_fragments, BPM, EDM, JNM, SBM, OJM, OLR, CNV_onco, and CNV_im are extracted.
		 
   Output: A feature directory corresponding to each feature type is built for each input directory, and includes the feature profiles corresponding to the files under the input directory.
   
   Example:

          Rscript EccDNAFE.R BPM,EDM,JNM,SBM,OJM,OLR,CNV_onco,CNV_im /path/to/dir1 /path/to/dir2 
   
8. **runClassifiers *<Control_samples_dir>* *<DiseaseCase_samples_dir>***

    Function: Build Random Forest Classifier based on each type of cf-eccDNA feature, and build the assemble model (CFECC) based on the classifiers of eight types of cf-eccDNA features.
    
    Input parameter:
    
          (1)<Control_samples_dir>: A directory containing eight sub-directories corresponding to eight types of cf-eccDNA features of control samples.
   			
           Control_samples_dir/
           ├──BPM
           │   ├──sample1_BPM.txt
           │   ├──sample2_BPM.txt
           ├──EDM
           │   ├──sample1_EDM.txt
           │   ├──sample2_EDM.txt
           ├──JNM
           │   ├──sample1_JNM.txt
           │   ├──sample2_JNM.txt
           ├──SBM
           │   ├──sample1_SBM.txt
           │   ├──sample2_SBM.txt
           ├──OJM
           │   ├──sample1_SBM.txt
           │   ├──sample2_SBM.txt
           ├──OLR
           │   ├──sample1_SBM.txt
           │   ├──sample2_SBM.txt
           ├──CNV_onco
           │   ├──sample1_CNV_onco.txt
           │   ├──sample2_CNV_onco.txt
           

    
          (2)<DiseaseCase_samples_dir>: A directory containing eight sub-directories corresponding to eight types of cf-eccDNA features of case samples.
           
           Control_samples_dir/
           ├──BPM
           │   ├──sample1_BPM.txt
           │   ├──sample2_BPM.txt
           ├──EDM
           │   ├──sample1_EDM.txt
           │   ├──sample2_EDM.txt
           ├──JNM
           │   ├──sample1_JNM.txt
           │   ├──sample2_JNM.txt
           ├──SBM
           │   ├──sample1_SBM.txt
           │   ├──sample2_SBM.txt
           ├──OJM
           │   ├──sample1_SBM.txt
           │   ├──sample2_SBM.txt
           ├──OLR
           │   ├──sample1_SBM.txt
           │   ├──sample2_SBM.txt
           ├──CNV_onco
           │   ├──sample1_CNV_onco.txt
           │   ├──sample2_CNV_onco.txt
           
    Output:
    
          (1) GridSearch_All_Features_Results.csv: this file store the parameters of classifiers on the training dataset.
    
          (2) ClassifyPerformance.csv: this file store the performance of eight classifiers and the assemble model (CFECC) on the validation dataset.
    
    Example:
    
           python runClassifiers.py /path/to/HCC /path/to/Control
    
