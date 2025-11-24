<img width="432" height="12" alt="image" src="https://github.com/user-attachments/assets/4e80f464-2e53-4d20-ae28-4c0b2e3c8843" /># cfeccDNACodes
# Introduction
Cell free eccDNA codes (cfeccDNACodes) contains scripts for linearized cell-free eccDNA (cf-eccDNA) fragments analysis, cf-eccDNA feature extraction and cf-eccDNA feature-based cancer detection models.
The functions provided in cf-eccDNACodes are illustrated as follows.
<br/>
<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="https://github.com/user-attachments/assets/a340c24a-5bd7-4a71-ba5e-acc5f64c6feb">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">cfDNAFE Function</div>
</center>
<br/>


## Getting Started
## Hardware requirements
`cf-eccDNACodes` requires enough RAM to support the in-memory operations.

### OS Requirements
This package is supported for *Linux* and *macOS*. The package has been tested on the following systems:
+ Linux: Ubuntu 16.04.7
+ macOS: ventura 13.2.1

### Prerequisites
    Python --- 3.7.13
    numpy --- 1.21.5
    pandas --- 1.3.5
    scikit-learn --- 1.0.2
    Perl --- 5.32
    R --- 4.3
 
### Installation
1. **Necessary Step:** Download from Github:
```
git clone https://github.com/xqpeng/cfeccDNACodes.git
cd cfeccDNACodes
```
2. **Recommended Step:** (Conda users, Conda version: 4.12.0) Create your environment and activate it:
```
conda env create -f environment.yml
conda activate cfeccDNACodes
```
## Usage
1. **sh map.sh *<input_dir>* *<output_dir>* *<genome_index_dir>***
   
   Function: Align the paired fastq files from the file directory *input_dir*, and output the sorted bam files to the file directory *output_dir*.
   Make sure you have built indexed reference genome and provided it in *genome_index_dir*.
   
2. **sh runSplit_read.sh *<input_dir>* *<output_dir>***
   
   Function: Extract split-aligned fragments from sorted BAM files.
   
3. **runDiscordant_read.sh *<input_dir>* *<output_dir>***

   Function: Extract discordantly-aligned fragments from sorted BAM files.
   
4. **perl SplitAnalysis.pl *<list_file>***

   Function:
   
   (1) Extract two ends of each fragment, following the sequencing direction 5'-3'.
   
   (2) Perform statistic analysis on the split-aligned fragments from each sample.
   
	 (3) Perform length analysis on the split-aligned fragments from a group of samples.

   (4) Perform lengths, endpoint distances, and identical 6-mer end motifs analysis on the split-aligned fragments from each group.

   Parameter *list_file* is a file contain the split-aligned SAM filenames of a group of samples.
   
5. **perl DiscordantAnalysis.pl  *<list_file>***
   
   Function:
   
   (1) Perform statistic analysis on the discordant-aligned fragments from each sample.
   
   (2) Perform  identical 6-mer end motifs analysis on the discordant-aligned fragments from each group

   Parameter *list_file* is a file contain the discordantly-aligned SAM filenames of a group of samples.
	
   
7. **Rscript EccDNAFE.R *<feature1,feature2,...>* *<dir1>* *<dir2>* ...**

   Function: Extract eight types of feature profile, including BPM, EDM, JNM, SBM, OJM, OLR, CNV_onco, and CNV_im for each sample in the input directories, such as *dir1* *dir2* ...
   
   Example: Rscript EccDNAFE.R JNM,SBM,OLR, /path/to/dir1 /path/to/dir2 /path/to/dir3

   Available features: BPM, EDM, JNM, SBM, OJM, OLR, CNV_onco, CNV_im
  
10. **runClassifiers *<Control_samples_dir>* *<DiseaseCase_samples_dir>* **

    Function: Build Random Forest Classifier based on each type of cf-eccDNA feature, and build assemble model (CFECC) based on the classifiers of eight types of cf-eccDNA features.
    
    
