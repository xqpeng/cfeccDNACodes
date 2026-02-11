Introduction
This is a pipeline for identification and analysis of H-DNA (triplex DNA) sites
from Same-Strand 5′ Ends fragments with endpoint distance less than 5 bp. 
The pipeline covers the entire workflow from format conversion, H-DNA structural searching, and multi-dimensional statistical visualization.

Getting Started
Hardware Requirements
CPU: Multi-core processor recommended (The pipeline is optimized for 18-thread parallel execution).
RAM: 16GB or higher recommended to support large-scale genomic in-memory vectorized operations.
Prerequisites
Ensure the following dependencies are installed:
Python (>=3.7):
numpy, pandas (for CSV data cleaning and filtering).
R (>=4.0):
Bioconductor: TriplexR (core search engine), Biostrings (sequence manipulation).
CRAN: tidyverse, ggplot2, ggpubr, data.table (high-speed processing), stringi (fast pattern matching), parallel (cluster management).
 

Usage & Pipeline Steps

Phase 1: Preprocessing
Step 1: Format Conversion 
Script: sequences_to_fasta.py
Input: SAM files in the Target_sequences directory.
Function: Converts SAM files in the Target_sequences directory to standardized FASTA format.

Step 2: Proximity Filtering 
Script: sequence_ls5bp_matching.R
Input: FASTA files (from Step 1) and LS_5bp reference list.
Function: Filters sequences with gaps less than 5bp using the LS_5bp reference list to generate refined FASTA files.

Phase 2: Triplex Identification
Step 3: Structural Search 
Script: triplex_parallel_search_pipeline.R
Input: LS5bp_FASTA files (from Step 2).
Function: Executes TriplexR to identify triplex structures. This script features checkpointing (resume capability) and 6-core parallelization.
Output:
all_triplex_sequences_detailed.csv: Detailed structural components including Stem1, Loop, and Stem2.
analysis_summary.csv: Global sequence metrics and H-DNA detection status (AnalysisStatus).

Phase 3: Feature Extraction & Optimization
Step 4: Best-Score Filtering 
Script: triplex_filter_max.py
Input: all_triplex_sequences_detailed.csv (from Step 3).
Function: Parallelly filters the highest-scoring triplex for each sequence to generate the best_triplex_sequences_detailed.csv.

Phase 4: Visualization & Statistics
Step 5: Detection Rate Analysis
Script: triplex_detection_stat_viz.R
Input Files: analysis_summary.csv and 5093_README.txt (Grouping File).
Description: Calculates and compares the H-DNA detection rate (percentage of reads containing triplexes) across different biological groups (Control, HCC, HBV, etc.).
Step 6: Score Distribution Analysis
Script: triplex_score_density_viz.R
Input Files: best_triplex_sequences_detailed.csv, 5093_README.txt(Grouping File).
Description: Generates density plots to analyze the stability scores of identified triplex structures, highlighting thermodynamic differences between groups.
Step 7: Length Distribution Analysis
Script: triplex_length_distribution_viz.R
Input Files: best_triplex_sequences_detailed.csv, 5093_README.txt(Grouping File).
Description: Analyzes the physical length characteristics of triplexes to determine if specific length motifs are associated with specific conditions.
Step 8: Positional Mapping (Distance to 5' End)
Script: triplex_vectorized_mapping_fast.R
Input Files: analysis_summary.csv, best_triplex_sequences_detailed.csv, 5093_README.txt (Grouping File), and Refined FASTA files (from Step 2).
Description: Uses high-speed vectorized mapping to calculate the precise physical distance from each H-DNA site to the nearest genomic 5' end, revealing fragmentomic positioning patterns.
