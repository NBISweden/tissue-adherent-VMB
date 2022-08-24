# How to run pipe

### 00_dada2_script
input: 

      Downloaded from ENA:
      - 109 fastq files Boston_run1,   # 108 Luminal samples 1 Tissue sample
      -  92 fastq files Boston_run2,   # 92 Tissue samples
      ../data/mapping file
      ../resources/RDP_database/rdp_species_assignment_16.fa.gz
      ../resources/RDP_database/rdp_train_set_16.fa.gz
      
Output:

      "../../results/00_dada2_script.Rmd"
      "../../results/filtered/samplename.filt.fastq"  # filtered fastq files
      "../results/dada2_output/Error_rates_fwd_2019_02_25_v1.pdf"
      "../results/dada2_output/Read_quality_aggregate_fwd_2019_02_25_v1.pdf"
      "../results/dada2_output/MiSeq_2019_02_25_v1_fwd_read_tracking_log.txt"
      "../results/dada2_output/read_track_fwd_L10_R225_2019_02_25_v1.pdf"
      "../results/dada2_output/MiSeq_2019_02_25_v1_preprocessing_single_nochim.RDS"
      "../../data/phyloseq_Boston_run1.RDS"         

### 01_taxonomy_improved
input: 

      "../../data/phyloseq_boston_r1.RDS",   
      "../../data/phyloseq_boston_r2.RDS",   
      "https://raw.githubusercontent.com/ctmrbio/BVAB-and-Lac-sequences/master/BVAB_rRNA_database.fa",
      "https://github.com/ctmrbio/optivag/raw/master/database/db/16S/v0.1/optivag_db.fasta.gz",
      "https://github.com/ctmrbio/optivag/raw/master/database/db/16S/v0.1/optivag_seqinfo.csv"
      
Output:

      "../results/01_taxonomy_improved_output/Seq_ID_key.csv"
      "../results/01_taxonomy_improved_output/BVAB_Query.fasta"
      "../results/01_taxonomy_improved_output/optivag_Query.fasta"
      "../results/01_taxonomy_improved_output/BVAB_hits.tsv"
      "../results/01_taxonomy_improved_output/optivag_hits.tsv"
      "../results/01_taxonomy_improved_output/SeqID_to_taxa.csv"
      "../results/01_taxonomy_improved_output/Unnasigned_tax_report.csv"
      
      "../results/01_taxonomy_improved_output/"ASV_CVL_V3.csv"
      "../results/01_taxonomy_improved_output/ASV_tissue_B1.csv"
      "../results/01_taxonomy_improved_output/ASV_tissue_B2.csv"

### 02_data_preprocessing
input:
      
      "../results/01_taxonomy_improved_output/data/"ASV_CVL_V3.csv"
      "../results/01_taxonomy_improved_output/data/ASV_tissue_B1.csv"
      "../results/01_taxonomy_improved_output/data/ASV_tissue_B2.csv"
      
Output:
      
      "../results/02_data_preprocessing_output/ASV_Luminal_raw_counts.csv"
      "../results/02_data_preprocessing_output/ASV_Tissue_raw_counts.csv"

### 03_normalize_data
input:
      
      "../data/Raw_gene_counts_matrix.csv"
      "../results/02_data_preprocessing_output/ASV_Luminal_raw_counts.csv"
      "../results/02_data_preprocessing_output/ASV_Tissue_raw_counts.csv"
      
Output:
    
      "../results/03_normalize_data_output/datasets_all_samples.RDS"
      "../results/03_normalize_data_output/Tissue_RNAseq_normalized.csv"
      "../results/03_normalize_data_output/ASV_tissue_normalized.csv"
      "../results/03_normalize_data_output/ASV_luminal_normalized.csv"

### 04_clustering
input: 
      
      "../data/metadata.csv"
      "../results/03_normalize_data_output/datasets_all_samples.RDS"
      
Output:

      "../results/04_clustering_output/participant_SNN_graph.csv"
      "../results/04_clustering_output/bacterial_communities.csv"
      "../results/04_clustering_output/bacterial_SNN_graph.csv"

### 04_run_picrust
input: 

      "../results/01_taxonomy_improved_output/ASV_CVL_V3_B1.csv"
      
      "https://raw.githubusercontent.com/picrust/picrust2/master/picrust2/default_files/description_mapfiles/ko_info.tsv.gz"
      "https://raw.githubusercontent.com/picrust/picrust2/master/picrust2/default_files/pathway_mapfiles/KEGG_pathways_to_KO.tsv"
      "https://raw.githubusercontent.com/picrust/picrust2/master/picrust2/default_files/description_mapfiles/KEGG_pathways_info.tsv.gz"

Output:
      
      "../results/05_picrust_output/Picrust_input_Luminal.fasta"
      "../results/05_picrust_output/Picrust_input_Luminal.tsv"
      
      "../results/05_picrust_output/" # dirctory with all outfiles from picrust
      Files needed for further anlysis:
      "../results/05_picrust_output/picrust_outfiles/out_2021-01-20/pred_metagenome_unstrat.tsv"
      "../results/05_picrust_output/out_2021-01-20/KO_predicted.tsv"

### Figure1.Rmd
input: 

      "../results/03_normalize_data_output/datasets_all_samples.RDS"
      "../data/metadata.csv"
      "../../resources/Scematic figure.pdf"
      "../results/04_clustering_output/participant_SNN_graph.csv"
      "../results/04_clustering_output/bacterial_communities.csv"
      
Output:

      "./Figures/Figure 1.pdf"
      "./Suppl.Tbl/Suppl.Tbl.01.xlsx"
      "./Suppl.Tbl/alpha_diversity.xlsx"

### Figure2.Rmd
input: 

      "../results/03_normalize_data_output/datasets_all_samples.RDS"
      "../data/metadata.csv"
      "../results/04_clustering_output/participant_SNN_graph.csv"
      "../results/04_clustering_output/bacterial_communities.csv"
      "../results/05_picrust_output/picrust_outfiles/out_2021-01-20/pred_metagenome_unstrat.tsv"
      
Output:

      "./Suppl.Tbl/Suppl.Tbl.02"
      "./Figures/Figure 2.pdf"


### Figure3.Rmd
input: 

      "../results/03_normalize_data_output/datasets_all_samples.RDS"
      "../data/metadata.csv"
      "../results/02_data_preprocessing_output/ASV_Luminal_raw_counts.csv"
      "../results/02_data_preprocessing_output/ASV_Tissue_raw_counts.csv"
      
Output:

      "./Figures/Figure 3.pdf"

### Figure4-5.Rmd
input: 

      "../results/03_normalize_data_output/datasets_all_samples.RDS"
      "../data/metadata.csv"
      "../../resources/KEGG_GO_database/c5.bp.v6.2.symbols.gmt.txt"
      "../../resources/KEGG_GO_database/c2.cp.kegg.v6.2.symbols.gmt.txt"
      
Output:

      "./Suppl.Tbl/","Suppl.Tbl.04"
      "./Suppl.Tbl/","Suppl.Tbl.05"
      "./Suppl.Tbl/","Suppl.Tbl.06"
      "./Suppl.Tbl/","Suppl.Tbl.07"
      "./Suppl.Tbl/","Suppl.Tbl.08"
      "./Suppl.Tbl/","Suppl.Tbl.09"
      "./Suppl.Tbl/","Suppl.Tbl.10"
      "./Suppl.Tbl/","Suppl.Tbl.11"
      "./Suppl.Tbl/","Suppl.Tbl.12"
      "./Figures/Figure 4.pdf"
      "./Figures/Figure 5.pdf"

### Figure6.Rmd
input: 

      "../results/03_normalize_data_output/datasets_all_samples.RDS"
      "../data/metadata.csv"
      "../../data/Protein_norm_MFI_Bradley_et.al.xlsx"
      "../../data/210125_Cytokines_visit3CVL.csv"
      
Output:

      "./Suppl.Tbl/","Suppl.Tbl.15"
      "./Suppl.Tbl/","Suppl.Tbl 16"
      "./Suppl.Tbl/","Suppl.Tbl.17"
      "./Figures/Figure 6.pdf"

