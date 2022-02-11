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

      '../../data/phyloseq_boston_r1.RDS',   
      '../../data/phyloseq_boston_r2.RDS',   
      "https://raw.githubusercontent.com/ctmrbio/BVAB-and-Lac-sequences/master/BVAB_rRNA_database.fa",
      "https://github.com/ctmrbio/optivag/raw/master/database/db/16S/v0.1/optivag_db.fasta.gz",
      "https://github.com/ctmrbio/optivag/raw/master/database/db/16S/v0.1/optivag_seqinfo.csv"
      
Output:

      "../../results/Seq_ID_key.csv"
      "../../results/BVAB_Query.fasta"
      "../../results/optivag_Query.fasta"
      "../../results/BVAB_hits.tsv"
      "../../results/optivag_hits.tsv"
      "../../data/SeqID_to_taxa.csv"
      "../../results/Unnasigned_tax_report.csv"
      
      "../../data/"ASV_CVL_V3.csv"
      "../../data/ASV_tissue_B1.csv"
      "../../data/ASV_tissue_B2.csv"

### 02_data_preprocessing
input: 
      
      "../../data/Clinical_visit_2_3_updatesept28.csv"
      
      "../../data/"ASV_CVL_V3.csv"
      "../../data/ASV_tissue_B1.csv"
      "../../data/ASV_tissue_B2.csv"
      
Output:
      
      "../results/encoder.h5"
      "../results/decoder.h5"
      "../results/aen.h5"
      
      "../results/ASV_tissue_V3_normalized_batch_corrected.csv"
      "../results/ASV_CVL_V3_normalized_batch_corrected.csv"
      "../results/Tissue_RNAseq_V3_normalized.csv"
      
### 03_clustering
input: 
      
      "../../../results/datasets_all_samples.RDS"
      "../../../results/metadata_integration.csv"
      
Output:

      "../../results/patient_SNN_graph.csv"
      "../../results/bacteria_SNN_graph.csv"

### 04_run_picrust
input: 

    

Output:

    "../results/picrust_outfiles/" # dirctory with all outfiles from picrust
    Files needed for further anlysis:
    "../results/picrust_outfiles/pred_metagenome_unstrat.tsv"
    "../results/picrust_outfiles/KO_predicted.tsv"

### Figure1.Rmd
input: 

      "../../../results/datasets_all_samples.RDS"
      "../../../results/metadata_integration.csv"
      "../../../results/bacteria_SNN_graph.csv"
      "../../../results/bacterial_communities.csv"
      
Output:

      "./Suppl.Tbl/","Suppl.Tbl.01"

### Figure2.Rmd
input: 

      "../../../results/datasets_all_samples.RDS"
      "../../../results/metadata_integration.csv"
      "../../../results/bacteria_SNN_graph.csv"
      "../../../results/bacterial_communities.csv"
      "../results/picrust_outfiles/pred_metagenome_unstrat.tsv"
      
Output:

      "./Suppl.Tbl/","Suppl.Tbl.02"
      "./Figures/Figure 2.pdf"


### Figure3.Rmd
input: 

      "../../../results/datasets_all_samples.RDS"
      "../../../results/metadata_integration.csv"
      
Output:

      "./Figures/Figure 1.pdf"

### 02b_evaluating_batch_correction
input: 

      "../../results/Tissue_RNAseq_V3_normalized.csv"
      "../../results/ASV_tissue_V3_normalized_batch_corrected.csv"
      "../../results/ASV_CVL_V3_normalized_batch_corrected.csv"
      "../../results/ASV_CVL_V2_normalized_batch_corrected.csv"
      "../../results/ASV_CVL_V2_normalized_NOT_batch_corrected.csv"
      
      "../../results/batches_CVL2.csv" ??
      
Output:

### 03_integrated_datasets
input: 

      "../../results/Tissue_RNAseq_V3_normalized.csv"
      "../../results/ASV_tissue_V3_normalized_batch_corrected.csv"
      "../../results/ASV_CVL_V3_normalized_batch_corrected.csv"
      "../../results/ASV_CVL_V2_normalized_NOT_batch_corrected.csv"
      
      "../../supplementary_files/h.all.v6.2.symbols.gmt.txt" ??
      
Output:

      "../../results/datasets_all_samples.RDS"
      "../../results/patient_SNN_graph.csv"
      "../../results/metadata_integration.csv"
      
      "../../results/metadata_association_results_complete.csv"
      "../../results/metadata_association_results_filtered.csv"
      "../../results/DGE_ASV_CVL_V2_normalized_NOT_batch_corrected.csv"
      "../../results/bacterial_communities.csv"
      "../../results/taxonomy.csv"
      
      "../../results/Hallmark_NESse_list.csv"
      "../../results/Hallmark_pvalues_list.csv"

### 04_
input: 
Output:

      
Output:

### manuscript_template
input: 
Output: