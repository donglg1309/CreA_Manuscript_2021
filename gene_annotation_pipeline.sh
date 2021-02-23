#The pipeline is for summit annotation to 
for files in /Users/dongliguo/Documents/ANAlyses_folder/CreA_Summary/Novel_plots/Fig1/Summit_genes_v7/all_together_latest_liguo_v6.bed
do
cp ${files} input_summit_binding.bed
annotatePeaks.pl input_summit_binding.bed /Users/dongliguo/Documents/reference_sequences/A_nidulans/reference/A_nidulans_FGSC_A4_version_s10-m04-r03_chromosomes.fasta -gtf /Users/dongliguo/Documents/reference_sequences/A_nidulans/annotation/A_nidulans_FGSC_A4_version_s10-m04-r03_features_with_chromosome_sequences_without_tRNA.gtf -size 1 > input_summit_binding.bed_annotation_200bp
########

R CMD BATCH annotated_genes.R

########
sed 's/transcript_id    /transcript_id /g' other_genes_files.gtf > files_v1.gtf
sed 's/ ;       gene_id /; gene_id /g' files_v1.gtf > files_v2.gtf
sed 's/ ;       gene_name       /; gene_name /g' files_v2.gtf > files_v3.gtf
sed 's/ ;/;/g' files_v3.gtf > files_v4.gtf
annotatePeaks.pl input_summit_binding.bed /Users/dongliguo/Documents/reference_sequences/A_nidulans/reference/A_nidulans_FGSC_A4_version_s10-m04-r03_chromosomes.fasta -gtf files_v4.gtf -size 1 > input_summit_binding.bed_annotation_others_200bp
R CMD BATCH annotated_genes_v2.R
cat merge_part_one.bed merge_part_two.bed > merged_binding_values.bed
cat output_files_values_promoter_first.bed output_files_values_others.bed merged_binding_values.bed > merged_binding_values_output_files.bed
cp merged_binding_values_output_files.bed ${files}_annotated_files.bed
done

