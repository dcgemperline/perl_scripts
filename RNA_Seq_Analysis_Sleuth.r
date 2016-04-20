library("sleuth")
base_dir <- "/home/david/GS_MS_bootstraps"
sample_id <- dir(file.path(base_dir, "output"))
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "output", id))
s2c <- read.table(file.path(base_dir, "experimental_design.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
so <- sleuth_prep(s2c, ~condition)
so <- sleuth_fit(so)
so <- sleuth_wt(so, which_beta="conditionGS")
so <- sleuth_wt(so, which_beta="conditionMS")

writeresults <- function(condition)
{
  results_table <- sleuth_results(so, condition)
  Bdistachyon_314_v3.1.annotation_info <- read.delim("~/Bdistachyon_314_v3.1.annotation_info.txt", row.names=1)
  result_table_annotated <- merge(results_table, Bdistachyon_314_v3.1.annotation_info, by.x="target_id", by.y="transcriptName", all.x=TRUE)
  result_table_annotated_subset <- subset(result_table_annotated, !is.na(pval))
  result_table_annotated_subset <- result_table_annotated_subset[order(-result_table_annotated_subset$b,result_table_annotated_subset$qval),]
  result_table_annotated_subset <- result_table_annotated_subset[which(result_table_annotated_subset$qval<=0.01),]
  write.csv(result_table_annotated_subset, file=paste(condition, "_DE_genes.csv", sep=""))
}

writeresults("conditionGS")
writeresults("conditionMS")
#sleuth_live(so)