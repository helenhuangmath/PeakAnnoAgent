## Annotate peaks using GREAT(Bejerano lab Stanford)

library(rGREAT)
library(GenomicRanges)
library(dplyr)

bed_input <- GRanges(peaks)

job <- submitGreatJob(bed_input, rule = "oneClosest", adv_oneDistance = 2000, species="mm10")

res <- plotRegionGeneAssociationGraphs(job)
peaks$peak_name <- paste0(peaks$Chr,peaks$Start,peaks$End)
peak_annot_Great <- data.frame(res)
peak_annot_Great$peak_name <- paste0(peak_annot_Great$seqnames,peak_annot_Great$start,peak_annot_Great$end)

peak_annot_Great2 <- merge(peaks,peak_annot_Great,by.x="peak_name",by.y="peak_name") %>% dplyr::select(-peak_name)
peak_annot_Great2 <- peak_annot_Great2 %>% arrange(Chr, Start, End)
peak_annot_Great2 <- peak_annot_Great2 %>% mutate(peak_id=paste0(peakN, "_", Chr, "_", Start, "_", End))

# Save peak annotations
saveRDS(peak_annot_Great2,file='peak_annot_Great.Rds')
write.csv(peak_annot_Great2,file='peak_annot_Great.csv',quote=F)



