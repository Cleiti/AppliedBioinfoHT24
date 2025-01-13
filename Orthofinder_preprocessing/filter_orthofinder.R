test<-read.table("/proj/naiss2024-6-218/miguel_analysis/orthol/orthofinder_mech_240125/Results_Jan25/WorkingDirectory/OrthoFinder/Results_Feb20/Orthogroups/Orthogroups.GeneCount.tsv", sep="\t", header=TRUE, row.names=1)

test[]<-lapply(test, function(x) as.numeric(as.character(x)))

test$prop_species<-apply(test[,1:292], 1, function(x) sum(x!=0)/292)

test_filtered<-test[test$Total>=100 & test$prop_species>=0.5,]
write.table(test_filtered, "/proj/uppstore2019013/nobackup/private/1000spider_master_project/orthofinder_filtered/Filtered_orthofinder_new.tsv", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)