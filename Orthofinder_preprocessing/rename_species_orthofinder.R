orthofinder <- read.delim("Filtered_orthofinder.tsv", header = TRUE, sep = "\t")
species <- read.delim("../renamed_to_speciesnames/rename_log.csv", header = TRUE, sep = ",")

# fixing format in old fastafile names from orthofinder
old_colnames <- colnames(orthofinder)
colnames_to_rename <- old_colnames[1:292]
typical <- colnames_to_rename[1:291]
cleaned_typical <- substr(typical, 2, nchar(typical))
l.sclopetarius <- colnames_to_rename[292]
colnames_to_rename <- c(cleaned_typical, l.sclopetarius)

# fixing format in file connecting old fastafile names to new names (species names)
species <- species[, c(1, 3)]
filenames <- species[,1]
filenames <- as.character(filenames)
filenames_cleaned <- sapply(strsplit(filenames, "/"), function(x) tail(x, 1))
filenames_cleaned_again <- gsub("-", ".", filenames_cleaned)
filenames_cleaned_again_again <- sub("\\.fasta$", "", filenames_cleaned_again)
species[,1] <- filenames_cleaned_again_again
new_names <- species[,2]
new_names <- as.character(new_names)
new_names_cleaned <- sapply(strsplit(new_names, "/"), function(x) tail(x, 1))
species[,2] <- new_names_cleaned

# renaming old fastafile names in orthofinder to new filenames (species names)
colnames_renamed <- c()
for (i in seq_along(colnames_to_rename)) {
  match_row <- which(species$file == colnames_to_rename[i])
  if (length(match_row) > 0) {
    colnames_renamed[i] <- species$new_name[match_row]
  }
}
old_colnames[1:292] <- colnames_renamed
colnames(orthofinder) <- old_colnames

write.table(orthofinder, "Filtered_orthofinder_speciesrenamed.tsv", sep = "\t", row.names = FALSE, quote = FALSE)