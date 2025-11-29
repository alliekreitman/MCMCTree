library(Biostrings)


fasta_sequences <- list.files(path = "/Users/alliekreitman/Documents/github_repos/MCMCTree/fasta/complete_data_patients", full.names =T)
sequences <- readDNAStringSet(fasta_sequences[1])

# create output file paths that add _RT and _pol to the sample name
output_rt <- gsub("\\.fasta$", "_RT.fasta", fasta_sequences[1])
output_protease <- gsub("\\.fasta$", "_protease.fasta", fasta_sequences[1])

# Extract protease
# protease is the first 297 nt of the gene
protease <- subseq(sequences, 1, 297) 
names(protease) <- paste0(names(protease), "_protease")
writeXStringSet(protease, output_protease)

# Extract RT
# the first part of rt region is teh rest of the gene cds
rt <- subseq(sequences, 298, 984)
names(rt) <- paste0(names(rt), "_RT")
writeXStringSet(rt, output_rt)