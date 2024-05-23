##install.packages("devtools")
##devtools::install_github("benjjneb/dada2")
library(dada2)

# On traite les trois marqueurs (16SV4, 18SV9 et ITS2) ensemble, on sort des seqtab_nochim communes
# puis on sépare les marqueurs dans QIIME2 avec $ qiime feature-classifier extract-reads

path <- "~/Desktop/mangroves_guadeloupe_rawdata/genbank_data_guadeloupe_litterbags/"

list.files(path)
setwd(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="R1_final.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2_final.fastq", full.names = TRUE))
fnFs <- sort(fnFs)
fnRs <- sort(fnRs)

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)
sample.names
# plotQualityProfile(fnFs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen=200, matchIDs=TRUE,
                     maxN=0, maxEE=c(3,3), rm.phix=TRUE,# trimLeft=15, trimRight=15,
                     compress=TRUE, multithread=TRUE)
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# plotErrors(errF, nominalQ=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # 
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# plot(table(nchar(getSequences(seqtab))))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))
sum(seqtab.nochim)/sum(seqtab)
# plot(table(nchar(getSequences(seqtab.nochim))))
write.csv(t(seqtab.nochim),"~/sync/mangroves/data_sequencages/guadeloupe/2022/dada2/seqtabnochim_guadeloupe_LB_final_mai24_notrim.csv")
uniquesToFasta(seqtab.nochim,"~/sync/mangroves/data_sequencages/guadeloupe/2022/dada2/seqtabnochim_guadeloupe_LB_final_mai24_notrim.fasta")

write.table(t(seqtab.nochim), "~/sync/mangroves/data_sequencages/guadeloupe/2022/dada2/seqtabnochim_guadeloupe_LB_final_mai24_notrim.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
uniquesToFasta(seqtab.nochim, "~/sync/mangroves/data_sequencages/guadeloupe/2022/dada2/seqtabnochim_guadeloupe_LB_final_mai24_notrim.fna", ids=colnames(seqtab.nochim))


# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track,"~/sync/mangroves/data_sequencages/guadeloupe/2022/dada2/seqtabnochim_guadeloupe_LB_track_final_mai24_notrim.txt")




######## élimination des amorces (trim)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen=200, matchIDs=TRUE,
                     maxN=0, maxEE=c(3,3), rm.phix=TRUE, trimLeft=15, trimRight=15,
                     compress=TRUE, multithread=TRUE)
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
# plotErrors(errF, nominalQ=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # 
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# plot(table(nchar(getSequences(seqtab))))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))
sum(seqtab.nochim)/sum(seqtab)
# plot(table(nchar(getSequences(seqtab.nochim))))
write.csv(t(seqtab.nochim),"~/sync/mangroves/data_sequencages/guadeloupe/2022/dada2/seqtabnochim_guadeloupe_LB_final_mai24_trim15.csv")
uniquesToFasta(seqtab.nochim,"~/sync/mangroves/data_sequencages/guadeloupe/2022/dada2/seqtabnochim_guadeloupe_LB_final_mai24_trim15.fasta")

write.table(t(seqtab.nochim), "~/sync/mangroves/data_sequencages/guadeloupe/2022/dada2/seqtabnochim_guadeloupe_LB_final_mai24_trim15.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
uniquesToFasta(seqtab.nochim, "~/sync/mangroves/data_sequencages/guadeloupe/2022/dada2/seqtabnochim_guadeloupe_LB_final_mai24_trim15.fna", ids=colnames(seqtab.nochim))


# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track,"~/sync/mangroves/data_sequencages/guadeloupe/2022/dada2/seqtabnochim_guadeloupe_LB_track_final_mai24_trim15.txt")
