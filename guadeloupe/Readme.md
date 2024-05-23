

DNA extraction

Samples were sent freezed to metropolitan France, without thawing, then were freeze-dried and crushed to powder in an obsidian mortar, carefully cleaned with an alcoholic tissue between each sample processing. Total genomic DNA from 50mg of dried samples and a standard microbial community (zymoBIOMICS Microbial Community Standard D6300, by ZYMO RESEARCH), here named “Ze”, were extracted using the NucleoSpin Soil kit (Macherey-Nagel) with a final elution volume of 50 µl following the manufacturer instructions. After this DNA extraction of samples and Ze, nucleic acid yield and purity were checked using a Nanodrop spectrophotometer (Thermo Fisher Scientific) and the concentration of each sample was equalized to final concentration of 10ng.µl-1 on a PCR plate of 96 wells.


2.3	Illumina library

In order to limit PCR biases, the first round of PCR consisted in 3 PCR replicates per sample, targeting the DNA coding for the V4-V5 hypervariable region of 16S RNA ribosomal with degenerate primers (Parada et al. 2016): 515F (GTGYCAGCMGCCGCGGTAA) and 926R (CCGYCAATTYMTTTRAGTTT). Two other primer pairs (18SV9 Fw: CCCTGCCHTTTGTACACAC, Rv: CCTTCYGCAGGTTCACCTAC; and ITS2 Fw: GTGAATCATCGAATCTTTGAA, Rv: TCCTCCGCTTATTGATATGC) were amplified, added to libraries and sequenced together with 16S-V4V5. Each primer was flanked in its 5’-end by a nucleotide sequence used for indexing at a later step, according to a protocol proposed by Nag et al. (2017). At this stage, 2 additional PCR blanks were done with water instead of extracted DNA. Each 12,5 µl reaction mix contained 1 µl of DNA (~10ng.µl-1), 0,25 µl of forward primer, 0,25 µl of reverse primer (10nM), 6,25µl of 2✕ Promega Green Master mix G2, 4,25µl of milliQ water. The PCR cycles consisted of of initial denaturing for 2 min at 94°C, followed by 30 cycles (denaturation 30 s at 94°C, hybridization 30 s at 51°C, elongation 45 s at 72 °C) and a final elongation during 5 min at 72°C. First PCR products were verified by electrophoresis on 1% agarose gel, re-amplified if negative until they were positive. Each PCR triplicate was pooled into one before the indexing PCR. Indexation PCR was realized in a 27.5 µl reaction mix containing 2 µl of first PCR products, 5 µl of reverse and forward index, 12,5µl of NEB Q5 2X mix and 8µl of milliQ water. This second PCR consisted of a initial denaturing for 30s at 98°C, followed by 30 cycles (denaturation 20s at 98°C, hybridization 20s at 60 °C, elongation 10s at 72°C) and final elongation 10s at 72°C. At this stage, one PCR blank was added with water instead of first PCR products. All indexed samples were pooled into a single low-bind tube and purified with magnetic beads (Nucleomag, Macherey Nagel, 1:1 ratio). Size range of final PCR products was verified by electrophoresis (Agilent BioAnalyzer, High-sensitivity), with an waited size peak around 420bp, then pooled in a final library, and sequenced on an Illumina MiSeq (one Miseq Reagent v3 kit 600 cycles and one nano MiSeq Reagent kit v2 kit 500 cycles for resequencing) in the Concarneau marine station (MNHN) to output demultiplexed fastq files.


Raw sequences (format fastq) are available here :
https://www.ncbi.nlm.nih.gov/sra/PRJNA1113982

Scripts in this repository
- dada2_mangroves_guadeloupe_LB.R : dereplication and merging of Illumina paired reads, chimeras removal.
- qiime2_guadeloupe_LB_analyse_2022.ipynb : separation of the different genetic markers, production of ASV and OTU tables.
