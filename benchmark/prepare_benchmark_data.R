# Builds benchmark/data/pxd028735_qehfx_dda_precursors.tsv.gz from two
# files in the ProteoBench repository (Apache-2.0):
#
#   test/data/quant/quant_lfq_ion_DDA_QExactive/WOMBAT_stand_ion_quant_mergedproline.csv
#     precursor-level quantification (WOMBAT-P) of the Q Exactive HF-X DDA runs
#     of PXD028735, 3 x condition A and 3 x condition B
#   test/data/intermediate_files/result_performance_MaxQuant_20241216_120952.csv
#     ProteoBench's MaxQuant result for the same runs, which labels each
#     species-unique precursor HUMAN / YEAST / ECOLI
#
# Species are assigned to WOMBAT protein groups through their peptide
# sequences: a protein group is kept only if every one of its peptides that
# appears in the species-unique MaxQuant set points to the same species.
# Contaminants (Cont_) are removed.
#
# Usage: Rscript benchmark/prepare_benchmark_data.R /path/to/ProteoBench
# Tested with ProteoBench commit 2943465ef132315452991ae86e420cd8564676ce.

args <- commandArgs(trailingOnly = TRUE)
pb <- if (length(args)) args[1] else "ProteoBench"
wombat <- file.path(pb, "test/data/quant/quant_lfq_ion_DDA_QExactive",
                    "WOMBAT_stand_ion_quant_mergedproline.csv")
mq_path <- file.path(pb, "test/data/intermediate_files",
                     "result_performance_MaxQuant_20241216_120952.csv")
stopifnot(file.exists(wombat), file.exists(mq_path))

strip <- function(x) gsub("[^A-Z]", "", gsub("\\[[^]]*\\]", "", x))

mq <- utils::read.csv(mq_path, stringsAsFactors = FALSE)
mq <- mq[mq$unique == 1, ]
mq$seq <- strip(sub("\\|.*$", "", mq$precursor.ion))
seq_species <- tapply(mq$species, mq$seq, function(s) paste(sort(unique(s)), collapse = "/"))
seq_species <- seq_species[!grepl("/", seq_species)]

w <- utils::read.csv(wombat, stringsAsFactors = FALSE, na.strings = c("NA", ""))
w <- w[!grepl("Cont_", w$protein_group), ]
w$seq <- strip(w$modified_peptide)
w$species <- unname(seq_species[w$seq])

pg_species <- tapply(w$species, w$protein_group, function(s) {
  s <- unique(s[!is.na(s)])
  if (length(s) == 1) s else NA_character_
})
w$species <- unname(pg_species[w$protein_group])
n_pg <- length(unique(w$protein_group))
w <- w[!is.na(w$species), ]

ab <- grep("^abundance_", names(w), value = TRUE)
out <- data.frame(precursor = paste0(w$modified_peptide, "/", w$charge),
                  protein_group = w$protein_group, species = w$species,
                  w[ab], check.names = FALSE)
names(out) <- sub("^abundance_", "", names(out))
for (k in c("A_1", "A_2", "A_3", "B_1", "B_2", "B_3")) out[[k]][out[[k]] <= 0] <- NA

dir.create("benchmark/data", recursive = TRUE, showWarnings = FALSE)
con <- gzfile("benchmark/data/pxd028735_qehfx_dda_precursors.tsv.gz", "w")
utils::write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
close(con)

ss <- data.frame(sample_id = c("A_1", "A_2", "A_3", "B_1", "B_2", "B_3"),
                 condition = rep(c("A", "B"), each = 3))
utils::write.csv(ss, "benchmark/data/samples.csv", row.names = FALSE)

pgs <- unique(out[c("protein_group", "species")])
cat(sprintf("Kept %d of %d non-contaminant protein groups (%d precursors).\n",
            nrow(pgs), n_pg, nrow(out)))
print(table(pgs$species))
