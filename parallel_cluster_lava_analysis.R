#!/usr/bin/env Rscript
# This line MUST be at the top of your R script
.libPaths(c(file.path(Sys.getenv("HOME"), "R", "libs"), .libPaths()))

# ----------------------------------------
# Local SNP-heritability and local r_g (PARALLEL VERSION)
# ----------------------------------------
library(LAVA)
library(foreach)
library(doParallel)

# 1) Set up directories and files
WORKDIR <- "path//ukbiobank"
BLOCKS  <- file.path(WORKDIR, "lava_blocks_hg19.txt")

# --- CHANGE 1: Corrected LD reference paths as you requested ---
ldref   <- list(
  AFR = file.path(WORKDIR, "ld_ref", "AFR"),
  EUR = file.path(WORKDIR, "ld_ref", "EUR")
)
traits    <- c("SystolicBP", "DiastolicBP", "PulseRate", "MaxHeartRate")
ancestry  <- c("AFR", "EUR")
outdir    <- file.path(WORKDIR, "lava_out")
dir.create(outdir, showWarnings = FALSE)

# 2) Loop over ancestries
for (anc in ancestry) {
  cat("=== LAVA for", anc, "===\n")

  ss <- file.path(WORKDIR,
                  sprintf("munged_%s_%s.sumstats.gz", traits, anc))
  if (!all(file.exists(ss))) {
    stop("Missing sumstats for ", anc)
  }

  # --- CHANGE 2: Added the required 'cases' and 'controls' columns ---
  info_df <- data.frame(
    phenotype = traits,
    cases     = 1, # Convention for continuous traits in this script
    controls  = 0, # Convention for continuous traits in this script
    filename  = ss
  )
  info_file <- file.path(WORKDIR, sprintf("lava_input_info_%s.txt", anc))
  write.table(info_df, info_file, sep = "\t", quote = FALSE, row.names = FALSE)

  lava_in <- process.input(
               input.info.file     = info_file,
               ref.prefix          = ldref[[anc]],
               sample.overlap.file = NULL
             )

  # 3) Set up the parallel environment
  n_cores <- as.numeric(Sys.getenv("SLURM_NTASKS", unset = 1))
  cat(sprintf("... Registering %d cores for parallel execution\n", n_cores))
  registerDoParallel(cores = n_cores)

  # 4) Read locus definitions and process each one in parallel
  loc_df <- read.loci(BLOCKS)

  results_list <- foreach(i = seq_len(nrow(loc_df))) %dopar% {
    locus_raw <- loc_df[i, ]
    locus     <- process.locus(locus_raw, lava_in)
    if (is.null(locus)) return(NULL)

    res <- run.univ.bivar(locus)

    meta <- data.frame(locus=locus$id, chr=locus$chr, start=locus$start, stop=locus$stop, n.snps=locus$n.snps)
    list(
      univ = cbind(meta, res$univ),
      bivar = if (!is.null(res$bivar)) cbind(meta, res$bivar) else NULL
    )
  }

  # 5) Combine results from all loci
  cat("... Combining results from all loci\n")
  results_list <- results_list[!sapply(results_list, is.null)]

  univ_df  <- do.call(rbind, lapply(results_list, function(x) x$univ))
  bivar_df <- do.call(rbind, lapply(results_list, function(x) x$bivar))

  # 6) Write out full tables
  write.table(univ_df,  file = file.path(outdir, paste0(anc, ".univ.lava")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(bivar_df, file = file.path(outdir, paste0(anc, ".bivar.lava")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  cat("-> Written:", basename(outdir), "/", anc, ".{univ,bivar}.lava\n\n")
}

cat("✓  Local h² / r_g finished\n")