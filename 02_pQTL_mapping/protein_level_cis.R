# Protein-level cis association for every measured protein, three models:
#   ASD_only                 abundance ~ genotype               (n = 48)
#   TD_only                  abundance ~ genotype               (n = 42)
#   combined_n90_dxadjusted  abundance ~ genotype + diagnosis   (n = 90)
# Abundance is the mean of peptide-wise standardized log intensities (ASD and TD
# scaled together); windows are gene span +/- 1 Mb (GRCh37). No further covariates.
# Inputs (relative to ASD_ROOT): pQTL/pQTL_re/pep_unique_ASD.txt, pep_unique_UHC.txt,
#   ASD_pQTL_imputated_sample90_rm29903_n_alt.tsv; derived/11_protein_maps/measured_proteins.tsv
# Outputs (relative to ASD_ROOT): derived/12_cis_summary_all_proteins/protein_level_cis_<model>_ALL.csv.gz
# Usage: Rscript protein_level_cis.R
suppressMessages({library(data.table)})
ROOT <- Sys.getenv("ASD_ROOT", ".")
RE   <- file.path(ROOT, "pQTL/pQTL_re")
PK   <- file.path(ROOT, "derived")
MAP  <- file.path(PK, "11_protein_maps")
OUT  <- file.path(PK, "12_cis_summary_all_proteins")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
say  <- function(...) { cat(format(Sys.time()), "-", ..., "\n"); flush.console() }
norm_id <- function(x) gsub("-", "_", x, fixed = TRUE)
WIN <- 1e6

## ---- regions -------------------------------------------------------------
M <- fread(file.path(MAP, "measured_proteins.tsv"))
M <- M[!is.na(chr) & chr %in% as.character(1:22)]
M[, `:=`(chr = as.integer(chr), lo = pmax(1, start - WIN), hi = end + WIN)]
M[, label := fifelse(is.na(gene) | gene == "", uniprot, gene)]
say("proteins with an autosomal window:", nrow(M), "on", uniqueN(M$chr), "chromosomes")

## ---- protein-level abundance, ASD and TD on a common scale ---------------
ea <- fread(file.path(RE, "pep_unique_ASD.txt")); setnames(ea, c("pid", norm_id(names(ea)[-1])))
eb <- fread(file.path(RE, "pep_unique_UHC.txt")); setnames(eb, c("pid", norm_id(names(eb)[-1])))
stopifnot(identical(ea$pid, eb$pid))
ea[, up := sub("_.*", "", pid)]
sa <- setdiff(names(ea), c("pid", "up")); sb <- setdiff(names(eb), "pid")
PROT <- list()
for (u in unique(M$uniprot)) {
  r <- which(ea$up == u); if (!length(r)) next
  Za <- NULL; Zb <- NULL
  for (j in r) {
    ya <- as.numeric(ea[j, ..sa]); yb <- as.numeric(eb[j, ..sb])
    v <- c(ya, yb); mu <- mean(v, na.rm = TRUE); s <- sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) next
    Za <- rbind(Za, (ya - mu)/s); Zb <- rbind(Zb, (yb - mu)/s)
  }
  if (is.null(Za)) next
  v <- c(colMeans(Za, na.rm = TRUE), colMeans(Zb, na.rm = TRUE))
  names(v) <- c(sa, sb)
  PROT[[u]] <- v
}
say("protein-level values built for", length(PROT), "proteins")
M <- M[uniprot %in% names(PROT)]
SAMPLES <- c(sa, sb)
DX <- c(rep(1L, length(sa)), rep(0L, length(sb)))      # 1 = ASD

## ---- genotypes: one pass, keep every variant in any window ---------------
gf  <- file.path(RE, "ASD_pQTL_imputated_sample90_rm29903_n_alt.tsv")
con <- file(gf, "r")
hdr <- strsplit(readLines(con, 1), "\t")[[1]]
gid <- norm_id(hdr[-1])
idx <- match(SAMPLES, gid); stopifnot(!anyNA(idx))
chunks <- list(); ids <- list(); nline <- 0L
repeat {
  l <- readLines(con, 20000); if (!length(l)) break
  nline <- nline + length(l)
  parts <- tstrsplit(l, "\t", fixed = TRUE); vid <- parts[[1]]
  sp <- tstrsplit(vid, ":", fixed = TRUE)
  ch <- suppressWarnings(as.integer(sp[[1]])); po <- suppressWarnings(as.integer(sp[[2]]))
  hit <- rep(FALSE, length(vid))
  for (i in seq_len(nrow(M)))
    hit <- hit | (!is.na(ch) & ch == M$chr[i] & !is.na(po) & po >= M$lo[i] & po <= M$hi[i])
  w <- which(hit); if (!length(w)) next
  m <- suppressWarnings(matrix(as.numeric(unlist(lapply(parts[-1], `[`, w))), nrow = length(w)))
  chunks[[length(chunks) + 1L]] <- m[, idx, drop = FALSE]
  ids[[length(ids) + 1L]] <- vid[w]
  if (nline %% 1000000 == 0) say("  scanned", nline, "variants, kept", sum(lengths(ids)))
}
close(con)
GT <- do.call(rbind, chunks); rownames(GT) <- unlist(ids)
rm(chunks, ids); gc(verbose = FALSE)
VP <- data.table(SNP = rownames(GT),
                 chr = as.integer(tstrsplit(rownames(GT), ":")[[1]]),
                 pos = as.integer(tstrsplit(rownames(GT), ":")[[2]]),
                 other_allele = tstrsplit(rownames(GT), ":")[[3]],
                 effect_allele = tstrsplit(rownames(GT), ":")[[4]])
say(sprintf("genotype matrix: %d variants x %d samples (%.1f GB)",
            nrow(GT), ncol(GT), as.numeric(object.size(GT))/2^30))

## ---- association, chromosome by chromosome ------------------------------
fit_one <- function(y, x, grp) {
  use <- is.finite(y) & is.finite(x)
  if (sum(use) < 20) return(NULL)
  xv <- x[use]; if (sd(xv) == 0) return(NULL)
  f <- if (is.null(grp)) summary(lm(y[use] ~ xv))$coefficients
       else              summary(lm(y[use] ~ xv + grp[use]))$coefficients
  if (nrow(f) < 2) return(NULL)
  list(beta = f[2,1], se = f[2,2], p = f[2,4], n = sum(use), eaf = mean(xv)/2)
}
TAGS <- c("ASD_only", "TD_only", "combined_n90_dxadjusted")
COLS <- list(ASD_only = which(DX == 1L), TD_only = which(DX == 0L),
             combined_n90_dxadjusted = seq_along(DX))
con_out <- lapply(TAGS, function(tg) {
  f <- gzfile(file.path(OUT, sprintf("protein_level_cis_%s_ALL.csv.gz", tg)), "w")
  writeLines("protein,uniprot,SNP,chr,pos,effect_allele,other_allele,beta,se,p,eaf,maf,n", f)
  f
})
names(con_out) <- TAGS
ntest <- setNames(integer(length(TAGS)), TAGS)

for (cc in sort(unique(M$chr))) {
  sel_v <- which(VP$chr == cc)
  if (!length(sel_v)) next
  Gc <- GT[sel_v, , drop = FALSE]; Vc <- VP[sel_v]
  Mc <- M[chr == cc]
  for (i in seq_len(nrow(Mc))) {
    u <- Mc$uniprot[i]; y_all <- PROT[[u]][SAMPLES]
    k <- which(Vc$pos >= Mc$lo[i] & Vc$pos <= Mc$hi[i])
    if (!length(k)) next
    for (tg in TAGS) {
      cs  <- COLS[[tg]]
      grp <- if (tg == "combined_n90_dxadjusted") DX[cs] else NULL
      y   <- y_all[cs]
      res <- lapply(k, function(kk) {
        r <- fit_one(y, Gc[kk, cs], grp); if (is.null(r)) return(NULL)
        sprintf("%s,%s,%s,%d,%d,%s,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%d",
                Mc$label[i], u, Vc$SNP[kk], Vc$chr[kk], Vc$pos[kk],
                Vc$effect_allele[kk], Vc$other_allele[kk],
                r$beta, r$se, r$p, r$eaf, min(r$eaf, 1 - r$eaf), r$n)
      })
      res <- unlist(res)
      if (length(res)) { writeLines(res, con_out[[tg]]); ntest[tg] <- ntest[tg] + length(res) }
    }
  }
  rm(Gc, Vc); gc(verbose = FALSE)
  say(sprintf("chr%-2d done | %d proteins | cumulative tests %s", cc, nrow(Mc),
              paste(sprintf("%s=%d", TAGS, ntest), collapse = "  ")))
}
for (f in con_out) close(f)

writeLines(c(
  "Protein-level cis association for every measured protein.",
  sprintf("Proteins: %d   Genotype matrix: %d variants x %d samples", nrow(M), nrow(GT), ncol(GT)),
  "",
  "Models:",
  "  ASD_only                 abundance ~ genotype               (n = 48)",
  "  TD_only                  abundance ~ genotype               (n = 42)",
  "  combined_n90_dxadjusted  abundance ~ genotype + diagnosis   (n = 90)",
  "",
  "Abundance: peptide-wise common z-score (ASD and TD scaled together), averaged",
  "per sample. effect_allele is the allele the dosage counts (the alt allele of",
  "the chr:pos:ref:alt identifier). Windows: gene span +/- 1 Mb, GRCh37.",
  "No further covariates are included."),
  file.path(OUT, "MODEL_NOTES.txt"))
say("done |", paste(sprintf("%s=%d", TAGS, ntest), collapse = "  "))
