# =====================================================================
# MR instrument strength and pleiotropy sensitivity (cis instruments)
# Per protein (exposure): nsnp, F-statistic (mean/min/%>10), IVW b/se/p,
# Cochran's Q, MR-Egger intercept, Steiger directionality, and MR-PRESSO
# (global test, number of outliers, outlier-corrected estimate).
#
# Input: harmonised exposure-outcome data (cis-pQTL instruments vs ASD GWAS),
#        taken from the saved two-sample MR workspace.
# =====================================================================
# Paths are relative to the project data root.
# Set the ASD_ROOT environment variable, or run R with that directory as the working directory.
setwd(Sys.getenv("ASD_ROOT", "."))

suppressMessages({library(data.table)})
OUT <- "output/mr_sensitivity"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# ---- load harmonised exposure/outcome data -------------------------------
# Preferred: a previously exported harmonised table. Otherwise it is extracted
# from the saved two-sample MR workspace (object `harmonized_asd`).
harm_rds <- file.path(OUT, "harmonised_cis.rds")
if (file.exists(harm_rds)) {
  dat <- as.data.frame(readRDS(harm_rds))
} else {
  mr_ws <- "pQTL/MR_cis_up.r"                 # saved MR workspace (RData)
  stopifnot(file.exists(mr_ws))
  e <- new.env(); load(mr_ws, envir = e)
  stopifnot("harmonized_asd" %in% ls(e))
  dat <- as.data.frame(get("harmonized_asd", envir = e))
  saveRDS(dat, harm_rds)                      # cache for re-runs
}
dat <- dat[is.na(dat$mr_keep) | dat$mr_keep, ]            # keep MR-usable SNPs
# coerce possibly comma-formatted numeric columns
num <- function(x) suppressWarnings(as.numeric(gsub(",","",as.character(x))))
for(cc in c("beta.exposure","se.exposure","beta.outcome","se.outcome",
            "pval.exposure","pval.outcome","samplesize.exposure","samplesize.outcome",
            "eaf.exposure","eaf.outcome")) if(cc %in% names(dat)) dat[[cc]] <- num(dat[[cc]])

have_presso <- requireNamespace("MRPRESSO", quietly=TRUE)
have_tsmr   <- requireNamespace("TwoSampleMR", quietly=TRUE)
cat("MRPRESSO:", have_presso, " TwoSampleMR:", have_tsmr, "\n")

ivw_fixed <- function(bx, by, sy){
  w <- 1/sy^2; b <- sum(w*bx*by)/sum(w*bx^2); se <- sqrt(1/sum(w*bx^2))
  list(b=b, se=se, p=2*pnorm(-abs(b/se)))
}
cochran_Q <- function(bx, by, sy, b){
  Q <- sum((by - b*bx)^2 / sy^2); df <- length(bx)-1
  list(Q=Q, df=df, p=pchisq(Q, df, lower.tail=FALSE))
}
egger_int <- function(bx, by, sy){
  s <- sign(bx); bx2 <- bx*s; by2 <- by*s                 # orient exposure positive
  w <- 1/sy^2
  fit <- lm(by2 ~ bx2, weights=w)
  co <- summary(fit)$coefficients
  list(intercept=co[1,1], se=co[1,2], p=co[1,4],
       slope=co[2,1], slope_se=co[2,2], slope_p=co[2,4])
}
# manual r^2 from p,n (TwoSampleMR get_r_from_pn style); robust to NA
r2_from_pn <- function(p, n){
  ok <- is.finite(p) & is.finite(n) & n>2
  out <- rep(NA_real_, length(p))
  p2 <- pmin(pmax(p[ok],1e-300),1)
  Fv <- qf(p2, 1, n[ok]-2, lower.tail=FALSE)
  out[ok] <- Fv/(Fv + n[ok] - 2)
  out
}
steiger_manual <- function(d){
  ne <- d$samplesize.exposure; ne[!is.finite(ne)] <- 90      # pQTL n~90
  no <- d$samplesize.outcome;  no[!is.finite(no)] <- 46351   # PGC ASD 2019
  r2e <- sum(r2_from_pn(d$pval.exposure, ne), na.rm=TRUE)
  r2o <- sum(r2_from_pn(d$pval.outcome,  no), na.rm=TRUE)
  list(r2_exp=r2e, r2_out=r2o, correct=r2e>r2o)
}

prot <- sort(unique(dat$exposure))
rows <- list()
for (g in prot){
  d <- dat[dat$exposure==g, ]
  d <- d[is.finite(d$beta.exposure) & is.finite(d$se.exposure) &
         is.finite(d$beta.outcome)  & is.finite(d$se.outcome) & d$se.exposure>0 & d$se.outcome>0, ]
  n <- nrow(d); if (n<1) next
  Fv <- (d$beta.exposure/d$se.exposure)^2
  r <- list(exposure=g, nsnp=n,
            F_mean=mean(Fv), F_min=min(Fv), pct_F_gt10=mean(Fv>10)*100)
  if (n>=2){
    iv <- ivw_fixed(d$beta.exposure, d$beta.outcome, d$se.outcome)
    q  <- cochran_Q(d$beta.exposure, d$beta.outcome, d$se.outcome, iv$b)
    r <- c(r, list(IVW_b=iv$b, IVW_se=iv$se, IVW_p=iv$p,
                   Q=q$Q, Q_df=q$df, Q_p=q$p))
  } else r <- c(r, list(IVW_b=d$beta.outcome/d$beta.exposure, IVW_se=NA, IVW_p=NA,
                        Q=NA, Q_df=NA, Q_p=NA))
  if (n>=3){
    eg <- tryCatch(egger_int(d$beta.exposure,d$beta.outcome,d$se.outcome), error=function(e) NULL)
    if(!is.null(eg)) r <- c(r, list(Egger_intercept=eg$intercept, Egger_int_se=eg$se, Egger_int_p=eg$p))
    else r <- c(r, list(Egger_intercept=NA,Egger_int_se=NA,Egger_int_p=NA))
  } else r <- c(r, list(Egger_intercept=NA,Egger_int_se=NA,Egger_int_p=NA))
  # Steiger (manual r2-from-p,n; robust for binary outcome approximation)
  st <- tryCatch(steiger_manual(d), error=function(e) list(r2_exp=NA,r2_out=NA,correct=NA))
  r <- c(r, list(Steiger_r2_exp=st$r2_exp, Steiger_r2_out=st$r2_out, Steiger_correct=st$correct))
  # MR-PRESSO
  pr <- list(PRESSO_global_p=NA, PRESSO_n_outlier=NA, PRESSO_b_corrected=NA)
  if (have_presso && n>=4){
    pres <- tryCatch({
      df <- data.frame(by=d$beta.outcome, bx=d$beta.exposure, sy=d$se.outcome, sx=d$se.exposure)
      m <- MRPRESSO::mr_presso(BetaOutcome="by", BetaExposure="bx",
              SdOutcome="sy", SdExposure="sx", data=df,
              OUTLIERtest=TRUE, DISTORTIONtest=TRUE, NbDistribution=1000,
              SignifThreshold=0.05, seed=1)
      gp <- m$`MR-PRESSO results`$`Global Test`$Pvalue
      outl <- m$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
      bc <- m$`Main MR results`$`Causal Estimate`[2]   # outlier-corrected
      list(PRESSO_global_p=gp,
           PRESSO_n_outlier=ifelse(is.null(outl)||identical(outl,"No significant outliers"),0,length(outl)),
           PRESSO_b_corrected=bc)
    }, error=function(e){ list(PRESSO_global_p=NA, PRESSO_n_outlier=NA, PRESSO_b_corrected=NA) })
    pr <- pres
  }
  r <- c(r, pr)
  rows[[g]] <- as.data.frame(r, stringsAsFactors=FALSE)
}
tab <- data.table::rbindlist(rows, fill=TRUE)
setorder(tab, -nsnp)
fwrite(tab, file.path(OUT, "MR_sensitivity_cis.csv"))
cat("\nwrote MR_sensitivity_cis.csv with", nrow(tab), "proteins\n")

# focused view: 9 prioritized + leading complement/MR proteins
focus <- c("C3","AHSG","MBL2","VASN","COL6A1","POSTN","VNN1","CFHR3","C1RL",
           "CFH","C6","F11","F2")
ft <- tab[exposure %in% focus]
fwrite(ft, file.path(OUT, "MR_sensitivity_cis_focus.csv"))
print(ft[, .(exposure,nsnp,F_mean,F_min,IVW_b,IVW_p,Q_p,Egger_int_p,Steiger_correct,PRESSO_global_p,PRESSO_n_outlier)])
cat("[done]\n")
