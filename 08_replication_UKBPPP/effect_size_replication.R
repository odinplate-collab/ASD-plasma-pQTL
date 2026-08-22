# cis-pQTL replication replacement: cis-pQTL replication, our discovery (N=48) vs UKB-PPP (N=33,493)
# Effect-size (z) concordance across shared cis SNPs for risk proteins.
# READ-ONLY on source; writes ONLY to output/replication.
suppressMessages({library(data.table)})
ROOT <- Sys.getenv("ASD_ROOT", ".")
OUT<-file.path(ROOT,"output/replication")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
UKB<-file.path(ROOT,"Co_localization/UKB_PPP")
log<-function(...)cat(format(Sys.time()),"-",...,"\n")
WIN<-1e6

prot<-data.table(
 gene=c("MBL2","AHSG","C3","C1RL"),
 uni=c("P11226","P02765","P01024","Q9NZP8"),
 chr=c(10,3,19,12),
 start=c(54531235,186329052,6677704,7695532),
 end=c(54535093,186340389,6720650,7741689),
 folder=c("MBL2/MBL2_P11226_OID30759_v1_Inflammation_II",
          "AHSG/AHSG_P02765_OID30706_v1_Inflammation_II",
          "C3/C3_P01024_OID30776_v1_Inflammation_II",
          "C1RL/C1RL_Q9NZP8_OID30721_v1_Inflammation_II"))

log("loading our cis pQTL ...")
pq<-fread(file.path(ROOT,"pQTL/pQTL_re/Results/cis_pQTL_ASD_result.tsv"))
pq[,uni:=tstrsplit(gene,"_",keep=1)]
pq<-pq[uni %in% prot$uni]
pq[,c("pchr","ppos","pref","palt"):=tstrsplit(SNP,":")]
pq[,`:=`(pchr=as.integer(pchr),ppos=as.integer(ppos),
         obeta=as.numeric(beta),ose=as.numeric(beta)/as.numeric(`t-stat`))]
pq[,oz:=obeta/ose]

load_ukb<-function(i){
  fp<-file.path(UKB,prot$folder[i],sprintf("discovery_chr%d_%s.gz",prot$chr[i],basename(prot$folder[i])))
  if(!file.exists(fp))return(NULL)
  d<-fread(fp); d[,c("ic","ip","a0","a1"):=tstrsplit(ID,":",keep=1:4)][,ip:=as.integer(ip)]
  d<-d[ip>=prot$start[i]-WIN & ip<=prot$end[i]+WIN & is.finite(BETA)&is.finite(SE)&SE>0]
  data.table(ppos=d$ip,uref=toupper(d$ALLELE0),ualt=toupper(d$ALLELE1),ubeta=d$BETA,use=d$SE,uz=d$BETA/d$SE)
}

all<-list(); lead<-list()
for(i in 1:nrow(prot)){
  g<-prot$gene[i]
  po<-pq[gene %in% pq[uni==prot$uni[i],.(mp=min(`p-value`)),by=gene][order(mp)][1]$gene & pchr==prot$chr[i]]
  po<-po[`FDR`<0.05]                      # our significant cis SNPs
  uk<-load_ukb(i); if(is.null(uk)||nrow(po)==0){log("skip",g);next}
  m<-merge(po[,.(ppos,pref,palt,obeta,ose,oz,op=`p-value`)],uk,by="ppos")
  # SNPs only (no indels) and exclude strand-ambiguous A/T, C/G pairs
  amb<-function(a,b){p<-paste0(pmin(a,b),pmax(a,b)); p %in% c("AT","CG")}
  m<-m[nchar(pref)==1 & nchar(palt)==1 & nchar(uref)==1 & nchar(ualt)==1]
  m<-m[!amb(pref,palt)]
  m[,flip:=fifelse(ualt==palt & uref==pref,1,fifelse(ualt==pref & uref==palt,-1,NA_real_))]
  m<-m[!is.na(flip)][,`:=`(uz=uz*flip, ubeta=ubeta*flip)]
  if(nrow(m)<3){log("few",g,nrow(m));next}
  m[,gene:=g]; all[[g]]<-m
  ld<-m[which.min(op)]
  lead[[g]]<-data.table(gene=g,nSNP=nrow(m),our_z=ld$oz,ukb_z=ld$uz,
                        r=cor(m$oz,m$uz),sign_concord=mean(sign(m$oz)==sign(m$uz)))
  log(sprintf("%s: nSNP=%d r=%.2f sign-concord=%.0f%%",g,nrow(m),cor(m$oz,m$uz),100*mean(sign(m$oz)==sign(m$uz))))
}
A<-rbindlist(all); LD<-rbindlist(lead)
fwrite(A,file.path(OUT,"pqtl_replication_snps.csv")); fwrite(LD,file.path(OUT,"pqtl_replication_summary.csv"))
overall_r<-cor(A$oz,A$uz); overall_sc<-mean(sign(A$oz)==sign(A$uz))
log(sprintf("OVERALL: %d SNPs, r=%.2f, sign-concordance=%.1f%%",nrow(A),overall_r,100*overall_sc))
print(LD)

# ---- plot: effect-direction concordance scatter (UKB z compressed) ----
cols<-c(MBL2="#B2182B",AHSG="#1B7837",C3="#762A83",C1RL="#2166AC")
A[,uy:=sign(uz)*log10(1+abs(uz))]      # compress huge UKB z for display
sc<-LD[,.(gene,sign_concord,r,nSNP)]
lab<-sapply(names(cols),function(g){x<-sc[gene==g]; if(nrow(x))sprintf("%s (%d SNPs, %.0f%% concordant)",g,x$nSNP,100*x$sign_concord) else paste0(g," (n/a)")})
pdf(file.path(OUT,"cis_pqtl_effectsize_replication.pdf"),width=5.6,height=5.2)
par(mar=c(4.5,4.8,3.2,1))
xl<-max(abs(A$oz))*1.05; yl<-max(abs(A$uy))*1.1
plot(NA,xlim=c(-xl,xl),ylim=c(-yl,yl),xlab="Discovery cis-pQTL z (this study, N=48)",
     ylab=expression(UKB-PPP~cis-pQTL~sign%*%log[10](1+abs(z))),
     main=sprintf("cis-pQTL replication in UKB-PPP (N=33,493)\noverall sign-concordance %.0f%% (%d SNPs)",100*overall_sc,nrow(A)))
rect(-xl,0,0,yl,col="#F5F5F5",border=NA); rect(0,-yl,xl,0,col="#F5F5F5",border=NA)  # discordant quadrants shaded grey
points(A$oz,A$uy,col=cols[A$gene],pch=16,cex=.7)
abline(h=0,v=0,col="grey60")
legend("topleft",legend=lab,col=cols,pch=16,bty="n",cex=.72)
dev.off()
cat("wrote cis_pqtl_effectsize_replication.pdf\n")

# ---- lead-SNP forest (direction + UKB significance) ----
LDf<-copy(LD)[,ukb_neglog10p:= (abs(ukb_z))^2/(2*log(10)) ]  # approx -log10p from z
pdf(file.path(OUT,"cis_pqtl_leadSNP_forest.pdf"),width=5,height=3)
par(mar=c(4.5,5,2,1))
LDf<-LDf[order(match(gene,c("MBL2","AHSG","C3","C1RL")))]
plot(NA,xlim=c(-1,1),ylim=c(.5,nrow(LDf)+.5),yaxt="n",xlab="lead cis-pQTL direction (sign of z)",ylab="",main="Lead cis-pQTL direction concordance")
axis(2,at=seq_len(nrow(LDf)),labels=rev(LDf$gene),las=1)
for(i in seq_len(nrow(LDf))){ y<-nrow(LDf)-i+1
  points(sign(LDf$our_z[i])*0.5,y,pch=17,col="black",cex=1.3)
  points(sign(LDf$ukb_z[i])*0.5,y,pch=16,col="#B2182B",cex=1.3)
}
abline(v=0,lty=2,col="grey60")
legend("bottom",legend=c("this study","UKB-PPP"),pch=c(17,16),col=c("black","#B2182B"),bty="n",horiz=TRUE,cex=.8)
dev.off()
cat("wrote cis_pqtl_leadSNP_forest.pdf\n")
writeLines("done",file.path(OUT,".repl_done"))
