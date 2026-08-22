# Regional cis-pQTL signal: OUR (N=48) vs UKB-PPP (N=33,493), same locus.
# Shows how similar our data is to UKB across each risk-protein cis region.
suppressMessages({library(data.table)})
ROOT <- Sys.getenv("ASD_ROOT", ".")
OUT<-file.path(ROOT,"output/replication")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
UKB<-file.path(ROOT,"Co_localization/UKB_PPP")
WIN<-5e5   # widened to capture C1RL lead (~370kb from gene) and C3
prot<-data.table(
 gene=c("MBL2","AHSG","C3","C1RL"),uni=c("P11226","P02765","P01024","Q9NZP8"),
 chr=c(10,3,19,12),start=c(54531235,186329052,6677704,7695532),end=c(54535093,186340389,6720650,7741689),
 folder=c("MBL2/MBL2_P11226_OID30759_v1_Inflammation_II","AHSG/AHSG_P02765_OID30706_v1_Inflammation_II",
          "C3/C3_P01024_OID30776_v1_Inflammation_II","C1RL/C1RL_Q9NZP8_OID30721_v1_Inflammation_II"))

pq<-fread(file.path(ROOT,"pQTL/pQTL_re/Results/cis_pQTL_ASD_result.tsv"))
pq[,uni:=tstrsplit(gene,"_",keep=1)]; pq<-pq[uni %in% prot$uni]
pq[,c("pc","pp","pr","pa"):=tstrsplit(SNP,":")][,`:=`(pc=as.integer(pc),pp=as.integer(pp))]

res<-list()
for(i in 1:nrow(prot)){
  g<-prot$gene[i]
  # our best peptide in region
  sub<-pq[uni==prot$uni[i] & pc==prot$chr[i] & pp>=prot$start[i]-WIN & pp<=prot$end[i]+WIN]
  if(nrow(sub)){bp<-sub[,.(mp=min(`p-value`)),by=gene][order(mp)][1]$gene; sub<-sub[gene==bp]}
  our<-sub[,.(pos=pp, nlp=-log10(`p-value`), src="This study (N=48)")]
  # UKB
  fp<-file.path(UKB,prot$folder[i],sprintf("discovery_chr%d_%s.gz",prot$chr[i],basename(prot$folder[i])))
  uk<-fread(fp); uk[,c("ic","ip"):=tstrsplit(ID,":",keep=1:2)][,ip:=as.integer(ip)]
  uk<-uk[ip>=prot$start[i]-WIN & ip<=prot$end[i]+WIN,.(pos=ip, nlp=LOG10P, src="UKB-PPP (N=33,493)")]
  d<-rbind(our,uk); d[,gene:=g]; res[[g]]<-d
  log_our<-if(nrow(our))max(our$nlp) else NA;
  cat(sprintf("%s: our nSNP=%d (peak %.1f), UKB nSNP=%d (peak %.1f)\n",g,nrow(our),max(c(our$nlp,0)),nrow(uk),max(c(uk$nlp,0))))
}
D<-rbindlist(res)
fwrite(D,file.path(OUT,"regional_compare_data.csv"))
cat("saved regional_compare_data.csv\n")
