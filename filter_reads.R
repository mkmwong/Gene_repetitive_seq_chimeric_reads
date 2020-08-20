library(plyr)
library(GenomicRanges)
library(biomaRt)

### defined constants to calculate potential ranges
INS_SIZE = 800
READ_LEN = 150

### read in command line arguments to save the output in the 
### correct directory
args <- commandArgs(trailingOnly = TRUE)
out_dir = args[1]
install_dir = args[2]

### get annotation of repeats of the type desired(Alu, L1, SVA)
### return a genomic range of those repeats
get_repeat_anno = function() {
  rep = read.table(paste0(install_dir,"/data/repeatmasker_hg19.bed"))
  alu = rep[grepl("Alu", rep$V4),]
  alu$V4 = "Alu"
  l1 = rep[grepl("L1", rep$V4),]
  l1$V4 = "L1"
  sva = rep[grepl("SVA", rep$V4),]
  sva$V4 = "SVA"
  summ = rbind(alu,l1,sva)
  summ = summ[-which(summ$V1 %in% unique(summ$V1)[25:length(unique(summ$V1))]),]
  summGR = GRanges(seqnames = summ$V1, ranges = IRanges(start = summ$V2, end = summ$V3), rep_type = summ$V4, strand = summ$V6)
  seqlevels(summGR) <- seqlevelsInUse(summGR)
  return(summGR)
}

### get gene name conversion df
### return a df with 2 columns(entrezid, and hgnc_symbol)
get_conv_names = function() {
  grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  conv_names = getBM(attributes=c("entrezgene_id","hgnc_symbol"), mart = grch37)
  colnames(conv_names)[2] = "gene_id"
  return(conv_names)
}

### filter mapped reads that have repeats at the potential insert window
### this depends on 1. whether the read mapped to genome is R1 or R2 and
### 2. whether the read is forward or reverse.
### return a GRanges of filtered reads
filter_repeat_ins = function(f_name, repGR) {
  tab = read.table(f_name, sep="")
  tab = tab[tab$V2 %in% unique(sort(tab$V2))[1:23],]
  tab$V6[which(tab$V5=="AluYa5")] = "Alu"
  tab$V6[which(tab$V5=="L1HS_3end")] = "L1"
  tab$V6[which(tab$V5=="L1HS_5end")] = "L1"
  tab$V6[which(tab$V5=="SVA_F")] = "SVA"
  if (grepl("for_R1", f_name) | grepl("for_R2", f_name)) {
    r_start = tab$V3+READ_LEN +1
    r_end = tab$V3 + READ_LEN*2 + INS_SIZE
  } else {
    r_start = tab$V3-INS_SIZE - READ_LEN +1
    r_end = tab$V3
  }
  tabGR = GRanges(seqnames = tab$V2, ranges = IRanges(start = r_start, end = r_end), readname = tab$V1, rep_name = tab$V6)
  idx = list()
  for( i in 1:length(tabGR)) {
    ol = findOverlaps(tabGR[i],repGR, type = 'any')
    if( length(ol) == 0) {
      idx = c(idx, i)
    }
    else{
      ol = ol[which(repGR[ol@to]$rep_type== tabGR[i]$rep_name)]
      if( length(ol) == 0) {
        idx = c(idx, i)
      }
    }
  }
  return(tabGR[unlist(idx)])
}

### get genic annotations
get_genic_anno = function() {
  df = read.table(paste0(install_dir,"/data/genic_annotation.bed"), header = TRUE)
  annoGR = GRanges(seqnames = df$seqnames, ranges = IRanges(start=df$start, end = df$end), 
                   gene_id = df$gene_id, type = df$type)
  seqlevels(annoGR) <- seqlevelsInUse(annoGR)
  return(annoGR)
}

### annotate reads with genic features
annotate_df = function(subtabGR, cov_names, annoGR) {
  ol = findOverlaps(subtabGR,annoGR, type = 'any')
  summ = cbind(as.data.frame(subtabGR[ol@from]),as.data.frame(annoGR[ol@to]))
  colnames(summ)[8:12] = paste0("feature_", colnames(summ)[8:12])
  summ_na = summ[is.na(summ$gene_id),]
  summ = na.omit(summ)
  summ = merge(summ, cov_names, by="gene_id")
  summ = rbind.fill(summ,summ_na)
  return(summ)
}

main = function() {
  print("Getting repeat annotations...")
  repGR = get_repeat_anno()
  print("Getting name conversions...")
  cov_names = get_conv_names()
  print("Getting genic annotations...")
  annoGR = get_genic_anno()
  print("Filtering reads...")
  for_R1 = filter_repeat_ins("coor_for_R1.out",repGR)
  for_R2 = filter_repeat_ins("coor_for_R2.out",repGR)
  rev_R1 = filter_repeat_ins("coor_rev_R1.out",repGR)
  rev_R2 = filter_repeat_ins("coor_rev_R2.out",repGR)
  print("Annotating reads...")
  for_R1_anno = annotate_df(for_R1, cov_names, annoGR)
  for_R2_anno = annotate_df(for_R2, cov_names, annoGR)
  rev_R1_anno = annotate_df(rev_R1, cov_names, annoGR)
  rev_R2_anno = annotate_df(rev_R2, cov_names, annoGR)
  print(length(for_R1_anno$seqnames))
  print(length(for_R2_anno$seqnames))
  print(length(rev_R1_anno$seqnames))
  print(length(rev_R2_anno$seqnames))
  lst = as.data.frame(matrix(nrow = 0 , ncol = length(colnames(for_R1_anno))))
  colnames(lst) = colnames(for_R1_anno)
  if (length(for_R1_anno$seqnames) > 0 ) { 
    lst = rbind(lst, for_R1_anno)
  }
  if (length(for_R2_anno$seqnames) > 0 ) { 
    lst = rbind(lst, for_R2_anno)
  }
  if (length(rev_R1_anno$seqnames) > 0 ) { 
    lst = rbind(lst, rev_R1_anno)
  }
  if (length(rev_R2_anno$seqnames) > 0 ) { 
    lst = rbind(lst, rev_R2_anno)
  }
  print("Returning final dataframe...")
  write.csv(lst,paste0(out_dir,"/fin.csv"),quote = FALSE, row.names = FALSE, col.names = FALSE)
}

main()
