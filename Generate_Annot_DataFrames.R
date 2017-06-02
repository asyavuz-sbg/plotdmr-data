refGenomes <- c("hg19", "hg38", "mm9", "mm10", "dm3", "dm6", "rn4", "rn5", "rn6")
for (i in 1:length(refGenomes)) {
  refGenome <- refGenomes[i]
  if (refGenome=="hg38"){
    library(org.Hs.eg.db)
    orgdb <- org.Hs.eg.db
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    
  } else if (refGenome=="hg19") {
    library(org.Hs.eg.db)
    orgdb <- org.Hs.eg.db
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  } else if (refGenome=="mm9") {
    library(org.Mm.eg.db)
    orgdb <- org.Mm.eg.db
    library(TxDb.Mmusculus.UCSC.mm9.knownGene)
    txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
  } else if (refGenome=="mm10") {
    library(org.Mm.eg.db)
    orgdb <- org.Mm.eg.db
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  } else if (refGenome=="dm3") {
    library(org.Dm.eg.db)
    orgdb <- org.Dm.eg.db
    library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
    txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
  } else if (refGenome=="dm6") {
    library(org.Dm.eg.db)
    orgdb <- org.Dm.eg.db
    library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
    txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
  } else if (refGenome=="rn4") {
    library(org.Rn.eg.db)
    orgdb <- org.Rn.eg.db
    library(TxDb.Rnorvegicus.UCSC.rn4.ensGene)
    txdb <- TxDb.Rnorvegicus.UCSC.rn4.ensGene
  } else if (refGenome=="rn5") {
    library(org.Rn.eg.db)
    orgdb <- org.Rn.eg.db
    library(TxDb.Rnorvegicus.UCSC.rn5.refGene)
    txdb <- TxDb.Rnorvegicus.UCSC.rn5.refGene
  } else if (refGenome=="rn6") {
    library(org.Rn.eg.db)
    orgdb <- org.Rn.eg.db
    library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
    txdb <- TxDb.Rnorvegicus.UCSC.rn6.refGene
  } else {
    stop("Unknown genome! Available genomes are hg19, hg38, mm9, mm10, dm3, dm6, rn4, rn5, rn6.")
  }
  
  txs <- transcriptsBy(txdb, by="gene")
  exs <- exonsBy(txdb, by="tx", use.names=TRUE)
  
  adhoc_df <- data.frame(chr=character(0), start=integer(0), end=integer(0), gene_ID=character(0),
                         exon_number=integer(0), strand = character(0), gene_name=character(0),
                         isoforms=integer(0))
  
  genes_of_interest <- names(txs)
  # Really dirty ad hoc dataframe creation for geneTrack parameter of plotRegions (bsseq)
  for (i in 1:length(genes_of_interest)){
    gene_id <- as.character(genes_of_interest[i])
    if (i%%100 == 1) {
      message(paste("We're at ",i," of ", length(genes_of_interest)," in transcript calculation ", "(genome: ", refGenome, ").", sep=""))
    }
    if (!is.null(txs[[gene_id]])){
      tx_order <- order(width(txs[[gene_id]]), decreasing = TRUE)
      tx_df <- data.frame()
      
      tx_name <- txs[[gene_id]][which(tx_order==1)]$tx_name
      exdf <- as.data.frame(exs[[tx_name]])
      exdf$isoforms <- rep(1, dim(exdf)[1])
      exdf$gene_ID <- rep(gene_id, dim(exdf)[1])
      exdf$gene_name <- NA
      colnames(exdf)[1] <- "chr"
      colnames(exdf)[8] <- "exon_number"
      exdf$width <- NULL
      tx_df <- rbind(tx_df, exdf[order(exdf$exon_number), colnames(adhoc_df)])
      
      adhoc_df <- rbind(adhoc_df, tx_df)
    }
  }
  if (refGenome %in% c("dm3", "dm6", "rn4")){
    adhoc_df$gene_name <- mapIds(orgdb, adhoc_df$gene_ID, "SYMBOL", "ENSEMBL", multiVals = "first")
  } else {
    adhoc_df$gene_name <- mapIds(orgdb, adhoc_df$gene_ID, "SYMBOL", "ENTREZID", multiVals = "first")
  }
  write.csv(adhoc_df, file=paste(refGenome, "_Exons.csv", sep=""), row.names=FALSE, quote=FALSE)
}