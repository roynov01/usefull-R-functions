library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(fgsea)
library(ggplot2)
library(Seurat)
library(utils)
library(dplyr)

########### convert ENSMBL genes to gene names #############################

ensmbl = function() {
  library(biomaRt)
  mart = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version=100)
  # check version in Ensemble archive (GRC38, MM8 etc.)
  attributes = listAttributes(mart)
  # choose which attributes you want in the table
  gene_ids = getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id","gene_biotype",
                                  "description"), mart = mart)
  return(gene_ids)
}


########### convert RDS dgecounts to CSV ###################################
dgeconts2csv = function(ensmbl_filepath="X:\\roy\\resources\\Ensemble\\ensemble_entrez_geneId_conversion.csv") {
  dge_files = choose.files(filters=matrix(c("RDS dgecounts file", "*.rds"),1, 2, byrow = TRUE))
  ensmbl = read.csv(ensmbl_filepath)
  
  for (file in dge_files) {
    df = readRDS(file)
    df = as.data.frame(as.matrix(df$exons$umicounts))
    df$gene = ensmbl$external_gene_name[match(rownames(df), ensmbl$ensembl_gene_id)]
    df = na.omit(distinct(df, gene, .keep_all=T))
    rownames(df) = df$gene
    df$gene = NULL
    write.csv(df, paste(file,".csv",sep=""))
  }
}

########### Mat-norm ##########################
matnorm = function(df) {
  sum_col = unname(colSums(df))
  normilized_df =  as.data.frame(t(t(df) / sum_col))
  return (normilized_df)
}

########### Pseudobulk #######################
sc2pseudobulk = function(seurat_object, clusters ="seurat_clusters" ,assay="RNA") {
  #' get mat-norm (RC normilized) table of the means of each cluster
  #' 
  #' @param seurat_object seurat object
  #' @param clusters name of seurat metadata that contains the cluster identity
  #' @param assay what assay to use (example: "RNA")
  #' @return Returns a dataframe with mat-norm (RC normilized) table of the means of each cluster
  df = as.data.frame(AverageExpression(seurat_object, group.by=clusters),assays=assay,slot="counts")
  sum_col = unname(colSums(df))
  normilized_df =  as.data.frame(t(t(df) / sum_col))
  normilized_df = normilized_df[,grepl(assay,colnames(normilized_df))]
  colnames(normilized_df) = gsub(paste(assay,".",sep=""),"", colnames(normilized_df))
  return (normilized_df)
}

########### DGE + volcano plot ###############
wilcox = function(row, seurat_meta, ident, celltype1, celltype2){
  w = wilcox.test(row[seurat_meta[[ident]]==celltype1], row[seurat_meta[[ident]]==celltype2], alternative="two.sided")
  p_val = w$p.value
  return(p_val)
}

DGE_seurat = function(seurat_object,celltype1,celltype2,pseudobulk,outpath=NA,ident="cluster",Thresh_GSEA=10^-4){
  # celltype2 is the control
  Idents(seurat_object) = seurat_object[[ident]]
  seurat_object = subset(seurat_object, idents=c(celltype1,celltype2))
  message(paste("DGE: ",celltype1,"_",celltype2,"\n",sep=""))
  df = as.data.frame(seurat_object@assays$RNA@data)
  pseudnum = 10^-7
  mean.1 = pseudobulk[,celltype1]
  mean.2 = pseudobulk[,celltype2]
  max_expression = pmax(mean.1, mean.2)
  ratio = (mean.1+pseudnum)/(mean.2+pseudnum)
  seurat_meta = seurat_object@meta.data
  message(" [calculating P val]")
  pval = apply(df, 1, wilcox,seurat_meta=seurat_object,ident=ident,celltype1=celltype1,celltype2=celltype2)
  qval = p.adjust(pval,method = "fdr")
  plotframe = data.frame(ratio=ratio,log2ratio=log2(ratio),max_expression=max_expression , pval=pval, qval=qval)
  plotframe=na.omit(plotframe)
  if (!is.na(outpath)){
    write.csv(x = plotframe,file = paste(outpath,"/dge_",celltype1,"_",celltype2,".csv",sep=""))
  }
  plotframe$gene = rownames(plotframe)
  mini = plotframe[plotframe$max_expression >= Thresh_GSEA,c("gene","ratio")]
  sortedRNK = mini[order(mini$ratio,decreasing = T),]
  if (!is.na(outpath)){
    savepath = paste(outpath,"/GSEA_",celltype1,"_",celltype2,".txt",sep="")
    write.table(sortedRNK, savepath, quote = F, row.names = F, col.names = T, sep = "\t")
    message(" [data saved] ")
  }
  return (list(dge=plotframe,GSEA=sortedRNK))
}

volcano_dge = function(dge, group1, group2, l2fc_threshold_up = 2, l2fc_threshold_down = -2, alpha=0.001,exp_thresh=10^-7, title="",ylimit=320) {
  dge = dge[dge$max_expression>exp_thresh,]
  dge$log10expression = log(dge$max_expression,10)
  dge$gene = rownames(dge)
  
  dge$expression = "NO"
  dge$expression[dge$log2ratio > l2fc_threshold_up & dge$qval < alpha] = "UP"
  dge$expression[dge$log2ratio < l2fc_threshold_down & dge$qval < alpha] = "DOWN"
  dge$delabel = NA
  dge$delabel[dge$expression!="NO"] = dge$gene[dge$expression!="NO"]
  pseudo_num = min(dge$qval[dge$qval>0])
  dge$qval[dge$qval==0] = pseudo_num
  plot = ggplot(data=dge, aes(x=log2ratio, y=-log10(qval), label=delabel,color=log10expression)) +
    geom_point() + 
    theme_minimal() +
    scale_color_continuous(low = "black", high = "red") +
    geom_text_repel(max.overlaps = Inf) +
    ggtitle(paste(title," - ",group1," vs ",group2,sep="")) +
    ylim(c(0,ylimit))
  return (plot)
}

########### Enricher + GSEA ##################
# tutorial at: https://www.youtube.com/watch?v=Bzu4_yDcBLY&t=428s&ab_channel=KimDill-McFarland 
Enricher <- function(gene_list, geneset="hallmark",organism="Mus musculus", log2fc_threshold=2, qval_threshold=0.05) {
  #' get enriched pathways based on a list of genes
  #' 
  #' @param gene_list a vector containing genes that are significanly changed between conditions
  #' @param geneset either "hallmark" or "KEGG", for more check ?msigdbr
  #' @param organism either "Mus musculus" or "Homo sapiens"
  #' @return Returns a dataframe with all relevent results. enriched in Hallmark gene sets (FDR < 0.05)
  
  if (geneset=="hallmark"){
    hallmark = msigdbr(species=organism,category="H" )
    GENESET = "HALLMARK_"
  }
  else if (geneset=="KEGG"){
    hallmark = msigdbr(species=organism,category="C2", subcategory = "KEGG" )
    GENESET = "KEGG_"
  }
  else {stop ("geneset need to be hallmark or KEGG")}
  
  hallmark.genes = hallmark[,c("gs_name", "gene_symbol")]
  
  enrich.hallmark = enricher(gene = gene_list, TERM2GENE = hallmark.genes)
  enrich.hallmark = enrich.hallmark@result %>%   #separate ratios into 2 columns of data
    separate(BgRatio, into=c("size.term","size.category"), sep="/") %>% 
    separate(GeneRatio, into=c("size.overlap.term", "size.overlap.category"),
             sep="/") %>% 
    #convert to numeric
    mutate_at(vars("size.term","size.category","size.overlap.term","size.overlap.category"),
              as.numeric) %>% 
    #Calculate k/K - how many genes are in the input divided by the size of the geneset (overlap)
    mutate("k.K"=size.overlap.term/size.term)
  
  enrich.hallmark$log10qval = -log(enrich.hallmark$qval,10)
  
  enrich.hallmark = enrich.hallmark %>% 
    filter(p.adjust <= qval_threshold) %>% 
    #Beautify descriptions by removing _ and HALLMARK
    mutate(Description = gsub(GENESET,"", Description),Description = gsub("_"," ", Description))
  return(enrich.hallmark)
}

plot_Enricher <- function(enrich.hallmark, qval_color="red",title="") {
  ggplot(enrich.hallmark,aes(x=reorder(Description, k.K), y=k.K,fill=log10qval)) +
    geom_col() +
    theme_classic() +
    coord_flip() +
    scale_fill_gradient(low="black", high=qval_color) +
    labs(y="Significant genes in set / Total genes in set",x="",title =title)+
    guides(fill=guide_legend(title="-log(Qval)",reverse=T))
}

GSEA_function <- function(gene_list,log2fc_ratio, geneset="hallmark",organism="Mus musculus",score_type="std", qval_threshold=0.05, nperm=1000) {
  if (geneset=="hallmark"){
    hallmark = msigdbr(species=organism,category="H" )
    GENESET = "HALLMARK_"
  }
  else if (geneset=="KEGG"){
    hallmark = msigdbr(species=organism,category="C2", subcategory = "KEGG" )
    GENESET = "KEGG_"
  }
  else {stop ("geneset need to be hallmark or KEGG")}
  names(log2fc_ratio) = gene_list
  hallmark.genes = hallmark[,c("gs_name", "gene_symbol")]
  
  hallmark.list <- hallmark.genes %>% 
    group_by(gs_name) %>% 
    summarise(all.genes = list(unique(gene_symbol))) %>% 
    deframe()
  
  enrich.hallmark <- fgseaSimple(pathways = hallmark.list,
                                 stats = log2fc_ratio,
                                 scoreType = score_type,
                                 nperm=1000)
  
  enrich.hallmark = enrich.hallmark %>% 
    filter(padj <= qval_threshold) %>% 
    mutate(pathway = gsub(GENESET,"", pathway),
           pathway = gsub("_"," ", pathway)) 
  enrich.hallmark$log10qval = -log(enrich.hallmark$padj,10)
  
  
  enrich.hallmark$num_of_genes = NA
  for (row in 1:nrow(enrich.hallmark)){enrich.hallmark$num_of_genes[row] = length(enrich.hallmark$leadingEdge[row][[1]])}
  
  enrich.hallmark$direction="upregulated"
  enrich.hallmark$direction[enrich.hallmark$NES<0] = "downregulated"
  enrich.hallmark$gene_ratio = enrich.hallmark$num_of_genes/enrich.hallmark$size
  return(enrich.hallmark)
}

plot_GSEA_bar = function(enrich.hallmark,title="GSEA") {
  ggplot(enrich.hallmark,aes(x=reorder(pathway, NES), y=NES,fill=direction)) +
    geom_col() +
    theme_classic() +
    coord_flip() +
    scale_fill_manual(values=c("blue","red")) +
    labs(y="Normalized enrichment score (NES)",x="",title =title) +
    theme(legend.position = "none")
}

plot_GSEA_dotplot = function(enrich.hallmark,title="GSEA") {
  # input is the output of "GSEA_function"
  ggplot(enrich.hallmark,aes(y=reorder(pathway, NES), x=NES)) +
    geom_point(aes(size=gene_ratio, color=direction, alpha=padj)) +
    theme_bw() +
    scale_color_manual(values=c("red","blue"))+
    scale_alpha(range = c(1, 0.2))+
    guides(size=guide_legend(title="Gene ratio"),color=F, alpha=guide_legend(title="Adjusted P-value"))+
    ylab(NULL) + xlab("Normalized enrichment score (NES)")+
    ggtitle(title) +
    facet_grid(.~direction) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

export_GSEA_json = function(GSEA, output_file) {
  # output file is path+name
  sink(file =  output_file)
  l = unname(split(GSEA, seq(nrow(GSEA)))) # convert data frame rows to list
  toJSON(l, pretty=T, auto_unbox = T) # convert to JSON
  sink(file = NULL)
}