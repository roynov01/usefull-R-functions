library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(fgsea)
library(ggplot2)


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