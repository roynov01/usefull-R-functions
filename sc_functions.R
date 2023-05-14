library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(fgsea)
library(ggplot2)
library(Seurat)
library(utils)
library(dplyr)
library(data.table)


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
  ind = unlist(lapply(df, is.numeric), use.names = FALSE)  
  numeric_df = df[,ind]
  sum_col = unname(colSums(numeric_df,na.rm = T))
  normilized_df =  as.data.frame(t(t(numeric_df) / sum_col))
  final_df = data.frame(df[,!ind],normilized_df)
  colnames(final_df) = c(colnames(df)[!ind],colnames(df)[ind])
  return (final_df)
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

GSEA_function <- function(gene_list,log2fc_ratio, geneset="hallmark",
                          organism="Mus musculus",score_type="std", qval_threshold=1,
                          nperm=1000) {
  # https://bioinformatics-core-shared-training.github.io/RNAseq_May_2020_remote/html/06_Gene_set_testing.html
  # https://bioconductor.org/packages/devel/bioc/manuals/fgsea/man/fgsea.pdf 
  # 
  #' get enriched pathways based on a list of genes
  #' 
  #' @param gene_list a vector containing genes that are significanly changed between conditions
  #' @param log2fc_ratio a vector of scores / FoldChange / Log2 Fold change / ranks et.
  #' @param geneset either "hallmark" or "KEGG", for more check ?msigdbr or go to: https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=CP:WIKIPATHWAYS
  #' Also can be a vector of category and subcategory, such as c("C2","KEGG"). 
  #' or a dataframe of two columns - first is the genesets, second is the gene
  #' @param organism either "Mus musculus" or "Homo sapiens"
  #' @param changes 
  #' @return Returns a dataframe with all relevent results. enriched in  gene set
  if (geneset=="hallmark"){
    hallmark = msigdbr(species=organism,category="H" )
    GENESET = "HALLMARK_"
  } else if (geneset=="KEGG"){
    hallmark = msigdbr(species=organism,category="C2", subcategory = "KEGG" )
    GENESET = "KEGG_"
  } else if (class(geneset)=="data.frame") {
    hallmark = geneset
    colnames(hallmark) = c("gs_name", "gene_symbol")
    GENESET = "customGeneset_"
  } else if (length(geneset) == 2){
    hallmark = msigdbr(species=organism,category=geneset[1], subcategory = geneset[2] )
    GENESET = paste(geneset[2],"_",sep="")
  } else {
    stop ("geneset problematic")
  }
  names(log2fc_ratio) = gene_list
  hallmark.genes = hallmark[,c("gs_name", "gene_symbol")]
  
  hallmark.list = hallmark.genes %>% 
    group_by(gs_name) %>% 
    summarise(all.genes = list(unique(gene_symbol))) %>% 
    deframe()
  
  enrich.hallmark = fgseaSimple(pathways = hallmark.list,
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
  
  enrich.hallmark$direction = "upregulated"
  enrich.hallmark$direction[enrich.hallmark$NES<0] = "downregulated"
  enrich.hallmark$gene_ratio = enrich.hallmark$num_of_genes/enrich.hallmark$size
  return(list(GSEA=enrich.hallmark,pathways=hallmark.list))
  
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
    # facet_grid(.~direction,scales = "free") +
    facet_wrap(.~direction,scales = "free_x") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

export_GSEA_csv = function(GSEA,output_file=NULL) {
  temp = GSEA
  temp$genes = ""
  for (i in 1:nrow(temp)){
    temp$genes[i] = paste(unlist(temp$leadingEdge[i]),collapse = '_')
  }
  temp$leadingEdge = NULL
  if (!is.null(output_file)){write.csv(temp,output_file)}
  return(temp)
}

export_GSEA_json = function(GSEA, output_file) {
  library(jsonlite)
  # output file is path+name
  sink(file =  output_file)
  l = unname(split(GSEA, seq(nrow(GSEA)))) # convert data frame rows to list
  toJSON(l, pretty=T, auto_unbox = T) # convert to JSON
  sink(file = NULL)
}

########### beautiful pie charts ##################
export_for_pie = function(genes, values, path="") {
  #' gets genes and values, transormes it into format accaptable by:
  #' https://bionic-vis.biologie.uni-greifswald.de/ 
  #' 
  #' @param genes a vector containing genes 
  #' @param values any value represent the size/portion of each gene
  #' @param path if provided - saves the data into file
  #' 
  translation = read.csv("X:\\roy\\resources\\Ensemble\\refseq_peptide_genes.csv")
  # translation = read.csv("X:\\roy\\resources\\Ensemble\\ensemble_entrez_geneId_conversion.csv")
  df = data.frame(genes,values)
  colnames(df) = c("external_gene_name","value")
   
  df_merged = merge(df,translation[,c("external_gene_name","refseq_peptide")],by="external_gene_name")

  df_merged = df_merged[df_merged$refseq_peptide != "",]
  
  dt = data.table(df_merged)
  dt = unique(dt, by = "external_gene_name")
  dt$refseq_peptide = gsub("NP_","",dt$refseq_peptide)
  dt$external_gene_name=NULL
  dt = as.data.frame(dt)
  dt = dt[,c(2,1)]
  
  if (path != ""){write.table(dt, path, sep="\t",row.names=F,col.names=F, quote = FALSE)}
  return(dt)
}

########### General R ############################
remove_duplicated_rows = function(df, column_containing_duplicates) {
  dt = data.table(df)
  dt = unique(dt, by = column_containing_duplicates)
  return (as.data.frame(dt))
}

averege_duplicated_rows = function(df, column_containing_duplicates) {
  return (aggregate( . ~ eval(parse(text=column_containing_duplicates)), df, mean, na.rm=T, na.action=NULL))
}

max_duplicated_rows = function(df, column_containing_duplicates) {
  new_df = data.frame(matrix(NA, nrow=0, ncol=ncol(df)))
  colnames(new_df) = colnames(df)
  for(i in 1:nrow(df)){
    gene = df[i,column_containing_duplicates]
    if (!(gene %in% new_df[,column_containing_duplicates])){
      new_df[nrow(new_df) + 1,] = df[i,]
      rownames(new_df)[nrow(new_df)] = gene
      next
    }
    cur_sum = sum(new_df[gene,colnames(new_df) != column_containing_duplicates])
    new_sum = sum(df[i,colnames(df) != column_containing_duplicates])
    if (new_sum > cur_sum){
      new_df[gene,] = df[i,]
    }
  }
  return(new_df)
}

########### intestinal plots and markers ############################

plot_marker = function(expression_mat, genes, celltypes=NA) {
  #' creates barplot of selected genes based on signature matrix
  #' 
  #' @param expression_mat signature matrix (rows are genes, columns are celltypes)
  #' @param genes vector of gene names
  #' @param celltypes vector of gene celltypes which would be highlighted by color
  plots = list()
  expression = expression_mat[genes,]
  expression$gene = rownames(expression)
  expression = pivot_longer(expression,cols=colnames(expression)[colnames(expression) != "gene"])
  expression$color = F
  expression$color[expression$name %in% celltypes] = T
  
  for (i in 1:length(genes)) {
    plot = expression[expression$gene==genes[i],]
    p = ggplot(plot,aes(name,value,fill=color)) + geom_bar(stat="identity")+
      coord_flip() +
      theme_bw() + 
      xlab("") + ylab("") +
      scale_fill_manual(values=c("blue","red")) +
      theme(legend.position="none") + 
      ggtitle(genes[i])
    plots[[i]] = p
  }
  return (plots)
}

plot_inna = function(genes, organism="mouse",celltypes=NA){
  #' plots expression in celltypes
  #' #' @param genes gene or genes to plot
  #' #' @param organism "human" or "mouse".
  #' #' @param celltypes vector of gene celltypes which would be highlighted by color
  genes = toupper(genes)
  if (organism=="human") {filename = "X:\\roy\\apicome\\visium_export_from_yotam\\P13cell_type_signature_matrix.csv"}
  if (organism=="mouse") {filename = "X:\\roy\\resources\\datasets\\innas_data.csv"}
  expression = read.csv(filename,row.names=1)
  # expression = expression[!grepl("^AC\\d|orf|^RP[SL]|^AP\\d|AL\\d|^LINC",rownames(expression)),]
  # expression = expression[!grepl("^Gm\\d|^BC|^Rp[sl]|Rik|^AI\\d|^AA",rownames(expression)),]
  rownames(expression) = toupper(rownames(expression))
  barplots = plot_marker(expression, genes,celltypes=celltypes)
  
  x = sapply(barplots, function(p, add){return(p+ggtitle(paste(p$lebels$title," (",organism,")",sep="")))}, add=organism)
  
  for (i in 1:length(barplots)){barplots[[i]] = barplots[[i]] + ggtitle(paste(barplots[[i]]$labels$title," (",organism,")",sep=""))}
  return(barplots)
}

find_markers = function(celltypes=NA,organism="human",celltype_name=NA,RATIO_THRESH=2,EXP_THRESH=10^-6,barplots=F) {
  #' Find and plot markers of intestine celltypes.
  #' to see available cell-types, run the function with organism variable only. then run it with celltype(s)
  #' 
  #' @param celltypes character vector of cell types. to see available cell-types, run the function with organism variable only.
  #' @param organism "human" or "mouse".
  #' @param celltype_name str. not necessary, used for title plot etc.
  #' @param RATIO_THRESH float. threshold of expression, above which markers will be returned.
  #' @param EXP_THRESH float. ratio threshold, above which markers will be returned.
  #' @param barplots bool. return a list of barplots for each marker?
  #' @return Returns vector of markers, expression matrix of the markers and plots
  if (organism=="human") {
    filename = "X:\\roy\\apicome\\visium_export_from_yotam\\P13cell_type_signature_matrix.csv"}
  if (organism=="mouse") {filename = "X:\\roy\\resources\\datasets\\innas_data.csv"}
  expression = read.csv(filename,row.names=1)
  if (typeof(celltypes)=="logical") {return(colnames(expression))}
  if (is.na(celltype_name)){celltype_name = celltypes[1]}
  
  expression = expression[!grepl("^AC\\d|orf|^RP[SL]|^AP\\d|AL\\d|^LINC",rownames(expression)),]
  expression = expression[!grepl("^Gm\\d|^BC|^Rp[sl]|Rik|^AI\\d|^AA",rownames(expression)),]
  
  chosen_cell = data.frame(expression[,celltypes])
  rownames(chosen_cell) = rownames(expression)
  colnames(chosen_cell) = celltypes
  
  other = expression[,!(colnames(expression) %in% celltypes)]
  
  chosen_avg = apply(chosen_cell,1,mean)
  other_max = apply(other,1,max)
  
  comparison = data.frame(chosen_avg,other_max)
  comparison$other_max[comparison$other_max==0] = 1e-6 # pseudonum
  
  comparison$ratio = comparison$chosen_avg/comparison$other_max
  comparison$fc_signif = comparison$ratio > RATIO_THRESH
  comparison$expressed = comparison$chosen_avg > EXP_THRESH
  comparison$signif = comparison$fc_signif & comparison$expressed
  
  plot = comparison
  plot$gene = rownames(plot)
  plot[c("chosen_avg","other_max")] = log(plot[c("chosen_avg","other_max")],10)
  p = ggplot(plot,aes(x=chosen_avg, y=other_max,colour=signif)) + geom_point() + 
    ggtitle(paste("Genes that are unique to ",celltype_name,sep=""))+
    geom_text_repel(plot[plot$signif==T,], mapping=aes(x=chosen_avg, y=other_max,label=gene),max.overlaps = 30)+
    theme_bw() + ylab("log10(max expression in other cells)") + xlab("log10(averege expression in celltype)")+
    theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_color_manual(values = c("grey", "red"))
  
  res = comparison[rownames(comparison)[comparison$signif==T],]
  res = res[order(res$ratio*res$chosen_avg,decreasing=T),c("chosen_avg","other_max","ratio")]
  celltype_markers = rownames(res)
  
  if (barplots==T) {
    barplots = plot_marker(expression, celltype_markers, celltypes)
  }
  
  return (list(markers=celltype_markers,markers_tbl=res,celltypes=celltypes,expression_mat=expression,plot=p, barplots=barplots))
}

plot_human_sections = function(gene){
  #' plots expression along the intestinal segments in human data from Rachel Zwick data in mouse or human
  #' @param gene string, gene to plot
  gene = toupper(gene)
  path = "X:\\Shalevi\\Rachel_Zwick"
  files = c("Enteroendocrine cells.xlsx",
            "Goblet cells.xlsx",
            "Tuft cells.xlsx",
            "Paneth cells.xlsx",
            "Crypt cells.xlsx",
            "Villus cells.xlsx")
  
  for(i in 1:length(files)){
    cur_path = paste(path,files[i],sep="/")
    name = sub(" cells.xlsx","",files[i])
    d = read_excel(cur_path)
    d = d[1:6]
    d_long = pivot_longer(d,cols=colnames(d[2:6]),names_to="section",values_to=name)
    if (i == 1){
      sections = d_long
    } else {
      sections[,name] = d_long[,name]
    }
  }
  sections$section = as.integer(gsub("mean ","",sections$section))
  plot = sections[sections["Gene name"] == gene,]
  plot = pivot_longer(plot,colnames(plot)[3:ncol(plot)], names_to = "celltype", values_to = "expression")
  ggplot(data=plot,aes(x=section,y=expression,color=celltype)) +
    geom_line(size=2) +
    theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    scale_color_manual(values=c("purple","lightblue","orange","darkgreen","blue","red")) +
    theme(legend.position="right") + 
    # theme(legend.position=c(0.77,0.9),legend.direction="horizontal",legend.title=element_blank()) +
    ggtitle("intestine axis") + xlab("") + ylab("")
}

plot_tip_base = function(gene,organism="mouse"){
  #' plots tip/base expression in mouse or human
  #' @param gene string, gene to plot
  #' @param organism "human" or "mouse".
  gene = toupper(gene)
  if(organism=="mouse"){
    apicome = read.csv("X:\\roy\\apicome\\analysis\\LCM_bulk\\outputs_apr23\\non_filtered/medians.csv")
  }
  if (organism=="human"){
    apicome = read.csv("X:\\roy\\apicome\\analysis\\LCM_human\\outputs_apr23\\non_filtered/medians.csv")
  }
  apicome$X=NULL
  apicome$gene = toupper(apicome$gene)
  apicome$tip = (apicome$apical_tip+apicome$basal_tip)/2
  apicome$base = (apicome$basal_base+apicome$basal_base)/2
  
  plot = apicome[apicome$gene %in% gene,c("gene","tip","base")]
  plot = pivot_longer(plot,colnames(plot)[2:3], names_to = "villus_zone", values_to = "expression")
  ggplot(data=plot,aes(x=villus_zone,y=expression)) +
    geom_bar(stat="identity",fill=c("pink","blue"),color="black",size=1) +
    theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    xlab("") + ggtitle("Villus zone") + ylab("")
}

plot_apicome = function(gene,organism="mouse"){
  #' plots apicome in mouse or human
  #' @param gene string, gene to plot
  #' @param organism "human" or "mouse".
  gene = toupper(gene)
  if(organism=="mouse"){
    apicome = read.csv("X:\\roy\\apicome\\analysis\\LCM_bulk\\outputs_apr23\\non_filtered/medians.csv")}
  if (organism=="human"){
    apicome = read.csv("X:\\roy\\apicome\\analysis\\LCM_human\\outputs_apr23\\non_filtered/medians.csv")}
  apicome$X=NULL
  apicome$gene = toupper(apicome$gene)
  apicome$apical = (apicome$apical_tip+apicome$apical_base)/2
  apicome$basal = (apicome$basal_tip+apicome$basal_base)/2
  plot = apicome[apicome$gene %in% gene,c("gene","apical","basal")]
  plot = pivot_longer(plot,colnames(plot)[2:3], names_to = "apicome", values_to = "expression")
  ggplot(data=plot,aes(x=apicome,y=expression)) +
    geom_bar(stat="identity",fill=c("lightgreen","yellow"),color="black",size=1) +
    theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    xlab("") + ggtitle("Apicome") + ylab("")
}

plot_intestines = function(gene,organism="mouse"){
  #' plots apicome, villus axis, intestinal axis and expression
  #' @param gene string, gene to plot
  #' @param organism "human" or "mouse".
  p1 = plot_inna(gene, organism=organism)
  p2 = plot_human_sections(gene)
  p3 = plot_tip_base(gene,organism=organism)
  p4 = plot_apicome(gene,organism=organism)
  p = ggarrange(p1[[1]],ggarrange(p2,p3,p4, nrow = 3),ncol = 2) 
  return(p)
}

