# to import this file, add: source("X:\\roy\\resources\\sc_functions.R")

########### convert ENSMBL genes to gene names #############################

ensmbl = function() {
  # Ensmbl API
  library(biomaRt)
  mart = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl", version=100)
  # check version in Ensemble archive (GRC38, MM8 etc.)
  attributes = listAttributes(mart)
  # choose which attributes you want in the table
  gene_ids = getBM(attributes = c("ensembl_gene_id","external_gene_name","entrezgene_id","gene_biotype",
                                  "description"), mart = mart)
  return(gene_ids)
}

########### SCRBseq ###################################
dgeconts2csv = function(organism) {
  library(dplyr)
  if(organism=="mouse"){ensmbl_filepath="X:\\roy\\resources\\Ensemble\\ensemble_entrez_geneId_conversion.csv"
  } else if (organism=="human"){"X:\\roy\\resources\\Ensemble\\human_ensemble_conversion.csv"}
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

scrb_zumi_to_umitable = function(
  input_dir = NULL,
  meta_path = NULL,
  output_dir = NULL,
  organism = "mouse",
  name_col = "Sample_name",
  pool_col = "pool",
  barcode_col = "BC seq"
){
  #' get umi_table from zUMI output (dgecounts.rds)
  #' @param input_dir folder containing X.dgecounts.rds, or NULL to get filedialog
  #' @param meta_path name of seurat metadata that contains the cluster identity
  #' @param output_dir optional, where to save the umi table
  #' @param organism human or mouse
  #' @param name_col,pool_col,barcode_col how are the columns in the metadata called (sample name, pool name, barcode)
  #' @return umi table, barplots of QC, and metadata table
  library(dplyr)
  library(readxl)
  library(tidyr)
  library(ggplot2)
  if (is.null(input_dir)){input_dir = choose.dir(caption="Choose input directory (zUMI output with dgecounts.rds files)")}
  if (is.null(meta_path)){meta_path = choose.files(caption="Choose metadata excel file",multi=F,filters=matrix(c("Metadata EXCEL", "*.xlsx"),1, 2, byrow = TRUE))}
  if(organism=="mouse"){ensmbl_filepath="X:\\roy\\resources\\Ensemble\\ensemble_entrez_geneId_conversion.csv"
  } else if (organism=="human"){ensmbl_filepath="X:\\roy\\resources\\Ensemble\\human_ensemble_conversion.csv"}
  meta = read_excel(meta_path,col_types="text")
  ensmbl = read.csv(ensmbl_filepath)
  dge_files = list.files(path = input_dir)
  dge_files = dge_files[grepl(".dgecounts.rds$",dge_files)]
  
  pools = gsub(".dgecounts.rds","",sapply(strsplit(dge_files, "\\\\"), function(x) tail(x, 1)))
  dataframes = vector("list", length = length(dge_files))
  plots = vector("list", length = length(dge_files)+1) 
  for (i in 1:length(dge_files)) { # create dataframes from RDS files, plot barplots of each pool
    file = paste(input_dir,dge_files[i],sep="/")
    df = readRDS(file)
    df = as.data.frame(as.matrix(df$exons$umicounts))
    df$gene = ensmbl$external_gene_name[match(rownames(df), ensmbl$ensembl_gene_id)]
    df = na.omit(distinct(df, gene, .keep_all=T))
    rownames(df) = df$gene
    df_long = pivot_longer(df, cols=-gene,names_to="barcodes")
    df_long = aggregate(value ~ barcodes, data = df_long, FUN = sum)
    df_long$value = log(df_long$value,10)
    p = ggplot(df_long,aes(barcodes,value,fill=value)) + 
      geom_bar(stat="identity") + 
      theme_bw() + labs(x="",y="log10(sum reads)",title=pools[i]) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(legend.position="none") +
      scale_fill_gradient(low = "blue", high = "red")
    plots[[i]] = p
    dataframes[[i]] = df
  }
  for (i in 1:length(dataframes)){ # for each pool, select only relevent barcodes, merge to one final table
    df = dataframes[[i]]
    cur_pool = meta[meta[,pool_col]==pools[i],c(name_col,barcode_col)]
    colnames(cur_pool) = c("name","barcode")
    results_cur = df[,cur_pool$barcode]
    colnames(results_cur) = cur_pool$name
    results_cur$gene = rownames(results_cur)
    if (i == 1){umi_table = results_cur
    } else {umi_table = merge(umi_table,results_cur,by="gene",all=T)}
  }
  umi_table[is.na(umi_table)] = 0
  umi_table = select(umi_table,all_of(meta[[name_col]])) # rearrange the order of columns based on metadata
  
  # create barplot of samples:
  bar = as.data.frame(colSums(umi_table[,colnames(umi_table)!="gene"]))
  colnames(bar) = "sum_umi"
  bar$name = factor(rownames(bar))
  p = ggplot(bar,aes(name,sum_umi)) + 
    geom_bar(stat="identity") + 
    theme_bw() + labs(x="",y="sum reads",title="final_umi_table") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position="none") +
    scale_x_discrete(limits = bar$name)
  plots[[length(plots)]] = p
  
  if ('gene' %in% colnames(umi_table)){
    rownames(umi_table) = umi_table$gene
    umi_table$gene = NULL
  }
  umi_table = as.data.frame(lapply(umi_table, as.integer),row.names=rownames(umi_table))
  if (!is.null(output_dir)){ # save
    write.csv(umi_table,paste0(output_dir,"/umi_table.csv"))
    write.csv(meta,paste0(output_dir,"/metadata.csv"))
  }
  return(list(umi_table=umi_table,plots=plots,meta=meta))
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

########### Seurat #######################
sc2pseudobulk = function(seurat_object, clusters ="seurat_clusters" ,assay="RNA") {
  #' get mat-norm (RC normilized) table of the means of each cluster
  #' 
  #' @param seurat_object seurat object
  #' @param clusters name of seurat metadata that contains the cluster identity
  #' @param assay what assay to use (example: "RNA")
  #' @return Returns a dataframe with mat-norm (RC normilized) table of the means of each cluster
  library(Seurat)
  
  df = as.data.frame(AverageExpression(seurat_object, group.by=clusters),assays=assay,slot="counts")
  sum_col = unname(colSums(df))
  normilized_df =  as.data.frame(t(t(df) / sum_col))
  normilized_df = normilized_df[,grepl(assay,colnames(normilized_df))]
  colnames(normilized_df) = gsub(paste(assay,".",sep=""),"", colnames(normilized_df))
  return (normilized_df)
}

wilcox = function(row, seurat_meta, ident, celltype1, celltype2){
  w = wilcox.test(row[seurat_meta[[ident]]==celltype1], row[seurat_meta[[ident]]==celltype2], alternative="two.sided")
  p_val = w$p.value
  return(p_val)
}

DGE_seurat = function(seurat_object,celltype1,celltype2,pseudobulk,outpath=NA,ident="cluster",Thresh_GSEA=10^-4){
  # celltype2 is the control
  library(Seurat)
  
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
  # input dge is the output of DGE_seurat()
  library(ggrepel)
  library(ggplot2)
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

VlnPlot1 = function(seuratObj, gene, shape="line", funct=mean,color="black",colors=NULL){
  #' plots a modified violin (with mean/median and transparent points)
  #' shape - "line" or "dot"
  #' colors = vector of colors for the violin
  library(Seurat)
  library(ggplot2)
  if (shape == "line") {size=25
  shape=95}
  if (shape == "dot") {size=4
  shape=16}
  
  p = VlnPlot(seuratObj, features=gene, pt.size = 0.1, cols=colors) 
  p$layers[[2]]$aes_params$alpha = 0.2
  p = p + theme(axis.title.x = element_blank()) +
    stat_summary(fun.y = funct, geom='point', size = size, colour = color,shape = shape)+
    NoLegend()
  return(p)
}

VlnPlot2 = function(seuratObj, genes, shape="line", funct=mean,color="black",colors=NULL){
  library(Seurat)
  library(gridExtra)
  library(ggplot2)
  
  if (shape == "line") {size=25
  shape=95}
  if (shape == "dot") {size=4
  shape=16}
  plots = vector(mode = "list", length = length(genes))
  for (i in 1:length(genes)){
    gene = genes[i]
    p = VlnPlot(seuratObj, features=gene, pt.size = 0.1, cols=colors) 
    p$layers[[2]]$aes_params$alpha = 0.2
    p = p + theme(axis.title.x = element_blank()) +
      stat_summary(fun.y = funct, geom='point', size = size, colour = color,shape = shape)+
      NoLegend()    
    if (i>1){p = p + theme(axis.title.y = element_blank())}
    plots[[i]] = p
  }
  nRow = floor(sqrt(length(plots)))
  combined = do.call("grid.arrange", c(plots, nrow=nRow))
  return(combined)
}
########### Enricher + GSEA ##################
# tutorial at: https://www.youtube.com/watch?v=Bzu4_yDcBLY&t=428s&ab_channel=KimDill-McFarland 
Enricher = function(gene_list, geneset="hallmark",organism="Mus musculus", log2fc_threshold=2, qval_threshold=0.05) {
  #' get enriched pathways based on a list of genes
  #' 
  #' @param gene_list a vector containing genes that are significanly changed between conditions
  #' @param geneset either "hallmark" or "KEGG", for more check ?msigdbr
  #' @param organism either "Mus musculus" or "Homo sapiens"
  #' @return Returns a dataframe with all relevent results. enriched in Hallmark gene sets (FDR < 0.05)
  library(msigdbr)
  library(clusterProfiler)
  library(fgsea)
  library(dplyr)
  
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

plot_Enricher = function(enrich.hallmark, qval_color="red",title="") {
  ggplot(enrich.hallmark,aes(x=reorder(Description, k.K), y=k.K,fill=log10qval)) +
    geom_col() +
    theme_classic() +
    coord_flip() +
    scale_fill_gradient(low="black", high=qval_color) +
    labs(y="Significant genes in set / Total genes in set",x="",title =title)+
    guides(fill=guide_legend(title="-log(Qval)",reverse=T))
}

GSEA_function = function(gene_list,log2fc_ratio, geneset="hallmark",
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
  library(msigdbr)
  library(fgsea)
  library(clusterProfiler)
  library(dplyr)

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
  library(ggplot2)
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
  library(ggplot2)
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
time_it = function(expr) {
#' executes the expression, but also measures how much time the execusion took
  start_time = Sys.time()  
  result = eval(expr)  
  end_time = Sys.time()  
  execution_time = end_time - start_time  
  cat("Execution Time:", execution_time, "\n")  
  return(result)  
}

size_of = function(object){
  #' returns size of object, in mb
  size = object.size(object)
  size_mb = size / (1024^2)
  return (as.integer(size_mb))
}

remove_duplicated_rows = function(df, column_containing_duplicates) {
  library(data.table)
  dt = data.table(df)
  dt = unique(dt, by = column_containing_duplicates)
  return (as.data.frame(dt))
}

averege_duplicated_rows = function(df, column_containing_duplicates) {
  return (aggregate( . ~ eval(parse(text=column_containing_duplicates)), df, mean, na.rm=T, na.action=NULL))
}

max_duplicated_rows = function(df, column_name) {
  df$medians = apply(df[sapply(df, is.numeric)],1,median)
  df = df[order(df$medians, decreasing = TRUE), ]
  res = df[!duplicated(df[,column_name]), ]
  res$medians = NULL
  return (res)
}

scatter = function(data, x_var, y_var, gene_var,title="") {
  #' data - dataframe with at least 3 columns
  #' x_var,y_var = values to scatter (column names)
  #' gene_var - the column name that contains the gene names
  #' title - title for the plot
  library(ggplot2)
  library(plotly)
  library(shiny)
  
  # Create ggplot2 scatter plot
  plot = ggplot(data, aes_string(x = x_var, y = y_var, text = gene_var)) +
    geom_point() + 
    ggtitle(title) + 
    theme_bw()
  # Convert ggplot2 plot to plotly
  plotly_plot = ggplotly(plot, tooltip = "text")
  # Define a Shiny app
  ui = shinyUI(fluidPage(plotlyOutput("scatter_plot")))
  server = shinyServer(function(input, output) {
    output$scatter_plot <- renderPlotly({
      plotly_plot %>%
        event_register("plotly_click") %>%
        onRender("function(el, x) {
                    el.on('plotly_click', function(d) {
                      var selectedGene = d.points[0].data.text[d.points[0].pointNumber];
                      var data = x[0].data;
                      for (var i = 0; i < data.length; i++) {
                        if (data[i].text === selectedGene) {
                          Plotly.restyle(el.id, 'marker.color', 'rgba(255, 0, 0, 0.8)', [i]);
                        } else {
                          Plotly.restyle(el.id, 'marker.color', 'rgba(31, 119, 180, 0.8)', [i]);}} });}")})})
  # Run the Shiny app
  shinyApp(ui = ui, server = server)
}

########### intestinal plots and markers ############################

#load apicome data:
load_human_apicome = function() {
  if (!exists("human_apicome")){
    human_apicome <<- read.csv("X:\\roy\\apicome\\analysis\\LCM_human\\outputs_june23\\both_batches\\filtered_protein_coding/medians.csv",row.names=1)
  } else {return(human_apicome)}
}

load_mouse_apicome = function() {
  if (!exists("mouse_apicome")){
    mouse_apicome = read.csv("X:\\roy\\apicome\\analysis\\LCM_bulk\\outputs_june23/non_filtered/medians.csv",row.names=1)
    mouse_apicome$gene = toupper(mouse_apicome$gene)
    mouse_apicome <<- mouse_apicome
  } else {return(mouse_apicome)}
}


plot_roy_marker = function(expression_mat, genes, celltypes=NA) {
  #' creates barplot of selected genes based on signature matrix
  #' 
  #' @param expression_mat signature matrix (rows are genes, columns are celltypes)
  #' @param genes vector of gene names
  #' @param celltypes vector of gene celltypes which would be highlighted by color
  library(ggplot2)
  library(tidyr)
  
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


load_human_expression = function(){
  # loads Yotams human intestines expression matrix
  if (!exists("expression_human")){
  expression_human <<- read.csv("X:\\roy\\apicome\\visium_export_from_yotam\\P13cell_type_signature_matrix.csv",row.names=1)
  } else {return(expression_human)}
}

load_mouse_expression = function(){
  # loads Innas expression matrix (mouse intestine)
  if (!exists("expression_mouse")){
    expression_mouse = read.csv("X:\\roy\\resources\\datasets\\innas_data.csv",row.names=1)
    rownames(expression_mouse) = toupper(rownames(expression_mouse))
    expression_mouse <<- expression_mouse
  } else {return(expression_mouse)}
}

plot_roy_inna = function(genes, organism="mouse",celltypes=NA){
  #' plots expression in celltypes
  #' #' @param genes gene or genes to plot
  #' #' @param organism "human" or "mouse".
  #' #' @param celltypes vector of gene celltypes which would be highlighted by color
  library(ggplot2)
  
  genes = toupper(genes)
  if (organism=="human") {expression = load_human_expression()}
  if (organism=="mouse") {expression = load_mouse_expression()}
  barplots = plot_roy_marker(expression, genes,celltypes=celltypes)
  x = sapply(barplots, function(p, add){return(p+ggtitle(paste(p$lebels$title," (",organism,")",sep="")))}, add=organism)
  for (i in 1:length(barplots)){barplots[[i]] = barplots[[i]] + ggtitle(paste(barplots[[i]]$labels$title," (",organism,")",sep=""))}
  return(barplots)
}

find_markers_intestines = function(celltypes=NA,organism="human",celltype_name=NA,RATIO_THRESH=2,EXP_THRESH=10^-6,barplots=F) {
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
  library(ggplot2)
  library(ggrepel)
  
  if (organism=="human") {expression = load_human_expression()}
  if (organism=="mouse") {expression = load_mouse_expression()}
  if (typeof(celltypes)=="logical") {return(colnames(expression))}
  if (is.na(celltype_name)){celltype_name = celltypes[1]}
  
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
  
  if (barplots==T) {barplots = plot_roy_marker(expression, celltype_markers, celltypes)}
  
  return (list(markers=celltype_markers,markers_tbl=res,celltypes=celltypes,expression_mat=expression,plot=p, barplots=barplots))
}

find_markers = function(signature_matrix,celltypes,RATIO_THRESH=2,EXP_THRESH=10^-6,barplots=F) {
  #' Find and plot markers of intestine celltypes.
  #' to see available cell-types, run the function with organism variable only. then run it with celltype(s)
  #' 
  #' @param signature_matrix colnames are celltypes, rows are genes 
  #' @param celltypes str. not necessary, used for title plot etc.
  #' @param RATIO_THRESH float. threshold of expression, above which markers will be returned.
  #' @param EXP_THRESH float. ratio threshold, above which markers will be returned.
  #' @param barplots bool. return a list of barplots for each marker?
  #' @return Returns vector of markers, expression matrix of the markers and plots
  library(ggplot2)
  library(ggrepel)
  chosen_cell = data.frame(signature_matrix[,celltypes])
  rownames(chosen_cell) = rownames(signature_matrix)
  colnames(chosen_cell) = celltypes
  other = signature_matrix[,!(colnames(signature_matrix) %in% celltypes)]
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
    ggtitle(paste("Genes that are unique to ",celltypes,sep=""))+
    geom_text_repel(plot[plot$signif==T,], mapping=aes(x=chosen_avg, y=other_max,label=gene),max.overlaps = 30)+
    theme_bw() + ylab("log10(max expression in other cells)") + xlab("log10(averege expression in celltype)")+
    theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    scale_color_manual(values = c("grey", "red"))
  res = comparison[rownames(comparison)[comparison$signif==T],]
  res = res[order(res$ratio*res$chosen_avg,decreasing=T),c("chosen_avg","other_max","ratio")]
  celltype_markers = rownames(res)
  if (barplots==T) {barplots = plot_roy_marker(expression, celltype_markers, celltypes)}
  return (list(markers=celltype_markers,markers_tbl=res,celltypes=celltypes,plot=p, barplots=barplots))
}

load_segments_human = function(){
#' loads human intestinal sections data from Rachel Zwick (unpublished)
  library(readxl)
  library(tidyr)
  if (!exists("human_sections")){
    path = "X:\\Shalevi\\Rachel_Zwick"
    files = c("Enteroendocrine cells.xlsx","Goblet cells.xlsx","Tuft cells.xlsx","Paneth cells.xlsx","Crypt cells.xlsx","Villus cells.xlsx")
    for(i in 1:length(files)){
      cur_path = paste(path,files[i],sep="/")
      name = sub(" cells.xlsx","",files[i])
      d = read_excel(cur_path)
      d = d[1:6]
      d_long = pivot_longer(d,cols=colnames(d[2:6]),names_to="section",values_to=name)
      if (i == 1){
        human_sections = d_long
      } else {human_sections[,name] = d_long[,name]}
    }
    human_sections$section = as.integer(gsub("mean ","",human_sections$section))
    colnames(human_sections)[colnames(human_sections)=="Gene name"] = "gene"
    human_sections <<- human_sections
  } else {return(human_sections)}
}

load_segments_mouse = function(){
#' loads mouse intestinal sections data from Rachel Zwick (unpublished)
  library(tidyr)
  if (!exists("mouse_sections")){
    path = "X:\\Shalevi\\Rachel_Zwick/mouse_human_comparison"
    files = c("Mouse Enteroendocrine cells","Mouse Goblet cells","Mouse Tuft cells","Mouse Crypt cells","Mouse Villus cells")
    for(i in 1:length(files)){
      cur_path = paste(path,files[i],sep="/")
      name = sub(" cells","",files[i])
      name = sub("Mouse ","",name)
      d = read.table(cur_path,sep="\t",header=T)
      d = d[1:6]
      d_long = pivot_longer(d,cols=colnames(d[2:6]),names_to="section",values_to=name)
      if (i == 1){
        mouse_sections = d_long
      } else {mouse_sections[,name] = d_long[,name]}
    }
    mouse_sections$section = as.integer(gsub("mean.segment.","",mouse_sections$section))
    mouse_sections$gene = toupper(mouse_sections$gene)
    mouse_sections <<- mouse_sections
  } else {return(mouse_sections)}
}
plot_roy_intestinal_sections = function(gene, organism="mouse"){
  #' plots expression along the intestinal segments in human data from Rachel Zwick data
  #' @param gene string, gene to plot
  library(ggplot2)
  library(tidyr)
  gene = toupper(gene)
  if(organism=="mouse"){
    sections = load_segments_mouse()
    colors = c("purple","lightblue","orange","blue","red")    
    }
  if (organism=="human"){
    sections = load_segments_human()
    colors = c("purple","lightblue","orange","darkgreen","blue","red")
    }
  plot = sections[sections["gene"] == gene,]
  plot = pivot_longer(plot,colnames(plot)[3:ncol(plot)], names_to = "celltype", values_to = "expression")
  ggplot(data=plot,aes(x=section,y=expression,color=celltype)) +
    geom_line(size=2) +
    theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    scale_color_manual(values=colors) +
    theme(legend.position="right") + 
    # theme(legend.position=c(0.77,0.9),legend.direction="horizontal",legend.title=element_blank()) +
    ggtitle(paste("intestine axis",organism,sep=" - ")) + xlab("") + ylab("")
}

load_human_visium = function(){
  #' loads Yotams Visium data
  if (!exists("human_visium_zonation")){
    human_visium_zonation = read.csv("X:\\roy\\apicome\\visium_export_from_yotam\\Supplementary Table XX zonation.csv")
    human_visium_zonation = human_visium_zonation[grepl("_median|gene_name",colnames(human_visium_zonation))]
    colnames(human_visium_zonation) = gsub("_median","",colnames(human_visium_zonation))
    colnames(human_visium_zonation)[colnames(human_visium_zonation)=="gene_name"] = "gene"
    human_visium_zonation <<- human_visium_zonation[,!(colnames(human_visium_zonation) %in% c("MusMucosa","SubMucosa"))]
  } else {return(human_visium_zonation)}
}

plot_roy_zonation = function(gene,organism="mouse"){
  #' plots expression levels in intestinal enterocytes along the villus axis.
  #' human is based on Visium data of Yotam, mouse is based on innas data
  library(ggplot2)
  library(tidyr)
  gene = toupper(gene)
  if(organism=="mouse"){
    expression = load_mouse_expression()
    zonation = expression[,grepl("nterocyte",colnames(expression))]
    colnames(zonation) = gsub("enterocyte_","",colnames(zonation))
    plot = zonation[gene,]
    plot = pivot_longer(plot,colnames(plot), names_to="villus_zone", values_to="expression")
    }
  if (organism=="human"){
    zonation = load_human_visium()
    plot = zonation[zonation$gene==gene,]
    plot = pivot_longer(plot,colnames(plot)[2:ncol(plot)], names_to="villus_zone", values_to="expression")
    }
  ggplot(data=plot,aes(x=villus_zone,y=expression)) +
    geom_bar(stat="identity",fill="lightgreen",color="black",size=1) +
    theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    xlab("") + ggtitle(paste("Villus zone - ",gene,sep="")) + ylab("")
}

plot_roy_apicome = function(gene,organism="mouse",exp_thresh=1e-8){
  #' plots the log2(apical/basal) from laser capture samples (median) in tip and base of villi
  library(ggplot2)
  library(ggrepel)
  
  gene = toupper(gene)
  if(organism=="mouse"){apicome = load_mouse_apicome()}
  if (organism=="human"){apicome = load_human_apicome()}
  plot = apicome[apicome$expression_min >= exp_thresh,]
  plot$signif = "no"
  plot$signif[plot$gene == gene] = "yes"
  p = ggplot(plot, aes(tip_a_b, base_a_b, color=signif)) +
    geom_point(aes(alpha=signif)) +
    geom_text_repel(data=plot[plot$signif != "no",],
                    aes(tip_a_b, base_a_b,label=gene),max.overlaps = 30, size=2) +
    ggtitle(paste("Apical/Basal score - ",gene,sep="")) +
    theme_bw() +
    theme(axis.line = element_line(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    xlab("Tip of villus log2(apical/basal)") + ylab("Base villus log2(apical/basal)")+
    scale_color_manual(values = c("black","red")) +
    theme(legend.position="none") + 
    geom_hline(yintercept=0, linetype="dashed",alpha = 0.3)+
    geom_vline(xintercept=0, linetype="dashed",alpha = 0.3)
  return(p)
}

plot_roy_intestines = function(gene,organism="mouse"){
  #' plots apicome, villus axis, intestinal axis and expression
  #' @param gene string, gene to plot
  #' @param organism "human" or "mouse".
  library(ggpubr)
  p1 = plot_roy_inna(gene,organism=organism)
  p2 = plot_roy_intestinal_sections(gene,organism=organism)
  p3 = plot_roy_zonation(gene,organism=organism)
  p4 = plot_roy_apicome(gene,organism=organism)
  p = ggarrange(p1[[1]],ggarrange(p2,p3,p4, nrow = 3),ncol = 2) 
  return(p)
}

plot_roy_apicome_intestines = function(human_gene,mouse_gene=NA) {
  #' plots the expression in mouse and human intestines, and the apicome log2(a/b) expression
  #' @param human_gene gene to plot
  #' @param mouse_gene in case human and mouse genes are different, put the mouse here
  library(gridExtra)
  if (typeof(mouse_gene)=="logical") {mouse_gene=human_gene}
  mouse_gene = toupper(mouse_gene)
  human_gene = toupper(human_gene)
  p1 = plot_roy_inna(mouse_gene,"mouse",celltypes=c("enterocyte_V1","enterocyte_V2","enterocyte_V3","enterocyte_V4","enterocyte_V5","enterocyte_V6"))
  p2 = plot_roy_inna(human_gene,"human",celltypes="Enterocyte_Epi")
  p3 = plot_roy_apicome(mouse_gene,organism="mouse") + ggtitle("mouse")
  p4 = plot_roy_apicome(human_gene,organism="human") + ggtitle("human")
  plot = grid.arrange(p1[[1]], grid.arrange(p3, p4, nrow=2),p2[[1]],ncol=3,widths=c(4, 2, 4))
  return(invisible(plot))
}

############ Tabula Muris and Sapiens ##################################

load_tabula_muris_pancreas = function(){
  library(tidyr)
  if (!exists("pancreas_tm")){
    pancreas_tm = read.csv("X:\\Common\\useful_datasets\\TabuleMuris\\tabulamuris_facs_pancreas.csv",row.names=1)
    pancreas_tm = matnorm(pancreas_tm)
    pancreas_tm$gene = toupper(rownames(pancreas_tm))
    pancreas_tm <<- pivot_longer(pancreas_tm, cols=-gene, names_to="celltype",values_to="expression")
  } else {return(pancreas_tm)}
}

plot_roy_pancreas = function(gene){
  #' based on Tabula Muris
  library(ggplot2)
  load_tabula_muris_pancreas()
  gene = toupper(gene)
  p = ggplot(pancreas_tm[pancreas_tm$gene==gene,],aes(celltype, expression)) +
    geom_bar(stat="identity",fill="cyan",color="black") +
    theme_bw() + xlab("") + ylab("log10(exp)") + ggtitle(gene) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(p)
}

load_tabula_sapiens_pancreas = function(){
  if (!exists("pancreas_ts")){
    pancreas_ts = read.csv("X:\\Common\\Tabula_Sapiens\\TabulaSapiens_pancreas.csv")
    pn = min(pancreas_ts$expression[pancreas_ts$expression>0])
    pancreas_ts <<- pancreas_ts
  } else {return (pancreas_ts)}
}

load_tabula_sapiens_liver = function(){
  if (!exists("liver_ts")){
    liver_ts = read.csv("X:\\Common\\Tabula_Sapiens\\TabulaSapiens_liver.csv")
    pn = min(liver_ts$expression[liver_ts$expression>0])
    liver_ts <<- liver_ts
  } else {return (liver_ts)}
}

load_tabula_sapiens_intestines = function(){
  if (!exists("intestines_ts")){
    intestines_ts = read.csv("X:\\Common\\Tabula_Sapiens\\TabulaSapiens_intestine.csv")
    pn = min(intestines_ts$expression[intestines_ts$expression>0])
    intestines_ts <<- intestines_ts
  } else {return (intestines_ts)}
}

plot_roy_tabula_sapiens = function(gene,organ){
  #' organ = "pancreas","intestine","liver"
  library(ggplot2)
  if (organ=="intestine"){plot = load_tabula_sapiens_intestines()
  } else if (organ=="liver") {plot = load_tabula_sapiens_liver()
  } else if (organ=="pancreas") {plot = load_tabula_sapiens_pancreas()}
  gene = toupper(gene)
  p = ggplot(plot[plot$gene_symbol==gene,],aes(celltype, expression)) +
    geom_bar(stat="identity",fill="cyan",color="black") +
    theme_bw() + xlab("") + ylab("log10(exp)") + ggtitle(gene) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1))
  return(p)
}

########### Human protein atlas ############################

load_translation_human = function(){
  #' loads ENSMBL id and gene symbol
  if (!exists("translation_human")){
    translation_human = read.csv("X:\\roy\\resources\\Ensemble\\human_ensemble_conversion.csv")
    translation_human = translation_human[!duplicated(translation_human$external_gene_name),c("ensembl_gene_id","external_gene_name")]
    translation_human <<- translation_human
  } else {return(translation_human)}
}

human_atlas_image_download = function(gene_name, dir_output,tissue="Small intestine",number_of_images=3){
  #' downloads IHC images of a selected gene from the human protein atlas
  #' notice that each protein can have few antibodies, all will be downloaded
  #' @param gene_name gene to download IHC images
  #' @param dir_output string of path - where to save the images in? a new folder will be created inside with the gene name
  #' @param tissue images of which tissue you want to download
  #' @param number_of_images how many images from each antibody to download? (will randomly be selected)
  library(HPAanalyze)
  library(utils)
  library(plyr)
  gene_name = toupper(gene_name)
  translation_human = load_translation_human()
  gene = translation_human$ensembl_gene_id[translation_human$external_gene_name==gene_name]
  tryCatch(expr={gene_xml=hpaXmlGet(gene)},error = function(err) {
    print(paste0("[ERROR] ",err)) 
    return(NULL)})
  if (!exists("gene_xml")){return(NULL)}
  antibody = hpaXmlAntibody(gene_xml)
  gene_exp = hpaXmlTissueExpr(gene_xml)
  if(length(gene_exp)==0) return()
  all_images = data.frame()
  for (i in 1:length(gene_exp)){
    cur_ab = gene_exp[[i]]
    if (nrow(cur_ab) == 0 && ncol(cur_ab) == 0){next()} # empty
    cur_ab_1 = cur_ab[cur_ab$tissueDescription1 == "Normal tissue, NOS" &
                        cur_ab$tissueDescription2 == tissue ,]
    # (cur_ab$intensity != "Negative" | is.na(cur_ab$intensity))
    cur_ab_1 = cur_ab_1[rowSums(is.na(cur_ab_1)) < ncol(cur_ab_1), ]
    
    if(nrow(cur_ab_1)==0){ # if no healthy tissue, use diseased
      cur_ab_1 = cur_ab[cur_ab$tissueDescription2 == tissue,]}
    # & (cur_ab$intensity != "Negative" | is.na(cur_ab$intensity))
    if(nrow(cur_ab_1)==0){next()} # no images of the tissue
    if (nrow(cur_ab_1) > number_of_images){
      cur_ab_1 = sample_n(cur_ab_1,number_of_images)
    }
    if(i==1) {
      all_images = cur_ab_1
    } else {all_images = rbind.fill(all_images,cur_ab_1)} 
  }
  all_images = all_images[rowSums(is.na(all_images)) < ncol(all_images), ]
  if (nrow(all_images)>0){
    path = paste0(dir_output,"/",gene_name)
    if (!dir.exists(path)){dir.create(path)} 
    for (j in 1:nrow(all_images)) {
      image_name = paste0(path,"/",gene_name,"_",all_images$tissueDescription2[j],"_",j,".jpg")
      download.file(all_images$imageUrl[j],destfile=image_name,mode="wb")
    }
  }
}

load_HPA = function(){
  if (!exists("hpa")){
    hpa = read.csv("X:\\Common\\useful_datasets\\human_protein_atlas_expression.csv")
    hpa <<- hpa
}}

plot_roy_hpa = function(gene){
  library(ggplot2)
  load_HPA()
  gene = toupper(gene)
  plot = hpa[hpa$gene == gene,]
  plot$tissue = factor(plot$tissue, levels = plot$tissue[order(plot$organ)])
  p = ggplot(plot, aes(tissue,nTPM,fill=organ)) + 
        geom_bar(stat="identity") +
        coord_flip() +
        theme_bw() + 
        xlab("") +
        theme(legend.pos="none",axis.line = element_line(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        ggtitle(gene)
  return(p)
}
