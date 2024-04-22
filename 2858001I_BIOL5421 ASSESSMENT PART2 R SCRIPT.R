
####Assessment Part 2####

####Preparing the data from the first part####

#Renaming and ordering the em by group
em=normalized_correctedcounts_slide 
normalized_correctedcounts_slide = normalized_correctedcounts_slide[,row.names(ss)]

#Sorting the ss by group
ss = ss[order(ss$sample_group, decreasing=TRUE),]

#Adding a significance column TRUE OR FALSE
de_GB1vsHC$sig = as.factor(de_GB1vsHC$p.adj < 0.05 & abs(de_GB1vsHC$log2fold) > 1)
de_GB2vsHC$sig = as.factor(de_GB2vsHC$p.adj < 0.05 & abs(de_GB2vsHC$log2fold) > 1)
de_GB1vsGB2$sig = as.factor(de_GB1vsGB2$p.adj < 0.05 & abs(de_GB1vsGB2$log2fold) > 1)


####Loading all libraries####

library(amap)
library(reshape2)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)       
library(eulerr)


####Plotting necessary functions####

#PCA
make_pc1_pc2 = function(colour_groups, e_data)
{
  # scale data
  e_data_scaled = na.omit(data.frame(t(scale(t(e_data)))))
  
  # run PCA
  xx = prcomp(t(e_data_scaled))
  pca_coordinates = data.frame(xx$x)
  
  # get % variation
  vars = apply(xx$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars),4) * 100
  prop_y = round(vars["PC2"] / sum(vars),4) * 100
  x_axis_label = paste("PC1 (" ,prop_x, "%)", sep="")
  y_axis_label = paste("PC2 (" ,prop_y, "%)", sep="")
  
  # plot  
  ggpFNC = ggplot(pca_coordinates, aes(x=PC1, y= PC2, colour = colour_groups)) +
    geom_point(size = 3) +
    labs(title = "PCA", x= x_axis_label, y= y_axis_label) +
    theme_bw()
  
  return(ggpFNC)
}

#Volcano
make_volcano = function(master,master_non_sig,master_sig_down_top5,master_sig_up_top5){
  
  master_sig_up=subset(master,p.adj < 0.05 & log2fold > 1)
  master_sig_down=subset(master,p.adj< 0.05 & log2fold < - 1)
  
  #Creating tables for top 5 significantly up- and top 5 down-regulated genes
  
  master_sig_up_top5 = master_sig_up[1:5,]
  master_sig_down_top5 = master_sig_down[1:5,]
  
  #Creating a table for non-significant genes
  master_non_sig = subset(master, p.adj > 0.05  | (log2fold < 1 & log2fold > -1) )
  
  # Adding new shared column and assigning fixed values of 3 groups
  master_non_sig $direction = "a"
  master_sig_down $direction = "b"
  master_sig_up$direction = "c"
  
  # Rejoining the tables
  master_main = rbind(master_non_sig, master_sig_down, master_sig_up)
  
  volcano_plot = ggplot(master, aes(x=log2fold, y=mlog10p , colour = direction)) +
    # adds the dots
    geom_point(data=master_non_sig,size=0.75,alpha=0.7,fill="black")+
    geom_point(data=master_sig_down,size=1.25,alpha=0.6,fill="deepskyblue")+
    geom_point(data=master_sig_up,size=1.25,alpha=0.6,fill="indianred")+ 
    #labels the top genes
    geom_label_repel (data=master_sig_down_top5, aes(label=row.names(master_sig_down_top5)),color="deepskyblue")+
    geom_label_repel (data=master_sig_up_top5, aes(label=row.names(master_sig_up_top5)),color="indianred")+
    #adds lines
    geom_vline(xintercept=-1, linetype="dashed") +
    geom_vline(xintercept=1, linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    #editing the legend
    scale_color_manual(values=c("black","deepskyblue", "indianred"), 
                       labels = c("No change", "Downregulated", "Upregulated"), name="")+
    guides(color = guide_legend( 
      override.aes=list(shape = 19)))+ 
    #editing the theme and adding labels
    theme(plot.title=element_text(size=18),axis.text.x = element_text(face="bold",size=9, hjust = 1), axis.text.y = element_text(face="bold",size=9),axis.title.x = element_text(size=12),axis.title.y= element_text(size=12))+labs(title = "Volcano Plot", x= "Log2 fold change", y= "-log10p") 
  
  return(volcano_plot)
}


#Heatmap
#Supplying the expression table and  vector of genes
make_heatmap = function(e_data, candidate_genes)
{
  
  # parses
  e_data_scaled_candidates = na.omit(data.frame(t(scale(t(e_data[candidate_genes,])))))
  hm.matrix = as.matrix(e_data_scaled_candidates)
  
  # does the y clustering
  y.dist = Dist(hm.matrix, method="spearman")
  y.cluster = hclust(y.dist, method="average")
  y.dd = as.dendrogram(y.cluster)
  y.dd.reorder = reorder(y.dd,0,FUN="average")
  y.order = order.dendrogram(y.dd.reorder)
  hm.matrix_clustered = hm.matrix[y.order,]
  
  # melt and plot
  hm.matrix_clustered = melt(hm.matrix_clustered)
  
  # plot
  ggp = ggplot(hm.matrix_clustered, aes(x=Var2, y=Var1, fill=value)) + 
    geom_tile() + 
    scale_fill_gradientn(colours = colorRampPalette(c("black","deepskyblue","white"))(100)) + 
    labs(x="", y="") + 
    theme(legend.title = element_blank(), legend.spacing.x = unit(0.25, 'cm'), axis.text.y=element_blank(),axis.text.x = element_blank(), axis.ticks=element_blank())
  scale_fill_gradientn(colours = colorRampPalette(c("black","deepskyblue","white"))(100), limits=c(-2,2)) + 
    
    return(ggp)
}

#Plotting the mega heatmap for all significant genes
sig_any = row.names(subset(master_main, sig.hc_gb1 == TRUE | sig.hc_gb2 == TRUE  | sig==TRUE ))
make_heatmap(em_symbols.s, sig_any)


#Rug
make_rug = function(ss)
{
  #rug colour theme
  rug_colours=c('red','orange','darkolivegreen3')
  
  #rug table for the different sample groups
  hm_groups.melt=as.matrix(as.numeric(as.factor(ss$sample_group)))
  hm_groups.melt=melt(hm_groups.melt)
  
  #making the rug
  rug=ggplot(hm_groups.melt,aes(x=Var1,y=Var2,fill=value))+
    geom_tile()+
    scale_fill_gradientn(colours=rug_colours)+
    
    #trimming off the headers, legends, everything
    theme(plot.margin =unit(c(0,1,1,1), "cm"), axis.line= element_blank(), axis.text.x= element_blank(), axis.title.x= element_blank(), axis.text.y= element_blank(), axis.ticks= element_blank(), axis.title.y= element_blank(), legend.position='none', panel.background= element_blank(), panel.border=element_blank(), panel.grid.major= element_blank(), panel.grid.minor= element_blank(), plot.background= element_blank())
  rug
  
  return(rug)
}

#Boxplot
make_boxplot=function(master, ss, em_scaled){
  
  lowest_3_genes=row.names(master)[1:3]
  lowest_3_genes
  
  top_gene_data = em_scaled[lowest_3_genes,]
  top_gene_data = data.frame(t(top_gene_data))
  top_gene_data$sample_groups = ss$sample_group
  
  # melting our table choosing the column we want to group by
  top_gene_data.m = melt(top_gene_data, id.vars="sample_groups")
  
  # making the plot 
  top_plot = ggplot(top_gene_data.m,aes(x=variable,y=value, fill=sample_groups))+geom_boxplot()+
    theme(plot.title=element_text(size=18),axis.text.x = element_text(face="bold", family="Arial",size=12,angle = 45, hjust = 1), axis.text.y = element_text(size=9,face="bold"),axis.title.x = element_text(size=17),axis.title.y= element_text(size=17))+
    labs(title="Expression Levels of 3 Most Significant Genes",y="Expression",x="Genes")+ylim(c(-1,2.5))+scale_fill_brewer(palette="Dark2")
  return(top_plot)
}


#Pathway
make_pathway = function(candidate_genes, id_type, org.db)
{
  # converts from ensembl Symbols to Entrez - needed for cluster profiler
  candidate_genes = bitr(candidate_genes, fromType = id_type, toType = "ENTREZID", OrgDb = org.db)
  
  # gets the enrichment
  enriched_gene_sets = enrichGO(gene = candidate_genes$ENTREZID,OrgDb = org.db, readable = T,ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)
  
  # Plot the enriched ontologies in various different ways:
  ggp = barplot(enriched_gene_sets, showCategory=10)
  
  # create a list to return
  results = list("top10_barplot" = ggp, "enriched_gene_sets" = enriched_gene_sets)
  
  return(results)
}

make_pathway(sig_genes,"SYMBOL",org.Hs.eg.db)


####Making the DE workflow function####
run_de_workflow=function(em, de, ss)
{
  #create master
  master = merge(em, de, by.x = 0, by.y = 0)
  
  # Change the names
  row.names(master) = master[,1]
  
  # remove nas
  master = na.omit(master)
  master= master[,-1]
  
  # sort by p
  master = master[order(master[,"p"]),]
  
  # create new columns
  master$mean = rowMeans(master[,1:84])
  master$sig = as.factor(master$p.adj < 0.05 & abs(master$log2fold) > 1)
  master$mlog10p = -log10(master$p.adj)
  
  # create sig_genes
  master_sig = subset(master, sig==TRUE)
  
  sig_genes = row.names(master_sig)
  
  # create em symbols
  em_symbols = master[,row.names(ss)]
  
  # create em symbols scaled
  em_scaled = na.omit(data.frame(t(scale(t(em_symbols)))))
  
  # create em sig
  em_symbols_sig = em_symbols[sig_genes,]
  em_scaled_sig = em_scaled[sig_genes,]
  
  # make the plots
  ggp.pca = make_pc1_pc2(ss$sample_group, em)
  ggp.hm = make_heatmap(em_symbols, sig_genes)
  path_results = make_pathway(sig_genes, "SYMBOL", org.Hs.eg.db)
  ggp_volcano = make_volcano(master,master_non_sig,master_sig_down_top5,master_sig_up_top5)
  ggp.boxplot = make_boxplot(master,ss,em_scaled)
  rug=make_rug(ss)
  results = list("master" = master_main, "master_sig" = master_sig , "sig_genes" =sig_genes, " em_symbols " = em_symbols, "pca_plot" = ggp.pca, "heatmap_plot"= ggp.hm,"rug"=rug, "pathway_results" = path_results, "ggp_volcano"= ggp_volcano, "boxplot"=ggp.boxplot)
  
  return(results)
}


####Running the workflow for each comparison####

#Getting the results
results_gb1vgb2=run_de_workflow(em, de_GB1vsGB2, ss)
results_gb1vhc=run_de_workflow(em, de_GB1vsHC, ss)
results_gb2vhc=run_de_workflow(em, de_GB2vsHC, ss)

#Plotting results on Volcano plots
results_gb1vgb2$ggp_volcano
results_gb1vhc$ggp_volcano
results_gb2vhc$ggp_volcano

#Plotting results Boxplots
results_gb1vgb2$boxplot
results_gb1vhc$boxplot
results_gb2vhc$boxplot


####Exploring Signatures####

#Merging the 3 DE tables
master_main = merge(de_GB1vsHC, de_GB2vsHC, by.x=0, by.y=0, suffixes=c(".hc_gb1",".hc_gb2"))
master_main = merge(master_main, de_GB1vsGB2, by.x=1, by.y=0)
master_main=merge(em,master_main, by.x=0, by.y=1)
master_main = transform(master_main, row.names = Row.names, Row.names = NULL)

#Getting em symbols and scaled
em_symbols = master_main[,row.names(ss)]
em_symbols.s = na.omit(data.frame(t(scale(t(em_symbols)))))

#Getting the significant genes for each comparison
sig_HCvsGB1 = row.names(subset(master_main, sig.hc_gb1 == TRUE ))
sig_HCvsGB2 = row.names(subset(master_main, sig.hc_gb2 == TRUE ))
sig_GB1vsGB2 = row.names(subset(master_main, sig == TRUE )) 

#Setting the levels for the ssample sheet
ss$sample_group = factor(ss$sample_group, levels=c("GB_1", "GB_2", "HC"))

#Venn diagram 
venn_data = list("HCvsGB1" = sig_HCvsGB1, "HCvsGB2" = sig_HCvsGB2, "GB1vsGB2" = sig_GB1vsGB2)
#plot Venn
plot(euler(venn_data, shape = "circle"), quantities = TRUE)


####Making heatmaps for signatures#### 

#Signature1
signature_1 = row.names(subset(master_main,(sig.hc_gb1 == TRUE & log2fold.hc_gb1 < -1 ) & (sig.hc_gb2 == TRUE & log2fold.hc_gb2 < -1) & (sig == FALSE )))
heatmap_sig1=make_heatmap(em_symbols.s,signature_1)
heatmap_sig1

#Signature2
signature_2 = row.names(subset(master_main,(sig.hc_gb1 == TRUE & log2fold.hc_gb1 >1 ) & (sig.hc_gb2 == FALSE) & (sig == TRUE & log2fold>1 )))
heatmap_sig2=make_heatmap(em_symbols.s,signature_2)
heatmap_sig2

#Signature3
signature_3= row.names(subset(master_main,(sig.hc_gb1 == FALSE ) & (sig.hc_gb2 == TRUE & log2fold.hc_gb2 >1 ) & (sig == FALSE )))
heatmap_sig3=make_heatmap(em_symbols.s,signature_3)
heatmap_sig3

#Signature4
signature_4 = row.names(subset(master_main,(sig.hc_gb1 == TRUE & log2fold.hc_gb1 >1 ) & (sig.hc_gb2 == TRUE & log2fold.hc_gb2 >1) & (sig == FALSE)))
heatmap_sig4=make_heatmap(em_symbols.s,signature_5)
heatmap_sig4


####Plotting signature metagenes####

#Boxplots

#Signature 1
#getting the em for a cluster 
signature_1_em = em_symbols.s[signature_1,]

#getting the metagene for each cluster
signature_1_metagene = data.frame(colMeans(signature_1_em))
names(signature_1_metagene) = "meta_expression"
signature_1_metagene$group = ss$sample_group

#plot
ggplot(signature_1_metagene, aes(x=group, y=meta_expression, fill=group)) + geom_boxplot()


#Signature2
signature_2_em = em_symbols.s[signature_2,]

#getting the metagene for each cluster
signature_2_metagene = data.frame(colMeans(signature_2_em))
names(signature_2_metagene) = "meta_expression"
signature_2_metagene$group = ss$sample_group

#plot
ggplot(signature_2_metagene, aes(x=group, y=meta_expression, fill=group)) + geom_boxplot()


#Signature3
signature_3_em = em_symbols.s[signature_3,]

#getting the metagene for each cluster
signature_3_metagene = data.frame(colMeans(signature_3_em))
names(signature_3_metagene) = "meta_expression"
signature_3_metagene$group = ss$sample_group

#plot
ggplot(signature_3_metagene, aes(x=group, y=meta_expression, fill=group)) + geom_boxplot()


#Signature4
signature_4_em = em_symbols.s[signature_4,]

#getting the metagene for each cluster
signature_4_metagene = data.frame(colMeans(signature_4_em))
names(signature_4_metagene) = "meta_expression"
signature_4_metagene$group = ss$sample_group

#plot
ggplot(signature_4_metagene, aes(x=group, y=meta_expression, fill=group)) + geom_boxplot()


#Pathway barplots
#Signature1
signature_1_path = bitr(signature_1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# gets the enrichment
enriched_gene_sets = enrichGO(gene = (signature_1_path)$ENTREZID,OrgDb = org.Hs.eg.db, readable = T,ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)

#Plot the enriched ontologies in various different ways:
ggp_path1 = barplot(enriched_gene_sets, showCategory=10)
ggp_path1

#Signature2
signature_2_path = bitr(signature_2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# gets the enrichment
enriched_gene_sets = enrichGO(gene = (signature_2_path)$ENTREZID,OrgDb = org.Hs.eg.db, readable = T,ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)

#Plot the enriched ontologies in various different ways:
ggp_path2 = barplot(enriched_gene_sets, showCategory=10)
ggp_path2


#Signature3
signature_3_path = bitr(signature_3, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# gets the enrichment
enriched_gene_sets = enrichGO(gene = (signature_3_path)$ENTREZID,OrgDb = org.Hs.eg.db, readable = T,ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)

#Plot the enriched ontologies in various different ways:
ggp_path3 = barplot(enriched_gene_sets, showCategory=10)
ggp_path3


#Signature4
signature_4_path = bitr(signature_4, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# gets the enrichment
enriched_gene_sets = enrichGO(gene = (signature_4_path)$ENTREZID,OrgDb = org.Hs.eg.db, readable = T,ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.10)

#Plot the enriched ontologies in various different ways:
ggp_path4 = barplot(enriched_gene_sets, showCategory=10)
ggp_path4




















