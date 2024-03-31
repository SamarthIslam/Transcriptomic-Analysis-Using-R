                          ####Part 1:Loading and Parsing Tables####

#Loading data using column headers and row names

em=read.table("C:/Users/samar/OneDrive/Desktop/Septic_Arthritis/em.csv",header=TRUE,row.names=1,sep="\t")
de_controls_v_arthritis=read.table("C:/Users/samar/OneDrive/Desktop/Septic_Arthritis/DE_hc_vs_sa.csv",header=TRUE,row.names=1,sep="\t")
annotations=read.table("C:/Users/samar/OneDrive/Desktop/Septic_Arthritis/annotations.csv",header=TRUE,row.names=1, sep= "\t")
sample_sheet=read.table("C:/Users/samar/OneDrive/Desktop/Septic_Arthritis/sample_sheet.csv",header=TRUE,row.names=1, sep= "\t") 

#Providing more sensible names to 'annotations'
names(annotations)=c("Gene_Symbol","Chromosome","Stop Position", "Start Position")

#Joining tables to give a Master table
master_temp=merge(em,annotations,by.x=0,by.y = 0)
master=merge(master_temp,de_controls_v_arthritis,by.x=1,by.y = 0)

#Editing the column names of the Master table
row.names(master)=master[,30]
names(master)[1]="Gene_ID"

#Making an expression table with only gene symbols 
em_symbols=master[,-c(1,30:36)]


#Removing NA's from Master table
master=na.omit(master)

#Sorting Master by increasing p-value
sorted_order=order(master[,36],decreasing=FALSE)
master=master[sorted_order,]
print(master)

#Making a new column for mean expression
value=rowMeans(master[,2:29])
master$Expression_means=value

#Adding a column for '-10logp' value
master$mlog10p=-log10(master$p)

#Adding a column for significant genes
master$sig=as.factor(master$p.adj< 0.05 & abs(master$log2fold) > 1.0)

#Making a scaled expression matrix for heatmaps and PCA
em_scaled=data.frame(t(scale(t(em_symbols))))
#Removing NA's
em_scaled=na.omit(em_scaled)

#Making a list of significant genes 

master_sig=subset(master,p.adj < 0.05 & abs(log2fold) >1)

#Creating an independent list of significant genes 
sig_genes=master_sig[,"Gene_Symbol"]
sig_genes

#Making expression tables of only significant genes

#adding significant genes to master sig table
sig_genes=rownames(master_sig)
em_symbols_sig=em_symbols[sig_genes,]

#Making table for scaled expression values
em_scaled_sig=em_scaled[sig_genes,]


 ####Part 2: Creating an auto-colored Volcano Plot to observe significant genes####

#Loading necessary packages
library(ggplot2)
library(ggrepel)

#Creating tables for significantly up- and down-regulated genes

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
master = rbind(master_non_sig, master_sig_down, master_sig_up)

#Making the plot
volcano_plot = ggplot(master, aes(x=log2fold, y=mlog10p , colour = direction)) +
  # adds the dots
  geom_point(data=master_non_sig,size=0.75,size=0.75,alpha=0.7,fill="black")+
  geom_point(data=master_sig_down,size=1.25,alpha=0.6,fill="deepskyblue")+
  geom_point(data=master_sig_up,size=1.25,alpha=0.6,fill="indianred")+ 
  #labels the top genes
  geom_label_repel (data=master_sig_down_top5, aes(label=Gene_Symbol),color="deepskyblue")+
  geom_label_repel (data=master_sig_up_top5, aes(label=Gene_Symbol),color="indianred")+
  #adds lines
  geom_vline(xintercept=-1, linetype="dashed") +
  geom_vline(xintercept=1, linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  #adds a y-axis limit
  ylim(c(0, 25))+
  #editing the legend
  scale_color_manual(values=c("black","deepskyblue", "indianred"), 
                     labels = c("No change", "Downregulated", "Upregulated"), name="")+
  guides(color = guide_legend( 
    override.aes=list(shape = 19)))+ 
  #editing the theme and adding labels
  theme(plot.title=element_text(size=18),axis.text.x = element_text(face="bold", family="Arial",size=9, hjust = 1), axis.text.y = element_text(face="bold",size=9),axis.title.x = element_text(size=12),axis.title.y= element_text(size=12))+labs(title = "Volcano Plot", x= "Log2 fold change", y= "-log10p") 

volcano_plot

#Saving the plot to disk
png("C:/Users/samar/OneDrive/Desktop/VolPlot.png",width=500,height=450)
print(volcano_plot)
dev.off()


                            ####Part 3:Creating a PCA Plot###

#Storing the sample sheet table in a new variable
ss = sample_sheet

#Converting em_scaled into a matrix of numerics
numeric_matrix=as.matrix(sapply(em_scaled,as.numeric))

#Transpose the matrix
pca=prcomp(t(numeric_matrix))

#Extracting the component data
pca_coordinates=data.frame(pca$x)

#Getting the % variance for  axis labels
vars=apply(pca$x,2,var)
prop_x = round(vars["PC1"] / sum(vars),4) * 100
prop_y = round(vars["PC2"] / sum(vars),4) * 100

#Storing the % variance data in 2 new vaiables
x_axis_label = paste("PC1 ", " (",prop_x, "%)",sep="")
y_axis_label = paste("PC2 ", " (",prop_y, "%)",sep="")


#Making the plot and coloring dots by sample groups 

pca_plot=ggplot(pca_coordinates,aes(x=PC1,y=PC2,color=ss$sample_group))+
  #adds the dots
   geom_point(size=3.5,shape=16)+
  #applying the colors and labels
  scale_color_manual(values=c("chartreuse4","indianred"),
                     labels=c("Healthy Controls", "Septic Arthritis",name=""))+
 
  #editing the theme and adding the labels
  theme(plot.title=element_text(size=18),axis.text.x = element_text(face="bold", family="Arial",size=9, hjust = 1), axis.text.y = element_text(face="bold",size=9),axis.title.x = element_text(size=12),axis.title.y= element_text(size=12))+
labs(title="PCA of Sample Groups", x=x_axis_label, y=y_axis_label)+guides(color = guide_legend(title = "Sample Groups"))
  
pca_plot

#Saving the plot to disk
png("C:/Users/samar/OneDrive/Desktop/PCA.png",width=500,height=350)
print(pca_plot)
dev.off()



     ####Part 4:Creating a Clustered Heatmap of Differentially Expressed Genes####
#Loading library
library(amap)

#Omiting any NA's form the scaled data
em_scaled_sig=na.omit(em_scaled_sig)

#Making a subset of 30 genes from the scaled data
cluster_30genes=em_scaled_sig [c(1:30),]

hm.matrix2 = as.matrix(cluster_30genes)

#getting distances
y.dist=Dist(hm.matrix, method="spearman")
#clustering using the distances (build the dendogram)
y.cluster=hclust(y.dist, method="average")
# extracting the dendrogram 
y.dd = as.dendrogram(y.cluster)
# untangling the dendrogram
y.dd.reorder = reorder(y.dd,0,FUN="average")
# getting the untangled gene order from the dendogram
y.order = order.dendrogram(y.dd.reorder)
# reordering the original matrix in the new order
hm.matrix_clustered = hm.matrix[y.order,]

#making the colour palette
colours = c("deepskyblue","orange","indianred")
palette = colorRampPalette(colours)(10)

#melting and plotting
hm.matrix_clustered = melt(hm.matrix_clustered)

heatmap_30genes = ggplot(hm.matrix_clustered, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile() + 
  scale_fill_gradientn(colours = palette)+
  labs(title="Clustered Heatmap",x="",y="Differentially Expressed Genes") + 
  theme(plot.title=element_text(size=18),axis.text.y=element_blank(),axis.text.x =element_blank(), axis.ticks = element_blank(),axis.title.y= element_text(size=14))
heatmap_30genes

#mention that value is expression z score

#Saving the plot to disk
png("C:/Users/samar/OneDrive/Desktop/Heatmap.png",width=850,height=1000)
print(heatmap_30genes)
dev.off()

# Making the rug for discrete  variable
groups_data = as.matrix(as.numeric(as.factor(ss$sample_group)))
groups_data = melt(groups_data)

# Making the heatmap
rug_colours = c("chartreuse4","orange","indianred")
rug = ggplot(groups_data, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + scale_fill_gradientn(colours = rug_colours)
rug

#Removing everything except plot
only_rug=ggplot(groups_data, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colours = rug_colours)+
  theme(plot.margin=unit(c(0,1,1,1), "cm"), axis.line=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.y=element_blank(),legend.position="none",panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.background=element_blank())
only_rug

#Saving the plot to disk
png("C:/Users/samar/OneDrive/Desktop/Rug.png",width=850,height=100)
print(only_rug)
dev.off()

                            ####Part 5: Creating boxplots####

#Making single boxplot for 10 most significant genes
lowest_10_genes=row.names(master_sig)[1:10]
lowest_10_genes

library(reshape2)
top_gene_data = em_scaled[lowest_10_genes,]
top_gene_data = data.frame(t(top_gene_data))
top_gene_data$sample_groups = ss$sample_group

# melting our table choosing the column we want to group by
top_gene_data.m = melt(top_gene_data, id.vars="sample_groups")

# making the plot 
top_plot = ggplot(top_gene_data.m,aes(x=variable,y=value, fill=sample_groups))+geom_boxplot()+
  theme(plot.title=element_text(size=18),axis.text.x = element_text(face="bold", family="Arial",size=12,angle = 45, hjust = 1), axis.text.y = element_text(size=9,face="bold"),axis.title.x = element_text(size=17),axis.title.y= element_text(size=17))+
  labs(title="Expression Levels of 10 Most Significant Genes",y="Expression",x="Genes")+ylim(c(-1,2.5))+scale_fill_brewer(palette="Dark2")
top_plot 
png("C:/Users/samar/OneDrive/Desktop/Boxplot.png",width=900,height=650)
print(top_plot)
dev.off()

#Making single boxplot for 5 most significantly upregulated genes
five_up_genes=row.names(master_sig_up_top5)[1:5]
five_up_genes

five_up_gene_data = em_scaled[five_up_genes,]
five_up_gene_data = data.frame(t(five_up_gene_data))
five_up_gene_data$sample_groups = ss$sample_group


# melting table choosing the column we want to group by
library(reshape2)
five_up_gene_data.m = melt(five_up_gene_data, id.vars="sample_groups")

# making the plot 
five_up_plot = ggplot(five_up_gene_data.m,aes(x=variable,y=value, fill=sample_groups))+geom_boxplot()+
#adding the theme and editing the labels
  theme(plot.title=element_text(size=18),axis.text.x = element_text(face="bold", family="Arial",size=12,angle = 45, hjust = 1), axis.text.y = element_text(size=9,face="bold"),axis.title.x = element_text(size=17),axis.title.y= element_text(size=17))+
  labs(title="Top 5 Upregulated Genes",y="Expression",x="Genes")+ylim(c(-1,2.5))+scale_fill_brewer(palette="Dark2")

five_up_plot

#Saving the plot to disk
png("C:/Users/samar/OneDrive/Desktop/BoxplotUp5.png",width=900,height=650)
print(five_up_plot)
dev.off()

#Making single boxplot for 5 most significantly downregulated genes
five_down_genes=row.names(master_sig_down_top5)[1:5]
five_down_genes

five_down_gene_data = em_scaled[five_down_genes,]
five_down_gene_data = data.frame(t(five_down_gene_data))
five_down_gene_data$sample_groups = ss$sample_group

# melting table choosing the column we want to group by
five_down_gene_data.m = melt(five_down_gene_data, id.vars="sample_groups")

# making the plot 
five_down_plot = ggplot(five_down_gene_data.m,aes(x=variable,y=value, fill=sample_groups))+geom_boxplot()+
#adding the theme and editing the labels
  theme(plot.title=element_text(size=18),axis.text.x = element_text(face="bold", family="Arial",size=12,angle = 45, hjust = 1), axis.text.y = element_text(size=9,face="bold"),axis.title.x = element_text(size=17),axis.title.y= element_text(size=17))+
  labs(title="Top 5 Downregulated Genes",y="Expression",x="Genes")+ylim(c(-1,2.5))+scale_fill_brewer(palette="Dark2")

five_down_plot

#Saving the plot to disk
png("C:/Users/samar/OneDrive/Desktop/BoxplotDwn5.png",width=900,height=650)
print(five_down_plot)
dev.off()

#Making a box plot of 3 genes of interest
three_genes=row.names(master_sig_up_top5)[2:4]
three_genes
three_gene_data = em_scaled[three_genes,]
three_gene_data = data.frame(t(three_gene_data))
three_gene_data$sample_groups = ss$sample_group

# melting our table choosing the column we want to group by
three_gene_data.m = melt(three_gene_data, id.vars="sample_groups")

three_genes_faceted=ggplot(three_gene_data.m,aes(x=sample_groups,y=value, facet=variable,fill=value)) + theme(plot.title=element_text(size=18),axis.text.x = element_text(face="bold", family="Arial",size=9,angle = 45, hjust = 1), axis.text.y = element_text(size=9),axis.title.x = element_text(size=15),axis.title.y= element_text(size=15))+
  geom_boxplot(aes(fill=sample_groups), show.legend=FALSE)+
  facet_wrap(~variable, ncol=3)+labs(y="Expression Level",x="Sample Groups",title="Faceted Boxplot of 3 Genes of Interest", size=12)+scale_fill_brewer(palette="Dark2")+ylim(c(-1,2.5))
three_genes_faceted

png("C:/Users/samar/OneDrive/Desktop/Boxplot3gn.png",width=900,height=650)
print(three_genes_faceted)
dev.off()


                        ####Part 6: Pathway Analyses####

#Loading the human gene database
library(org.Hs.eg.db)

#Creating variables for the 5 upregulated genes
five_sig_up_genes=row.names(master_sig_up_top5)
five_sig_genes_up_entrez = bitr(five_sig_up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ora_results_top5 = enrichGO(gene = five_sig_genes_up_entrez$ENTREZID, OrgDb = org.Hs.eg.db, readable = T,  ont = "BP",  pvalueCutoff = 0.05,  qvalueCutoff = 0.10)

#Making the cnet plot
cnet_up5= cnetplot(ora_results_top5, categorySize="pvalue")
cnet_up5
png("C:/Users/samar/OneDrive/Desktop/CnetUP.png",width=900,height=750)
print(cnet_up5)
dev.off()

png("C:/Users/samar/OneDrive/Desktop/CnetDOWN.png",width=900,height=750)
print(cnet_down5)
dev.off()

#Running STRING to observe interaction network

#Loading STRING 
library(STRINGdb)

#naming variables to store 10 genes with lowest p-values
lowest_10_genes_table = data.frame(lowest_10_genes)
names(lowest_10_genes_table) = "gene"

#loading database
string_db = STRINGdb$new( version="11.5", species=9606, score_threshold=200, network_type="full", input_directory="")

# mapping  genes to database
string_mapped = string_db$map(lowest_10_genes_table, "gene", removeUnmappedRows = TRUE )

#plotting
string_db$plot_network(string_mapped)
string_mapped

#saving the plot to disk
png("C:/Users/samar/OneDrive/Desktop/STIRNG.png",width=900,height=850)
print(string_mapped)
dev.off()
