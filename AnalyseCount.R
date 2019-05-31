require(tidyverse)
require(edgeR)
require(ggrepel)
mytheme = theme_bw()
data  = read.csv("Data/Raw_Counts_RNA-Seq_CetreSossah.txt",sep=",",header=T,row.names = 1)
data %>% dim
glimpse(data)
# transform as a matrix
mdata = as.matrix(data)
mdata %>%  head
# make some graph with raw counts
data %>% gather(sample,count,1:36) %>% ggplot + geom_histogram(aes(x=log(count))) + facet_wrap(~sample) + mytheme
ggsave("Figures/rawcount_histo.pdf")


# Convert to count per million, plot as function to raw count as a way to find cpm value correponding to 10 reads c ounts

mdatacpm = cpm(mdata)
plot(mdatacpm[,4],mdata[,4],xlim=c(0,3),ylim=c(0,50))

# 0.25 cpm is an okay value to filter. We will keep gene with at leats  0.25 cpm in three samples
abovecpm = mdatacpm > 0.75 
table(rowSums(abovecpm))
keep = rowSums(abovecpm) >= 3 

# number of retained genes :
summary(keep)

# 9305 / 17342 , I'm okay with that

filtmdata = mdata[keep,]

DG = DGEList(counts = filtmdata)

ggplot(data.frame(name = colnames(DG),libsize = DG$samples$lib.size,type = infos$type,time=infos$time,sample=infos$sample) %>% arrange(.,sample,time) ) + geom_bar(aes(x=name,y=libsize,fill=sample),stat="identity") + facet_wrap(~ type,scale="free_x")  + scale_fill_brewer(palette ="Dark2") + mytheme + theme(axis.text.x = element_text(angle=45,hjust =1 ))
ggsave("Figures/Unnormalized_lib_size.pdf")

logcount = cpm(DG$counts,log=T)

infos$name=as.factor(infos$name)


datalogcpm = data.frame(logcount) %>% gather(name,count,1:36) %>% left_join(infos,by = "name")

ggplot(datalogcpm %>% arrange(.,sample,time)) + geom_violin(aes(x=name,y=count,fill=sample))  + facet_wrap(~ type,scale="free_x")  + scale_fill_brewer(palette ="Dark2") + mytheme + theme(axis.text.x = element_text(angle=45,hjust =1 ))
ggsave("Figures/LogCPM_violin_count.pdf")
ggsave("Figures/LogCPM_violin_count.png")

ggplot(datalogcpm %>% arrange(.,sample,time)) + geom_jitter(aes(x=name,y=count,color=sample))  + facet_wrap(~ type,scale="free_x")  + scale_color_brewer(palette ="Dark2") + mytheme + theme(axis.text.x = element_text(angle=45,hjust =1 ))
ggsave("Figures/LogCPM_jitter.pdf")


ggplot(data.frame(name = colnames(DG),libsize = DG$samples$lib.size,type = infos$type,time=infos$time,sample=infos$sample) %>% arrange(.,sample,time) ) + geom_bar(aes(x=name,y=libsize,fill=sample),stat="identity") + facet_wrap(~ type,scale="free_x")  + scale_fill_brewer(palette ="Dark2") + mytheme + theme(axis.text.x = element_text(angle=45,hjust =1 ))



# Apply tmm normalization for lib size
DG = calcNormFactors(DG)



# MDS with all data

mdata = plotMDS(DG)
dfmdf=data.frame(x=mdata$x,y=mdata$y)
ggplot(dfmdf %>% mutate(name = rownames(dfmdf)) %>% left_join(infos,by="name" )) + geom_point(aes(x=x,y=y,color=type,shape=as.factor(time)),size=3) + geom_text_repel(aes(x=x,y=y,color = type,label= name)) + mytheme + ggtitle("MDS plot")
ggsave("Figures/MDS_All_DATA.pdf")


# MDS with 50  highly expressed genes
mdata50 = plotMDS(DG,top =50)
dfmdf=data.frame(x=mdata50$x,y=mdata50$y)
ggplot(dfmdf %>% mutate(name = rownames(dfmdf)) %>% left_join(infos,by="name" )) + geom_point(aes(x=x,y=y,color=type,shape=as.factor(time)),size=3) + geom_text_repel(aes(x=x,y=y,color = type,label= name)) + mytheme
ggsave("Figures/MDS_All_DATA_50_high_expressed_genes.pdf")
ggsave("Figures/MDS_All_DATA_50_high_expressed_genes.png")

# # Parsing GO annotation of genes for further analysis
# GO = read.csv("Data/GO_annot.txt",sep="|")
# GO$geneID = as.character(GO$NCBI.gene.ID) 

# Parse gene function description from vector base + GOslim
Desc = read.csv("Data/Gene_description.txt",sep="\t",header=T)
Long_description = Desc %>% group_by(NCBI.gene.ID,Gene.description) %>% summarize(GOslims = toString(GOSlim.GOA.Description)) %>% ungroup
colnames(Long_description)[1] = "geneID"
Long_description$geneID =  as.character(Long_description$geneID)




# Analysis with only Dengue and some mock data
# MDS with dengue and mock data

#dengue_data = data %>% select(contains("Dengue"),contains("MOCKB")) %>% select(-contains("MOCKB24hc")) %>% select(-contains("Dengue6ja"))
dengue_data = data %>% select(contains("Dengue"),contains("MOCKB")) #  %>% select(-contains("MOCKB24hc")) %>% select(-contains("Dengue6ja"))
# transform as a matrix
mddata = as.matrix(dengue_data)
mddatacpm = cpm(mddata)
abovedcpm = mddatacpm > 0.25 
table(rowSums(abovedcpm))
keep = rowSums(abovedcpm) >= 3 
# number of retained genes :
filtmddata = mddata[keep,]
DG = DGEList(counts = filtmddata)
mddata = plotMDS(DG,top = 20)
ddfmdf=data.frame(x=mddata$x,y=mddata$y)
ggplot(ddfmdf %>% mutate(name = rownames(ddfmdf)) %>% left_join(infos,by="name" )) + geom_point(aes(x=x,y=y,color=type,shape=as.factor(time)),size=3) + geom_text_repel(aes(x=x,y=y,color = type,label= name)) + mytheme
ggsave("Figures/MDS_Dengue_DATA.png")
# differentially expresssed genes:
# First remove crappy mock sample:
dengue_data = data %>% select(contains("Dengue"),contains("MOCKB"))  %>% select(-contains("MOCKB24hc")) #%>% select(-contains("Dengue6ja"))
mdengue_data = mdatacpm %>% as.data.frame %>% select(contains("Dengue"),contains("MOCKB"))  %>% select(-contains("MOCKB24hc")) #%>% select(-contains("Dengue6ja"))

# transform as a matrix
mddata = as.matrix(dengue_data)
mddatacpm = cpm(mddata)
abovedcpm = mddatacpm > 0.25 
table(rowSums(abovedcpm))
keep = rowSums(abovedcpm) >= 3 
# number of retained genes :
filtmddata = mddata[keep,]
DG = DGEList(counts = filtmddata)
DG = calcNormFactors(DG)
DG$samples$group = factor(c(1,1,1,2,2,2,3,3,4,4,4))
design = model.matrix(~ DG$samples$group)
DG = estimateDisp(DG, design)
fit = glmQLFit(DG,design)

# Test early response
qlf_dengue_24_vs_mockB24 =  glmTreat(fit, coef=3,lfc = 1)
DE_dengue_24_vs_mockB24  =  topTags(qlf_dengue_24_vs_mockB24,n=1000)
DE_dengue_24_vs_mockB24$table$geneID=rownames(DE_dengue_24_vs_mockB24$table)
df_dengue_24_vs_mockB24  = DE_dengue_24_vs_mockB24 $table %>% filter(FDR<=0.05)
df_dengue_24_vs_mockB24  = df_dengue_24_vs_mockB24   %>% left_join(Long_description,by="geneID")
# To plot the gene expression profile:
mdengue_data = mdatacpm %>% as.data.frame %>% select(contains("Dengue"),contains("MOCKB24"))#  %>% select(-contains("MOCKB24hc")) #%>% select(-contains("Dengue6ja"))
mdengue_data = mdatacpm %>% as.data.frame %>% select(contains("Dengue"),contains("MOCK"))  %>% select(-contains("MOCKB24hc")) %>% select(-contains("Dengue6ja"))
dengue_DE_vals = mdengue_data[df_dengue_24_vs_mockB24$geneID,] %>% t %>% as.data.frame 
dengue_DE_vals$sample = rownames(dengue_DE_vals)
dengue_DE_vals = dengue_DE_vals %>% gather(geneID,count,1:(length(colnames(dengue_DE_vals))-1))
dengue_DE_vals = dengue_DE_vals %>% full_join(df_dengue_24_vs_mockB24, by = "geneID")
dengue_DE_vals = dengue_DE_vals %>% mutate(FCsign = ifelse(logFC < 0,"up","down")) 
ggplot(dengue_DE_vals) + geom_line(aes(x=sample,y=count,group = geneID, color = FCsign)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the dengue infected samples") + scale_color_discrete(name="",labels = c("Down-regulated","Up-regulated"))
ggsave("Figures/dengue_24_expression_profile.png")


# Test late response
qlf_dengue_6j_vs_mockB6j =  glmTreat(fit, contrast=c(0,-1,0,1),lfc = 1)
DE_dengue_6j_vs_mockB6j  =  topTags(qlf_dengue_6j_vs_mockB6j,n=1000)
DE_dengue_6j_vs_mockB6j$table$geneID=rownames(DE_dengue_6j_vs_mockB6j$table)
df_dengue_6j_vs_mockB6j  = DE_dengue_6j_vs_mockB6j $table %>% filter(FDR<=0.05)
df_dengue_6j_vs_mockB6j  = df_dengue_6j_vs_mockB6j   %>% left_join(Long_description,by="geneID")
# To plot the gene expression profile:
mdengue_data = mdatacpm %>% as.data.frame %>% select(contains("Dengue"),contains("MOCKB24"))#  %>% select(-contains("MOCKB24hc")) #%>% select(-contains("Dengue6ja"))
mdengue_data = mdatacpm %>% as.data.frame %>% select(contains("Dengue"),contains("MOCK"))  %>% select(-contains("MOCKB24hc")) %>% select(-contains("Dengue6ja"))
dengue_DE_vals = mdengue_data[df_dengue_6j_vs_mockB6j$geneID,] %>% t %>% as.data.frame 
dengue_DE_vals$sample = rownames(dengue_DE_vals)
dengue_DE_vals = dengue_DE_vals %>% gather(geneID,count,1:(length(colnames(dengue_DE_vals))-1))
dengue_DE_vals = dengue_DE_vals %>% full_join(df_dengue_6j_vs_mockB6j, by = "geneID")
dengue_DE_vals = dengue_DE_vals %>% mutate(FCsign = ifelse(logFC < 0,"up","down")) 
ggplot(dengue_DE_vals) + geom_line(aes(x=sample,y=count,group = geneID, color = FCsign)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the dengue infected samples") + scale_color_discrete(name="",labels = c("Down-regulated","Up-regulated"))
ggsave("Figures/dengue_6j_expression_profile.png")

# TEST WITH RVF_data

#RVF_data = data %>% select(contains("RVF"),contains("MOCKB")) %>% select(-contains("MOCKB24hc")) %>% select(-contains("RVF6ja"))
RVF_data = data %>% select(contains("RVF"),contains("MOCKB"))   %>% select(-contains("MOCKB24hc")) #%>% select(-contains("RVF6ja"))
# transform as a matrix
mddata = as.matrix(RVF_data)
mddatacpm = cpm(mddata)
abovedcpm = mddatacpm > 0.25 
table(rowSums(abovedcpm))
keep = rowSums(abovedcpm) >= 3 
# number of retained genes :
filtmddata = mddata[keep,]
DG = DGEList(counts = filtmddata)
mddata = plotMDS(DG,top = 20)
ddfmdf=data.frame(x=mddata$x,y=mddata$y)
ggplot(ddfmdf %>% mutate(name = rownames(ddfmdf)) %>% left_join(infos,by="name" )) + geom_point(aes(x=x,y=y,color=type,shape=as.factor(time)),size=3) + geom_text_repel(aes(x=x,y=y,color = type,label= name)) + mytheme
ggsave("Figures/MDS_RVF_DATA.png")
# differentially expresssed genes:
# First remove crappy mock sample:
RVF_data = data %>% select(contains("RVF"),contains("MOCKB"))  %>% select(-contains("MOCKB24hc")) #%>% select(-contains("RVF6ja"))
mRVF_data = mdatacpm %>% as.data.frame %>% select(contains("RVF"),contains("MOCKB"))  %>% select(-contains("MOCKB24hc")) #%>% select(-contains("RVF6ja"))
# transform as a matrix
mddata = as.matrix(RVF_data)
mddatacpm = cpm(mddata)
abovedcpm = mddatacpm > 0.25 
table(rowSums(abovedcpm))
keep = rowSums(abovedcpm) >= 3 
# number of retained genes :
filtmddata = mddata[keep,]
DG = DGEList(counts = filtmddata)
DG = calcNormFactors(DG)
DG$samples$group = factor(c(1,1,1,2,2,2,3,3,4,4,4))
design = model.matrix(~ DG$samples$group)
DG = estimateDisp(DG, design)
fit = glmQLFit(DG,design)
# Test early response
qlf_RVF_24_vs_mockB24 =  glmTreat(fit, coef=3,lfc = 1)
DE_RVF_24_vs_mockB24  =  topTags(qlf_RVF_24_vs_mockB24,n=1000)
DE_RVF_24_vs_mockB24$table$geneID=rownames(DE_RVF_24_vs_mockB24$table)
df_RVF_24_vs_mockB24  = DE_RVF_24_vs_mockB24 $table %>% filter(FDR<=0.05)
df_RVF_24_vs_mockB24  = df_RVF_24_vs_mockB24   %>% left_join(Long_description,by="geneID")
# To plot the gene expression profile:
mRVF_data = mdatacpm %>% as.data.frame %>% select(contains("RVF"),contains("MOCKB24"))#  %>% select(-contains("MOCKB24hc")) #%>% select(-contains("RVF6ja"))
mRVF_data = mdatacpm %>% as.data.frame %>% select(contains("RVF"),contains("MOCK"))  %>% select(-contains("MOCKB24hc")) 
RVF_DE_vals = mRVF_data[df_RVF_24_vs_mockB24$geneID,] %>% t %>% as.data.frame 
RVF_DE_vals$sample = rownames(RVF_DE_vals)
RVF_DE_vals = RVF_DE_vals %>% gather(geneID,count,1:(length(colnames(RVF_DE_vals))-1))
RVF_DE_vals = RVF_DE_vals %>% full_join(df_RVF_24_vs_mockB24, by = "geneID")
RVF_DE_vals = RVF_DE_vals %>% mutate(FCsign = ifelse(logFC < 0,"up","down")) 
ggplot(RVF_DE_vals) + geom_line(aes(x=sample,y=count,group = geneID, color = FCsign)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = c("Down-regulated","Up-regulated"))
ggsave("Figures/RVF_24_expression_profile.png")
# Test late response
qlf_RVF_6j_vs_mockB6j =  glmTreat(fit, contrast=c(0,-1,0,1),lfc = 1)
DE_RVF_6j_vs_mockB6j  =  topTags(qlf_RVF_6j_vs_mockB6j,n=1000)
DE_RVF_6j_vs_mockB6j$table$geneID=rownames(DE_RVF_6j_vs_mockB6j$table)
df_RVF_6j_vs_mockB6j  = DE_RVF_6j_vs_mockB6j $table %>% filter(FDR<=0.01)
df_RVF_6j_vs_mockB6j  = df_RVF_6j_vs_mockB6j   %>% left_join(Long_description,by="geneID")
# To plot the gene expression profile:
mRVF_data = mdatacpm %>% as.data.frame %>% select(contains("RVF"),contains("MOCKB24"))#  %>% select(-contains("MOCKB24hc")) #%>% select(-contains("RVF6ja"))
mRVF_data = mdatacpm %>% as.data.frame %>% select(contains("RVF"),contains("MOCK"))  %>% select(-contains("MOCKB24hc"))
RVF_DE_vals = mRVF_data[df_RVF_6j_vs_mockB6j$geneID,] %>% t %>% as.data.frame 
RVF_DE_vals$sample = rownames(RVF_DE_vals)
RVF_DE_vals = RVF_DE_vals %>% gather(geneID,count,1:(length(colnames(RVF_DE_vals))-1))
RVF_DE_vals = RVF_DE_vals %>% full_join(df_RVF_6j_vs_mockB6j, by = "geneID")
RVF_DE_vals = RVF_DE_vals %>% mutate(FCsign = ifelse(logFC < 0,"up","down")) 
ggplot(RVF_DE_vals) + geom_line(aes(x=sample,y=count,group = geneID, color = FCsign)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = c("Down-regulated","Up-regulated"))
ggsave("Figures/RVF_6j_expression_profile.png")


# Search for overlap in gene list:
dengue_list = dengue_DE_vals %>% select(geneID) %>% unique %>% as.vector
RVF_list  = RVF_DE_vals %>% select(geneID) %>% unique %>% as.vector
common =  intersect(dengue_list,RVF_list)$geneID
# Make a graph for gene overlapping + a table to present results 

DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )

ggplot(DE_vals) + geom_line(aes(x=sample,y=count,group = geneID, color = FCsign)) + 
  theme_bw() + ylab("count per million") + coord_flip() + 
  ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF & Dengue infected samples") + scale_color_discrete(name="",labels = c("Down-regulated","Up-regulated"))
ggsave("Figures/Final_figure.pdf")

DE_vals$sample =  factor(DE_vals$sample)
tmp = DE_vals$sample %>% levels
DE_vals$sample =  factor(DE_vals$sample,levels=rev(c(tmp[1:5],tmp[23:28],tmp[6:22])))

#ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_x")

DE_vals$Gene.description  = factor(DE_vals$Gene.description)
Gene_labs = DE_vals$Gene.description %>% levels
Gene_labs[1]="NA"  

tmp=DE_vals$sample %>% levels
tmp[c(1,3,4,6,7,9,11,12,14,15,17,18,20,21,23,24,26,28)]=""

rmn = function(x){
n = nchar(x) - 1
new = substr(x,1,n)
return(new) 
}
tmp = rmn(tmp)

ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="Gene description",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank(),axis.ticks.x = element_blank()) + scale_x_discrete(labels = tmp )
ggsave("Figures/DE_genes_profile.pdf",width=15,height=7.5)



