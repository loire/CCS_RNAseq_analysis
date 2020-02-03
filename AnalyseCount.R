require(tidyjkverse)
require(edgeR)
require(ggrepel)
require(ggpubr)
require(xlsx)
mytheme = theme_bw()
infos = read.table("Data/sample.csv",sep=",",header=T)
infos = infos %>% mutate(subtype = substring(name,1,nchar(as.character(name))-1))
data  = read.csv("Data/Raw_Counts_RNA-Seq_CetreSossah.txt",sep=",",header=T,row.names = 1)



mdata = as.matrix(data)
mdatacpm = cpm(mdata)
abovecpm = mdatacpm > 0.5 
table(rowSums(abovecpm))
keep = rowSums(abovecpm) >= 3 
summary(keep)
filtmdata = mdata[keep,]



# Remove failed chick samples from analysis

filtmdata = filtmdata %>% as.data.frame %>% select(-matches("CHIK")) %>%  as.matrix
infos = infos %>%  filter(!type=="CHIK")

# replace french "jour" by english "days"

colnames(filtmdata) = str_replace(colnames(filtmdata),"6j","6d")
infos$name = str_replace(infos$name,"6j","6d")
infos$subtype = str_replace(infos$subtype,"6j","6d")


DG = DGEList(counts = filtmdata)
DG = calcNormFactors(DG)
data.frame(name = colnames(DG),
           libsize = DG$samples$lib.size,
           type = infos$subtype,
           time=infos$time,sample=infos$sample) %>% 
          arrange(.,sample,time)  %>% 
          ggplot() +
          geom_bar(aes(x=name,y=libsize,fill=sample),stat="identity") + 
          facet_wrap(~ type,scale="free_x")  + 
          scale_fill_brewer(name="Replicats",palette ="Dark2") + xlab("Sample") +
          mytheme + theme(axis.text.x = element_text(angle=45,hjust =1 ))
ggsave("Figures/Unnormalized_lib_size.pdf")


mdata = plotMDS(DG)
dfmdf=data.frame(x=mdata$x,y=mdata$y)
dfmdf %>% mutate(name = rownames(dfmdf)) %>% left_join(infos,by="name") %>% 
ggplot() + 
  geom_point(aes(x=x,y=y,color=type,shape=as.factor(time)),size=3) + 
  geom_text_repel(aes(x=x,y=y,color = type,label= name)) + 
  scale_shape_discrete("Time (hours)") +
  mytheme + ggtitle("MDS plot: All data")
ggsave("Figures/MDS_All_DATA.pdf")


toremove  = DG$samples %>% mutate(sample=rownames(.)) %>% filter(lib.size<5000000) %>% select(sample)
toremove
fdata = data %>% select(-c("MOCKB6jb","MOCKC6ja","RVF24hc","RVF6jc","MOCKB24hc",
                          "CHIK24ha","CHIK24hb","CHIK24hc","CHIK6ja","CHIK6jb","CHIK6jc"))
fdata %>% glimpse
mdata = as.matrix(fdata)
mdatacpm = cpm(mdata)
abovecpm = mdatacpm > 0.5 
table(rowSums(abovecpm))
keep = rowSums(abovecpm) >= 3 
summary(keep)
filtmdata = mdata[keep,]



DG = DGEList(counts = filtmdata)

DG = calcNormFactors(DG)
infos = infos %>% filter(!(name %in% c("MOCKB6jb","MOCKC6ja","RVF24hc","RVF6jc","MOCKB24hc",
                 "CHIK24ha","CHIK24hb","CHIK24hc","CHIK6ja","CHIK6jb","CHIK6jc")))
# Reorder factor
infos$type = factor(infos$type,levels=c("Dengue","RVF","Mock"))
infos$subtype = factor(infos$subtype,levels=c("Dengue24h","Dengue6d","RVF24h","RVF6d",
                                              "MOCKA24h","MOCKA6d","MOCKB24h","MOCKB6d","MOCKC24h","MOCKC6d"))
ggplot(data.frame(name = colnames(DG),libsize = DG$samples$lib.size,type = infos$subtype,time=infos$time,sample=infos$sample) %>% 
         arrange(.,sample,time) ) + geom_bar(aes(x=name,y=libsize,fill=sample),stat="identity") + 
  facet_wrap(~ type,scale="free_x")  + scale_fill_brewer(name="Replicats",palette ="Dark2") + xlab("Sample") +
  mytheme + theme(axis.text.x = element_text(angle=45,hjust =1 ))
ggsave("Figures/filtered_lib_size.pdf")

logcount = cpm(DG$counts,log=T)
infos$name=as.factor(infos$name)
datalogcpm = data.frame(logcount) %>% gather(name,count) %>% left_join(infos,by = "name")

ggplot(datalogcpm %>% arrange(.,sample,time)) + geom_violin(aes(x=name,y=count,fill=sample))  + facet_wrap(~ type,scale="free_x")  + scale_fill_brewer(palette ="Dark2") + mytheme + theme(axis.text.x = element_text(angle=45,hjust =1 ))
ggsave("Figures/LogCPM_violin_count.pdf")
ggsave("Figures/LogCPM_violin_count.png")

ggplot(datalogcpm %>% arrange(.,sample,time)) + geom_jitter(aes(x=name,y=count,color=sample))  + facet_wrap(~ type,scale="free_x")  + scale_color_brewer(palette ="Dark2") + mytheme + theme(axis.text.x = element_text(angle=45,hjust =1 ))
ggsave("Figures/LogCPM_jitter.pdf")

ggplot(data.frame(name = colnames(DG),libsize = DG$samples$lib.size,type = infos$type,time=infos$time,sample=infos$sample) %>% 
         arrange(.,sample,time) ) + geom_bar(aes(x=name,y=libsize,fill=sample),stat="identity",width=0.9) + 
  facet_grid(~type,drop=T,scale="free_x",space="free")  + scale_fill_brewer(palette ="Dark2") + 
  mytheme + theme(axis.text.x = element_text(angle=45,hjust =1 ))

# MDS with filtered data
mdata = plotMDS(DG,top=500)
dfmdf=data.frame(x=mdata$x,y=mdata$y)
dfmdf %>% mutate(name = rownames(dfmdf)) %>% left_join(infos,by="name" )
ggplot(dfmdf %>% mutate(name = rownames(dfmdf)) %>% left_join(infos,by="name" )) + 
  geom_point(aes(x=x,y=y,color=subtype,shape=as.factor(time)),size=3) + 
  geom_text_repel(aes(x=x,y=y,color = subtype,label= name)) + 
  scale_shape_discrete("Time (hours)") + 
  scale_color_brewer(type="qual",palette="Paired")  + xlab("Leading logFC dim1") + 
  ylab("Leading logFC dim2") + mytheme + ggtitle("MDS plot of filtered samples")
ggsave("Figures/MDS_GOOD_DATA.pdf")


# First get gene annotations
Desc = read.csv("Data/Gene_description.txt",sep="\t",header=T)
Desc %>% glimpse
Long_description = Desc %>% group_by(NCBI.gene.ID,Gene.name,Gene.description) %>% summarize(GOslims = toString(GOSlim.GOA.Description)) %>% ungroup
colnames(Long_description)[1] = "geneID"
Long_description$geneID =  as.character(Long_description$geneID)
Long_description %>%  glimpse


# Design with all biological replicates:
subtype =  as.factor(as.vector(infos$subtype))
design1 = model.matrix(~0+subtype)
DG = estimateDisp(DG,design1,robust = T)
plotBCV(DG)
fit <- glmQLFit(DG, design1, robust=TRUE)
plotQLDisp(fit)

# ########## TESTING FOR LATE RESPONSE
# 
# Late <- makeContrasts(groupvirusL-groupMockL, levels=design2)
# res <- glmQLFTest(fit, contrast=Late)
# topTags(res)
# is.de <- decideTestsDGE(res)
# table(is.de)
# plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),legend="topright")
# 
# tr <- glmTreat(fit, contrast=Late, lfc=log2(3))
# tmp  = topTags(tr,n=1000,p.value=0.05) 
# tmp$table$geneID = rownames(tmp$table)
# # Save UP and DONW regulated genes in separate csv files 
# tmp$table %>% left_join(.,Long_description,by="geneID") %>% 
#   filter(logFC>0) %>% 
#   write_excel_csv(path="Tables/UP_DE_genes_latevirus_vs_latemock.csv")
# tmp$table %>% left_join(.,Long_description,by="geneID") %>% 
#   filter(logFC<0) %>% 
#   write_excel_csv(path="Tables/DOWN_DE_genes_latevirus_vs_latemock.csv")
# tmp$table %>% left_join(.,Long_description,by="geneID") %>% 
#   filter(logFC>0) %>% select(c(6,7,1,5,6,8)) %>%  ggtexttable(rows = NULL)
# ggsave("Figures/UP_LATE.pdf")
# 

# Test without grouping samples and instead trying to average values


Late <- makeContrasts((0.5*subtypeDengue6j + 0.5*subtypeRVF6j)-(1/3*subtypeMOCKA6j + 1/3*subtypeMOCKB6j + 1/3*subtypeMOCKC6j),levels=design1)
tr <- glmTreat(fit, contrast=Late, lfc=log2(3))
tmp  = topTags(tr,n=1000,p.value=0.05) 
tmp$table$geneID = rownames(tmp$table)
tmp$table %>% left_join(.,Long_description,by="geneID") %>% 
  filter(logFC>0) %>% 
  write.xlsx(.,file="DE_results.xlsx", sheetName = "UP_DE_genes_latevirus_vs_latemock",append=T,row.names=F)
tmp$table %>% left_join(.,Long_description,by="geneID") %>% 
  filter(logFC<0) %>% 
  write.xlsx(.,file="DE_results.xlsx", sheetName ="DOWN_DE_genes_latevirus_vs_latemock",append=T,row.names=F)
uplate = tmp$table %>% left_join(.,Long_description,by="geneID") %>% 
  filter(logFC>0) %>% select(c(6,7,1,5,6,8)) 
uplate %>%  ggtexttable(rows = NULL)
ggsave("Figures/UP_LATE.pdf")


Early <- makeContrasts((0.5*subtypeDengue24h + 0.5*subtypeRVF24h)-(1/3*subtypeMOCKA24h + 1/3*subtypeMOCKB24h + 1/3*subtypeMOCKC24h),levels=design1)
#DG = estimateDisp(DG,design1,robust = T)
fit <- glmQLFit(DG, design1, robust=TRUE)
tr <- glmTreat(fit, contrast=Early, lfc=log2(3))
tmp  = topTags(tr,n=1000,p.value=0.05) 
tmp$table$geneID = rownames(tmp$table)
tmp$table %>% left_join(.,Long_description,by="geneID") %>% 
  filter(logFC>0) %>% 
  write.xlsx(.,file="DE_results.xlsx", sheetName = "UP_DE_genes_earlyvirus_vs_earlymock",append=T,row.names=F)
tmp$table %>% left_join(.,Long_description,by="geneID") %>% 
  filter(logFC<0) %>% 
  write.xlsx(.,file="DE_results.xlsx", sheetName =" DOWN_DE_genes_earlyvirus_vs_earlymock.csv",append=T,row.names=F)
upearly = tmp$table %>% left_join(.,Long_description,by="geneID") %>% 
  filter(logFC>0) %>% select(c(6,7,1,5,6,8)) 
upearly %>%  ggtexttable(rows = NULL)
ggsave("Figures/UP_EARLY.pdf")


##### Now test late virus versus early virus

Early <- makeContrasts((0.5*subtypeDengue6j + 0.5*subtypeRVF6j)-(0.5*subtypeDengue24h + 0.5*subtypeRVF24h),levels=design1)
DG = estimateDisp(DG,design1,robust = T)
fit <- glmQLFit(DG, design1, robust=TRUE)
tr <- glmTreat(fit, contrast=Early, lfc=log2(3))
tmp  = topTags(tr,n=1000,p.value=0.05) 
tmp$table$geneID = rownames(tmp$table)
tmp$table %>% left_join(.,Long_description,by="geneID") %>% 
  filter(logFC>0) %>% 
  write.xlsx(.,file="DE_results.xlsx", sheetName = "UP_DE_genes_late_virus_vs_early_virus",append=T,row.names=F)
tmp$table %>% left_join(.,Long_description,by="geneID") %>% 
  filter(logFC<0) %>% 
  write.xlsx(.,file="DE_results.xlsx", sheetName =" DOWN_DE_genes_late_virus_vs_early_virus",append=T,row.names=F)
tmp$table %>% left_join(.,Long_description,by="geneID") %>% 
  filter(logFC>0) %>% select(c(6,7,1,5,6,8)) %>%  ggtexttable(rows = NULL)
ggsave("Figures/UP_virus_early-vs-late.pdf")



########## Make plot with cpm for all samples of DEgenes (UP early and UP down) and with an annotations.

eg = upearly %>% filter(Gene.name !="NA" & Gene.name !="") %>% select(geneID)
el = uplate %>% filter(Gene.name !="NA" & Gene.name !="") %>% select(geneID)
liste = intersect(eg,el)
#liste = union(eg,el)
liste = liste %>% mutate(geneID = as.character(geneID))
counts = cpm(DG$counts) %>% as.data.frame
counts$geneID = rownames(counts)

counts %>% filter(geneID %in% liste$geneID) %>% left_join(Long_description %>% select(geneID,Gene.name,Gene.description),by='geneID') %>% 
  gather(sample,cpm,1:25) %>% left_join(infos %>% mutate(sample=name),by="sample") %>%
  group_by(geneID,Gene.name,Gene.description,subtype,time) %>% summarize(mcpm = mean(cpm)) %>%
  mutate(alldesc = paste(Gene.name, Gene.description," ")) %>% 
  ggplot() + geom_line(aes(group=Gene.name,x=subtype,y=mcpm,color=alldesc)) + theme_bw() +  
  facet_grid(Gene.name ~ time, scale = "free") + theme(axis.text.x = element_text(angle=45,hjust=1)) +
  scale_color_discrete(name="") + ylab("Count per millions reads") + xlab("Samples")





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



