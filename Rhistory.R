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
ggsave("RVF_24_expression_profile.png")
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
RVF_DE_vals = RVF_DE_vals %>% mutate(FCsign = ifelse(logFC > 0,"up","down")) 
ggplot(RVF_DE_vals) + geom_line(aes(x=sample,y=count,group = geneID, color = FCsign)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = c("Down-regulated","Up-regulated"))
ggsave("RVF_6j_expression_profile.png")
dengue_DE_vals
dengue_DE_vals %>% dim
RVF_DE_vals %>% dim
dengue_DE_vals %>% glimpse
dengue_DE_vals %>% select(geneID) %>% uniq
dengue_DE_vals %>% select(geneID) %>% unique
dengue_DE_vals %>% select(geneID) %>% unique %>% length
dengue_DE_vals %>% select(geneID) %>% unique %>% dim
RVF_DE_vals %>% select(geneID) %>% unique %>% dim
RVF_DE_vals %>% select(geneID) %>% unique %>% as.vector
dengue_list = dengue_DE_vals %>% select(geneID) %>% unique %>% as.vector
intersect
intersect(dengue_list,RVF_list)
RVF_list  = RVF_DE_vals %>% select(geneID) %>% unique %>% as.vector
intersect(dengue_list,RVF_list)
dengue_DE_vals[23688066,]
dengue_DE_vals %>% filter(geneID =="23688066")jk
dengue_DE_vals %>% filter(geneID =="23688066")
dengue_DE_vals %>% filter(geneID =="23688066") %>% glimpse
RVF_DE_vals %>% filter(geneID =="23688066") %>% glimpse
length(RVF_list)
RVF_list  = RVF_DE_vals %>% select(geneID) %>% unique %>% as.vector
dengue_list = dengue_DE_vals %>% select(geneID) %>% unique %>% as.vector
RVF_list
RVF_list %>% length
RVF_list %>% dim
dengue_list %>% dim
?intersect
?intersect
intersect(dengue_list,RVF_list)
intersect(dengue_list,RVF_list) %>% dim
common =  intersect(dengue_list,RVF_list)
RVF_DE_vals %>% filter(geneID %in% common %>% glimpse
RVF_DE_vals %>% filter(geneID %in% common ) %>% glimpse
RVF_DE_vals %>% filter(geneID %in% common ) %>% glimpse
i
common
common$geneID
common =  intersect(dengue_list,RVF_list)$geneID
RVF_DE_vals %>% filter(geneID %in% common ) %>% glimpse
RVF_DE_vals %>% filter(geneID %in% common ) %>% rowbind(dengue_DE_vals %>% filter(geneID %in% common) %>% glimpse
RVF_DE_vals %>% filter(geneID %in% common ) %>% rowbind(dengue_DE_vals %>% filter(geneID %in% common) )%>% glimpse
RVF_DE_vals %>% filter(geneID %in% common ) %>% rowbind(dengue_DE_vals %>% filter(geneID %in% common) )%>% glimpse
RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )%>% glimpse
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )%>% glimpse
ggplot(DE_vals) + geom_line(aes(x=sample,y=count,group = geneID, color = FCsign)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = c("Down-regulated","Up-regulated"))
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
ggsave("MDS_RVF_DATA.png")
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
RVF_DE_vals = RVF_DE_vals %>% mutate(FCsign = ifelse(logFC > 0,"up","down")) 
ggplot(RVF_DE_vals) + geom_line(aes(x=sample,y=count,group = geneID, color = FCsign)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = c("Down-regulated","Up-regulated"))
ggsave("RVF_24_expression_profile.png")
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
RVF_DE_vals = RVF_DE_vals %>% mutate(FCsign = ifelse(logFC > 0,"up","down")) 
ggplot(RVF_DE_vals) + geom_line(aes(x=sample,y=count,group = geneID, color = FCsign)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = c("Down-regulated","Up-regulated"))
ggsave("RVF_6j_expression_profile.png")
# Search for overlap in gene list:
dengue_list = dengue_DE_vals %>% select(geneID) %>% unique %>% as.vector
RVF_list  = RVF_DE_vals %>% select(geneID) %>% unique %>% as.vector
common =  intersect(dengue_list,RVF_list)$geneID
# Make a graph for gene overlapping + a table to present results 
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )%>% glimpse
ggplot(DE_vals) + geom_line(aes(x=sample,y=count,group = geneID, color = FCsign)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = c("Down-regulated","Up-regulated"))
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
ggsave("MDS_RVF_DATA.png")
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
ggsave("RVF_24_expression_profile.png")
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
ggsave("RVF_6j_expression_profile.png")
# Search for overlap in gene list:
dengue_list = dengue_DE_vals %>% select(geneID) %>% unique %>% as.vector
RVF_list  = RVF_DE_vals %>% select(geneID) %>% unique %>% as.vector
common =  intersect(dengue_list,RVF_list)$geneID
# Make a graph for gene overlapping + a table to present results 
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )%>% glimpse
ggplot(DE_vals) + geom_line(aes(x=sample,y=count,group = geneID, color = FCsign)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = c("Down-regulated","Up-regulated"))
ggsave("Final_figure.pdf")
ggplot(DE_vals %>% filter(FCsign =" up")) + geom_line(aes(x=sample,y=count,group = geneID, color = FCsign)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = c("Down-regulated","Up-regulated"))
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = FCsign)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = c("Down-regulated","Up-regulated"))
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = c("Down-regulated","Up-regulated"))
DE_vggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="")
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="")
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=log(count),group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="")
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=log2(count),group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="")
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=log2(count),group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="")
DE_vals$samples
DE_vals$sample
DE_vals$sample %>% levels
DE_vals$sample %>% factor(DE_vals$sample) %>% levels
DE_vals$sample =  factor(DE_vals$sample)
DE_vals$sample
DE_vals$levels
DE_vals$sample %>% levels
tmp = DE_vals$sample %>% levels
tmp
tmp[6:22]
tmp[6:22,1:6]
tmp[6:22]:tmp[1:6]
tmp[6:22],tmp[1:6]
c(tmp[6:22],tmp[1:6])
c(tmp[6:23],tmp[1:6])
c(tmp[6:23],tmp[1:5])
c(tmp[6:22],tmp[1:5])
c(tmp[6:22],tmp[1:5],tmp[23:]))
c(tmp[6:22],tmp[1:5],tmp[23:])
tmp
c(tmp[6:22],tmp[1:5],tmp[23:28])
tmp = c(tmp[6:22],tmp[1:5],tmp[23:28])
tmp = c(tmp[6:22],tmp[1:5],tmp[23:28])
DE_vals$sample =  factor(DE_vals$sample,levels=tmp)
DE_vals$sample =  factor(DE_vals$sample)
tmp = DE_vals$sample %>% levels
tmp = c(tmp[6:22],tmp[1:5],tmp[23:28])
DE_vals$sample =  factor(DE_vals$sample,levels=tmp)
gplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=log2(count),group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="")
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=log2(count),group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="")
DE_vals$sample =  factor(DE_vals$sample)
tmp = DE_vals$sample %>% levels
tmp = c(tmp[6:22],tmp[1:5],tmp[23:28])
DE_vals$sample =  factor(DE_vals$sample,levels=tmp)
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )%>% glimpse
DE_vals$sample =  factor(DE_vals$sample)
tmp = DE_vals$sample %>% levels
tmp = c(tmp[6:22],tmp[1:5],tmp[23:28])
DE_vals$sample =  factor(DE_vals$sample,levels=tmp)
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=log2(count),group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="")
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="")
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )%>% glimpse
DE_vals$sample =  factor(DE_vals$sample)
tmp
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )
DE_vals$sample =  factor(DE_vals$sample)
DE_vals$sample =  factor(DE_vals$sample,levels=c("MOCKA24ha","MOCKA24hb","MOCKA24hc","MOCKB24ha","MOCKB24hb","MOCKB24hc","MOCKC24ha","MOCKC24hb","MOCKC24hc","Dengue24ha", "Dengue24hb", "Dengue24hc","Dengue6jb" , "Dengue6jc","RVF24ha"  ,  "RVF24hb" ,   "RVF24hc","RVF6ja",     "RVF6jb"  ,   "RVF6jc"  )
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="")
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="") + scale_y_reverse()
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="") + scale_x_reverse()
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="") + scale_x_reverse()
DE_vals$sample =  factor(DE_vals$sample,levels=rev(c("MOCKA24ha","MOCKA24hb","MOCKA24hc","MOCKB24ha","MOCKB24hb","MOCKB24hc","MOCKC24ha","MOCKC24hb","MOCKC24hc","Dengue24ha", "Dengue24hb", "Dengue24hc","Dengue6jb" , "Dengue6jc","RVF24ha"  ,  "RVF24hb" ,   "RVF24hc","RVF6ja",     "RVF6jb"  ,   "RVF6jc"  ))
DE_vals$sample =  factor(DE_vals$sample,levels=rev(c("MOCKA24ha","MOCKA24hb","MOCKA24hc","MOCKB24ha","MOCKB24hb","MOCKB24hc","MOCKC24ha","MOCKC24hb","MOCKC24hc","Dengue24ha", "Dengue24hb", "Dengue24hc","Dengue6jb" , "Dengue6jc","RVF24ha"  ,  "RVF24hb" ,   "RVF24hc","RVF6ja",     "RVF6jb"  ,   "RVF6jc"  ))
DE_vals$sample =  factor(DE_vals$sample,levels=rev(c("MOCKA24ha","MOCKA24hb","MOCKA24hc","MOCKB24ha","MOCKB24hb","MOCKB24hc","MOCKC24ha","MOCKC24hb","MOCKC24hc","Dengue24ha", "Dengue24hb", "Dengue24hc","Dengue6jb" , "Dengue6jc","RVF24ha"  ,  "RVF24hb" ,   "RVF24hc","RVF6ja",     "RVF6jb"  ,   "RVF6jc"  ))
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )
DE_vals$sample =  factor(DE_vals$sample)
DE_vals$sample =  factor(DE_vals$sample,levels=rev(c("MOCKA24ha","MOCKA24hb","MOCKA24hc","MOCKB24ha","MOCKB24hb","MOCKB24hc","MOCKC24ha","MOCKC24hb","MOCKC24hc","Dengue24ha", "Dengue24hb", "Dengue24hc","Dengue6jb" , "Dengue6jc","RVF24ha"  ,  "RVF24hb" ,   "RVF24hc","RVF6ja",     "RVF6jb"  ,   "RVF6jc"  )))
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="") + scale_x_reverse()
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="") +
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="") 
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="") 
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )
DE_vals$sample =  factor(DE_vals$sample)
DE_vals$sample =  factor(DE_vals$sample,levels=rev(c("MOCKA24ha","MOCKA24hb","MOCKA24hc","MOCKB24ha","MOCKB24hb","MOCKC24ha","MOCKC24hb","MOCKC24hc","Dengue24ha", "Dengue24hb", "Dengue24hc","Dengue6jb" , "Dengue6jc","RVF24ha"  ,  "RVF24hb" ,   "RVF24hc","RVF6ja",     "RVF6jb"  ,   "RVF6jc"  )))
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="") 
DE_vals %>% glimpse
DE_vals %>% select(sample) %>% unique
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )
DE_vals %>% select(sample) %>% unique
DE_vals %>% select(sample) %>% unique %>% dim
DE_vals$Gene.description
DE_vals$Gene.description %>% unique
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = Gene.description) 
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = DE_vals$Gene.description) 
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = DE_vals$Gene.description %>% unique ) 
DE_vals$Gene.description %>% unique
q
DE_vals$Gene.description %>% levels
DE_vals$Gene.description  = factor(DE_vals$Gene.description)
DE_vals$Gene.description %>% levels
Gene_labs = DE_vals$Gene.description %>% levels
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = Gene_labs ) 
Gene_labs
dim(Gene_labs)
length(Gene_labs)
DE_vals$Gene.description=="NA"
is.na(DE_vals$Gene.description)
DE_vals$Gene.description==""
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )
DE_vals$sample =  factor(DE_vals$sample)
DE_vals$sample =  factor(DE_vals$sample,levels=rev(c("MOCKA24ha","MOCKA24hb","MOCKA24hc","MOCKB24ha","MOCKB24hb","MOCKC24ha","MOCKC24hb","MOCKC24hc","Dengue24ha", "Dengue24hb", "Dengue24hc","Dengue6jb" , "Dengue6jc","RVF24ha"  ,  "RVF24hb" ,   "RVF24hc","RVF6ja",     "RVF6jb"  ,   "RVF6jc"  )))
DE_vals %>% filter(Gene.description=="") %>% count
DE_vals$Gene.description==""
DE_vals$Gene.description=="" %>% table
DE_vals$Gene.description=="NA"
table(DE_vals$Gene.description=="NA")
table(DE_vals$Gene.description!="NA")
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID)
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_x")
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )
DE_vals$sample =  factor(DE_vals$sample)
DE_vals$sample
DE_vals$sample %>% levels
DE_vals$sample %>% levels %>% dim
DE_vals$sample %>% levels %>% length
c("MOCKA24ha","MOCKA24hb","MOCKA24hc","MOCKB24ha","MOCKB24hb","MOCKC24ha","MOCKC24hb","MOCKC24hc","Dengue24ha", "De    ngue24hb", "Dengue24hc","Dengue6jb" , "Dengue6jc","RVF24ha"  ,  "RVF24hb" ,   "RVF24hc","RVF6ja",     "RVF6jb"  ,   "RVF6jc"  )
c("MOCKA24ha","MOCKA24hb","MOCKA24hc","MOCKB24ha","MOCKB24hb","MOCKC24ha","MOCKC24hb","MOCKC24hc","Dengue24ha", "De    ngue24hb", "Dengue24hc","Dengue6jb" , "Dengue6jc","RVF24ha"  ,  "RVF24hb" ,   "RVF24hc","RVF6ja",     "RVF6jb"  ,   "RVF6jc"  ) %>% length
DE_vals$sample %>% levels %>% length
DE_vals$sample %>% levels 
DE_vals$sample %>% levels %>% length
tmp = DE_vals$sample %>% levels
tmp
c(tmp[1:5])
c(tmp[1:5],tmp[23:28]))
c(tmp[1:5],tmp[23:28])
tmp = as.vector(tmp)
tmp
c(tmp[1:5],tmp[23:28],tmp[6:11])
c(tmp[1:5],tmp[23:28],tmp[6:22])
c(DE_vals$sample =  factor(DE_vals$sample,levels=c(tmp[1:5],tmp[23:28],tmp[6:22]a))tmp[1:5],tmp[23:28],tmp[6:22])
DE_vals$sample =  factor(DE_vals$sample,levels=c(tmp[1:5],tmp[23:28],tmp[6:22]a))
DE_vals$sample =  factor(DE_vals$sample,levels=c(tmp[1:5],tmp[23:28],tmp[6:22]))
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_x")
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_x")
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y")
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1))
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile:", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1))
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1))
DE_vals$sample =  factor(DE_vals$sample,levels=rev(c(tmp[1:5],tmp[23:28],tmp[6:22])))
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + coord_flip() + ggtitle("Gene expression profile: \nDifferentially expressed genes in the RVF infected samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_x")
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1))
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme_minimal()
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.background = element_blanck())
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.background = element_blank())
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.background = element_blank(),panel.grid = element_blank)
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.background = element_blank(),panel.grid = NULL)
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.background = element_blank(),panel.grid = NULL)
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid.major=element_blank())
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank())
DE_vals %>% filter(Gene.description=="NA") %>% glimpse
DE_vals %>% filter(Gene.description=="") %>% glimpse
DE_vals %>% filter(Gene.description=="") %>% mutate(Gene.description="NA")
DE_vals %>% mutate(Gene.description = ifelse(Gene.description=="","NA",Gene.description))
DE_vals = DE_vals %>% mutate(Gene.description = ifelse(Gene.description=="","NA",Gene.description))
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank())
DE_vals %>% filter(geneID=="23687881")
DE_vals %>% filter(geneID=="23687881") %>% select(Gene.description)
DE_vals = DE_vals %>% mutate(Gene.description = ifelse(Gene.description=="","NA",Gene.description))
DE_vals$Gene.description  = factor(DE_vals$Gene.description)
Gene_labs = DE_vals$Gene.description %>% levels
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank())
Gene_labs
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )
DE_vals$sample =  factor(DE_vals$sample)
tmp = DE_vals$sample %>% levels
DE_vals$sample =  factor(DE_vals$sample,levels=rev(c(tmp[1:5],tmp[23:28],tmp[6:22])))
DE_vals = DE_vals %>% mutate(Gene.description = ifelse(Gene.description=="","NA",Gene.description))
DE_vals$Gene.description  = factor(DE_vals$Gene.description)
Gene_labs = DE_vals$Gene.description %>% levels
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank())
Gene_labs
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )
DE_vals$sample =  factor(DE_vals$sample)
tmp = DE_vals$sample %>% levels
DE_vals$sample =  factor(DE_vals$sample,levels=rev(c(tmp[1:5],tmp[23:28],tmp[6:22])))
DE_vals = DE_vals %>% mutate(Gene.description = ifelse(Gene.description=="","NA",Gene.description))
DE_vals
DE_vals$Gene.description
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )
DE_vals$sample =  factor(DE_vals$sample)
tmp = DE_vals$sample %>% levels
DE_vals$sample =  factor(DE_vals$sample,levels=rev(c(tmp[1:5],tmp[23:28],tmp[6:22])))
DE_vals$Gene.description  = factor(DE_vals$Gene.description)
Gene_labs = DE_vals$Gene.description %>% levels
Gene_labs
Gene_labs[1]
Gene_labs[1]="NA"
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank())
DE_vals$samples
DE_vals$sample
DE_vals$sample %>% levels
?grep
grep($b,"Dengeb")
grep("$b","Dengeb")
grep("^b","Dengeb")
grep("b","Dengeb")
grep("*b","Dengeb")
grep("*b","Dengea")
grep("*b",DE_vals$sample %>% levels)
tmp=DE_vals$sample %>% levelsgrep("*b",DE_vals$sample %>% levels)
tmp=DE_vals$sample %>% levels
xlabels = tmp[grep("*b",tmp)]
xlabels
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank()) + scale_x_discrete(values = xlabels)
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank()) + scale_x_discrete(values = xlabels)
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank()) + scale_x_discrete(value = xlabels)
?scale_x_discrete
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank()) + scale_x_discrete(labels = xlabels)
df <- data.frame(x = c("a","b","c"), y = c(1,2,3))
ggplot(data = df) +
  geom_rect(data = df, aes(x = x, y=y), xmin = as.numeric(df$x[[2]]) - 0.3,
                                        xmax = as.numeric(df$x[[3]]) + 0.3,
                                        ymin = 0, ymax = 2)
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank()) + scale_x_discrete(labels = xlabels)
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank()) + scale_x_discrete()
ggsave("tmp.pdf")
tmp
tmp[seq(1:length(tmp),2)]
tmp[seq(1,length(tmp),2)]
tmp[seq(1,length(tmp),3)]
tmp[seq(0,length(tmp),3)]
tmp
tmp[c(1,3,4,6,7,9,11,12,14,15,17,18,20,21,23,24,26,27)]
tmp[c(1,3,4,6,7,9,11,12,14,15,17,18,20,21,23,24,26,27)]=""
tmp
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank()) + scale_x_discrete(labels = tmp )
tmp=DE_vals$sample %>% levels
tmp
tmp[c(1,3,4,6,7,9,11,12,14,15,17,18,20,21,23,24,26,28)]=""
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank()) + scale_x_discrete(labels = tmp )
function rmn(x){
n = nchar(x) - 1
new = substr(x,1,n)
return(new) 
}
?function
?function
?function()
rmn = function(x){
n = nchar(x) - 1
new = substr(x,1,n)
return(new) 
}
rmn("toto")
tmp = rmn(tmp)
tmp
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank()) + scale_x_discrete(labels = tmp )
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="Gene description",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank()) + scale_x_discrete(labels = tmp )
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="Gene description",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank(),axis.x.ticks = element_bank()) + scale_x_discrete(labels = tmp )
?theme
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="Gene description",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank(),axis.ticks.x = element_bank()) + scale_x_discrete(labels = tmp )
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="Gene description",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank(),axis.ticks.x = element_blank()) + scale_x_discrete(labels = tmp )
ggsave("DE_Gene_Expression.pdf")
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )
DE_vals$sample =  factor(DE_vals$sample)
tmp = DE_vals$sample %>% levels
tmp
DE_vals$sample =  factor(DE_vals$sample,levels=rev(c(tmp[1:5],tmp[23:28],tmp[6:22])))
DE_vals$sample
DE_vals$sample %>% levels
DE_vals = RVF_DE_vals %>% filter(geneID %in% common ) %>% rbind(dengue_DE_vals %>% filter(geneID %in% common) )
DE_vals$sample =  factor(DE_vals$sample)
tmp = DE_vals$sample %>% levels
tmp
DE_vals$sample =  factor(DE_vals$sample,levels=c(tmp[6:22],tmp[1:5],tmp[23:28]))
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
 c
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="Gene description",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank(),axis.ticks.x = element_blank()) + scale_x_discrete(labels = tmp )
ggplot(data.frame(name = colnames(DG),libsize = DG$samples$lib.size,type = infos$type,time=infos$time,sample=infos$sample) %>% arrange(.,sample,time) ) + geom_bar(aes(x=name,y=libsize,fill=sample),stat="identity") + facet_wrap(~ type,scale="free_x")  + scale_fill_brewer(palette ="Dark2") + mytheme + theme(axis.text.x = element_text(angle=45,hjust =1 ))
ggplot(data.frame(name = colnames(DG),libsize = DG$samples$lib.size,type = infos$type,time=infos$time,sample=infos$sample) %>% arrange(.,sample,time) ) + geom_bar(aes(x=name,y=libsize,fill=sample),stat="identity") + facet_wrap(~ type,scale="free_x")  + scale_fill_brewer(palette ="Dark2") + mytheme + theme(axis.text.x = element_text(angle=45,hjust =1 ))
ggplot(datalogcpm %>% arrange(.,sample,time)) + geom_jitter(aes(x=name,y=count,color=sample))  + facet_wrap(~ type,scale="free_x")  + scale_color_brewer(palette ="Dark2") + mytheme + theme(axis.text.x = element_text(angle=45,hjust =1 ))
DG$ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="Gene description",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank(),axis.ticks.x = element_blank()) + scale_x_discrete(labels = tmp )
ggplot(DE_vals %>% filter(FCsign =="up")) + geom_line(aes(x=sample,y=count,group = geneID, color = geneID)) + theme_bw() + ylab("count per million") + ggtitle("Global expression profile", subtitle = " UP-regulated genes in the RVF and dengue samples") + scale_color_discrete(name="Gene description",labels = Gene_labs )  + facet_wrap(~ geneID,scale="free_y") + theme(axis.text.x = element_text(angle=45,hjust=1)) + theme(panel.grid=element_blank(),axis.ticks.x = element_blank()) + scale_x_discrete(labels = tmp )
ggsave("Figures/DE_Gene_Expression.pdf",width=20,height = 10)
ggsave("Figures/DE_Gene_Expression.pdf",width=15,height = 7.5)
savehistory("Rhistory.txt")
