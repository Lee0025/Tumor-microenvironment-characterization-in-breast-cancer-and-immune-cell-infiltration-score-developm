###主成份分析
exp <- DEG_gene_exp$DEG_expression
s <- w$`Correlative Analysis DEG_to ImmuneScore`[!w$`Correlative Analysis DEG_to ImmuneScore`$Correlation=='Not-Positive',]$Gene
s1 <- w$`Correlative Analysis DEG_to ImmuneScore`[w$`Correlative Analysis DEG_to ImmuneScore`$Correlation=='Positive',]$Gene
s2 <- w$`Correlative Analysis DEG_to ImmuneScore`[w$`Correlative Analysis DEG_to ImmuneScore`$Correlation=='Negative',]$Gene
exp1 <- exp[s1,]
exp2 <- exp[s2,]
exp1 <- apply(exp1,2,as.numeric)
rownames(exp1) <- s1
exp2 <- apply(exp2,2,as.numeric)
rownames(exp2) <- s2
exp1 <- t(exp1)
exp2 <- t(exp2)
exp1 <- as.data.frame(exp1)
exp2 <- as.data.frame(exp2)
###PCA
pca1 <- prcomp(exp1,scale. = TRUE)
biplot(pca1)
pca2 <- prcomp(exp2,scale. = TRUE)
biplot(pca2)
###scree plot
screeplot(pca1,npcs = 30,type = 'line')
screeplot(pca2,npcs = 17,type = 'line')
factoextra::fviz_eig(pca1,width=0.8,addlabels = TRUE)
factoextra::fviz_eig(pca2,width=0.8,addlabels = TRUE)
###scores
x1 <- as.data.frame(pca1$x)
x2 <- as.data.frame(pca2$x)
x1$Sum <- apply(x1,1,sum)
x2$Sum <- apply(x2,1,sum)
Score <- x1$Sum-x2$Sum
Score <- as.data.frame(Score)
rownames(Score) <- rownames(x2)
#####生存数据最佳ICI_Score截断值选择
library(pROC)
# label: 金标准，0 1 变量
# pred: 模型预测值，连续变量

cal_metrics <- function(label, pred){
  roc.p=pROC::roc(label, pred)
  if (roc.p$auc>0.5){
    cutoff=roc.p$thresholds[which.max(roc.p$sensitivities+roc.p$specificities)]
    sensitivity=roc.p$sensitivities[which.max(roc.p$sensitivities+roc.p$specificities)]
    specificity=roc.p$specificities[which.max(roc.p$sensitivities+roc.p$specificities)]
    df=data.frame(type='positive classification',
                  auc=round(roc.p$auc,3),cutoff=cutoff,
                  sensitivity=sensitivity,specificity=specificity)
    return(df)
  }
  else{
    cutoff=roc.p$thresholds[which.min(roc.p$sensitivities+roc.p$specificities)]
    sensitivity=roc.p$sensitivities[which.min(roc.p$sensitivities+roc.p$specificities)]
    specificity=roc.p$specificities[which.min(roc.p$sensitivities+roc.p$specificities)]
    df=data.frame(type='negative classification',
                  auc=1-round(roc.p$auc,3),cutoff=cutoff,
                  sensitivity=1-sensitivity,specificity=1-specificity)
    return(df)
  }
}
###
s <- OV_data$`Summary data`[,c(1,10,11)]
rownames(s) <- s$ID
ID <- rownames(Score)
Score <- cbind(ID,Score)
Score <- Score[,-1]
s <- s[ID,]
cuttoff_select <- merge(Score,s,by.x = 'ID',by.y = 'ID')
cuttoff_select$Statu <- ifelse(cuttoff_select$Statu=='Alive',0,1)
cal_metrics(cuttoff_select$Statu,cuttoff_select$Score)
cuttoff_select$Day <- as.numeric(cuttoff_select$Day)
cuttoff_select$Cut <- ifelse(cuttoff_select$Score > -0.3755676,'A','B')
save(cuttoff_select,file = 'cutoff_select.Rdata')
cuttoff_select$ICI_Score <- ifelse(cuttoff_select$Cut=='A','High','Low')
library(survival)
fit <- survfit(Surv(Day,statu)~ICI_score,KM)
library(survminer)
ggsurvplot(fit,risk.table=TRUE,#生存统计统计表
           
           conf.int=FALSE,#添加置信区间带
           
           palette = c("red","green"),#颜色设置
           
           pval=TRUE,#log-rank检验
           
           pval.method=TRUE,#添加检验text
           surv.median.line = "hv")
###
s <- OV_data$`Summary data`[,c(2,41)]
s$ImmuneScore <- as.numeric(s$ImmuneScore)
p <- ggplot(s,aes(x=ICI_Cluster,y=ImmuneScore))+geom_boxplot(aes(fill=ICI_Cluster))+theme_classic()+theme(legend.position = 'top')+
     scale_fill_manual(values = c('red','blue','green'))
p
###
D <- OV_data$`Summary data`[,c(10:41)]
D$Statu <- ifelse(D$Statu=='Alive',0,1) 
D <- apply(D,2,as.numeric)
fun <- function(x){
  a <- 10*x
  return(a)
}
D[,3:30] <- apply(D[,3:30],2,fun)
D <- as.data.frame(D)
x <- colnames(D)[3:ncol(D)]
library(survival)
outTab <- data.frame()
pre_survive_data4 <- D
for(i in x){
  expr=pre_survive_data4[,i]
  cox=coxph(Surv(Day,Statu)~ expr,pre_survive_data4)
  coxsummary=summary(cox)
  outTab=rbind(outTab,cbind(gene=i,HR=round(coxsummary$coefficients[,'exp(coef)'],2),
                            ###HR表示风险比
                            z=round(coxsummary$coefficients[,'z'],2),##z值
                            '95%CI'=paste(round(coxsummary$conf.int[,3],2),
                                          round(coxsummary$conf.int[,4],2),sep = '-'),##95%CI
                            pvalue=round(coxsummary$coefficients[,'Pr(>|z|)'],8)))
  ##p值
}
outTab$pvalue <- as.numeric(outTab$pvalue)
outTab <- outTab[order(outTab$pvalue),]
colnames(outTab)[1] <- 'Cell Type'
###
s <- OV_data$`Summary data`[,c(1,2,10)]
rownames(s) <- s$ID
s <- s[cuttoff_select$ID,]
s <- merge(s,cuttoff_select[,c(1,6)],by.x = 'ID',by.y = 'ID')
library(ggalluvial)
ggplot(data = s,
       aes(axis1 = ICI_Score, axis2 = ICI_Cluster, axis3 = Statu)) +
       scale_x_discrete(limits = c("ICI_Score", "ICI_Cluster", "Statu"), expand = c(.1, .05)) +
       geom_alluvium(aes(fill = Statu),show.legend = TRUE)+
       geom_stratum()+geom_text(stat = "stratum", label.strata = TRUE)+
       theme_minimal()+scale_fill_manual(values = c('green','red'))+
       theme(axis.text.x = element_text(colour="grey20",size=10,face="bold"))
 
# 查看教程
vignette(topic = "ggalluvial", package = "ggalluvial")
######Immune Check-point: CD274, CTLA4, HAVCR2,IDO1, LAG3, and PDCD1 as immune-checkpoint-relevant signatures
##CD8A, CXCL10, CXCL9, GZMA, GZMB, IFNG, PRF1,TBX2, and TNF as immune-activity-related signatures.
exp <- PcRNA_ICI_PD1$exp
exp <- as.data.frame(exp)
a <- c('CD274','CTLA4','HAVCR2','IDO1','LAG3','PDCD1','CD8A','CXCL10','CXCL9','GZMA','GZMB','IFNG','PRF1','TBX2','TNF')
expa <- exp[a,]
a <- rownames(expa)
expa <- apply(expa,2,as.numeric)
expa <- as.data.frame(expa)
rownames(expa) <- a
expa <- t(expa)
ID <- rownames(expa)
expa <- cbind(ID,expa)
expa <- as.data.frame(expa)
s2 <- s[,c(1,4)]
expa <- merge(s2,expa,by.x = 'ID',by.y = 'ID')
expa[,3:17] <- apply(expa[,3:17],2,as.numeric)
fun2 <- function(x){
  z <- log2(x+1)
  return(z)
}
fun3 <- function(x){
        z <- log2(x+1)
        return(z)
}
expa[,3:17] <- apply(expa[,3:17],2,fun3)
expa[,3:17] <- apply(expa[,3:17],2,scale)
expa <- expa[,-1]
library(reshape2)
library(ggsignif)

expa <- melt(expa)
p <- ggplot(expa,aes(variable,value,fill=ICI_Score))+geom_boxplot(outlier.size = 0.1)+theme_classic()+
     theme(legend.position = 'top')+theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.ticks.x = element_line(colour = "black", size = 0.5))+scale_fill_manual(values = c('red','green'))+
     labs(x='',y='Relative Normalized Expression Value')+theme(axis.text.x = element_text(colour="grey20",size=10,face="bold"))+
     geom_signif(comparisons = list(c('High','Low')),inherit.aes = TRUE,test = "wilcox.test")
     
     
p
####
com <- vector()
expb <- exp[a,]
expb <- t(expb)
ID <- rownames(expb)
expb <- cbind(ID,expb)
expb <- expb[s2$ID,]
expb <- merge(s2,expb,by.x = 'ID',by.y = 'ID')
expb <- as.data.frame(expb)
expb[,3:ncol(expb)]<- apply(expb[,3:ncol(expb)],2,as.numeric)
x <- colnames(expb)[3:ncol(expb)]
for(i in x){
  a <- expb[expb$ICI_Score=='Low',][,i]
  b <- expb[expb$ICI_Score=='High',][,i]
  d <- kruskal.test(list(a,b))$p.value
  if(length(com)==0){
    com=d
  }else{
    com=c(com,d)
  }
}
com <- as.data.frame(com)
rownames(com) <- x
colnames(com) <- 'pvalue'
com <- cbind(x,com)
colnames(com) <- c('gene','pvalue')
getwd()
setwd("/Users/llls2012163.com/TCGA/TCGA_BRCA")
com_data <- list(expb,com,annotation='First Exp; Second:Kruskal_result;Second:Annotation')
names(com_data) <- c('expression','Kruskal_testResult','Annotation')
save(com_data,file = 'Gene_expression_KW.test.Rdata')
setwd("/Users/llls2012163.com/R-scrpit/ImmuneFilither")
####
gene <- w$`Correlative Analysis DEG_to ImmuneScore`[!w$`Correlative Analysis DEG_to ImmuneScore`$Correlation=='Not-Positive',]$Gene
exp <- DEG_gene_exp$DEG_expression[gene,]
exp <- apply(exp,2,as.numeric)
rownames(exp) <- gene
exp <- PcRNA_ICI_PD1$exp
exp <- as.data.frame(exp)
s <- cuttoff_select[,c(1,5)]
exp <- t(exp)
exp <- as.data.frame(exp)
ID <- rownames(exp)
exp <- cbind(ID,exp)
exp <- merge(s,exp,by.x = 'ID',by.y = 'ID')
exp <- exp[order(exp$Cut),]
rownames(exp) <- exp$ID
exp <- exp[,-2]
exp <- exp[,-1]
exp <- t(exp)
exp <- as.data.frame(exp)
ID <- rownames(exp)
exp <- apply(exp,2,as.numeric)
rownames(exp) <- ID
library(edgeR)
group <- c(rep(1,553),rep(2,427))
y <- DGEList(counts = exp1,group = group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.size=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
et <- as.data.frame(et)
library(GSEABase)
colnames(clinical) <- clinical[1,]
clinical <- clinical[-1,]
clinical <- clinical[!duplicated(clinical$case_submitter_id),]
clin <- clinical[26:nrow(clinical),]
clin <- clin[,c(2,154)]
ID <- clin$case_submitter_id
library(limma)
ID <- strsplit2(ID,'-')
x <- paste(ID[,1],'_',ID[,2],sep = '')
ID <- paste(x,'_',ID[,3],sep = '')
clin <- cbind(ID,clin)
rownames(clin) <- ID
s <- OV_data$`Summary data`
x <- s$ID
clin <- clin[x,]
clin <- clin[,-2]
s <- merge(clin,s,by.x='ID',by.y='ID')
OV_data$`Summary data` <- s
exp <- PcRNA_ICI_PD1$exp
exp <- t(exp)
exp <- as.data.frame(exp)
ID <- rownames(exp)
exp <- cbind(ID,exp)
a <- cuttoff_select[,c(1,5)]
exp <- merge(a,exp,by.x='ID',by.y='ID')
exp <- exp[order(exp$Cut),]
rownames(exp) <- exp$ID
exp <- exp[,c(-1,-2)]
exp <- t(exp)
exp <- as.data.frame(exp)
ID <- rownames(exp)
exp <- apply(exp,2,as.numeric)
rownames(exp) <- ID
exp <- as.data.frame(exp)
exp1 <- cbind(exp[,428:ncol(exp)],exp[,1:427])
ID <- rownames(et)
et <- cbind(ID,et)
b <- et[,1:2]
b <- b[order(b$logFC),]
library(clusterProfiler)
library(org.Hs.eg.db)
ID <- rownames(et)
c <- bitr(ID,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
b <- b[c$SYMBOL,]
SYMBOL <- rownames(b)
b <- cbind(SYMBOL,b)
d <- merge(c,b,by.x='SYMBOL',by.y='SYMBOL')
d1 <- d[order(d$logFC),]
d1 <- d1[!duplicated(d1$ENTREZID),]
ID <- d1$ENTREZID
rownames(d1) <- d1$ENTREZID
d1 <- d1[,c(-1)]
geneList <- d1$logFC
names(geneList) <- d1$ENTREZID
geneList=sort(geneList,decreasing = T)
kegmt<-read.gmt("c2.cp.kegg.v7.4.entrez.gmt") #读gmt文件
KEGG1<-GSEA(geneList = geneList,TERM2GENE = kegmt) #GSEA分析
######
KEGG<- gseKEGG(geneList = geneList, organism = "hsa")
library(ggplot2)
library(enrichplot)
dotplot(KEGG) #出点图 
dotplot(KEGG,color="pvalue")  #按p值出点图 
gseaplot2(KEGG, row.names(KEGG@result)[c(1,2,10,33,34,39,41,30,44,45)],pvalue_table = FALSE,ES_geom = 'line',subplots = 1:2)
save(KEGG,file = 'KEGG_ICI_Score.Rdata')
save(KEGG1,file = 'KEGG_GSEA.Rdata')
dev.off()
#######
s <- OV_data$`Summary data`[,c(1,2,11,12)]
rownames(s) <- s$ID
s1 <- cuttoff_select[,c(1,5)]
s <- s[s1$ID,]
s2 <- merge(s1,s,by.x='ID',by.y='ID')
s2 <- s2[,-1]
library(reshape2)
s2 <- melt(s2)
s2$Statu <- ifelse(s2$Statu=='Alive',0,1)
s2$Day <- as.numeric(s2$Day)
'High_Score+PT '<- s2[s2$Cut=='A' & s2$treatment_type=='PT',]
'High_Score+RT' <- s2[s2$Cut=='A' & s2$treatment_type=='RT',]
'Low_Score+PT' <- s2[s2$Cut=='B' & s2$treatment_type=='PT',]
'Low_Score+RT' <- s2[s2$Cut=='B' & s2$treatment_type=='RT',]
Group <- c(rep('High_Score+PT',208),rep('High_Score+RT',215),rep('Low_Score+PT',275),rep('Low_Score+RT',269))
s3 <- rbind(`High_Score+PT `,`High_Score+RT`,`Low_Score+PT`,`Low_Score+RT`)
s3 <- cbind(Group,s3)
Strategy <- s3
save(Strategy,file = 'StrategyChoice.Rdata')
#######
s4 <- OV_data$TMB[,c(1,4)]
s5 <- cuttoff_select[,c(1,5)]
s4 <- s4[s5$ID,]
s6 <- merge(s4,s5,by.x='ID',by.y='ID')
s6$total_perMB <- as.numeric(s6$total_perMB)
library(ggplot2)
p <- ggplot(s6,aes(x=Cut,y=total_perMB,color=Cut))+geom_bar(stat = '')+theme_classic2()+
     labs(x='',y='TMB')
p
a <- s6[s6$Cut=='A',][,2]
b <- s6[s6$Cut=='B',][,2]
kruskal.test(list(a,b))
######
s7 <- cuttoff_select[,c(1,2)]
s8 <- OV_data$`Summary data`[,c(1,3,45)]
rownames(s7) <- s7$ID
s7 <- s7[s8$ID,]
s9 <- merge(s7,s8,by.x='ID',by.y='ID')
s9 <- s9[!s9$total_perMB==91.44,]
s9 <- s9[!s9$total_perMB==74.2,]
model.lm<-lm(formula = Score ~ total_perMB, data = s9)
summary(model.lm)
p <- ggplot(s9,aes(x=Score,y=total_perMB))+geom_point(aes(color=ICI_Cluster))+
     stat_smooth(method='lm',formula = y~x,colour='black')+
     scale_color_manual(values = c('red','blue','green'))+theme_classic()+theme(legend.position = c(0.9,0.8))+
     labs(x='ICI_Score',y='Tumor Burden Mutation')
     
p
####
s <- OV_data$`Summary data`[,c(1,11,12,45)]
s$Statu <- ifelse(s$Statu=='Alive',0,1)
s$Day <- as.numeric(s$Day)
s$TMB <- ifelse(s$total_perMB > median(s$total_perMB),'High','Low')
s$TMB1 <- ifelse(s$total_perMB>0.37,'High','Low')
s <- s[,-5]
names(s$TMB1) <- 'TMB'
colnames(s) <- c('ID','Statu','Day','total_perMB','TMB')
s1 <- cuttoff_select[,-2]
s2 <- s[,-c(2:4)]
s3 <- merge(s1,s2,by.x = 'ID',by.y = 'ID')
'H_ICI Score+H_TMB' <- s3[s3$Cut=='A'&s3$TMB=='High',]
'H_ICI Score+L_TMB' <- s3[s3$Cut=='A'&s3$TMB=='Low',]
'L_ICI Score+H_TMB' <- s3[s3$Cut=='B'&s3$TMB=='High',]
'L_ICI Score+L_TMB' <- s3[s3$Cut=='B'&s3$TMB=='Low',]
Subtype <- c(rep('H_ICI Score+H_TMB',319),rep('H_ICI Score+L_TMB',104),rep('L_ICI Score+H_TMB',423),rep('L_ICI Score+L_TMB',121))
s4 <- rbind(`H_ICI Score+H_TMB`,`H_ICI Score+L_TMB`,`L_ICI Score+H_TMB`,`L_ICI Score+L_TMB`)
s4 <- cbind(Subtype,s4)
b <- cuttoff_select[,c(1,2)]
s5 <- merge(b,s,by.x = 'ID',by.y = 'ID')
#######
library(maftools)
id <- '/Users/llls2012163.com/TCGA/TCGA_BRCA/gdc_download_20210318_073529.382183/6c93f518-1956-4435-9806-37185266d248/TCGA_BRCA.maf.gz'
lam1 <- read.maf(maf =id)
oncoPrint(lam1)
oncoplot(lam1,top = 20,draw_titv = TRUE)
dev.off()
library(ComplexHeatmap)
oncoPrint(lam1@data)
write.mafSummary(maf = lam1,basename = 'laml')
d <- as.data.frame(lam1@data)
d1 <- d[,c(1,16,9)]
d1$Variant_Classification <- as.character(d1$Variant_Classification)
#####
d2 <- data.frame()
m <- d1$Hugo_Symbol[!duplicated(d1$Hugo_Symbol)]
n <- d1$Tumor_Sample_Barcode[!duplicated(d1$Tumor_Sample_Barcode)]
######
library(limma)
s <- d$Tumor_Sample_Barcode
s <- strsplit2(s,'-')[,1:3]
Tumor_Sample_Barcode <- paste(s[,1],'_',s[,2],sep='')
Tumor_Sample_Barcode <- paste(Tumor_Sample_Barcode,'_',s[,3],sep = '')
d$Tumor_Sample_Barcode <- Tumor_Sample_Barcode
rm(Tumor_Sample_Barcode)
n <- cuttoff_select$ID
d <- lam1@data 
d1 <- d[,c(1,16,9)]
m <- d1$Hugo_Symbol[!duplicated(d1$Hugo_Symbol)]
n <- d1$Tumor_Sample_Barcode[!duplicated(d1$Tumor_Sample_Barcode)]
m <- d1$Hugo_Symbol
n <- d1$Tumor_Sample_Barcode
t <- rep(NA,14873460)
d2 <- matrix(data = t,nrow = 15177,ncol = 980)
colnames(d2) <- n
rownames(d2) <- m
d2 <- as.data.frame(d2)
d2 <- as.matrix(d2)
d1$Variant_Classification <- as.character(d1$Variant_Classification)
library(reshape2)
d2  <- dcast(d1,Hugo_Symbol~Tumor_Sample_Barcode)
d1$Variant_Classification <- ifelse(d1$Variant_Classification=='Frame_Shift_Del',1,
                                    ifelse(d1$Variant_Classification=='Frame_Shift_Ins',2,
                                    ifelse(d1$Variant_Classification=='In_Frame_Del',3,
                                    ifelse(d1$Variant_Classification=='In_Frame_Ins',4,
                                    ifelse(d1$Variant_Classification=='Missense_Mutation',5,
                                    ifelse(d1$Variant_Classification=='Nonsense_Mutation',6,
                                    ifelse(d1$Variant_Classification=='Nonstop_Mutation',7,
                                    ifelse(d1$Variant_Classification=='Splice_Site',8,
                                    ifelse(d1$Variant_Classification=='Translation_Start_Site',9,0)))))))))
oncoplot(maf = lam1,top = 30)
Clin <- lam1@clinical.data
library(limma)
ID <- strsplit2(Clin$Tumor_Sample_Barcode,'-')[,1:3]    
x <- paste(ID[,1],'_',ID[,2],sep = '')
ID <- paste(x,'_',ID[,3],sep = '')
Clin <- cbind(Clin,ID)
ID1 <- Clinical_BRCA[,1]
ID1 <- strsplit2(ID1,'-')
ID2 <- paste(ID1[,1],'_',ID1[,2],sep = '')
ID <- paste(ID2,'_',ID1[,3],sep = '')
Clinical_BRCA <- cbind(ID,Clinical_BRCA)
Clinical_BRCA <- Clinical_BRCA[,-2]
Clinical_BRCA <- Clinical_BRCA[!duplicated(Clinical_BRCA[,1]),]
rownames(Clinical_BRCA) <- Clinical_BRCA[,1]
Clinical_BRCA <- Clinical_BRCA[Clin$ID,]
Clinical_BRCA <- merge(Clin,Clinical_BRCA,by.x = 'ID',by.y = 'ID')
s <- cuttoff_select[,c(1,5)]
s$Cut <- ifelse(s$Cut=='A','High','Low')
colnames(s) <- c('ID','ICI_Score')
Clinical_BRCA <- merge(Clinical_BRCA,s,by.x = 'ID',by.y = 'ID')
a <- `ICI_clusterCons.k=3.consensusClass`
colnames(a) <- c('ID','ICI_Cluster')
a$ICI_Cluster <- ifelse(a$ICI_Cluster==1,'Cluster_A',ifelse(a$ICI_Cluster==2,'Cluster_B','Cluster_C'))
Clinical_BRCA <- merge(Clinical_BRCA,a,by.x = 'ID',by.y = 'ID')
CLINICAL <- Clinical_BRCA[,-1]
CLINICAL$Stage <- ifelse(CLINICAL$Stage=='Stage III' | CLINICAL$Stage=='Stage IIIA' |CLINICAL$Stage=='Stage IIIB'|CLINICAL$Stage=='Stage IIIC'|CLINICAL$Stage=='Stage IV','Stage_III_IV','Stage_I_II')
CLINICAL <- CLINICAL[order(CLINICAL$ICI_Score,decreasing = TRUE),]
write.csv(CLINICAL,file = 'CLINICAL.clin')
id1 <- "/Users/llls2012163.com/TCGA/TCGA_BRCA/CLINICAL.clin"
lam1 <- read.maf(maf =id,clinicalData = id1)
lam1@clinical.data <- lam1@clinical.data[,-1]
oncoplot(maf = lam1,top = 30,clinicalFeatures = c('ICI_Score','ICI_Cluster','Gender','Statu','Stage'),sortByAnnotation = TRUE,annotationColor = annocor,showTitle = FALSE)
color <- list(ICI_Score=c(High='red',Low='green'),ICI_Cluster=c(Cluster_A='red',Cluter_B='blue',Cluter_C='green'),Gender=c(male='red',female='green'),Statu=c(Alive='#00BFFF',Dead='#DC143C'))
annocor <- list(ICI_Score=c(High='red',Low='green'),
                ICI_Cluster=c(Cluster_A='red',Cluster_B='blue',Cluster_C='green'),
                Statu=c(Alive='#00BFFF',Dead='#DC143C'),
                Stage=c(Stage_I_II='#00FF7F',Stage_III_IV='#FF4500'),
                Gender=c(female='#FFFF00',male='#FF0000'))
oncoplot(maf = lam1,top = 30,clinicalFeatures = c('ICI_Score','ICI_Cluster','Gender','Statu','Stage'),
                 sortByAnnotation = TRUE,annotationColor = annocor,showTitle = FALSE,
         legend_height = 3,anno_height = 2,legendFontSize = 1,fill = TRUE,drawRowBar = TRUE,drawBox = FALSE,colbar_pathway = TRUE)
###
s <- `PD1&PDL1`
s <- s[,-1]
fun <- function(x){
    a <- log2(x+1)
    return(a)
}
s[,c(3,4)] <- apply(s[,c(3,4)],2,fun)
library(ggplot2)
library(ggpubr)
com <- list(c('Cluster_A','Cluster_B'),c('Cluster_B','Cluster_C'),c('Cluster_A','Cluster_C'))
p <- ggplot(s)+geom_violin(aes(x=ICI_Cluster,y=PDCD1,fill=ICI_Cluster))+
     theme_classic()+theme(legend.position = 'none')+scale_fill_manual(values = c('red','blue','green'))+
     labs(x='',y='The Relative Expression of PD1')+geom_boxplot(aes(x=ICI_Cluster,y=PDCD1),outlier.size = 0.1)
p
p1 <- ggplot(s)+geom_boxplot(aes(x=ICI_Cluster,y=PDCD1,fill=ICI_Cluster),outlier.size = 0.1)+
      theme_classic()+theme(legend.position = 'none',axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 15,colour = 'black'),axis.title.y = element_text(size = 15,colour = 'black'))+scale_fill_manual(values = c('red','blue','green'))+
      labs(x='',y='The Relative Expression of PD1')
p1
p2 <- ggplot(s)+geom_boxplot(aes(x=ICI_Cluster,y=PDCD1LG2,fill=ICI_Cluster),outlier.size = 0.1)+
      theme_classic()+theme(legend.position = 'none',axis.text.x = element_text(size = 14,color="black",angle = 0),axis.text.y = element_text(size = 15,colour = 'black'),axis.title.y = element_text(size = 15,colour = 'black'))+scale_fill_manual(values = c('red','blue','green'))+
      labs(x='',y='The Relative Expression of PDL1')
p2
write.csv(s,file = 'PD1&PDL1.csv')
#########
s <- OV_data$`Summary data`[,c(3,13:42)]
s[,c(2:31)] <- apply(s[,c(2:31)],2,as.numeric)
fun <- function(x){
   a <- scale(x)
   return(a)
}
s[,c(2:31)] <- apply(s[,c(2:31)],2,fun)
library(reshape2)
s <- melt(s)
library(ggplot2)
p <- ggplot(s,aes(x=variable,y=value,fill=ICI_Cluster))+geom_boxplot(outlier.size = 0.1)+theme_classic2()+
     labs(x='',y='Relative Normalized Value')+theme(legend.position = 'top',axis.text.x = element_text(size = 9,color="black",angle =45,vjust = 1,hjust = 1,face = 'bold'))+
     scale_fill_manual(values = c('red','blue','green'))
     


p
######
library(clusterProfiler)
a <- clusterProfiler::read.gmt("genesets.gmt")
Glu_genset <- a$gene[!duplicated(a$gene)]
exp <- PcRNA_ICI_PD1$exp
q <- rownames(exp)
b <- intersect(q,Glu_genset)
exp <- exp[b,]
exp <- as.data.frame(exp)
b <- rownames(exp)
exp <- apply(exp,2,as.numeric)
rownames(exp) <- b
z <- prcomp(exp,scale. = TRUE)
w <- z$x
colnames(w) <- colnames(exp)
rm(z)
rm(w)
exp <- t(exp)
exp <- as.data.frame(exp)
exp <- apply(exp,1,fun)
z <- prcomp(exp,scale. = FALSE)
w <- z$x
w <- as.data.frame(w)
colnames(w) <- colnames(exp)
rownames(w) <- b
w <- t(w)
sum <- apply(w,1,sum)
sum <- as.data.frame(sum)
ID <- rownames(sum)
sum <- cbind(ID,sum)
h <- OV_data$`Summary data`[,c(1,3)]
u <- intersect(h$ID,sum$ID)
sum <- sum[u,]
sum <- merge(sum,h,by.x = 'ID',by.y = 'ID')
library(ggplot2)
p <- ggplot(sum,aes(x=ICI_Cluster,y=sum,fill=ICI_Cluster))+geom_boxplot()
p
sum$sum <- scale(sum$sum)
x1 <- sum1[sum1$ICI_Cluster=='Cluster_A',]$sum1
x2 <- sum1[sum1$ICI_Cluster=='Cluster_B',]$sum1
x3 <- sum1[sum1$ICI_Cluster=='Cluster_C',]$sum1
kruskal.test(list(x1,x2))
kruskal.test(list(x1,x3))
kruskal.test(list(x1,x2,x3))
kruskal.test(list(x2,x3))
exp <- PcRNA_ICI_PD1$exp[b,]
exp <- as.data.frame(exp)
exp <- apply(exp,2,as.numeric)
rownames(exp) <- b
fun1 <- function(x){
    a <- log2(x+1)
    return(a)
}
exp <- apply(exp,2,fun1)
rownames(exp) <- b
exp <- t(exp)
exp <- apply(exp,1,fun)
z <- prcomp(exp,scale. = FALSE)
z$x
rm(sum)
sum1 <- apply(z$x,1,sum)
sum1 <- as.data.frame(sum1)
ID <- rownames(sum1)
sum1 <- cbind(ID,sum1)
sum1 <- sum1[u,]
sum1 <- merge(sum1,h,by.x = 'ID',by.y = 'ID')
p <- ggplot(sum1,aes(x=ICI_Cluster,y=sum1,fill=ICI_Cluster))+geom_boxplot()+theme_classic()
p
colnames(sum1) <- c('ID','Glucose metabolism Score','ICI_Cluster')
sum1 <- sum1[order(sum1$ICI_Cluster),]
write.csv(sum1,'Glucose metabolism Score.csv')
h1 <- OV_data$`Summary data`[,c(1,42)]
sum1 <- merge(sum1,h1,by.x = 'ID',by.y = 'ID')
sum1$ImmuneScore <- as.numeric(sum1$ImmuneScore)
p1 <- ggplot(sum1,aes(x=`Glucose metabolism Score`,y=ImmuneScore))+geom_point()+geom_smooth()+theme_classic()
p1
h1 <- OV_data$`Summary data`[,c(1,13:41)]
sum2 <- sum1[,c(-1,-3,-33)]
StromalScore <- sum1$StromalScore
sum2 <- cbind(sum2[,c(1:2)],StromalScore,sum2[,c(3:ncol(sum2))])

sum2 <- apply(sum2,2,as.numeric)
sum2 <- as.data.frame(sum2)
sum2 <- apply(sum2,2,scale)
sum2 <- as.data.frame(sum2)
library(reshape2)
sum2 <- melt(sum2,id.vars = 'Glucose metabolism Score')
sum2$value <- as.numeric(sum2$value)
sum2$`Glucose metabolism Score` <- as.numeric(sum2$`Glucose metabolism Score`)
colnames(sum2) <- c('Glucose','Variable','Value')
library(ggpubr)
p2 <- ggplot(data = sum2,aes(x=Value,y=Glucose))+geom_point(size=0.05,color='red')+facet_wrap(~Variable,ncol = 6,scales = 'free_x',labeller = labeller(size=2))+
      labs(x='',y='')+geom_smooth(method=lm)+
      stat_cor(method = "pearson",label.x = -1, label.y = 2,digits = 2,label.sep = '\n',label.x.npc = 'right')+
      theme_light()
p2
#########
s <- OV_data$`Summary data`[,c(1,3,45)]
h <- merge(s,cuttoff_select,by.x = 'ID',by.y = 'ID')
library(ggplot2)
library(ggpubr)
p <- ggplot(h,aes(x=Score,y=total_perMB))+geom_point(size=1.2)+scale_color_manual(values = c('red','blue','green'))+
     theme_classic()+labs(x='ICI_Score',y='Tumor Mutation Burden')+
     theme(legend.position = c(0.8,0.8),axis.title.x = element_text(face = 'bold'),axis.title.y = element_text(face = 'bold',size = 14))+
     geom_smooth(method = 'lm')+
     stat_cor(method = "pearson",label.x = 20, label.y = 15,digits = 2,label.sep = '',label.x.npc = 'right')+
     geom_point(aes(color=ICI_Cluster))
p
##########
library(maftools)
s <- read.maf(maf = id,clinicalData = id1)
Lowbar <- s@clinical.data[s@clinical.data$ICI_Score=='Low',]$Tumor_Sample_Barcode
Highbar <- s@clinical.data[s@clinical.data$ICI_Score=='High',]$Tumor_Sample_Barcode
Lmb <- data.frame()
for(i in Lowbar){
   a <- s@data[s@data$Tumor_Sample_Barcode==i,]
   print(i)
   if(dim(Lmb)[1]==0){
     Lmb=a
   }else{
     Lmb <- rbind(Lmb,a)
   }
}
Hmb <- data.frame()
for(i in Highbar){
  b <- s@data[s@data$Tumor_Sample_Barcode==i,]
  print(i)
  if(dim(Lmb)[1]==0){
    Hmb=a
  }else{
    Hmb <- rbind(Lmb,a)
  }
}
id
varpL <- data.frame()
for(i in a){
  c <- s@variants.per.sample[s@variants.per.sample$Tumor_Sample_Barcode==i,]
  if(dim(varpL)==0){
    varpL=c
  }else{
    varpL =rbind(varpL,c)
  }
}
varpH <- data.frame()
for(i in b){
   d <- s@variants.per.sample[s@variants.per.sample$Tumor_Sample_Barcode==i,]
   if(dim(varpH)==0){
     varpH=d
   }else{
     varpH=rbind(varpH,d)
   }
}
setClass('maftools::MAF',
        slots=list(data='data.frame',variants.per.sample='data.frame',variant.type.summary='data.frame',
                   variant.classification.summary='data.frame',gene.summary='data.frame',summary='data.frame',maf.silent='data.frame',
                   clinical.data='data.frame'))
High <- new('maftools::MAF')
Low <- new('maftools::MAF')
High@data=Hmb
Low@data=Lmb
High@variants.per.sample=varpH
Low@variants.per.sample=varpL
vtyL <- data.frame()
for(i in a){
  u <- s@variant.type.summary[s@variant.type.summary$Tumor_Sample_Barcode==i,]
  if(dim(vtyL)[1]==0){
    vtyL=u
  }else{
    vtyL=rbind(vtyL,u)
  }
}
Low@variant.type.summary=vtyL
vtyH <- data.frame()
for(i in b){
   u <- s@variant.type.summary[s@variant.type.summary$Tumor_Sample_Barcode==i,]
   if(dim(vtyH)[1]==0){
     vtyH=u
   }else{
     vtyH=rbind(vtyH,u)
   }
}
High@variant.type.summary=vtyH
vch <- data.frame()
for(i in b){
  c<- s@variant.classification.summary[s@variant.classification.summary$Tumor_Sample_Barcode==i,]
  if(dim(vch)[1]==0){
    vch=c
  }else{
    vch=rbind(vch,c)
  }
}
High@variant.classification.summary <- vch
vcl <- data.frame()
for(i in a){
  c <- s@variant.classification.summary[s@variant.classification.summary$Tumor_Sample_Barcode==i,]
  if(dim(vcl)==0){
    vcl=c
  }else{
    vcl=rbind(vcl,c)
  }
}
Low@variant.classification.summary <- vcl
oncoplot(maf = High,top = 10)
msh <- data.frame()
for(i in b){
  c <- s@maf.silent[s@maf.silent$Tumor_Sample_Barcode==i,]
  if(dim(msh)[1]==0){
    msh=c
  }else{
    msh=rbind(msh,c)
  }
}
High@maf.silent=msh
msL <- data.frame()
for(i in a){
  c <- s@maf.silent[s@maf.silent$Tumor_Sample_Barcode==i,]
  if(dim(msL)[1]==0){
    msL=c
  }else{
    msL=rbind(msL,c)
  }
}
Low@maf.silent=msL
oncoplot(maf = Low,top = 10)
######
library(wordclud2)
install.packages('wordcloud2')
data2 <- as.data.frame(table(s@data$Hugo_Symbol))
da2 <- subset(data2,Freq >=3)
wordcloud2(da2)
library(dplyr)
library(tidyverse)
library(reshape)
library(reshape2)
####
library(maftools)
s <- read.maf(maf = id,clinicalData = id1)
Lowbar <- subset(s@clinical.data,ICI_Score=='Low')$Tumor_Sample_Barcode
Highbar<- s@clinical.data[s@clinical.data$ICI_Score=='High',]$Tumor_Sample_Barcode
Lowmaf <- subsetMaf(maf = s,clinQuery = "ICI_Score %in% 'Low'")######high light code
Highmaf <- subsetMaf(maf = s,clinQuery = "ICI_Score %in% 'High'")#######high light code
coOncoplot(m1=Lowmaf,m2=Highmaf,m1Name = 'Low',m2Name = 'High')
oncoplot(Lowmaf)
oncoplot(Highmaf)
oncoplot(maf = Lowmaf,top = 30,clinicalFeatures = c('ICI_Score','ICI_Cluster','Gender','Statu','Stage'),
         sortByAnnotation = TRUE,annotationColor = annocor,showTitle = FALSE,
         legend_height = 3,anno_height = 2,legendFontSize = 1,fill = TRUE,drawRowBar = TRUE,drawBox = FALSE,colbar_pathway = TRUE)
oncoplot(maf = Highmaf,top = 30,clinicalFeatures = c('ICI_Score','ICI_Cluster','Gender','Statu','Stage'),
         sortByAnnotation = TRUE,annotationColor = annocor,showTitle = FALSE,
         legend_height = 3,anno_height = 2,legendFontSize = 1,fill = TRUE,drawRowBar = TRUE,drawBox = FALSE,colbar_pathway = TRUE)
lollipopPlot2(m1=Lowmaf,m2=Highmaf,m1_name = 'Low',m2_name = 'High',gene = 'TP53',AACol1 = 'HGVSp_Short',AACol2 = 'HGVSp_Short')
com <- mafCompare(m1=Lowmaf,m2=Highmaf,m1Name = 'Low',m2Name = 'High',minMut = 5)
forestPlot(mafCompareRes = com, pVal = 0.1, color = c('royalblue', 'maroon'), geneFontSize = 0.8)
write.csv(com$results[com$results$pval<0.05,],file = 'TMB_compare.csv')
#############Mtabolism
s <- OV_data$`Summary data`[,c(1,13:40)]
h <- merge(Glucose.metabolism.Score[,-1],s,by.x = 'ID',by.y = 'ID')
h1 <- h[,c(-1,-3)]
h1 <- apply(h1,2,as.numeric)
h1 <- as.data.frame(h1)
library(reshape2)
h2 <- melt(h1,id.vars = 'Glucose.metabolism.Score')
library(ggplot2)
library(ggpubr)
p <- ggplot(data = h2,aes(x=value,y=Glucose.metabolism.Score))+facet_wrap(~variable,ncol = 5,scales = 'free_x')+geom_point(size=0.05,color='red')+
     geom_smooth(method = 'lm')+scale_fill_manual(values = 'red')+stat_cor(method = "pearson",label.x = -0.1, label.y = 2,digits = 2,label.sep = '\n',label.x.npc = 'right')
p
######Glucose
w <- OV_data$`Summary data`[,c(1,2,5:12)]
z <- merge(Glucose.metabolism.Score,w,by.x = 'ID',by.y = 'ID')
z$Day <- as.numeric(z$Day)
z$Statu <- ifelse(z$Statu=='Alive',0,1)
cox=coxph(Surv(Day,Statu)~ Glucose.metabolism.Score,z)
summ <- summary(cox)
outTab=rbind(outTab,cbind(gene=i,HR=round(coxsummary$coefficients[,'exp(coef)'],2),
                          ###HR表示风险比
                          z=round(coxsummary$coefficients[,'z'],2),##z值
                          '95%CI'=paste(round(coxsummary$conf.int[,3],2),
                                        round(coxsummary$conf.int[,4],2),sep = '-'),##95%CI
                          pvalue=round(coxsummary$coefficients[,'Pr(>|z|)'],2)))
cal_metrics(z$Statu,z$Glucose.metabolism.Score)
z$Glucose.Score <- ifelse(z$Glucose.metabolism.Score> -2.518711,'High','Low')
library(survival)
library(survminer)
fit <- survfit(Surv(Day,Statu)~Glucose.Score,z)
ggsurvplot(fit,risk.table=TRUE,#生存统计统计表
           
           conf.int=FALSE,#添加置信区间带
           
           palette = c("red","green","red"),#颜色设置
           
           pval=TRUE,#log-rank检验
           
           pval.method=TRUE,#添加检验text
           surv.median.line = "hv")
##########
exp <- exp_seq[,c(5,8,9)]
library(tidyverse)
library(reshape2)
exp <- dcast(	
  data = exp,	
 gene_id~ submitted_sample_id,	
 value.var = 'raw_read_count',
 fun.aggregate = mean,
 na.rm=TRUE
)	
#############
library(maftools)
id <- '/Users/llls2012163.com/TCGA/TCGA_BRCA/gdc_download_20210318_073529.382183/6c93f518-1956-4435-9806-37185266d248/TCGA_BRCA.maf.gz'
id1 <- "/Users/llls2012163.com/TCGA/TCGA_BRCA/CLINICAL.clin"
tmb <- read.maf(maf = id,clinicalData = id1)
tmb@clinical.data <- tmb@clinical.data[,-1]
Lowmaf <- subsetMaf(maf = tmb,clinQuery = "ICI_score %in% 'Low'")######high light code
Highmaf <- subsetMaf(maf = tmb,clinQuery = "ICI_score %in% 'High'")#######high light code
tmb@clinical.data$Age <- as.numeric(tmb@clinical.data$Age)
s <- tmb@clinical.data$Age
Age=c(colorRampPalette(colors=c('#008000','#ADFF2F',"#FFA500"))(range(s)/4),colorRampPalette(color=c('#ADFF2F',"#FFA500","#FF0000"))(range(s)/4))
annocor <- list(ICI_Cluster=c(Cluster_A='red',Cluster_B='blue',Cluster_C='green'),
                Statu=c(Alive='#00BFFF',Dead='#DC143C'),
                Stage=c(Stage_I_II='#00FF7F',Stage_III_IV='#FF4500'),
                Gender=c(female='#FFFF00',male='#FF0000'),
                Age=Age,
                ICI_score=c(High='red',Low='green'))
oncoplot(maf = list(Lowmaf,Highmaf),top = 30,clinicalFeatures = c('ICI_score','ICI_Cluster','Gender','Statu','Stage'),
         sortByAnnotation = TRUE,annotationColor = annocor,showTitle = TRUE,
         legend_height = 3,anno_height = 2,legendFontSize = 1,fill = TRUE,drawRowBar = TRUE,drawBox = FALSE,colbar_pathway = TRUE)
oncoplot(maf = Highmaf,top = 30,clinicalFeatures = c('ICI_score','ICI_Cluster','Gender','Statu','Stage'),
         sortByAnnotation = TRUE,annotationColor = annocor,showTitle = TRUE,
         legend_height = 3,anno_height = 2,legendFontSize = 1,fill = TRUE,drawRowBar = TRUE,drawBox = FALSE,colbar_pathway = TRUE)
############ICGC-KR icohort
exp_KR <- exp
exp_TCGA <- PcRNA_ICI_PD1$exp
exp_TCGA <- as.data.frame(exp_TCGA)
exp_KR <- exp_KR[rownames(exp_TCGA),]
BiocManager::install('sva')
a
########Normal Sample Immune Cell iflater condition
library(GSVA)
exp <- pcExpdata$Exp[,1:113]
x <- rownames(exp)
exp <- apply(exp,2,as.numeric)
exp <- as.data.frame(exp)
exp <- as.matrix(exp)
rownames(exp) <- x
NormalsampleICI <- gsva(expr = exp,gset.idx.list = z,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
tumorsampleICI <- OV_data$SSgsva
tumorsampleICI <- t(tumorsampleICI)
tumorsampleICI <- as.data.frame(tumorsampleICI)
tumorsampleICI <- tumorsampleICI[-1,]
id <- rownames(tumorsampleICI)
tumorsampleICI <- apply(tumorsampleICI,2,as.numeric)
rownames(tumorsampleICI) <- id
s <- cbind(NormalsampleICI,tumorsampleICI)
colanno <- data.frame(Group=c(rep('Normal',113),rep('Tumor',981)))
rownames(colanno) <- colnames(s)
col <- list(Group=c(Normal='#FFFF00',Cluster_A='red',Cluster_B='blue',Cluster_C='green'))
w <- OV_data$`Summary data`[,c(1,3,13:40)]
w <- w[order(w$ICI_Cluster),]
w <- t(w)
colnames(w) <- w[1,]
w <- w[c(-1,-2),]
s <- cbind(NormalsampleICI,w)
s <- as.data.frame(s)
s <- apply(s,2,as.numeric)
rownames(s) <- id
s <- as.data.frame(s)
s2 <- cbind(NormalsampleICI,w)
s2 <- apply(s2,2,as.numeric)
rownames(s2) <- id
anno <- data.frame(c(rep('Normal',113),rep('Tumor',968)))
rownames(anno) <- colnames(s2)
colnames(anno) <- 'Group'
color <- list(Group=c(Normal='green',Tumor='red'))
colanno <- data.frame(c(rep('Normal',113),rep('Cluster_A',465),rep('Cluster_B',425),rep('Cluster_C',78)))
colnames(colanno) <- 'Group'
rownames(colanno) <- c(colnames(NormalsampleICI),h$ID)
library(pheatmap)
bk <- c(seq(-3,0,by=0.01),seq(0.01,3,by=0.01))
pheatmap(s,scale='row',show_colnames = FALSE,show_rownames = TRUE,cluster_rows = FALSE,cluster_cols = FALSE,
         annotation_col = anno,annotation_colors = color,
         color=c(colorRampPalette(colors=c("blue","white"))(length(bk)/2),colorRampPalette(color=c("white","red"))(length(bk)/2)),
         legend_breaks=seq(-3,3,1),breaks = bk,fontsize_row=10)
#######
s <- t(s)
ID <- rownames(s)
s <- cbind(ID,s)
s <- as.data.frame(s)
Group <- c(rep('Normal',113),rep('Tumor',968))
s <- cbind(Group,s)
s <- s[,-2]
library(reshape2)
s <- melt(s,id.vars = 'Group')
library(ggplot2)
s$value <- as.numeric(s$value)
p <- ggplot(s,aes(x=variable,y=value,fill=Group))+geom_boxplot(outlier.size = 0.05)+theme_classic()+theme(legend.position = 'top')+
     scale_fill_manual(values = c('green','red'))+theme(axis.text.x = element_text(angle = 45,hjust = 1,face = 'bold'))+
     labs(x='',y='')

p
########Comparision between Normal and Tumor Group ICI
result <- vector()
for(i in id){
  a <- kruskal.test(list(s[i,1:113],s[i,114:ncol(s)]))$p.value
  if(length(result)==0){
    result=a
  }else{
    result=c(result,a)
  }
}
result <- cbind(id,result)
colnames(result) <- c('Cell Type','Pvalue')
save(result,file = 'NormalvsTumor(ICI).Rdata')
######
s <- tmb@data
library(reshape2)
s1 <- s[,c(1,9,16)]
s1 <- dcast(s1,Hugo_Symbol~Tumor_Sample_Barcode+Variant_Classification)
s1$Variant_Classification <- as.character(s1$Variant_Classification)
s1 <- dcast(s1,Hugo_Symbol~Tumor_Sample_Barcode,value.var = 'Variant_Classification')
library(tidyr)
library(sva)
#########tidy the BRCA-KR Cohort
exp <- exp_seq[,c(1,8,9)]
library(reshape2)
exp <- dcast(exp,gene_id~icgc_donor_id,fun.aggregate = mean,value.var = "normalized_read_count")
rownames(exp) <- exp$gene_id
Clinical_BRCA_KR <- donor[,c(1,5,6,7,9,14,17)]
exp <- exp[rownames(PcRNA_ICI_PD1$exp),]
exp_BRCA_KR <- exp
BRCA_KR <- list(exp_BRCA_KR,Clinical_BRCA_KR)
names(BRCA_KR) <- c('exp_BRCA_KR','Clinical_BRCA_KR')
BRCA_KR$exp_BRCA_KR <- exp
save(BRCA_KR,file = 'BRCA_KR.Rdata')
exp <- BRCA_KR$exp_BRCA_KR
exp <- exp[,-1]
exp1 <- exp[ID,]
is.na(exp1)
exp <- exp[,-1]
scores <- as.data.frame(scores)
scores$Tumor_Purity <- cos(0.6049872018+0.0001467884*scores$ESTIMATEScore)
library(GSVA)
library(clusterProfiler)
exp <- as.matrix(exp)
ssgsva <- gsva(expr = exp,gset.idx.list = z,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
ssgsva <- as.data.frame(ssgsva)
BRCA_KR_OV <- list(BRCA_KR,scores,ssgsva)
names(BRCA_KR_OV) <- c('KR_exp&Clinical','estimateScore','ICI_abundance')
save(BRCA_KR_OV,file = 'BRCA_KR_OV.Rdata')
ssgsva <- t(ssgsva)
q <- cbind(ssgsva,scores)[,c(-31,-32)]
s <- OV_data$`Summary data`[,c(1,13:42)]
ID <- rownames(q)
q <- cbind(ID,q)
rownames(q) <- 1:50
w <- rbind(q,s)
w <- t(w)
w <- as.data.frame(w)
colnames(w) <- w[1,]
w <- w[-1,]
w <- apply(w,2,as.numeric)
w <- as.data.frame(w)

rownames(w) <- colnames(s)[-1]
w1 <- apply(t(w),2,scale)
rownames(w1) <- colnames(w)
w1 <- t(w1)
w1 <- as.matrix(w1)
w1 <- as.data.frame(w1)
w1 <- as.matrix(w1)
class(w1)
q[,2:ncol(q)] <- apply(q[,2:ncol(q)],2,scale)
s[,2:ncol(s)] <- apply(s[,2:ncol(s)],2,as.numeric)
s[,2:ncol(s)] <- apply(s[,2:ncol(s)],2,scale)
w2 <- rbind(q,s)
w2 <- t(w2)
colnames(w2) <- w2[1,]
w2 <- w2[-1,]
w2 <- apply(w2,2,as.numeric)
rownames(w2) <- rownames(w1)
h <- `TCGA&KRconsensus_Cluster.k=3.consensusClass`[1:50,]
h1 <- BRCA_KR_OV$`KR_exp&Clinical`$Clinical_BRCA_KR[,c(1,3,7)]
h4 <- merge(h1,h,by.x = 'icgc_donor_id',by.y = 'V1')
h4$ICI_Cluster <- ifelse(h4$V2==1,'Cluster_A',ifelse(h4$V2==2,'Cluster_B','Cluster_C'))
h4$donor_vital_status <- ifelse(h4$donor_vital_status=='Alive',0,1)
h3 <- rbind(BRCA_KR_OV$ICI_abundance,t(BRCA_KR_OV$estimateScore[,c(1,2)]))
ID <- rownames(h3)
h3 <- apply(h3,1,scale)
rownames(h3) <- h[1:50,1]
h3 <- t(h3)
eset <- getGEO('GSE138536',destdir = '.',
               AnnotGPL = T,
               getGPL =T) 
###DEG consensus cluster
exp <- DEG_gene_exp$DEG_expression
id <- rownames(exp)
exp1 <- t(exp)
exp1 <- as.data.frame(exp1)
exp1 <- apply(exp1,2,as.numeric)
exp1 <- apply(exp1,2,scale)
exp1 <- as.data.frame(exp1)
colnames(exp1) <- rownames(exp)
rownames(exp1) <- colnames(exp)
library(ConsensusClusterPlus)
exp1 <- as.matrix(exp1)
exp1 <- t(exp1)
p <- ConsensusClusterPlus(d=exp1,maxK = 8,reps =1000,
                          pFeature = 1,clusterAlg = 'hc',title = 'Gene_consensusBRCA',
                          innerLinkage = 'average',finalLinkage = 'average',distance = 'pearson',
                          ml=NULL,tmyPal = NULL,writeTable = TRUE,weightsItem = NULL,
                          weightsFeature = NULL,verbose = F,corUse = 'everything',plot = 'png')
s <- Gene_consensusBRCA_k_3_consensusClass
colnames(s) <- c('id','group')
a <- OV_data$`Summary data`[,c(1,11,12)]
colnames(s) <- c('ID','Gene_cluster')
s1 <- merge(a,s,by.x = 'ID',by.y = 'ID')
s1$statu <- ifelse(s1$Statu=='Alive',0,1)
s1$Day <- as.numeric(s1$Day)
s1$Gene_cluster1 <- ifelse(s1$Gene_cluster==1,"cluster_A",ifelse(s1$Gene_cluster==2,'cluster_B','cluster_C'))
#######compare immune checkpoint between low and high ici group
s <- ICI_Score$`KM data`[,c(1,6)]
expa <- merge(s,expa,by.x = 'ID',by.y = 'ID')
expa2 <- expa
expa2[,3:ncol(expa2)] <- apply(expa2[,3:ncol(expa2)],2,as.numeric)
FUN <- function(x){
      y=log2(x+1)
      return(y)
}
expa2[,3:ncol(expa2)] <- apply(expa2[,3:ncol(expa2)],2,FUN)
library(reshape2)
expa1 <- melt(expa1[,2:ncol(expa)],id.vars = 'ICI_score')
expa1$ICI_score <- factor(expa1$ICI_score,levels = c('Low','High'))
expa1$value <- as.numeric(expa1$value)
expa1$value1 <- scale(expa1$value)
library(ggplot2)
library(ggprism)
p <- ggplot(expa1,aes(x=variable,y=value,fill=ICI_score))+geom_boxplot(outlier.size = 0.02)+theme_prism(base_size = 10)+labs(x='',y='Relative expression value')+
     theme(legend.position = 'top',axis.text.x = element_text(size = 9,color="black",angle =45,vjust = 1,hjust = 1))+scale_fill_manual(values = c('#00FF00','#FF0000'))
p
outTab <- c()
id <- colnames(expa2)[3:ncol(expa2)]
for(i in id){
    x1 <- expa2[expa2$ICI_score=='Low',i]
    x2 <- expa2[expa2$ICI_score=='High',i]
    x3 <- wilcox.test(x1,x2)$p.value
    outTab <- if(length(outTab)==0){
      outTab <- x3
    }else{
      outTab <- c(outTab,x3)
    }
}
outTab <- cbind(id,outTab)
colnames(outTab) <- c('gene','pvalue')
#########治疗方案与ICI-score
s1 <- OV_data$`Summary data`[,c(1,2,11,12)]
s2 <- merge(s1,s,by.x = 'ID',by.y = 'ID')
a1 <- s2[s2$treatment_type=='PT' & s2$ICI_score=='High',]
a2 <- s2[s2$treatment_type=='PT' & s2$ICI_score=='Low',]
a3 <- s2[s2$treatment_type=='RT' & s2$ICI_score=='High',]
a4 <- s2[s2$treatment_type=='RT' & s2$ICI_score=='Low',]
a5 <- rbind(a1,a2,a3,a4)
a5$Group <- c(rep('Ct+High',136),rep('Ct+Low',347),rep('Rt+High',113),rep('Rt+Low',371))
a5$statu <- ifelse(a5$Statu=='Alive',0,1)
a5$Day <- as.numeric(a5$Day)
a6 <- a5[a5$Group=='Ct+High' | a5$Group=='Ct+Low',]
a7 <- a5[a5$Group=='Rt+High' | a5$Group=='Rt+Low',]
fig5 <- list(expa1,outTab,a5)
names(fig5) <- c('Expgene','wilcox','KMdata')
save(fig5,file = 'Fig5.Rdata')
#######TMB与ICI-score
s <- fig5$KMdata[,c(1,5)]
s1 <- OV_data$`Summary data`[,c(1,45)]
s2 <- merge(s,s1,by.x = 'ID',by.y = 'ID')
s2 <- s2[order(s2$ICI_score),]
write.csv(s2,file = 'TMB_ICIscore.csv')
id <- clin$Tumor_Sample_Barcode
library(limma)
id1 <- strsplit2(id,'-')[,1:3]
id2 <- paste(id1[,1],'_',id1[,2],sep = '')
id <- paste(id2,'_',id1[,3],sep = '')
clin <- cbind(id,clin)
clin <- clin[,-11]
clin <- merge(clin,s2,by.x = 'id',by.y = 'ID')
clin <- clin[,-1]
write.csv(clin,file = 'CLINICAL.clin')
s3 <- OV_data$`Summary data`[,c(1,3)]
s4 <- ICI_Score$`KM data`[,c(1,4)]
s5 <- merge(s3,s4,by.x = 'ID',by.y = 'ID')
s5 <- merge(s5,s2,by.x = 'ID',by.y = 'ID')
library(ggplot2)
install.packages('ggpmisc')
library(ggpmisc)
library(ggpubr)
library(ggprism)
p <- ggplot(s5,aes(x=Score1+1,y=log2(total_perMB+1)))+geom_point(aes(color=ICI_Cluster),size=1)+geom_smooth(method = 'lm',color='black')+
     scale_color_manual(values = c('red','blue','green'))+
     theme_prism(base_size = 12)+theme(legend.position = c(0.8,0.85))+
     labs(x='ICI_score',y='Normalized Tumor Mutation Burden')
     
     
     
     
p
####TMB
s5$TMB <- ifelse(s5$total_perMB > 0.37,'High','Low')
s6 <- OV_data$`Summary data`[,c(1,11,12)]
s5 <- merge(s5,s6,by.x = 'ID',by.y = 'ID')
s5$statu <- ifelse(s5$Statu=='Alive',0,1)
s5$Day <- as.numeric(s5$Day)
'High_ICI+High_TMB' <- s5[s5$ICI_score=='High' & s5$TMB=='High',]
'High_ICI+Low_TMB'  <- s5[s5$ICI_score=='High' & s5$TMB=='Low',]
'Low_ICI+Low_TMB'   <- s5[s5$ICI_score=='Low' & s5$TMB=='Low',]
'Low_ICI+High_TMB'  <- s5[s5$ICI_score=='Low' & s5$TMB=='High',]
s6 <- rbind(`High_ICI+Low_TMB`,`High_ICI+High_TMB`,`Low_ICI+High_TMB`,`Low_ICI+Low_TMB`)
s6$Group <- c(rep('H_ICI+L_TMB',58),rep('H_ICI+H_TMB',191),rep('L_ICI+H_TMB',551),rep('L_ICI+L_TMB',167))
save(s6,file = 'Fig6.Rdata')
library(maftools)
s7 <- mafCompare(m1=Lowmaf,m2=Highmaf,m1Name = 'Low',m2Name = 'High',minMut = 5)
data <- s7$results
data <- as.data.frame(data)
save(data,file = 'mutCom.Rdata')
c <- c('TP53','PIK3CA','TTN','CDH1','MUC16','MAP3K1','MUC4','KMT2C','NCOR1','NEB','PTEN','USH2A','HMCN1','SYNE1','FLG','MUC17','RYR2','RUNX1','SPTA1','FAT3',
       'RYR3','DMD','NF1','MAP2K4','GOLGA6L6','SYNE2','HUWE1','ARID1A','ZFHX4')
d <- c('TBX3','KMT2C','AHNAK','F8','VPS13C','MXRA5','PKHD1L1','CCDC168','DNAH10','PIK3R1','USH2A','ABCA13',
       'CACNA1E','TPR')
e <- c(c,d)
e <- e[!duplicated(e)]
rownames(data) <- data$Hugo_Symbol
f <- data[e,]
colnames(f)[1:4] <- c('Gene symbol','Low_ICI score','High_ICI score','Pvalue')
library(sjPlot)
f <- f[order(f$Pvalue),]
tab_df(f[,1:4],title = '')
s6 <- as.data.frame(s6)
z <- cor(as.numeric(s6$Score1),as.numeric(s6$total_perMB),method = 'pearson')
z <- cor.test(as.numeric(s6$Score1),as.numeric(s6$total_perMB),method = 'pearson')
#######Metabric_cohort
data_clinical_sample <- data_clinical_sample[,-2
                                             ]
clin <- merge(data_clinical_patient,data_clinical_sample,by.x = "PATIENT_ID",by.y = "PATIENT_ID" )
clin1 <- clin[,s]
MetaBric_CorhortClin <- list(clin,clin1)
names(MetaBric_CorhortClin) <- c('Clinical','clin')
save(MetaBric_CorhortClin,file = 'MetabricClin.Rdata')
rm(MetaBric_CorhortClin)
rownames(MetaBric_Corhort) <- MetaBric_Corhort$Hugo_Symbol
rownames(MetaBric_Corhort) <- MetaBric_Corhort$Group.1
expa <- MetaBric_Corhort[a,]
expb <- MetaBric_Corhort[b,]
expa <- na.omit(expa)
expb <- na.omit(expb)
expa <- expa[,-1]
expb <- expb[,-1]
expa <- t(expa)
expb <- t(expb)
pcaa <- prcomp(expa,scale. = TRUE,center = TRUE)
pcab <- prcomp(expb,scale. = TRUE,center = TRUE)
ISA <- as.data.frame(pcaa$x)
ISB <- as.data.frame(pcab$x)
ISA <- ISA$PC1
ISB <- ISB$PC1
score <- ISA-ISB
ID <- rownames(as.data.frame(pcaa$x))
score <- cbind(ID,score)
clin <- MetaBric_CorhortClin$clin[,c(1,6,8,9,10)]
rownames(clin) <- clin$PATIENT_ID
library(limma)
ID1 <- strsplit2(ID,'[.]')
ID <- paste(ID1[,1],'-',ID1[,2],sep = '')
score <- as.data.frame(score)
score$ID <- ID
clin <- clin[ID,]
KM <- merge(clin,score,by.x = 'PATIENT_ID',by.y = 'ID')
Statu <- KM$VITAL_STATUS
Statu <- strsplit2(Statu,' ')
Statu <- as.data.frame(Statu)
Statu$ID <- KM$PATIENT_ID
Statu <- Statu[!Statu$V4=='Causes',]
Statu <- Statu[,c(1,5)]
colnames(Statu) <- c('Statu','ID')
KM <- merge(KM,Statu,by.x = 'PATIENT_ID',by.y = 'ID')
KM$statu <- ifelse(KM$Statu=='Living',0,1)
KM$score <- as.numeric(KM$score)
cal_metrics(KM$statu,KM$score)
KM$ICI_score <- ifelse(KM$score > -12.12817,'Low','High')
KM$ICI_score <- ifelse(KM$score > median(KM$score),'High','Low')
KM <- KM[,-9]
KM$Day <- KM$OS_MONTHS*30
Metabric_ICI <- list(expa,expb,pcaa,pcab,ISA,ISB,score,KM)
names(Metabric_ICI) <- c('expa','expb','PCAA','PCAB','ISA','ISB','ICI_Score','KMdata')
save(Metabric_ICI,file = 'Metabric_ICI.Rdata')
KM <- ICI_Score$`KM data`
s <- Metabric_ICI$KMdata[,c(1,8,9)]
ici <- `ICI_consensus_metabric.k=3.consensusClass`
colnames(ici) <- c('ID','ICI_Cluster')
ID <- colnames(ssgsva)
library(limma)
ID <- strsplit2(ID,'[.]')
ID <- paste(ID[,1],'-',ID[,2],sep = '')
ici$ID <- ID
s <- merge(s,ici,by.x ="PATIENT_ID",by.y = 'ID')
s$ICI_Cluster <- ifelse(s$ICI_Cluster==1,'Cluster_A',
                        ifelse(s$ICI_Cluster==2,'Cluster_B','Cluster_C'))
colnames(ssgsva) <- ID
library(pheatmap)
ann <- cbind(s$ICI_Cluster)
ann <- as.data.frame(ann)
rownames(ann) <- s$PATIENT_ID
colnames(ann) <- 'ICI_Cluster'
ann$ID <- rownames(ann)
ann <- ann[order(ann$ICI_Cluster),]
ID <- rownames(ann)
ann <- ann[,-2]
ann <- as.data.frame(ann)
rownames(ann) <- ID
ssgsva2 <- ssgsva[,s$PATIENT_ID]
annocor <- list(ICI_Cluster=c(Cluster_A='red',Cluster_B='blue',Cluster_C='green'))
bk <- c(seq(-4,0,by=0.01),seq(0.01,4,by=0.01))
ID <- colnames(ssgsva2)
ssgsva2 <- t(ssgsva2)
ssgsva2 <- cbind(ID,ssgsva2)
ssgsva2 <- as.data.frame(ssgsva2)
ssgsva2 <- merge(s,ssgsva2,by.x = "PATIENT_ID",by.y = 'ID')
ssgsva2 <- ssgsva2[order(ssgsva2$ICI_Cluster),]
rownames(ssgsva2) <- ssgsva2$PATIENT_ID
ssgsva2 <- t(ssgsva2)
ssgsva2 <- as.data.frame(ssgsva2)
ssgsva2 <- ssgsva2[-c(1:4),]
ID <- rownames(ssgsva2)
ssgsva2 <- apply(ssgsva2,2,as.numeric)
rownames(ssgsva2) <- ID
pheatmap(ssgsva2,scale = 'row',show_colnames = FALSE,cluster_rows = FALSE,cluster_cols = FALSE,
         annotation_col = ann,annotation_colors = annocor,color=c(colorRampPalette(colors=c("blue","white"))(length(bk)/2),colorRampPalette(color=c("white","red"))(length(bk)/2)),
         legend_breaks=seq(-3,3,1),breaks = bk,fontsize_row=10
         )
wd <- "/Users/llls2012163.com/R-scrpit/ImmuneFilither/metabric"
setwd(wd)
clin <- MetaBric_CorhortClin$clin
library(limma)
id <- strsplit2(clin$VITAL_STATUS,' ')
id <- as.data.frame(id)
id$ID <- clin$PATIENT_ID
id <- id[!id$V4=='Causes',]
id <- id[,c(1,5)]
colnames(id) <- c('statu','ID')
clin <- merge(clin,id,by.x = 'PATIENT_ID',by.y = 'ID')
MetaBric_CorhortClin$clin1 <- clin
save(MetaBric_CorhortClin,file = '(MetaBric_CorhortClin.Rdata')
geneSymble <- MetaBric_Corhort$Group.1
MetaBric_Corhort <- as.data.frame(MetaBric_Corhort)
rownames(MetaBric_Corhort) <- MetaBric_Corhort$Group.1
MetaBric_Corhort <- MetaBric_Corhort[,-1]
cola <- strsplit2(colnames(scores),'[.]')
cola <- paste(cola[,1],'-',cola[,2],sep = '')
colnames(scores) <- cola
save(MetaBric_Corhort,file = 'metabricRNAseq.Rdata')
library(GSVA)
MetaBric_Corhort <- as.matrix(MetaBric_Corhort)
ssgsva <- gsva(expr = MetaBric_Corhort,gset.idx.list = z,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
scores <- as.data.frame(scores)
scores <- t(scores)
ssgsva <- t(ssgsva)
ssgsva <- as.data.frame(ssgsva)
scores <- as.data.frame(scores)
ID <- rownames(scores)
scores <- cbind(ID,scores)
ID <- rownames(ssgsva)
ssgsva <- cbind(ID, ssgsva)
ssI <- merge(ssgsva,scores,by.x = 'ID',by.y = "ID")
rownames(ssI) <- ssI$ID
ID <- intersect(clin$PATIENT_ID,ssI$ID)
ssI <- ssI[ID,]
save(ssI,file = 'metabricImEs.Rdata')
ssI <- ssI[,c(-1,-32)]
ssI <- t(ssI)
ss <- apply(ssI,2,scale)
ss <- as.data.frame(ss)
rownames(ss) <- rownames(ssI)
mads=apply(ssI,1,mad)
d <- ssI
d=d[rev(order(mads))[1:30],]
d = sweep(d,1, apply(d,1,median,na.rm=T))
s <- `ICI_consensus_metabric.k=3.consensusClass`
colnames(s) <- c('ID','ICI_Cluster')
s$ICI_Cluster <- ifelse(s$ICI_Cluster==1,'Cluster_A',
                        ifelse(s$ICI_Cluster==2,'Cluster_B',
                               'Cluster_C'))
f <- clin[,c(1,6,12)]
f <- merge(f,s,by.x = 'PATIENT_ID',by.y = 'ID')
f$Statu <- ifelse(f$statu=='Living',0,1)
f$Day <- 30*f$OS_MONTHS
#############################################################relationship between deg expression value and immune score
DEG <- deg$DEG
exp <- pcExpdata$Exp[,114:ncol(pcExpdata$Exp)]
exp <- exp[DEG,]
id <- colnames(exp)
library(limma)
id <- strsplit2(id,'-')[,1:3]
id1 <- paste(id[,1],'_',id[,2],sep = '')
id <- paste(id1,'_',id[,3],sep = '')
colnames(exp) <- id
barcode <- OV_data$`Summary data`$ID
id <- intersect(id,barcode)
exp <- exp[,id]
DEG <- rownames(exp)
ImmuneScore <- OV_data$estimate
ImmuneScore <- ImmuneScore[id,]
exp <- t(exp)
exp <- as.data.frame(exp)
barcode <- rownames(exp)
ImmuneScore <- ImmuneScore[barcode,]
ImmuneScore <- ImmuneScore[,-1]
ImmuneScore <- apply(ImmuneScore,2,as.numeric)
rownames(ImmuneScore) <- barcode
exp <- apply(exp,2,as.numeric)
rownames(exp) <- barcode
exp <- as.data.frame(exp)
ImmuneScore <- as.data.frame(ImmuneScore)
exp <- exp[,DEG]
outTab1 <- data.frame()
for(i in DEG){
  x <- exp[,i]
  y <- ImmuneScore$StromalScore
  z <- cor.test(x,y,method = 'pearson')
  m <- c(z$statistic,z$estimate,z$p.value)
  m <- t(as.data.frame(m))
  rownames(m) <- i
  if(dim(outTab)[1]== 0){
    outTab1 <- m
  }else{
      outTab1 <- rbind(outTab1,m)
    }
}
colnames(outTab1) <- c('t','cor','pvalue')
outTab1$relationship <- ifelse(outTab1$pvalue < 0.05 & outTab1$cor < -0.3,'Negative',
                              ifelse(outTab1$pvalue < 0.05 & outTab1$cor > 0.3,'Positive','Non'))
a <- intersect(rownames(outTab[outTab$relationship=='Positive',]),rownames(outTab1[outTab1$relationship=='Positive',]))
corr_deg_TCGA <- list(deg_immScore=outTab,deg_strScore=outTab1)
exp <- exp[,rownames(outTab[outTab$relationship=='Positive',])]
exp <- t(exp)
exp <- as.data.frame(exp)
colnames(s) <- c('ID','Gene_Cluster')
s$Gene_Cluster <- ifelse(s$Gene_Cluster==1,'Cluster_A',
                         ifelse(s$Gene_Cluster==2,'Cluster_B',
                                'Cluster_C'))
s1 <- OV_data$`Summary data`[,c(1,11,12)]
s3 <- merge(s1,s,by.x = 'ID',by.y = 'ID')
s2$Day <- as.numeric(s2$Day)
s2$Statu <- ifelse(s2$Statu=='Alive',0,1)
s2 <- s2[,-1]
GeneClusterTcga <- list(ImmCordegexp = exp,Gene_clusterOS =s3)
save(GeneClusterTcga,file = 'TCGA_GeneCluster.Rdata')
####################################relationship between gene_cluster and DEG value positively related to immune score
id <- s3[order(s3$Gene_Cluster),]$ID
exp <- exp[,id]
m <- c(rep(1,329),rep(2,605),rep(3,33))
m <- as.vector(m)
m <- as.numeric(m)
outTab3 <- data.frame()
n <- rownames(exp)
for(i in n){
  x <- as.numeric(exp[i,])
  y <- cor.test(m,x,method = 'pearson')
  z <- c(y$statistic,y$estimate,y$p.value)
  q <- t(as.data.frame(z))
  rownames(q) <- i
  if(dim(outTab3)[1]== 0){
    outTab3 <- q
  }else{
    outTab3 <- rbind(outTab3,q)
  }
  
}
colnames(outTab3) <- c('t','cor','pvalue')
outTab3 <- as.data.frame(outTab3)
###########################prognostic value with gene_cluster
s <- as.data.frame(`deg_consensus.k=4.consensusClass`)
colnames(s) <- c('ID','Gene_cluster')
which(s$Gene_cluster==4)
s <- s[-827,]
table(s$V2)
s$Gene_cluster <- ifelse(s$Gene_cluster==1,'Cluster_A',
                         ifelse(s$Gene_cluster==2,'Cluster_B',
                                'Cluster_C'))
s2 <- merge(s,s1,by.x = 'ID',by.y = 'ID')
s2$Day <- as.numeric(s2$Day)
s2$Statu <- ifelse(s2$Statu=='Alive',0,1)
s2$Gene_cluster <- ifelse(s2$Gene_cluster=='Cluster_A','A',
                          ifelse(s2$Gene_cluster=='Cluster_B','B','C'))
#########################
gene <- c(rownames(corr_deg_TCGA$deg_immScore[corr_deg_TCGA$deg_immScore$relationship=='Positive',]),
          rownames(corr_deg_TCGA$deg_strScore[!corr_deg_TCGA$deg_strScore$relationship=='Non',]))
gene <- gene[!duplicated(gene)]
exp1 <- exp1[,114:ncol(exp1)]
library(limma)
ID <- strsplit2(colnames(exp1),'-')[,1:3]
ID1 <- paste(ID[,1],'_',ID[,2],sep = '')
ID <- paste(ID1,'_',ID[,3],sep = '')
colnames(exp1) <- ID
ID <- intersect(OV_data$`Summary data`$ID,ID)
exp1 <- exp1[,ID]
save(exp1,file = 'TCGApatientsexp.Rdata')
exp1 <- exp1[gene,]
gene <- rownames(exp1)
exp1 <- apply(exp1,2,as.numeric)
rownames(exp1) <- gene
exp1 <- as.data.frame(exp1)
d1 <- d[,-which(s$V2==3)]
s <- as.data.frame(`deg_consensus.k=3.consensusClass`)
table(s$V2)
Gene_ClusterLogRank <- list(data=s2,fit=fit,plot=p)
save(Gene_ClusterLogRank,file = 'Gene_ClusterLogRank.Rdata')
save(corr_deg_TCGA,file = 'Corr_deg_TCGA.Rdata')
rm(list = ls())
cor <- corr_deg_TCGA$deg_immScore
cor <- cor[cor$relationship=='Positive',]
gene <- rownames(cor)
exp <- exp1[gene,]
gene <- rownames(exp)
exp <- apply(exp,2,as.numeric)
rownames(exp) <- gene
exp <- as.data.frame(exp)
exp <- t(exp)
Pca <- prcomp(exp1,center = TRUE,scale. = TRUE)
ICI_socre <- Pca$x[,1]
ICI_socre <- cbind(ID=rownames(Pca$x),Score=Pca$x[,1])
ICI_socre <- as.data.frame(ICI_socre)
ICI_socre$Score <- as.numeric(ICI_socre$Score)
ICI_socre <- merge(Gene_ClusterLogRank$data,ICI_socre,by.x = 'ID',by.y = 'ID')
cutoff <- cal_metrics(ICI_socre$Statu,ICI_socre$Score)
ICI_socre$ICI_group <- ifelse(ICI_socre$Score > 5.669764,'High','Low')
a <- ICI_socre[,c(-1,-2,-5)]
ICI_socre <- list(exp=exp,pca=Pca,ICI_socre=ICI_socre,fit=fit,plot=p,a=a)
save(ICI_socre,file = 'ICI_score.Rdata')
###################metabric ICI socre
rownames(MetaBric_Corhort) <- MetaBric_Corhort$Group.1
s <- rownames(MetaBric_Corhort)
MetaBric_Corhort <- MetaBric_Corhort[,-1]
x <- intersect(gene,s)
exp <- MetaBric_Corhort[x,]
exp <- t(exp)
exp <- as.data.frame(exp)
fun <- function(x){
       y <- ifelse(x=NA,0,x)
       return(x)
}
exp1 <- apply(exp,2,fun)
Clin <- MetaBric_CorhortClin$clin[,c(1,6,8)]
Statu <- Clin$VITAL_STATUS
library(limma)
Statu <- strsplit2(Statu,' ')
Clin <- cbind(Clin,Statu)
colnames(Clin) <- c('ID','OS','statu','Statu')
Clin <- Clin[!Clin$V4=='Causes',]
Clin <- Clin[,c(-5,-6,-7)]
Clin <- Clin[order(Clin$OS),]
rownames(Clin) <- 1:2012
Clin <- Clin[1:1484,]
ID <- Clin$ID
ID <- strsplit2(ID,'[.]')
ID <- paste(ID[,1],'_',ID[,2],sep = '')
Clin$ID <- ID
Clin <- Clin[,-3]
ID <- rownames(exp)
rownames(exp) <- ID
table(is.na(exp))
dim(exp)
exp <- exp[Clin$ID,]
ID <- intersect(ID,Clin$ID)
exp <- exp[ID,]
exp1 <- na.omit(exp)
ICI_socre <- merge(Clin,ICI_socre,by.x = 'ID',by.y = 'ID')
ICI_socre$Statu <- ifelse(ICI_socre$Statu=='Died',1,
                          ifelse(ICI_socre$Statu=='Living',0,2))
ICI_socre <- ICI_socre[!ICI_socre$Statu==2,]
ICI_socre$Score <- as.numeric(ICI_socre$Score)
ICI_socre$OS <- as.numeric(ICI_socre$OS)
a <- cal_metrics(ICI_socre$Statu,ICI_socre$Score)
ICI_socre$ICI_score <- ifelse(ICI_socre$Score > a$cutoff,'High','Low')
table(ICI_socre$ICI_score)
ICI_socreMetabric <- list(exp=exp,pca=Pca,ICI_socre=ICI_socre,a=a,fit=fit,plot=p)
ICI_socreMetabric$annotation <- 'The OS were not signigicant, and the DEGs was 274 less than 314'
Clin <- MetaBric_CorhortClin$clin[,c(1,9,10)]
FRS <- Clin$RFS_STATUS
FRS <- strsplit2(FRS,':')
Clin <- cbind(Clin,FRS)
Clin <- Clin[,c(-2,-5)]
colnames(Clin) <- c('ID','Months','Statu')
ID <- strsplit2(Clin$ID,'-')
ID <- paste(ID[,1],'_',ID[,2],sep = '')
Clin$ID <- ID
ICI_socre <- merge(ICI_socre,Clin,by.x = 'ID',by.y = 'ID')
ICI_socre$Statu.y <- as.numeric(ICI_socre$Statu.y)
ICI_socreMetabric$FPS <- list(ICI_socre=ICI_socre,fit=fit,plot=p)
ICI_socreMetabric$FPS$annotation <- 'Moths is the FPS, Statu.y is the statu of FPS'
save(ICI_socreMetabric,file = 'ICI_scoreMetabric.Rdata')
###############GSE78200 Validation corhort
GSE78220 <- GSE78220_PatientFPKM
rownames(GSE78220) <- GSE78220$Gene
clin <- phe[,c(1,13,14,22,16,17)]
library(limma)
s <- cbind(ID=clin$title,sex=strsplit2(clin$characteristics_ch1.3,':')[,2],
           age=strsplit2(clin$characteristics_ch1.4,':')[,2],
           biospyTime= strsplit2(clin$characteristics_ch1.12,':')[,2],
           Day=strsplit2(clin$characteristics_ch1.6,':')[,2],
           Statu=strsplit2(clin$characteristics_ch1.7,':')[,2])
clin <- s
gene <- colnames(ICI_socre$exp)
exp <- GSE78220[gene,]
exp <- as.data.frame(exp)
rownames(exp) <- exp$Gene
exp <- na.omit(exp)
exp <- exp[,-1]
s <- rownames(exp)
exp <- apply(exp,2,as.numeric)
rownames(exp) <- s
GSE78220 <- list(exp=GSE78220,clin=clin)
save(GSE78220,file = 'GES78200.Rdata')
exp <- t(exp)
pca <- prcomp(exp,scale. = F,center = TRUE)
ICI_socre <- list(exp=exp,pca=pca,Score=cbind(id=rownames(pca$x),ICI_socre=pca$x[,1]))
ICI_socre$Score[,1] <- strsplit2(ICI_socre$Score[,1],'[.]')[,1]
score <- ICI_socre$Score
score <- merge(clin,score,by.x='ID',by.y='id')
score$statu <- ifelse(score$Statu==' Dead',1,0)
score$ICI_score <- as.numeric(score$ICI_socre)
score <- score[,-7]
cutoff <- cal_metrics(score$statu,score$ICI_score)
colnames(score) <- c(colnames(score)[1:7],'score')
score$ICI_score <- ifelse(score$score > cutoff$cutoff,'High','Low')
score$Day <- as.numeric(score$Day)
b <- phe[,c(1,8)]
a <- score[score$biospyTime==' pre-treatment',]
colnames(b) <- c('ID','Response')
z <- strsplit2(b$Response,' ')
z <- as.data.frame(z)
q <- strsplit2(z$V2,'_')
q <- as.data.frame(q)
Response <- paste(q$V2,'_',z$V3,sep = '')
b$Response <- Response
score <- merge(score,b,by.x = 'ID',by.y = 'ID')
z <- score[!score$score=='8065.68673890251',]
z$cut <- ifelse(z$score > cut$cutoff,'High','Low')
GSE72800ICI_ScoreOS <- list(exp=exp,pca=pca,cutoff=cutoff,score=score,fit=fit,plot=p)
save(GSE72800ICI_ScoreOS,file = 'GSE72800ValidationICI.Rdata')
##########Imvigor210 corhort validation
library(IMvigor210CoreBiologies)
rm(list = ls())
p_load(DESeq, lsmeans, spatstat)
data(cds)
expMatrix <- counts(cds)
eff_length2 <- fData(cds)[,c("entrez_id","length","symbol")]
rownames(eff_length2) <- eff_length2$entrez_id
head(eff_length2)
feature_ids <- rownames(expMatrix)
expMatrix <- expMatrix[feature_ids %in% rownames(eff_length2),]
mm <- match(rownames(expMatrix),rownames(eff_length2))
eff_length2 <- eff_length2[mm,]

x <- expMatrix/eff_length2$length
eset <- t(t(x)/colSums(x))*1e6
summary(duplicated(rownames(eset)))
BiocManager::install('pacman')
library(pacman)
eset <- IOBR::anno_eset(eset = eset,
                        annotation = eff_length2,
                        symbol = "symbol",
                        probe = "entrez_id",
                        method = "mean")
tumor_type <- "blca"
if(max(eset)>100) eset <- log2(eset+1)
pdata <- pData(cds)
colnames(pdata) <- gsub(colnames(pdata),pattern = " ",replacement = "_")
pdata <- rownames_to_column(pdata[,c("binaryResponse",
                                     "FMOne_mutation_burden_per_MB",
                                     "Neoantigen_burden_per_MB",
                                     "censOS","os")],var = "ID")
colnames(pdata)<-c("ID","BOR_binary","TMB","TNB","status","time")
exp <- as.data.frame(x)
exp <- cbind(probes=rownames(exp),exp)
exp <- merge(eff_length2,exp,by.x = 'entrez_id',by.y = 'probes')
exp <- exp[,c(-1,-2)]
exp <- aggregate(exp,by=list(exp$symbol),FUN = mean)
exp <- exp[,-2]
exp <- exp[-1,]
rownames(exp) <- exp$Group.1
exp <- exp[,-1]
gene <- colnames(ICI_socre$exp)
exp1 <- exp[gene,]
exp1 <- na.omit(exp1)
exp1 <- t(exp1)
pca <- prcomp(exp1,center = TRUE,scale. = TRUE)
score <- pca$x[,1]
ICI_socre <- cbind(id=rownames(pca$x),score)
ICI_socre <- as.data.frame(ICI_socre)
ICI_socre$score <- as.numeric(ICI_socre$score)
clin <- cbind(id=rownames(pdata),pdata[,21:22])
ICI_socre <- merge(ICI_socre,clin,by.x = 'id',by.y = 'id')
colnames(ICI_socre) <- c('id','score','Day','Statu')
cutoff <- cal_metrics(ICI_socre$Statu,ICI_socre$Day)
ICI_socre$ICI_score <- ifelse(ICI_socre$score > cutoff$cutoff,'High','Low')
ICI_socre$ICI_score <- ifelse(ICI_socre$score > median(ICI_socre$score),'High','Low')
install.packages('randomForest')
gene <- rownames(corr_deg_TCGA$deg_immScore)
exp1 <- exp1[gene,]
gene <- rownames(exp1)
exp1 <- apply(exp1,2,as.numeric)
rownames(exp1) <- gene
exp1 <- t(exp1)
exp1 <- as.data.frame(exp1)
exp1 <- cbind(ID=rownames(exp1),exp1)
clin <- OV_data$`Summary data`[,c(1,11)]
exp1 <- merge(clin,exp1,by.x = 'ID',by.y = 'ID')
exp1$Statu <- ifelse(exp1$Statu=='Alive',0,1)
gene <- colnames(exp1)[3:ncol(exp1)]
library(limma)
gene <- colnames(exp2)[2:ncol(exp2)]
gene <- strsplit2(gene,'-')
gene <- as.data.frame(gene)
gene <- ifelse(gene$V2=='1' |gene$V2=='2' |gene$V2=='3' |gene$V2=='4' |gene$V2=='5'|
                 gene$V2=='685N3.1' | gene$V2=='DOA' | gene$V2=='DOB' |gene$V2=='DQA1'|
                 gene$V2=='DQA2' | gene$V2=='F' | gene$V2=='G', paste(gene[,1],'_',gene[,2],sep = ''),gene$V1)
colnames(exp2)[2:ncol(exp2)] <- gene 
a <- paste(colnames(exp1)[3:ncol(exp1)],collapse = ' + ')
colnames(exp1)[3:ncol(exp1)] <- gene
exp2 <- exp1[,-1]
library(randomForest)
forest <- randomForest(Statu ~., exp2)
exp2$Statu <- as.factor(exp2$Statu)
s <- forest$importance
s <- as.data.frame(s)
s <- cbind(gene=rownames(s),s)
s <- s[order(s$MeanDecreaseGini,decreasing = TRUE),]
b <- importance(forest,type = 1)
###################################################################
id1 <- "/Users/llls2012163.com/TCGA/TCGA_BRCA/OV_data.Rdata"
id2 <- "/Users/llls2012163.com/R-scrpit/ImmuneFilither/gene_cluster_TCGA/TCGApatientsexp.Rdata"
id3 <- "/Users/llls2012163.com/R-scrpit/ImmuneFilither/gene_cluster_TCGA/Corr_deg_TCGA.Rdata"
load(id1)
load(id2)
load(id3)
gene <- rownames(corr_deg_TCGA$deg_immScore)
exp <- exp1[gene,]
gene <- rownames(exp)
exp <- apply(exp,2,as.numeric)
rownames(exp) <- gene
if(max(exp)>100) exp <- log2(exp+1)############标准化数据
s <- as.data.frame(`geneClusterTcga.k=3.consensusClass`)
table(s$V2)
s[s$V2==3,]$V1
which(colnames(exp)== "TCGA_A2_A3XW")
exp <- exp[,-765]
clin <- OV_data$`Summary data`[,c(1,11,12)]
s <- merge(clin,s,by.x = 'ID',by.y = 'V1')
colnames(s)[4] <- 'Gene_cluster'
s$Gene_cluster <- ifelse(s$Gene_cluster==1,'A',
                         ifelse(s$Gene_cluster==2,'B','C'))
s$Day <- as.numeric(s$Day)
s$statu <- ifelse(s$Statu=='Alive',0,1)
Gene_ClusterLogRankNew <- list(exp=exp,d=d,data=s,fit=fit,plot=p)
save(Gene_ClusterLogRankNew,file = 'ConsensusClusterGeneNew.Rdata')
exp <- as.data.frame(exp)
exp <- t(exp)
exp <- cbind(ID=rownames(exp),exp)
exp <- merge(s[,c(1,4)],exp)
exp2 <- exp[,-1]
exp2$Gene_cluster <- as.factor(exp2$Gene_cluster)
exp3 <- exp1[gene,]
library(limma)
s <- strsplit2(gene,'_')
s <- ifelse(s$V2=='2' | s$V2=='DOB',paste(s[,1],'-',s[,2],sep = ''),s[,1])
exp3 <- exp1[s,]
gene <- rownames(exp3)
exp3 <- apply(exp3,2,as.numeric)
if(max(exp3)>100) exp3 <- log2(exp3+1)############标准化数据
rownames(exp3) <- gene
exp3 <- as.data.frame(exp3)
exp3 <- t(exp3)
PCA <- prcomp(exp3,center = TRUE,scale. = TRUE)
ICI_score <- cbind(ID=rownames(PCA$x),score=PCA$x[,1])
s <- OV_data$`Summary data`[,c(1,11,12)]
clin <- merge(s,ICI_score,by.x = 'ID',by.y = 'ID')
clin$statu <- ifelse(clin$Statu=='Alive',0,1)
clin$Day <- as.numeric(clin$Day)
clin$score <- as.numeric(clin$score)
cutoff <- cal_metrics(clin$statu,clin$score)
clin$ICI_score <- ifelse(clin$score > cutoff$cutoff,'High','Low')
TCGA_ICIscore <- list(gene=gene,data=exp3,PCA=PCA,clin=clin,fit=fit,plot=p,cutoff=cutoff)
TCGA_ICIscore$annotation <- 'The top 1/3 improtant gene from RandomForest were seleted to construt ICI Score'
save(TCGA_ICIscore,file = 'TCGA_ICIscore.Rdata')
##########################################Imvigor 210 validation
eset <- as.data.frame(eset)
exp4 <- cbind(id=rownames(eset),eset)
exp4 <- merge(eff_length2,exp4,by.x = 'entrez_id',by.y = 'id')
exp4 <- exp4[,c(-1,-2)]
rownames(exp4) <- exp4$symbol
colnames(exp4)[1] <- 'id'
exp4 <- as.data.frame(exp4)
id <- exp4$id
exp4[,2:ncol(exp4)] <- apply(exp4[,2:ncol(exp4)],2,as.numeric)
exp4 <- aggregate(exp4,by=list(exp4$id),FUN = mean,na.rm=TRUE)
exp4 <- exp4[-1,]
rownames(exp4) <- exp4$Group.1
exp4 <- exp4[,c(-1,-2)]
exp5 <- exp4[gene,]
exp5 <- na.omit(exp5)
exp5 <- t(exp5)
PCA <- prcomp(exp5,center = TRUE,scale. = TRUE)
ICI_score <- cbind(id=rownames(PCA$x),score=PCA$x[,1])
ICI_score <- as.data.frame(ICI_score)
clin <- cbind(id=rownames(pdata),statu=pdata$censOS,OS=pdata$os)
clin <- merge(clin,ICI_score,by.x='id',by.y='id')
clin[,2:4] <- apply(clin[,2:4],2,as.numeric)
cutoff <- cal_metrics(clin$statu,clin$score)
clin$ICI_score <- ifelse(clin$score > cutoff$cutoff,'High','Low')#################
IMvigor210ICI <- list(exp=exp5,PCA=PCA,cutoff=cutoff,clin=clin,fit=fit,plot=p)
save(IMvigor210ICI,file = 'IMvigor210CorhortValidate.Rdata')
pdata$binaryResponse <- as.character(pdata$binaryResponse)
s <- cbind(id=rownames(pdata),Response=pdata[,2])
clin <- merge(clin,s,by.x='id',by.y = 'id')
clin1 <- na.omit(clin)
library(dplyr)
a <- data.frame(table(clin1$ICI_score,clin1$Response))
a$percent <- a$percent*100
a$label <- paste0(sprintf('%.1f',a$percent),'%')
pvalue <- chisq.test(c(54,14,189,41,ncol=2))$p.value
library(ggplot2)
library(ggprism)
library(plyr)
p1 <- ggplot(a,aes(Var1,percent,fill=Var2))+
      geom_bar(stat = 'identity',position = position_stack())+
      scale_fill_manual(values = c('green','red'),label=c('SD/PD','CR/PR'))+
      scale_y_continuous(labels = scales::percent_format(scale=1))+
      labs(x='ICI_score',y='Percent Weight',fill='')+
      geom_text(aes(label=label),vjust=0,size=6,color='black')+
      annotate(geom = 'text',cex=6,x=1.5,y=90,label=paste0('Chi-squared test P',ifelse(pvalue<0.001,'<0.001',paste0('=',round(pvalue,3),color='black'))))+
      theme_prism(base_size = 12)+theme(legend.position = 'top')
      
      p1
IMvigor210ICI$Stack <- list(data=a,plot=p1)      
b <- clin[,c(4,6)]
b <- na.omit(b)
Pvalue <- wilcox.test(b[b$Response=='SD/PD',]$score,b[b$Response=='CR/PR',]$score)
p2 <- ggplot(b,aes(Response,score,fill=Response))+geom_boxplot()+
      scale_fill_manual(values = c('green','red'))+labs(x='Anti-PD1',y='ICI_score')+
      theme_prism(base_size = 12)+theme(legend.position = 'top')+
      annotate(geom = 'text',cex=6,x=1.2,y=-25,label=paste0('Wilcoxon test P','< 0.05'))
p2
IMvigor210ICI$Bar <- list(data=b,plot=p2)
save(IMvigor210ICI,file = 'IMvigor210Validate.Rdata')
########Metabric
gene <- TCGA_ICIscore$gene
exp <- MetaBric_Corhort[gene,]
rownames(MetaBric_Corhort) <- MetaBric_Corhort$Group.1
exp <- MetaBric_Corhort[gene,]
exp <- na.omit(exp)
exp <- exp[,-1]
exp <- t(exp)
PCA <- prcomp(exp,center = TRUE,scale. = TRUE)
ICI_score <- cbind(id=rownames(PCA$x),score=PCA$x[,1])
clin <- Metabric_ICI$KMdata[,c(1,2,8)]
ICI_score <- as.data.frame(ICI_score)
id <- ICI_score$id
library(limma)
id <- strsplit2(id,'[.]')
id <- paste(id[,1],'-',id[,2],sep = '')
ICI_score <- cbind(id,ICI_score)
ICI_score <- ICI_score[,-2]
ICI_score <- merge(ICI_score,clin,by.x='id',by.y='PATIENT_ID')
ICI_score$score <- as.numeric(ICI_score$score)
cutoff <- cal_metrics(ICI_score$statu,ICI_score$score)
ICI_score$ICI_score <- ifelse(ICI_score$score > cutoff$cutoff,'High','Low')
clin <- MetaBric_CorhortClin$clin[,c(1,9,10)]
clin$RFS <- strsplit2(clin$RFS_STATUS,':')[,1]
clin <- clin[,-2]
ICI_score <- merge(clin,ICI_score,by.x = 'PATIENT_ID',by.y = 'id')
ICI_score$RFS <- as.numeric(ICI_score$RFS)
Metabric_ICIScore <- list(exp=exp,PCA=PCA,ICI_socre=ICI_score,cutoff=cutoff,LogRankOS=list(fit=fit,plot=p),LogRankFPS=list(fit=fit1,plot=p1))
save(Metabric_ICIScore,file = 'Metabric_iciscore.Rdata')
rm(list = ls())
######GSE78200 PD1
exp <- GSE78220$exp
exp[,2:ncol(exp)] <- apply(exp[,2:ncol(exp)],2,as.numeric)
exp1 <- exp[,-1]
if(max(exp1)>100) exp1 <- log2(exp1+1)############标准化数据
rownames(exp1) <- exp$Gene
exp1 <- t(exp1)
exp1 <- as.data.frame(exp1)
gene <- TCGA_ICIscore$gene
exp1 <- t(exp1)
s <- intersect(gene,rownames(exp1))
exp1 <- exp1[s,]
exp1 <- t(exp1)
PCA <- prcomp(exp1,center = TRUE,scale. = FALSE)
ICI_Score <- cbind(id=rownames(PCA$x),score=PCA$x[,1])
ICI_Score <- as.data.frame(ICI_Score)
id <- strsplit2(ICI_Score$id,'[.]')[,1]
ICI_Score$id <- id
ICI_Score$score <- as.numeric(ICI_Score$score)
clin <- GSE78220$clin
ICI_Score <- merge(clin,ICI_Score,by.x='ID',by.y='id')
ICI_Score$statu <- ifelse(ICI_Score$Statu==' Alive',0,1)
cutoff <- cal_metrics(ICI_Score$statu,ICI_Score$score)
ICI_Score$ICI_score <- ifelse(ICI_Score$score > cutoff$cutoff,'High','Low')
ICI_Score$Day <- as.numeric(ICI_Score$Day)
ICI_score <- na.omit(ICI_Score)
GSE78220iciVa <- list(exp=exp1,PCA=PCA,cutoff=cutoff,ICI_socre=ICI_score,fit=fit,plot=p)
save(GSE78220iciVa,file = 'GSE78200ici_score.Rdata')
phe <- eset$GSE78220_series_matrix.txt.gz@phenoData
phe <- phe@data
s <- strsplit2(phe$characteristics_ch1.1,':')[,2]
s <- cbind(id=phe$title,s)
s <- as.data.frame(s)
phe <- s
phe$s <- ifelse(phe$s=='Complete Response','CR',
                ifelse(phe$s=='Partial Response','PR','PD'))
colnames(s) <- c('id','Response')
ICI_score <- merge(ICI_score,phe,by.x = 'ID',by.y = 'id')
ICI_score <- ICI_score[,-10]
ICI_score$Response <- ifelse(ICI_score$s=='PD','PD','CR/PR')
p1 <- ggplot(ICI_score,aes(Response,score))+geom_boxplot()
p1
rm(gse18728_ICI)
############GSE91061
exp <- read.csv("~/R-scrpit/ImmuneFilither/Immune therapy/GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv", header=TRUE)
exp <- as.data.frame(exp)
s <- eff_length2[,c(1,3)]
exp <- merge(s,exp,by.x = 'entrez_id',by.y = 'X')
exp <- exp[,-1]
exp <- aggregate(exp,by=list(exp$symbol),FUN = mean)
exp <- exp[,-2]
colnames(exp)[1] <- 'id'
if(max(exp[,2:ncol(exp)])>100) exp[,2:ncol(exp)] <- log2(exp[,2:ncol(exp)]+1)############标准化数据
rownames(exp) <- exp$id
gene <- TCGA_ICIscore$gene
exp1 <- exp[gene,]
exp1 <- na.omit(exp1)
clin <- eset$GSE91061_series_matrix.txt.gz@phenoData@data
m <- colnames(exp1)[2:ncol(exp1)]
m <- strsplit2(m,'[.]')
m <- paste(m[,1],'_',m[,2],sep = '')
colnames(exp1)[2:ncol(exp1)] <- m
exp1 <- exp1[,-1]
exp1 <- t(exp1)
PCA <- prcomp(exp1,center = TRUE,scale. = TRUE)
ICI_score <- cbind(id=rownames(PCA$x),score=PCA$x[,1])
library(XML)
id ='/Users/llls2012163.com/R-scrpit/ImmuneFilither/GSE35640_family.xml'
s <- xmlParse(file = id)
rootload <- xmlRoot(s)
data <- xmlSApply(rootload,function(x) xmlSApply(x,xmlValue))
data <- data.frame(t(data),row.names = NULL)
data$Sample.2
exp <- BRCA_KR$exp_BRCA_KR
gene <- TCGA_ICIscore$gene
exp <- exp[gene,]
exp <- na.omit(exp)
exp <- exp[,-1]
exp <- log2(exp+1)
exp <- t(exp)
pca <- prcomp(exp,center = TRUE,scale. = TRUE)
ICI_score <- cbind(ID=rownames(pca$x),score=pca$x[,1])
clin <- BRCA_KR$Clinical_BRCA_KR[,c(1,3,7)]
colnames(clin) <- c('ID','Statu','Day')
s <- merge(clin,ICI_score,by.x = 'ID',by.y = 'ID')
s$score <- as.numeric(s$score)
s$statu <- ifelse(s$Statu=='alive',0,1)
cutoff <- cal_metrics(s$statu,s$score)
s$ICI_score <- ifelse(s$score > cutoff$cutoff,'High','Low')
ICI_Score <- list(exp=exp,PCA=pca,cutoff=cutoff,ICI_score=ICI_score,data=s,plot=p)
save(ICI_Score,file = 'ICI_socreBRCA_KRcohort.Rdata')
##############
clin <- TCGA_ICIscore$clin
s <- OV_data$`Summary data`[,c(1,3)]
s <- merge(s,clin,by.x = 'ID',by.y = 'ID')
s1 <- s[s$ICI_Cluster=='Cluster_A',]
s2 <- s[s$ICI_Cluster=='Cluster_B',]
s3 <- s[s$ICI_Cluster=='Cluster_C',]
###########
s <- OV_data$`Summary data`
clin <- Gene_ClusterLogRankNew$data
clin <- clin[,c(1,4)]
s <- merge(clin,s,by.x = 'ID',by.y = 'ID')
s1 <- s[,c(1,2,14:43)]
colnames(s1)
m <- colnames(s1)[3:32]
output <- vector()
s1[,3:32] <- apply(s1[,3:32],2,as.numeric)
for(i in m){
   x =s1[s1$Gene_cluster=='A',i]
   y=s1[s1$Gene_cluster=='B',i]
   z= s1[s1$Gene_cluster=='C',i]
   a =kruskal.test(list(x,y,z))$p.value
   if(length(output)==0){
     output=a
   }else{
     output=c(output,a)
   }
   
}
output <- cbind(cellType=m,pvalue=output)
output <- as.data.frame(output)
output$pvalue <- as.numeric(output$pvalue)
write.csv(output,'KruskiltestGeneCluster.csv')
rm(list = ls())
outTab <- as.data.frame(outTab)
outTab$pvalue <- as.numeric(outTab$pvalue)
outTab$pvalue1 <- ifelse(outTab$pvalue < 0.0001,'****',
                         ifelse(outTab$pvalue < 0.001,'***',
                                ifelse(outTab$pvalue < 0.01,'**',
                                       ifelse(outTab$pvalue < 0.05,'*','ns'))))
outTab1 <- outTab[,-5]
colnames(outTab1)[5] <-'pvalue'
exp <- TCGA_ICIscore$data
clin <- TCGA_ICIscore$clin
s <- clin[,c(1,6)]
exp <- as.data.frame(exp)
exp <- cbind(ID=rownames(exp),exp)
exp <- merge(s,exp,by.x = 'ID',by.y = 'ID')
exp <- exp[order(exp$ICI_score),]
m <- colnames(exp[,3:ncol(exp)])
FC <- vector()
exp1 <- exp[exp$ICI_score=='High',][,c(-1,-2)]
exp2 <- exp[exp$ICI_score=='Low',][,c(-1,-2)]
a <- apply(exp1,2,mean)
b <- apply(exp2,2,mean)
FC <- a/b
FC <- cbind(Gene=m,FC=FC)
FC <- as.data.frame(FC)
FC <- FC[order(FC$FC),]
FC$FC <- as.numeric(FC$FC)
FC$FC <- log2(FC$FC)
colnames(FC)[2] <- 'logFC'
#####################
s <- OV_data$`Summary data`[,c(1,2)]
s1 <- TCGA_ICIscore$clin
s2 <- merge(s1,s,by.x='ID',by.y='ID')
s2$subgroup <- ifelse(s2$ICI_score=='High' & s2$treatment_type=='PT','PT+High',
                      ifelse(s2$ICI_score=='High' & s2$treatment_type=='RT','RT+High',
                             ifelse(s2$ICI_score=='Low'& s2$treatment_type=='PT','PT+Low',
                                    'RT+Low')))
ICI_socrePerSub <- list(Data=s2,fit=fit,plot=p)
save(ICI_socrePerSub,file = 'ICI_scorePerSub.Rdata')
rm(list = ls())
##########
gene <- TCGA_ICIscore$gene
exp <- exp1[gene,]
gene <- rownames(exp)
exp <- apply(exp,2,as.numeric)
exp <- log2(exp+1)
rownames(exp) <- gene
s <- OV_data$`Summary data`[,c(1,41,42)]
exp <- t(exp)
exp <- cbind('ID'=rownames(exp),exp)
exp <- as.data.frame(exp)
s <- merge(s,exp,by.x='ID',by.y='ID')
s <- s[,-1]
s <- apply(s,2,as.numeric)
s <- as.data.frame(s)
#############
pvalue <- vector()
coef <- vector()
for(i in gene){
    a <- s$StromalScore
    b <- s[,i]
    C <- cor.test(a,b,method = 'pearson')
    if(length(pvalue)==0){
      pvalue=c$p.value
    }else{
      pvalue=c(pvalue,C$p.value)
    }
    if(length(coef)==0){
      coef=C$estimate
    }else{
      coef=c(coef,C$estimate)
    }
}
ImmuCor <- cbind(gene=gene,cor=coef,pvalue=pvalue)
StraCor <- cbind(gene=gene,cor=coef,pvalue=pvalue)
ImmuCor[,2:3] <- apply(ImmuCor[,2:3],2,as.numeric)
ImmuCor <- as.data.frame(ImmuCor)
Imgene <- ImmuCor[ImmuCor$pvalue < 0.05,]$gene
StraCor[,2:3] <- apply(StraCor[,2:3],2,as.numeric)
StraCor <- as.data.frame(StraCor)
Stgene <- StraCor[StraCor$pvalue < 0.05,]$gene
###############BRCA-FR
exp <- exp_array.BRCA.FR[,c(1,8,9)]
library(reshape2)
exp1 <- dcast(exp,gene_id~icgc_donor_id)
rownames(exp1) <- exp1$gene_id
gene <- TCGA_ICIscore$gene
exp2 <- exp1[gene,]
exp2 <- na.omit(exp2)
exp2 <- exp2[,-1]
exp2 <- t(exp2)
s <- prcomp(exp2,center = TRUE,scale. = TRUE)
ICI_score <- cbind(ID=rownames(s$x),score=s$x[,1])
m <- donor.BRCA.FR[,c(1,6,9,11)]
m$Day <- m$donor_age_at_last_followup-m$donor_age_at_diagnosis
m$Day <- 365*m$Day
m$Statu <- ifelse(m$donor_vital_status=='alive',0,1)
m <- merge(m,ICI_score,by.x='icgc_donor_id',by.y='ID')
m$score <- as.numeric(m$score)
cutoff <- cal_metrics(m$Statu,m$score)
ICI_score <- ifelse(m$score > cutoff$cutoff,'High','Low')
clin <- TCGA_ICIscore$clin
tmb <- OV_data$`Summary data`[,c(1,45,46)]
s <- merge(clin,tmb,by.x='ID',by.y='ID')
s <- s[order(s$ICI_score),]
write.csv(s,file = 'TMB_ICIscore.csv')
wilcox.test(s[s$ICI_score=='Low',]$total_perMB,s[s$ICI_score=='High',]$total_perMB)
h <- Gene_ClusterLogRankNew$data
s <- merge(s,h,by.x='ID',by.y='ID')
library(ggplot2)
library(ggprism)
p <- ggplot(s,aes(score,total_perMB_log,color=Gene_cluster))+geom_point(size=0.8)+
     theme_prism(base_size = 10)+labs(x='ICI score',y='Normalized Tumor mutation burden')+
     geom_smooth(aes(x=score,y=total_perMB_log),method = "lm",color='black',se=TRUE,size=0.5,span = 0.4)+
     theme(legend.position = 'top')+scale_color_manual(values = c(A='red',B='blue',C='green'))

p
s <- s[,c(-9,-10,-12)]
colnames(s)[c(2,3,5)] <- c('Statu','Day','statu')
cutoff <- cal_metrics(s$statu,s$total_perMB)
s$TMB <- ifelse(s$total_perMB > cutoff$cutoff,'High','Low')
s$group <- ifelse(s$ICI_score=='High' & s$TMB=='High','High_ICI+High_TMB',
                  ifelse(s$ICI_score=='High' & s$TMB=='Low','High_ICI+Low_TMB',
                         ifelse(s$ICI_score=='Low' & s$TMB=='High','Low_ICI+High_TMB',
                                'Low_ICI+Low_TMB')))
ICIscore_TMB <- list(cutoff=cutoff,data=s,plot1=p,plot2=p1,plot3=p2,annotation='Figure4')
save(ICIscore_TMB,file = 'ICI_TMB.Rdata')
s <- ICIscore_TMB$data
s <- s[,c(1,6,9,10)]
a <- a[,-12]
library(limma)
ID <- a$Tumor_Sample_Barcode
ID <- strsplit2(ID,'-')[,1:3]
ID1 <- paste(ID[,1],'_',ID[,2],sep = '')
ID <- paste(ID1,'_',ID[,3],sep = '')
a$V1 <- ID
a <- merge(s,a,by.x = 'ID',by.y = 'V1')
a <- cbind(a[,c(-1,-2,-3,-4)],a[,1:4])
a$Gene_cluster <- ifelse(a$Gene_cluster=='A','Cluster_A',
                         ifelse(a$Gene_cluster=='B','Cluster_B','Cluster_C'))
write.csv(a,file = 'CLINICAL.clin')
