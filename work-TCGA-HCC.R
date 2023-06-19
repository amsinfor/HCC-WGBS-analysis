#-------validation of WGBS+RNAseq results in TCGA and GEO datasets
library(parallel)
library(data.table)
library(RColorBrewer)
library(pheatmap)
library(ggpubr);
library(ggfortify)#
library(limma) #
library(sva)#
library(bnstruct)
library(Rtsne)
library(minfi);
library(R.utils);
library(MASS)
library(glmnet);
library(pROC)
library(cowplot)
library(minfi);
library(Rtsne)
library(MASS)
library(bnstruct)

#---------------- dataset 
#:four target genes
c("CYP2C8","CYP2C19","CYP2A6","ADH1A")->target_genes;

############################################################################################################################
#-----------------------------------------------------------------------------------------------------------------------
#----------GEO dataset
GEO_merged_methy_objRef$getConditionGroup("SampleType",c("Normal","Tumor"))->GEO_merged_methy_group;
GEO_merged_methy_objRef$getGroupMatrix("SampleType",c("Normal","Tumor"))->GEO_merged_methy_matrix;
#--- do ranktest for DMProbes 
do_ranktest_diff(GEO_merged_methy_matrix,GEO_merged_methy_group,is_pair=FALSE)->GEO_merged_methy_ranktest;
#-----------------------------------------------------------------------------------------------------------------------
#-- TCGA HCC
TCGA_HCC_methy_objRef$getConditionGroup("SampleTypeCode",c("Normal","Tumor"))->TCGA_HCC_methy_TN_group;
TCGA_HCC_methy_objRef$getGroupMatrix("SampleTypeCode",c("Normal","Tumor"))->TCGA_HCC_methy_TN_matrix;
draw_by_MDsplot_by_matrix(TCGA_HCC_methy_TN_matrix,TCGA_HCC_methy_TN_group[-which(TCGA_HCC_methy_TN_group$SampleID%in%TCGA_outliers),],show_label=F)
draw_by_tsne_by_matrix(TCGA_HCC_methy_TN_matrix,TCGA_HCC_methy_TN_group[-which(TCGA_HCC_methy_TN_group$SampleID%in%TCGA_outliers),],show_label=F)
#- do ranktest for DMProbes
do_ranktest_diff(TCGA_HCC_methy_TN_matrix,TCGA_HCC_methy_TN_group[-which(TCGA_HCC_methy_TN_group$SampleID%in%TCGA_outliers),],is_pair=FALSE)->TCGA_HCC_methy_ranktest;
#---pheatmap for four target genes：
TCGA_HCC_methy_ranktest$TransID[which(TCGA_HCC_methy_ranktest$gName%in%target_genes & TCGA_HCC_methy_ranktest$FDR<=1)]->target_transIDs;
prepare_DMC_heatmap(TCGA_HCC_methy_TN_matrix,TCGA_HCC_methy_ranktest,TCGA_HCC_methy_TN_group,"Condition",target_transIDs)->TCGA_HCC_DMC.sort;#top 1%
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(250),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(250))->fill_colors;
#：
pheatmap(TCGA_HCC_DMC.sort$Methyd,show_rownames=T,cluster_rows=F,cluster_cols=F,show_colnames=F,annotation_row=TCGA_HCC_DMC.sort$DMC,col=fill_colors,fontsize_col=10,annotation_col=TCGA_HCC_DMC.sort$Condition)

############################################################################################################################
############################################################################################################################
#--------TCGA-HCC RNA-seq data：
TCGA_HCC_exp_objRef$getConditionGroup("SampleTypeCode",c("Normal","Tumor"))->TCGA_HCC_exp_TN_group;
TCGA_HCC_exp_objRef$getGroupMatrix("SampleTypeCode",c("Normal","Tumor"))->TCGA_HCC_exp_TN_matrix;
#- do ranktest for DEGs
do_ranktest_diff(TCGA_HCC_exp_TN_matrix,TCGA_HCC_exp_TN_group[-which(TCGA_HCC_exp_TN_group$SampleID%in%TCGA_exp_outliers),],is_pair=FALSE)->TCGA_HCC_exp_ranktest;
#：removed low expressed genes
TCGA_HCC_exp_ranktest[which(TCGA_HCC_exp_ranktest$AveExpr+TCGA_HCC_exp_ranktest$B>=10),]->TCGA_HCC_exp_ranktest_filter;
TCGA_HCC_exp_ranktest_filter$TransID[TCGA_HCC_exp_ranktest_filter$gName%in%target_genes]->target_exp_transIDs;
#---pheatmap for DEGs
prepare_DMC_heatmap(TCGA_HCC_exp_TN_matrix,TCGA_HCC_exp_ranktest,TCGA_HCC_exp_TN_group,"Condition",target_exp_transIDs)->TCGA_HCC_DEG.sort;#top 1%
pheatmap(log2(TCGA_HCC_DEG.sort$Methyd+1),show_rownames=T,cluster_rows=F,cluster_cols=F,show_colnames=F,annotation_row=TCGA_HCC_DEG.sort$DMC,col=fill_colors,fontsize_col=10,annotation_col=TCGA_HCC_DEG.sort$Condition)
#--------- circoheatmap for TCGA-HCC 
library(circlize)
library(ComplexHeatmap)
#
colorRamp2(c(-10,0,15), c("green", "white", "#E41A1C"))->col_fun1
circos.clear()
circos.heatmap(t(log2(TCGA_HCC_exp_TN_matrix[target_exp_transIDs,TCGA_HCC_exp_TN_group$SampleID[which(TCGA_HCC_exp_TN_group$Condition=="Tumor")]]+1)), col = col_fun1,dend.side = "inside",rownames.side = "none",dend.track.height = 0.6)
Legend(title = "TCGA", col_fun = col_fun1)->lgd;
grid.draw(lgd)
#--survival analysis for three subgroups
library(survival)
#--
prepare_pheatmap_kmeans_cluster(TCGA_HCC_exp_tumor_heatmap,3,is_row=FALSE)->test_;
"ExpCluster"->names(test_)[2];
merge(TCGA_HCC_exp_factor,test_,by.x="A0_Samples",by.y="A0_Samples")->test_TCGA_HCC_exp_factor;
draw_survial_curve_custom_v2(test_TCGA_HCC_exp_factor,ncol(test_TCGA_HCC_exp_factor),365,myd_colors,NULL)#pretty good
#--------------------------------------------------------------------------------------------------------
#------在TCGA数据集上，甲基化与表达的关系；------------------- 未使用，因为甲基化与表达差异非常大
HumanMethylation450_annotats[which(HumanMethylation450_annotats$UCSC_RefGene_Name%in%target_genes),c(2,12:13,22,24)]->target_genes_annotats;
paste(target_genes_annotats$UCSC_RefGene_Name,target_genes_annotats$Name,sep="|")->target_genes_annotats$TransID;
#-对每个探针进行计算：表达与甲基化的相关性
lapply(target_transIDs,function(cx){
	unlist(strsplit(cx,split="\\|"))->cx_splits;
	cx_splits[1]->cx_g;
	#cx_splits[2]->cx_p;
	#-
	intersect(TCGA_HCC_exp_factor$A0_Samples,colnames(TCGA_HCC_methy_TN_matrix))->shared_samples
	TCGA_HCC_exp_TN_matrix[grep(paste(cx_g,"\\|",sep=""),rownames(TCGA_HCC_exp_TN_matrix)),shared_samples]->cx_exp;
	TCGA_HCC_methy_TN_matrix[cx,shared_samples]->cx_methy;
	cor.test(cx_exp,cx_methy)->cx_cor1;
	c(cx_cor1$estimate,cx_cor1$p.value);
})->test_res;
data.frame("TransID"=target_transIDs,matrix(unlist(test_res),ncol=2,byrow=T))->test_res;
c("Cor","CorP")->names(test_res)[2:3]
merge(target_genes_annotats,test_res,by.x="TransID",by.y="TransID")->target_genes_annotats;
target_genes_annotats[order(target_genes_annotats$MAPINFO),]->target_genes_annotats;
#--- t-SNE 
TCGA_HCC_exp_TN_matrix[target_exp_transIDs,TCGA_HCC_exp_TN_group$SampleID[which(TCGA_HCC_exp_TN_group$Condition=="Tumor")]]->test_matrix;
test_TCGA_HCC_exp_factor[,c("A0_Samples","ExpCluster")]->test_group;
c("SampleID","Condition")->colnames(test_group);
test_group$SampleID->rownames(test_group);
test_group[colnames(test_matrix),]->test_group;
draw_by_tsne_by_matrix(test_matrix,test_group)
############################################################################################################################
############################################################################################################################
#：Clinical features of the three subgroups
ggboxplot(test_TCGA_HCC_exp_factor,x="ExpCluster",y="Retinol_score",fill="ExpCluster",yscale="log10")+stat_compare_means()+theme(legend.position="n")->p1;
#: Stage_T Stage_N Stage_M   Stage Grade
ggboxplot(subset(test_TCGA_HCC_exp_factor,Stage_T!="Un"),x="Stage_T",y="Retinol_score",fill="Stage_T",yscale="log10",order=c("T1","T2","T3","T4"))+stat_compare_means()+theme(legend.position="n")->p2;
#:
ggboxplot(subset(test_TCGA_HCC_exp_factor,Stage_N!="Un"),x="Stage_N",y="Retinol_score",fill="Stage_N",yscale="log10",order=c("N0","N1"))+stat_compare_means()+theme(legend.position="n")->p3;
#:
ggboxplot(test_TCGA_HCC_exp_factor,x="Stage_M",y="Retinol_score",fill="Stage_M",yscale="log10",order=c("M0","M1","MX"))+stat_compare_means()+theme(legend.position="n")->p4;
#: 
ggboxplot(subset(test_TCGA_HCC_exp_factor,Stage!="Un"),x="Stage",y="Retinol_score",fill="Stage",yscale="log10",order=c("Stage I","Stage II","Stage III","Stage IV"))+stat_compare_means()+theme(legend.position="n")->p5;
#: grade
ggboxplot(test_TCGA_HCC_exp_factor[!is.na(test_TCGA_HCC_exp_factor$Grade),],x="Grade",y="Retinol_score",fill="Grade",yscale="log10",order=c("G1","G2","G3","G4"))+stat_compare_means()+theme(legend.position="n")->p6;
#: virus
ggboxplot(test_TCGA_HCC_exp_factor[!is.na(test_TCGA_HCC_exp_factor$Virus),],x="Virus",y="Retinol_score",fill="Virus",yscale="log10",order=c("HBV","HCV"))+stat_compare_means()+theme(legend.position="n")->p7;
#：age: <50,50-60,60-70,>70
ggboxplot(test_TCGA_HCC_exp_factor,x="Age_group",y="Retinol_score",fill="Age_group",yscale="log10",order=c("<50","50-60","60-70",">=70"))+stat_compare_means()+theme(legend.position="n")->p8;
#------
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,nrow=2,byrow=T)

#---------------------------------------------
#：---------- Genomic characteristics of the three subgroups
#:1. mutation Signatures, TMB
subsetMaf(maf = HCC_SNV_maf, clinQuery = "ExpCluster=='C1'")->HCC_SNV_maf_C1;
plotmafSummary(maf = HCC_SNV_maf_C1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE,top=20)
subsetMaf(maf = HCC_SNV_maf, clinQuery = "ExpCluster=='C2'")->HCC_SNV_maf_C2;
plotmafSummary(maf = HCC_SNV_maf_C2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE,top=20)
subsetMaf(maf = HCC_SNV_maf, clinQuery = "ExpCluster=='C3'")->HCC_SNV_maf_C3;
plotmafSummary(maf = HCC_SNV_maf_C3, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE,top=20)
#：TMB
tcgaCompare(maf = HCC_SNV_maf_C1, cohortName = 'LIHC-C1', logscale = TRUE, capture_size = 35.8,tcga_cohorts="LIHC")->HCC_SNV_maf_C1_mutload;
tcgaCompare(maf = HCC_SNV_maf_C2, cohortName = 'LIHC-C2', logscale = TRUE, capture_size = 35.8,tcga_cohorts="LIHC")->HCC_SNV_maf_C2_mutload;
tcgaCompare(maf = HCC_SNV_maf_C3, cohortName = 'LIHC-C3', logscale = TRUE, capture_size = 35.8,tcga_cohorts="LIHC")->HCC_SNV_maf_C3_mutload;
#: boxplot
HCC_SNV_maf_C1_mutload$mutation_burden_perSample[HCC_SNV_maf_C1_mutload$mutation_burden_perSample$cohort=="LIHC-C1",]->test_1;
HCC_SNV_maf_C2_mutload$mutation_burden_perSample[HCC_SNV_maf_C2_mutload$mutation_burden_perSample$cohort=="LIHC-C2",]->test_2
HCC_SNV_maf_C3_mutload$mutation_burden_perSample[HCC_SNV_maf_C3_mutload$mutation_burden_perSample$cohort=="LIHC-C3",]->test_3;
rbind(test_1,test_2)->HCC_SNV_maf_C123_mutload
rbind(HCC_SNV_maf_C123_mutload,test_3)->HCC_SNV_maf_C123_mutload;

#：2. enriched mutations
clinicalEnrichment(maf = HCC_SNV_maf, clinicalFeature = 'ExpCluster')->HCC_SNV_expCluster_enrichment;
plotEnrichmentResults(enrich_res = HCC_SNV_expCluster_enrichment, pVal = 0.05, geneFontSize = 0.5, annoFontSize = 0.6)
#: oncoplot
unique(HCC_SNV_expCluster_enrichment$groupwise_comparision$Hugo_Symbol[which(HCC_SNV_expCluster_enrichment$groupwise_comparision$p_value<0.05)])->enriched_genes;
oncoplot(maf = HCC_SNV_maf, genes = enriched_genes, clinicalFeatures = 'ExpCluster',sortByAnnotation=T,fontSize=0.5,bgCol="white")
#：3. mutation Signatures
library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
library('NMF')
#---for subgroup C1
trinucleotideMatrix(maf = HCC_SNV_maf_C1,ref_genome = "BSgenome.Hsapiens.UCSC.hg38")->HCC_SNV_maf_C1_tnm
plotApobecDiff(tnm = HCC_SNV_maf_C1_tnm, maf = HCC_SNV_maf_C1, pVal = 0.2)
#：step 1：ntry
estimateSignatures(mat = HCC_SNV_maf_C1_tnm, nTry = 6)->test_;
plotCophenetic(res = test_)
#：step 2：
extractSignatures(mat = HCC_SNV_maf_C1_tnm, n = 4)->test_sig;
#：step 3：Compate against updated version3 60 signatures 
compareSignatures(nmfRes = test_sig, sig_db = "SBS")->test_sig_cosm;
#：step 4: heatmap 
pheatmap::pheatmap(mat = test_sig_cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures C1")
#----------------
#-- for subgroup C2 
trinucleotideMatrix(maf = HCC_SNV_maf_C2,ref_genome = "BSgenome.Hsapiens.UCSC.hg38")->HCC_SNV_maf_C2_tnm
plotApobecDiff(tnm = HCC_SNV_maf_C2_tnm, maf = HCC_SNV_maf_C2, pVal = 0.2)
estimateSignatures(mat = HCC_SNV_maf_C2_tnm, nTry = 6)->test_;
plotCophenetic(res = test_)
#：step 1：
extractSignatures(mat = HCC_SNV_maf_C2_tnm, n = 3)->test_sig;
#：step 2：Compate against updated version3 60 signatures 
compareSignatures(nmfRes = test_sig, sig_db = "SBS")->test_sig_cosm;
#：step 3：
pheatmap::pheatmap(mat = test_sig_cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures C2")
#----------------
#-- for subgroup C3
trinucleotideMatrix(maf = HCC_SNV_maf_C3,ref_genome = "BSgenome.Hsapiens.UCSC.hg38")->HCC_SNV_maf_C3_tnm
plotApobecDiff(tnm = HCC_SNV_maf_C3_tnm, maf = HCC_SNV_maf_C3, pVal = 0.2) #无显著富集的基因
estimateSignatures(mat = HCC_SNV_maf_C3_tnm, nTry = 6)->test_;
plotCophenetic(res = test_)
#：step 1：
extractSignatures(mat = HCC_SNV_maf_C3_tnm, n = 3)->test_sig;
#：step 2：Compate against updated version3 60 signatures 
compareSignatures(nmfRes = test_sig, sig_db = "SBS")->test_sig_cosm;
#：step 3：
pheatmap::pheatmap(mat = test_sig_cosm$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures C3")

#--------------------------------------------
#：4、tumor microenvironment
#: 
c("Leukocyte.Fraction","Stromal.Fraction","Intratumor.Heterogeneity","Proliferation","Wound.Healing","Macrophage.Regulation","Lymphocyte.Infiltration.Signature.Score","IFN.gamma.Response","TGF.beta.Response","SNV.Neoantigens","Indel.Neoantigens","Silent.Mutation.Rate","Nonsilent.Mutation.Rate","Number.of.Segments","Fraction.Altered","Aneuploidy.Score","Homologous.Recombination.Defects","BCR.Evenness","BCR.Shannon","BCR.Richness","TCR.Shannon","TCR.Richness","TCR.Evenness","CTA.Score")->immune_signatures;
c("Th1.Cells","Th2.Cells","Th17.Cells","B.Cells.Memory","B.Cells.Naive","Dendritic.Cells.Activated","Dendritic.Cells.Resting","Eosinophils","Macrophages.M0","Macrophages.M1","Macrophages.M2","Mast.Cells.Activated","Mast.Cells.Resting","Monocytes","Neutrophils","NK.Cells.Activated","NK.Cells.Resting","Plasma.Cells","T.Cells.CD4.Memory.Activated","T.Cells.CD4.Memory.Resting","T.Cells.CD4.Naive","T.Cells.CD8","T.Cells.Follicular.Helper","T.Cells.gamma.delta","T.Cells.Regulatory.Tregs","Lymphocytes")->LM22_immune_cells;
#: LM22_immune_cells
melt(TCGA_HCC_exp_immuneland,id.vars="ExpCluster",variable.name="Immune_cells",value.name="Score",measure.vars=LM22_immune_cells)->test_;
ggboxplot(test_,x="ExpCluster",y="Score",fill="ExpCluster",xlab="")+facet_wrap(~Immune_cells,scales="free")+stat_compare_means()+theme(legend.position="n")
#：immune_signatures
melt(TCGA_HCC_exp_immuneland,id.vars="ExpCluster",variable.name="immune_signatures",value.name="Score",measure.vars=immune_signatures)->test_;
ggboxplot(test_,x="ExpCluster",y="Score",fill="ExpCluster",xlab="")+facet_wrap(~immune_signatures,scales="free")+stat_compare_means()+theme(legend.position="n")
############################################################################################################################
############################################################################################################################
#:----------------single cells of HCC

#--------------
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
#---------------------------------------------------------------------------
#-------------- remove batch effects
# cells >= 200
names(which(table(GSE151530_seurat_obj$S_ID)>=200))->keep_samples;
subset(GSE151530_seurat_obj,S_ID%in%keep_samples)->GSE151530_seurat_obj;
#：split 
SplitObject(GSE151530_seurat_obj, split.by = "S_ID")->GSE151530_HCC_list;

#step 1. normalization and the top 2000 most variable genes
lapply(X = GSE151530_HCC_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})->GSE151530_HCC_list;
#step 2. Selecting highly variable genes across datasets
SelectIntegrationFeatures(object.list = GSE151530_HCC_list)->GSE151530_HCC_seleted_features;
#step 3. anchor identification
FindIntegrationAnchors(object.list = GSE151530_HCC_list, anchor.features = GSE151530_HCC_seleted_features)->GSE151530_HCC_anchors;
#step 4. integration
IntegrateData(anchorset = GSE151530_HCC_anchors)->GSE151530_HCC_combined
#----------
#：standard analysis for integrated dataset 
"integrated"->DefaultAssay(GSE151530_HCC_combined);
# Run the standard workflow for visualization and clustering
ScaleData(GSE151530_HCC_combined, verbose = FALSE)->GSE151530_HCC_combined
RunPCA(GSE151530_HCC_combined, npcs = 30, verbose = FALSE)->GSE151530_HCC_combined;
RunUMAP(GSE151530_HCC_combined, reduction = "pca", dims = 1:30)->GSE151530_HCC_combined
FindNeighbors(GSE151530_HCC_combined, reduction = "pca", dims = 1:30)->GSE151530_HCC_combined
FindClusters(GSE151530_HCC_combined, resolution = 0.5)->GSE151530_HCC_combined
#-
#: Retinol scores: the average expression values of four target genes
"RNA"->DefaultAssay(GSE151530_HCC_combined);
GetAssayData(GSE151530_HCC_combined,slot="data")->test_;
as.data.frame(test_)->test_;
apply(test_[keep_target_genes,],2,mean)->test_;
AddMetaData(GSE151530_HCC_combined,metadata=test_,col.name="Retinol_score")->GSE151530_HCC_combined;
#----
FeaturePlot(GSE151530_HCC_combined, features = c(keep_target_genes,"Retinol_score"),cols=fill_colors)& theme(plot.title = element_text(size=10))
#--
DimPlot(subset(GSE151530_HCC_combined,Type!="unclassified"), reduction = "umap", group.by = "Type",cols=myd_colors)+xlab("UMAP 1")+ylab("UMAP 2")+ggtitle("Cell type")
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(150),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(350))->fill_colors;
FeaturePlot(GSE151530_HCC_combined, features = "Retinol_score",cols=fill_colors)+xlab("UMAP 1")+ylab("UMAP 2")+ggtitle("Retinol score")
#:
VlnPlot(subset(GSE151530_HCC_combined,Type!="unclassified"),features = "Retinol_score",pt.size=0,sort="decreasing",group.by="Type")+theme(plot.title = element_text(size=10))+NoLegend()

#-----------------------------------------------------------------------------------------------------------
#--- re-clustring Malignant cells 
"integrated"->DefaultAssay(GSE151530_HCC_combined);
subset(GSE151530_HCC_combined,Type=="Malignant cells")->GSE151530_HCC_combined_subset;
#:
FindVariableFeatures(GSE151530_HCC_combined_subset, selection.method = "vst", nfeatures = 2000)->GSE151530_HCC_combined_subset;
ScaleData(GSE151530_HCC_combined_subset, verbose = FALSE)->GSE151530_HCC_combined_subset
RunPCA(GSE151530_HCC_combined_subset, npcs = 30, verbose = FALSE)->GSE151530_HCC_combined_subset;
RunUMAP(GSE151530_HCC_combined_subset, reduction = "pca", dims = 1:30)->GSE151530_HCC_combined_subset
FindNeighbors(GSE151530_HCC_combined_subset, reduction = "pca", dims = 1:30)->GSE151530_HCC_combined_subset
FindClusters(GSE151530_HCC_combined_subset, resolution = 0.5)->GSE151530_HCC_combined_subset
#:--plot 
SetIdent(GSE151530_HCC_combined_subset, value = "seurat_clusters")->GSE151530_HCC_combined_subset;
DimPlot(GSE151530_HCC_combined_subset, reduction = "umap", label = TRUE, repel = TRUE)+NoLegend()+xlab("UMAP 1")+ylab("UMAP 2")+ggtitle("Subclusters of malignant cells");
#:
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(250),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(310))->fill_colors;
FeaturePlot(GSE151530_HCC_combined_subset, features = "Retinol_score",label=T,cols=fill_colors)+theme(plot.title = element_text(size=10))+xlab("UMAP 1")+ylab("UMAP 2")+ggtitle("Retinol score")
VlnPlot(GSE151530_HCC_combined_subset,features = "Retinol_score",pt.size=0,sort="decreasing")+theme(plot.title = element_text(size=10))+NoLegend()+ylab("Retinol score")+ggtitle("")

#:------------------------------
#:------------------------------
#：8、differential analysis for Malignant cells
# find markers for every cluster compared to all remaining cells, report only the positive ones
FindAllMarkers(GSE151530_HCC_combined_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)->Malignant_cells_all_markers;
Malignant_cells_all_markers %>%group_by(cluster) %>%slice_max(n =10, order_by = avg_log2FC)->top_markers;
#------------------------#-----------------------------------------------------------------------------------------------
#---DEGs between the Subpopulations of Malignant cells
SetIdent(GSE151530_HCC_combined_subset, value = "seurat_clusters")->GSE151530_HCC_combined_subset;
FindMarkers(GSE151530_HCC_combined_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,ident.1=5,ident.2=3)->Malignant_cells_ident12_vs_ident3_markers;
rownames(Malignant_cells_ident12_vs_ident3_markers)->Malignant_cells_ident12_vs_ident3_markers$gene;
#-----
#: DEGs between the three subgroups in TCGA-HCC dataset
test_TCGA_HCC_exp_factor[,c("A0_Samples","ExpCluster")]->test_group;
c("SampleID","Condition")->colnames(test_group);
test_group$SampleID->rownames(test_group);
test_group[colnames(test_matrix),]->test_group;
subset(test_group,Condition!="C2")->test_group;
TCGA_HCC_exp_TN_matrix[,test_group$SampleID]->test_matrix;
#: comparing DEGs between scRNA-seq and TCGA-HCC 
do_ranktest_diff(test_matrix,test_group,is_pair=FALSE)->TCGA_HCC_exp_retinol_group_ranktest;
TCGA_HCC_exp_retinol_group_ranktest[which(TCGA_HCC_exp_retinol_group_ranktest$FDR<0.05),]->TCGA_HCC_exp_retinol_group_ranktest_filter;
intersect(TCGA_HCC_exp_retinol_group_ranktest_filter$gName,Malignant_cells_ident12_vs_ident3_markers$gene)->TCGA_scRNA_shared_genes;
#-GO enrichment
enrichGO(gene=TCGA_scRNA_shared_genes,OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05)->TCGA_scRNA_shared_genes.GO;
TCGA_scRNA_shared_genes.GO[order(TCGA_scRNA_shared_genes.GO$Count,decreasing=T),]->TCGA_scRNA_shared_genes.GO;
ggbarplot(TCGA_scRNA_shared_genes.GO[1:15,],x="Description",y="Count",orientation="horizontal",fill="qvalue")+gradient_fill(c("red","gray"))+ggtitle("Enriched GO terms for overlapped genes")

#####################################################################################################################################################
#:------------------------------
#：cell types annotation
celldex::MonacoImmuneData()->celldex_ref;
as.SingleCellExperiment(DietSeurat(GSE151530_HCC_combined_subset))->GSE151530_HCC_combined_sc_obj;
# 将注释的细胞类型添加到 seurat 对象中
GSE151530_HCC_sc_celldex_main$pruned.labels->pruned_labels;
rownames(GSE151530_HCC_sc_celldex_main)->names(pruned_labels);
AddMetaData(GSE151530_HCC_combined_subset,pruned_labels,col.name="Monaco_main_pruned_labels")->GSE151530_HCC_combined_subset;
#：
GSE151530_HCC_sc_celldex_fine$pruned.labels->GSE151530_HCC_combined_subset$Monaco_fine_pruned_labels
#:-----------------------
#：using HumanPrimaryCellAtlasData
HumanPrimaryCellAtlasData()->celldex_ref2;
#：
SingleR(test = GSE151530_HCC_combined_sc_obj,assay.type.test = 1,ref = celldex_ref2,labels = celldex_ref2$label.main)->GSE151530_HCC_sc_celldex_main2
SingleR(test = GSE151530_HCC_combined_sc_obj,assay.type.test = 1,ref = celldex_ref2,labels = celldex_ref2$label.fine)->GSE151530_HCC_sc_celldex_fine2;
GSE151530_HCC_sc_celldex_main2$pruned.labels->GSE151530_HCC_combined_subset$HPC_main_pruned_labels;
GSE151530_HCC_sc_celldex_fine2$pruned.labels->GSE151530_HCC_combined_subset$HPC_fine_pruned_labels;
#--heatmap showing the annotated results
library(superheat)
table(GSE151530_HCC_combined_subset$seurat_clusters,GSE151530_HCC_combined_subset$HPC_main_pruned_labels)->test_;
superheat(test_,bottom.label.text.angle = 90,X.text=test_,X.text.size=4,heat.pal=c("white","red"),X.text.col="black",heat.lim = c(0.2, 1),extreme.values.na = FALSE,row.title="Subpopulations",row.title.size=5,column.title="Annotated cell types",column.title.size=5,grid.hline.size=0.6,grid.vline.size=0.6,bottom.label.size=0.6,left.label.size=0.2,grid.vline.col="gray",grid.hline.col="gray")

#------------------------------
#--plot with annotated cell types
names(which(table(GSE151530_HCC_combined_subset$HPC_main_pruned_labels)<10))->other_cells;
"Other"->GSE151530_HCC_combined_subset$HPC_main_pruned_labels[which(GSE151530_HCC_combined_subset$HPC_main_pruned_labels%in%other_cells)]
#
#Idents(GSE151530_HCC_combined_subset);
SetIdent(GSE151530_HCC_combined_subset, value = "HPC_main_pruned_labels")->GSE151530_HCC_combined_subset;
DimPlot(GSE151530_HCC_combined_subset, label = T , repel = T, label.size = 3) + NoLegend()->p1_annotated;
#
SetIdent(GSE151530_HCC_combined_subset, value = "HPC_fine_pruned_labels")->GSE151530_HCC_combined_subset;
DimPlot(GSE151530_HCC_combined_subset, label = T , repel = T, label.size = 3) + NoLegend()->p2_annotated;
#
SetIdent(GSE151530_HCC_combined_subset, value = "HPC_main_pruned_labels")->GSE151530_HCC_combined_subset;
FeaturePlot(GSE151530_HCC_combined_subset, features = "Retinol_score",cols=fill_colors,label=T)& theme(plot.title = element_text(size=10))->p3_annotated
p1_annotated+p3_annotated
#: VlnPlot
SetIdent(GSE151530_HCC_combined_subset, value = "HPC_main_pruned_labels")->GSE151530_HCC_combined_subset;
VlnPlot(GSE151530_HCC_combined_subset, features = "Retinol_score",pt.size=0,sort="decreasing") & theme(plot.title = element_text(size=10))+NoLegend()->p4_annotated# data,scale.data


#:------------------------------------------------------------
#：GO and KEGG enrichment analysis for marker genes of the Subpopulations in Malignant cells
#-------------------------GO+KEGG: subpopulation = 14
unique(Malignant_cells_all_markers$gene[which(Malignant_cells_all_markers$cluster==5)])->filtered_marker_genes;
enrichGO(gene=filtered_marker_genes,OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05)->GO_KEGG.selected_genes.GO1;
#--- subpopulation = 12 
unique(Malignant_cells_all_markers$gene[which(Malignant_cells_all_markers$cluster==14)])->filtered_marker_genes;
enrichGO(gene=filtered_marker_genes,OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05)->GO_KEGG.selected_genes.GO2;
#--- subpopulation = 0
unique(Malignant_cells_all_markers$gene[which(Malignant_cells_all_markers$cluster==9)])->filtered_marker_genes;
enrichGO(gene=filtered_marker_genes,OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05)->GO_KEGG.selected_genes.GO3;
#--- subpopulation = 12
unique(Malignant_cells_all_markers$gene[which(Malignant_cells_all_markers$cluster==12)])->filtered_marker_genes;
enrichGO(gene=filtered_marker_genes,OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05)->GO_KEGG.selected_genes.GO4;
#:plot 
barplot(dropGO(GO_KEGG.selected_genes.GO1,level=c(1,2,3,4,6,7)),showCategory=5,label_format=100)+ ggtitle("GO terms for cluster 5 maker genes")->p1_GO;
barplot(GO_KEGG.selected_genes.GO2, showCategory=5,label_format = 100) + ggtitle("GO terms for cluster 14 maker genes")->p2_GO;
barplot(dropGO(GO_KEGG.selected_genes.GO3,level=c(1:4)), showCategory=5,label_format = 100) + ggtitle("GO terms for cluster 9 maker genes")->p3_GO;
barplot(dropGO(GO_KEGG.selected_genes.GO4,level=c(1:4,5)), showCategory=5,label_format = 100) + ggtitle("GO terms for cluster 12 maker genes")->p4_GO;
(p1_GO|p2_GO)/(p3_GO|p4_GO)
