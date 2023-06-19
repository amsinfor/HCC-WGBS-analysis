#------------------------------------HCC WGBS 
library(parallel)
library(data.table)
library(RColorBrewer)
library(pheatmap)
library(pROC);
library(glmnet)
library(ggpubr);
library(cowplot)#combine multiple plots
library(reshape2);
library(ggbeeswarm)
library(Rtsne)
library(limma)
library(MASS)

#--
c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(8,"Dark2"),brewer.pal(12,"Set3"),brewer.pal(8,"Accent"),brewer.pal(11,"Spectral")[c(1,4,10,11)])->myd_colors;
get_palette("npg",10)->myd_colors_npg;
#############################################################################################
#----------------------WGBS data for 97 samples : 24 paired + 66 paired + 7 (4 tumor and 3 normal)
#--
HCC_97samples_WGBS_objRef$getConditionGroup("SampleType",c("Normal","Tumor"))->HCC_97samples_group;
HCC_97samples_WGBS_objRef$getGroupMatrix("SampleType",c("Normal","Tumor"))->HCC_97samples_WGBS_matrix;
#--do ranktest for DMCpGs identification
do_ranktest_diff(HCC_97samples_WGBS_matrix[grep("chrY_",rownames(HCC_97samples_WGBS_matrix)),],HCC_97samples_group,is_pair=FALSE)->HCC_97samples_WGBS_diff_ranktest;#
for(cx in rev(shared_chrs)[-1]){
	paste(cx,"_",sep="")->cx;
	print(cx);flush.console();
	do_ranktest_diff(HCC_97samples_WGBS_matrix[grep(cx,rownames(HCC_97samples_WGBS_matrix)),],HCC_97samples_group,is_pair=FALSE)->cx_ttest;# 
	rbind(HCC_97samples_WGBS_diff_ranktest,cx_ttest)->HCC_97samples_WGBS_diff_ranktest;
	print(paste(cx,"is over!"));flush.console();
}

#-------------- MDSplot showing the structure of 97 samples,
as.numeric(HCC_97samples_WGBS_diff_ranktest$AveExpr-HCC_97samples_WGBS_diff_ranktest$B)->HCC_97samples_WGBS_diff_ranktest$DeltaB;
HCC_97samples_WGBS_diff_ranktest[which(HCC_97samples_WGBS_diff_ranktest$FDR<0.05 & abs(log2(HCC_97samples_WGBS_diff_ranktest$logFC))>=log2(2)),]->HCC_97samples_WGBS_diff_ranktest_filter;
draw_by_MDsplot_by_matrix(HCC_97samples_WGBS_matrix[HCC_97samples_WGBS_diff_ranktest_filter$TransID,],HCC_97samples_group)

#---------------------------
summary_CpG_values_by_chr<-function(expd,colData){
	unique(colData$Condition)->colData_groups;
	lapply(colData_groups,function(cx){
		colData$SampleID[which(colData$Condition==cx)]
	})->colData_group_samples;
	#----------------
	expd$Probe->expd_probes;
	expd[,-c(1,2)]->expd;
	#--
	detectCores()->no_cores;
	makeCluster(6)->c1;
	clusterExport(c1,c("expd","colData_group_samples","colData_groups"),envir=environment());
	parLapply(c1,1:nrow(expd),function(i){
		c()->i_group_values;
		for(gx in colData_group_samples){
			as.numeric(expd[i,gx])->gx_values;
			c(i_group_values,mean(gx_values,na.rm=T),sd(gx_values,na.rm=T),quantile(gx_values,probs=0.25,na.rm=T),quantile(gx_values,probs=0.75,na.rm=T))->i_group_values;
		}
		i_group_values;
	})->expd_group_values;
	stopCluster(c1);
	#--
	matrix(unlist(expd_group_values),ncol=8,byrow=T)->expd_group_values;
	paste(rep(colData_groups,each=4),rep(c("mean","sd","Q25","Q75"),length(colData_groups)),sep="_")->colnames(expd_group_values);
	#--
	as.numeric(unlist(lapply(expd_probes,function(cxp){unlist(strsplit(cxp,split="_"))[2]})))->cx_pos;
	unlist(lapply(expd_probes,function(cxp){unlist(strsplit(cxp,split="_"))[1]}))->cx_size;
	#---
	data.frame("Chr"=cx_size,"Position"=cx_pos,expd_group_values)->expd_group_values;
	expd_group_values[order(expd_group_values$Position,decreasing=F),]->expd_group_values;
	#--
	1:nrow(expd_group_values)->expd_group_values$Index;
	return(expd_group_values);
}
#---by chrs
data.frame()->HCC_97samples_CpG_value_summary;
for(cx in shared_chrs){
	HCC_97samples_keep_probes_DATA[grep(paste(cx,"_",sep=""),HCC_97samples_keep_probes_DATA$Probe),]->test_data;
	summary_CpG_values_by_chr(test_data,HCC_97samples_group[-which(HCC_97samples_group$SampleID%in%HCC_methy_outliers),])->test_;
	rbind(HCC_97samples_CpG_value_summary,test_)->HCC_97samples_CpG_value_summary;
	print(paste(cx,"is over!"));flush.console();
}
#--------by bins：
summary_CpG_values_by_bins<-function(expd_sum,window_size=1000,move_step=500){
	makeCluster(4)->c1;
	clusterExport(c1,c("expd_sum"),envir=environment());
	parLapply(c1,shared_chrs[-c(23,24)],function(cx){
		which(expd_sum$Chr==cx)->cx_index;
		seq(1,length(cx_index),move_step)->cx_steps;
		c()->cx_values;
		for(i in cx_steps){
			ifelse(i+window_size>length(cx_index),length(cx_index),i+window_size)->i_end;
			cx_index[i:(i_end-1)]->cx_i_index;
			mean(expd_sum$Normal_mean[cx_i_index])->cx_i_n_values;
			mean(expd_sum$Tumor_mean[cx_i_index])->cx_i_t_values;
			c(cx_values,cx_i_n_values,cx_i_t_values)->cx_values;
		}
		cx_values;
	})->expd_chr_values;
	stopCluster(c1);
	#----------
	unlist(lapply(expd_chr_values,length))/2->chr_size;
	rep(shared_chrs,chr_size)->chr_size;
	matrix(unlist(expd_chr_values),ncol=2,byrow=T)->expd_chr_values;
	c("Normal_mean","Tumor_mean")->colnames(expd_chr_values);
	data.frame("Chr"=chr_size,expd_chr_values)->expd_sum;
	1:nrow(expd_sum)->expd_sum$Index;
	return(expd_sum);
}
summary_CpG_values_by_bins(HCC_97samples_CpG_value_summary[,c("Chr","Normal_mean","Tumor_mean")],300,100)->HCC_97samples_CpG_value_summary_bins;
HCC_97samples_CpG_value_summary_bins->test_data;
#-----------------plot:
HCC_97samples_WGBS_diff_ranktest_filter$Probe->plot_probes;
unlist(lapply(plot_probes,function(px){unlist(strsplit(px,split="_"))->px_split;paste(px_split[1:2],collapse="_")}))->plot_probes;
#--
plot(test_data$Index,test_data$Normal_mean,type="p",cex=0.1,pch=20,ylim=c(-100,100),col=myd_colors_npg[1],xaxt='n',yaxt='n');
axis(side=2,at=seq(-100,100,20),las=1)
points(test_data$Index,-test_data$Tumor_mean,cex=0.1,pch=20,col=myd_colors_npg[2]);
unlist(lapply(shared_chrs,function(cx){which(test_data$Chr==cx)->cx_index;cx_index[length(cx_index)]}))->test_ablines;
abline(v=test_ablines,lty=2,lwd=0.5,col="black")
text(x=test_ablines,y=100,labels=shared_chrs[-c(23,24)]);
abline(h=c(mean(test_data$Normal_mean),-mean(test_data$Tumor_mean)),lty=1,lwd=1,col="gray")
#----------------------
#----density for the methylation levels of CpGs
melt(HCC_97samples_CpG_value_summary,value.name="Methylation",variable.name="SampleType",id.vars="Index",measure.vars=c("Normal_mean","Tumor_mean"))->test_;
ggdensity(data=test_,x="Methylation",fill="SampleType",palette="jco",size=1,linetype=2,color="SampleType")+scale_x_continuous(breaks=seq(0,100,10));
#-######################################################################################################################
#-- Using the CBS algorithm to identify breakpoints between HCC and NATs
library(PSCBS);
library(iClusterPlus);
#--- HCC/NAT delta values
data.frame()->HCC_97samples_delta_CBS;
for(ix in 1:length(shared_chrs)){
	#--
	shared_chrs[ix]->ix_chr;
	HCC_97samples_CpG_value_summary[,c("Chr","Position","Tumor_mean","Normal_mean")]->test_data;
	(test_data$Tumor_mean-test_data$Normal_mean)->test_data$delta_b;
	test_data[,c("Chr","Position","delta_b")]->test_data;
	test_data[which(test_data$Chr==ix_chr),]->test_data;
	#test_data[which(test_data$Position>58064006 & test_data$Position<58075006),]->test_data
	c("chromosome","x","y")->colnames(test_data);
	ix->test_data$chromosome;
	dropSegmentationOutliers(test_data)->test_data;
	findLargeGaps(test_data, minLength = 1e+04)->gaps;
	gapsToSegments(gaps)->knownSegments;
	#-CBS
	segmentByCBS(test_data, knownSegments = knownSegments,seed = 48879, verbose = -10,min.width=3,joinSegments=FALSE)->fit;
	#-- parse results
	getSegments(fit, simplify = TRUE)->fit_seg;
	"TN"->fit_seg$sampleName
	(fit_seg$end-fit_seg$start)->fit_seg$Region;
	fit_seg[!is.na(fit_seg$mean),]->fit_seg;
	rbind(HCC_97samples_delta_CBS,fit_seg)->HCC_97samples_delta_CBS;
	print(paste0(ix_chr," is over!"));flush.console();
}
#-
"HyperM"->HCC_97samples_delta_CBS$DMR_type;
"HypoM"->HCC_97samples_delta_CBS$DMR_type[which(HCC_97samples_delta_CBS$mean<0)];
#-----------------
c("Sample","Chr","DMR_start","DMR_end","nbrOfLoci","MeanDeltaB","Region","DMR_type")->colnames(HCC_97samples_delta_CBS);
#--- map the regions to genes
change_values(HCC_97samples_delta_CBS,2,shared_chrs)->HCC_97samples_delta_CBS;
map_DMR_to_genome_multiHit(HCC_97samples_delta_CBS,hg38_gtf)->HCC_97samples_delta_CBS_multiMapped;
#-
c("Sample","Chr","Start","End","nbrOfLoci","MeanDeltaB","Region","DMR_type","DMR_location","gName","gType","mappedHit")->colnames(HCC_97samples_delta_CBS_multiMapped);
HCC_97samples_delta_CBS_multiMapped[,c("Chr","Start","End","nbrOfLoci","MeanDeltaB","Region","DMR_type","DMR_location","gName","gType","mappedHit")]->HCC_97samples_delta_CBS_multiMapped;
paste(HCC_97samples_delta_CBS_multiMapped$Chr,HCC_97samples_delta_CBS_multiMapped$Start,HCC_97samples_delta_CBS_multiMapped$End,sep="_")->HCC_97samples_delta_CBS_multiMapped$DMR_id

#####################################################################################################
##############################
#-------------- ROC curve analysis for the four genes' DMRs (HCC vs NATs)
c("chr4_99275265_99275805","chr19_40843446_40843616","chr10_95034460_95061102","chr10_94759925_94760906")->target_genes_DMRs;
prepare_DMR_values(HCC_97samples_keep_probes_DATA, HCC_97samples_delta_CBS_multiMapped,target_genes_DMRs)->target_genes_values;
merge(HCC_97samples_WGBS_objRef$CliniInfo,target_genes_values,by.x="A0_Samples",by.y="SampleID")->target_genes_values
add_Y_labels(target_genes_values,"SampleType",c("Normal","Tumor"))->target_genes_values
#---
par(mfrow=c(2,2))
for(px in target_genes_DMRs){
	do_logistic_fit(target_genes_values,"Y_label",px)->px_logit;
}

#-------------------
"No"->HCC_97samples_delta_CBS_multiMapped$DMR_type;
"HyperM"->HCC_97samples_delta_CBS_multiMapped$DMR_type[which(HCC_97samples_delta_CBS_multiMapped$MeanDeltaB>=10)]
"HypoM"->HCC_97samples_delta_CBS_multiMapped$DMR_type[which(HCC_97samples_delta_CBS_multiMapped$MeanDeltaB<(-15))]
HCC_97samples_delta_CBS_multiMapped[which(HCC_97samples_delta_CBS_multiMapped$DMR_type=="HyperM"),]->HCC_97samples_TN_CBS_merged_mapped_hyperM;
HCC_97samples_delta_CBS_multiMapped[which(HCC_97samples_delta_CBS_multiMapped$DMR_type=="HypoM"),]->HCC_97samples_TN_CBS_merged_mapped_hypoM;

#-DMR length：hyper-DMR and hypo-DMR
ggboxplot(HCC_97samples_delta_CBS_multiMapped,x="DMR_type",y="Region",yscale="log10",color="DMR_type",palette="npg",outlier.shape=2)+stat_compare_means()
#-DMR count：hyper-DMR and hypo-DMR
table(HCC_97samples_delta_CBS_multiMapped$DMR_type)->test_pie;
pie(test_pie,clockwise = TRUE,col = myd_colors[1:3],labels=paste(paste(round(test_pie/sum(test_pie)*100,2),"%",sep=""),test_pie,sep="(n="))
legend("topleft",fill=myd_colors[1:3],legend=c("Hyper-DMR","Hypo-DMR","Not"))

#-----------#----------examples of CBS algorithm
library(trackViewer)
prepare_lolliplot_data<-function(expd,colData,DMR_testd,gene,show_names=TRUE,by_gene=FALSE){
	#---get gene range:
	DMR_testd[which(DMR_testd$gName==gene),]->DMR_testd;
	#merge_two_nearest_CpG_values(DMR_testd)->DMR_testd;
	if(by_gene){
		hg38_gtf$Strand[hg38_gtf$GeneName==gene]->gene_strand;
		if(gene_strand=="+"){
			(hg38_gtf$Start[hg38_gtf$GeneName==gene]-2000)->gene_start;
			(hg38_gtf$End[hg38_gtf$GeneName==gene]+500)->gene_end;
		}else{
			(hg38_gtf$Start[hg38_gtf$GeneName==gene]-5000)->gene_start;
			(hg38_gtf$End[hg38_gtf$GeneName==gene]+2000)->gene_end;
		}
	}else{
		(min(DMR_testd$Start)-100)->gene_start;
		max(DMR_testd$End)+100->gene_end;
	}
	
	#-----------
	expd[grep(paste(gene,"|",sep="",fixed=T),rownames(expd)),]->expd;
	#--排序：
	lapply(rownames(expd),function(px){
		unlist(strsplit(px,split="\\|"))[2]->px;
		as.numeric(unlist(strsplit(px,split="_"))[2])->px_start;
	})->keep_orders;
	order(unlist(keep_orders))->keep_orders;
	expd[keep_orders,]->expd;
	#--
	lapply(rownames(expd),function(px){
		px->px_copy;
		unlist(strsplit(px,split="\\|"))[2]->px;
		as.numeric(unlist(strsplit(px,split="_"))[2])->px_start;
		#---
		if(px_start>=gene_start & px_start<=gene_end){
			px_copy;
		}
	})->keep_rows;
	unlist(keep_rows)->keep_rows;
	expd[keep_rows,]->expd;
	print(dim(expd));flush.console();
	#----calculate mean values between normal and ESCC 
	unique(colData$Condition)->colData_groups;
	unlist(lapply(colData_groups,function(gx){
		as.character(colData$SampleID[which(colData$Condition==gx)])->gx_samples;
		apply(expd[,gx_samples],1,mean,na.rm=T)
	}))->gene_values;
	#matrix(gene_values,ncol=length(colData_groups),byrow=F)->gene_values;
	#colData_groups->colnames(gene_values);
	#--------prepare CpGs positoins;
	unlist(lapply(keep_rows,function(rx){
		unlist(strsplit(rx,split="\\|"))[2]->rx;
		unlist(strsplit(rx,split="_"))[1];
	}))->gene_chrs;
	unlist(lapply(keep_rows,function(rx){
		unlist(strsplit(rx,split="\\|"))[2]->rx;
		as.numeric(unlist(strsplit(rx,split="_"))[2]);
	}))->gene_starts;
	unlist(lapply(keep_rows,function(rx){
		unlist(strsplit(rx,split="\\|"))[2];
	}))->cpg_names;
	unlist(lapply(colData_groups,function(gx){paste(gx,cpg_names,sep="|")}))->cpg_names;
	#---
	if(show_names){
		GRanges(rep(gene_chrs,length(colData_groups)),IRanges(rep(gene_starts,length(colData_groups)),width=1, names=cpg_names))->trackViewer_gr1;
	}else{
		GRanges(rep(gene_chrs,length(colData_groups)),IRanges(rep(gene_starts,length(colData_groups)),width=1))->trackViewer_gr1;
	}
	gene_values->trackViewer_gr1$score;
	#GRanges(gene_chrs,IRanges(gene_starts,width=1, names=cpg_names))->trackViewer_gr2;
	#gene_values[,2]->trackViewer_gr2$score;
	#--prepare range:
	GRanges(gene_chrs[1],IRanges(min(gene_starts),max(gene_starts)))->trackViewer_range;
	#----------prepare DMR:
	unique(DMR_testd[,c("DMR_id","Start","End","Chr")])->trackViewer_features;
	GRanges(trackViewer_features$Chr, IRanges(start=trackViewer_features$Start,end=trackViewer_features$End,names=trackViewer_features$DMR_id))->trackViewer_features
	#---

	list(trackViewer_gr1,trackViewer_features,trackViewer_range)->res;
	c("gr","feature","range")->names(res);
	return(res);
}
#--
HCC_97samples_group[-which(HCC_97samples_group$SampleID%in%HCC_methy_outliers),]->test_group;
HCC_97samples_WGBS_matrix[,test_group$SampleID]->test_matrix;
test_matrix[grep("MGMT\\|",rownames(test_matrix)),]->test_matrix;
prepare_lolliplot_data(test_matrix,test_group,HCC_97samples_delta_CBS_multiMapped[which(HCC_97samples_delta_CBS_multiMapped$gName=="MGMT" & HCC_97samples_delta_CBS_multiMapped$DMR_type!="No"),],"MGMT",show_names=F,by_gene=F)->trackViewer_listd;
"#FF8833"->trackViewer_listd$feature$fill;
0.08->trackViewer_listd$feature$height;
seq(0,40,10)->yaxis;
yaxis->names(yaxis);
length(trackViewer_listd$gr)/2->gr_col_len;
1->trackViewer_listd$gr$alpha
rep(c( "#51C6E6","#FF8833"),c(gr_col_len,gr_col_len))->trackViewer_listd$gr$color;
rep(c("circle", "diamond"),c(gr_col_len,gr_col_len))->trackViewer_listd$gr$shape;
rep(c(0.3,0.3),c(gr_col_len,gr_col_len))->trackViewer_listd$gr$cex;
paste0("legend", as.numeric(factor(trackViewer_listd$gr$shape)))->trackViewer_listd$gr$legend;
"DMR"->trackViewer_listd$feature$feature;
"DMR"->trackViewer_listd$feature$featureLayerID;
#---prepare gene_features
prepare_gene_gtf_annotation<-function(gname){
	read.table("F:/dataset/hg/gencode.v40.annotation_parsed_gtf.txt",header=T,sep="\t")->hg38_gtf_for_annotation;
	hg38_gtf_for_annotation[which(hg38_gtf_for_annotation$Element%in%c("exon","UTR")),]->hg38_gtf_for_annotation;
	#-
	hg38_gtf_for_annotation[hg38_gtf_for_annotation$GeneName==gname,]->gtf_subset;
	gtf_subset[order(gtf_subset$Start,decreasing=F),]->gtf_subset;
	GRanges(seqnames=gtf_subset$Chr,ranges=IRanges(start=gtf_subset$Start,end=gtf_subset$End),strand=gtf_subset$Strand)->gene_feature;
	gtf_subset$Element->gene_feature$feature;
	gtf_subset$GeneID->gene_feature$id;
	paste("exon",1:nrow(gtf_subset),sep="_")->gene_feature$exon;
	gtf_subset$TransID->gene_feature$transcript;
	gtf_subset$GeneName->gene_feature$gene
	#-
	"Gene"->gene_feature$featureLayerID;
	#"#51C6E6"->gene_feature$fill;
	0.08->gene_feature$height;
	#--
	unique(gtf_subset$Element)->gene_feature_annots;
	myd_colors_npg[1:length(gene_feature_annots)]->gene_feature_colors;
	gene_feature_annots->names(gene_feature_colors);
	gene_feature_colors[gene_feature$feature]->gene_feature$fill;
	paste(gene_feature$feature,gene_feature$exon,sep="_")->names(gene_feature)
	#--
	return(gene_feature);
}
#-
prepare_gene_gtf_annotation("MGMT")->gene_feature;
#---
gene_feature$id[1]->trackViewer_listd$feature$id;
paste(gene_feature$id[1],length(gene_feature)+1,sep="_")->trackViewer_listd$feature$exon;
gene_feature$transcript[1]->trackViewer_listd$feature$transcript;
gene_feature$gene[1]->trackViewer_listd$feature$gene;
#--plot delta T/N
length(trackViewer_listd$gr)->gr_rows;
trackViewer_listd$gr[1:(gr_rows/2),]->TN_delta_gr;
(trackViewer_listd$gr$score[(gr_rows/2+1):gr_rows]-trackViewer_listd$gr$score[1:(gr_rows/2)])->TN_delta_gr$score;
unlist(lapply(TN_delta_gr$score,function(sx){ifelse(sx>0,"#FF8833","#51C6E6")}))->TN_delta_gr$color;
lolliplot(TN_delta_gr[sample(1:length(TN_delta_gr),100),],c(trackViewer_listd$feature,gene_feature,ranges=trackViewer_listd$range),yaxis=yaxis,cex=1.6)


######################################################################################################################
######################################################################################################################
#：Showing consistency between 450k and WGBS. four genes for example: ADH1A（AP002026.1），CYP2C8,CYP2A6(AC008537.1),CYP2C19(AL583836.1)
#：the coordinates of 450k probes should be liftover from hg19 to hg38 
#-----------------
#：CYP2C19
#
as.data.frame(trackViewer_listd$gr)->test_gr;
length(trackViewer_listd$gr)->gr_rows;
test_gr[1:(gr_rows/2),]->test_gr;
trackViewer_listd$gr$score[(gr_rows/2+1):gr_rows]->test_gr$score_t;
(test_gr$score_t-test_gr$score)->test_gr$delta_score;
TCGA_methyd_TN_ranktest[which(TCGA_methyd_TN_ranktest$gName=="CYP2C19"),]->test_probes;
HCC_97samples_delta_CBS_multiMapped[grep("CYP2C19",HCC_97samples_delta_CBS_multiMapped$gName),]->test_dmr
test_dmr[test_dmr$DMR_type!="No",]->test_dmr;
#-plot
ggplot(test_gr,aes(x=start))+geom_area(aes(y=delta_score,fill="delta_score"),alpha=0.5)+geom_point(data=test_probes,aes(x=Start,y=DeltaB*100))+geom_text(data=test_probes,aes(x=Start,y=DeltaB*100,label=Probe),angle = 90)+geom_tile(data=test_dmr,aes(x=Start,y=MeanDeltaB,width=Region,height=3,fill=DMR_type))+scale_fill_brewer(palette="Set1")
#-----------------
#：CYP2A6
as.data.frame(trackViewer_listd$gr)->test_gr;
length(trackViewer_listd$gr)->gr_rows;
test_gr[1:(gr_rows/2),]->test_gr;
trackViewer_listd$gr$score[(gr_rows/2+1):gr_rows]->test_gr$score_t;
(test_gr$score_t-test_gr$score)->test_gr$delta_score;
TCGA_methyd_TN_ranktest[which(TCGA_methyd_TN_ranktest$gName=="CYP2A6"),]->test_probes;
HCC_97samples_delta_CBS_multiMapped[grep("CYP2A6",HCC_97samples_delta_CBS_multiMapped$gName),]->test_dmr
test_dmr[test_dmr$DMR_type!="No",]->test_dmr;
#-plot 
ggplot(test_gr,aes(x=start))+geom_area(aes(y=delta_score,fill="delta_score"),alpha=0.5)+geom_point(data=test_probes,aes(x=Start,y=DeltaB*100))+geom_text(data=test_probes,aes(x=Start,y=DeltaB*100,label=Probe),angle = 90)+geom_tile(data=test_dmr,aes(x=Start,y=MeanDeltaB,width=Region,height=3,fill=DMR_type))+scale_fill_brewer(palette="Set1")
#-----------------
#：CYP2C8
as.data.frame(trackViewer_listd$gr)->test_gr;
length(trackViewer_listd$gr)->gr_rows;
test_gr[1:(gr_rows/2),]->test_gr;
trackViewer_listd$gr$score[(gr_rows/2+1):gr_rows]->test_gr$score_t;
(test_gr$score_t-test_gr$score)->test_gr$delta_score;
TCGA_methyd_TN_ranktest[which(TCGA_methyd_TN_ranktest$gName=="CYP2C8"),]->test_probes;
HCC_97samples_delta_CBS_multiMapped[grep("CYP2C8",HCC_97samples_delta_CBS_multiMapped$gName),]->test_dmr
test_dmr[test_dmr$DMR_type!="No",]->test_dmr;
#-plot 
ggplot(test_gr,aes(x=start))+geom_area(aes(y=delta_score,fill="delta_score"),alpha=0.5)+geom_point(data=test_probes,aes(x=Start,y=DeltaB*100))+geom_text(data=test_probes,aes(x=Start,y=DeltaB*100,label=Probe),angle = 90)+geom_tile(data=test_dmr,aes(x=Start+2643,y=MeanDeltaB,width=1905,height=3,fill=DMR_type))+scale_fill_brewer(palette="Set1")
#-----------------
#：ADH1A
as.data.frame(trackViewer_listd$gr)->test_gr;
length(trackViewer_listd$gr)->gr_rows;
test_gr[1:(gr_rows/2),]->test_gr;
trackViewer_listd$gr$score[(gr_rows/2+1):gr_rows]->test_gr$score_t;
(test_gr$score_t-test_gr$score)->test_gr$delta_score;
TCGA_methyd_TN_ranktest[which(TCGA_methyd_TN_ranktest$gName=="ADH1A"),]->test_probes;
HCC_97samples_delta_CBS_multiMapped[grep("ADH1A",HCC_97samples_delta_CBS_multiMapped$gName),]->test_dmr
test_dmr[test_dmr$DMR_type!="No",]->test_dmr;
# plot 
ggplot(test_gr,aes(x=start))+geom_area(aes(y=delta_score,fill="delta_score"),alpha=0.5)+geom_point(data=test_probes,aes(x=Start,y=DeltaB*100))+geom_text(data=test_probes,aes(x=Start,y=DeltaB*100,label=Probe),angle = 90)+geom_tile(data=test_dmr,aes(x=Start,y=MeanDeltaB,width=Region,height=3,fill=DMR_type))+scale_fill_brewer(palette="Set1")


######################################################################################################################
######################################################################################################################
#-- Analysis of transcriptome data
#--------------------------------prepare Normal and Tumor matrix
HCC_tpm_objRef$getConditionGroup("SampleType",c("Normal","Tumor"))->HCC_TN_group;
HCC_tpm_objRef$getGroupMatrix("SampleType",c("Normal","Tumor"))->HCC_TN_matrix;
#------do rank test for all genes:
do_ranktest_diff(HCC_TN_matrix,HCC_TN_group[-which(HCC_TN_group$SampleID%in%HCC_expr_outliers),],FALSE)->HCC_TN_tpm_ranktest;
#---------------------------------------------------------------------------------------------
#------DEGs identification
HCC_TN_tpm_ranktest[which(HCC_TN_tpm_ranktest$FDR<0.05),]->HCC_TN_tpm_ranktest_filter;
HCC_TN_tpm_ranktest_filter[which(HCC_TN_tpm_ranktest_filter$logFC!="Inf"),]->HCC_TN_tpm_ranktest_filter;
HCC_TN_tpm_ranktest_filter[abs(HCC_TN_tpm_ranktest_filter$AveExpr+HCC_TN_tpm_ranktest_filter$B)>=10,]->HCC_TN_tpm_ranktest_filter;
HCC_TN_tpm_ranktest_filter[which(HCC_TN_tpm_ranktest_filter$logFC<1),]->HCC_TN_tpm_ranktest_filter_down;
HCC_TN_tpm_ranktest_filter[which(HCC_TN_tpm_ranktest_filter$logFC>1),]->HCC_TN_tpm_ranktest_filter_up;
#--Volcano plot
"Not"->HCC_TN_tpm_ranktest$Filter
"Pass"->HCC_TN_tpm_ranktest$Filter[abs(HCC_TN_tpm_ranktest$AveExpr+HCC_TN_tpm_ranktest$B)>=10]
#:
(-log10(HCC_TN_tpm_ranktest$P.value))->HCC_TN_tpm_ranktest$log10_P;
log2(HCC_TN_tpm_ranktest$logFC)->HCC_TN_tpm_ranktest$log2_FC;
#:
"Not"->HCC_TN_tpm_ranktest$Type;
"Up"->HCC_TN_tpm_ranktest$Type[which(HCC_TN_tpm_ranktest$FDR<0.05 & HCC_TN_tpm_ranktest$logFC>1)];
"Down"->HCC_TN_tpm_ranktest$Type[which(HCC_TN_tpm_ranktest$FDR<0.05 & HCC_TN_tpm_ranktest$logFC<1)]
#:
ggscatter(subset(HCC_TN_tpm_ranktest,Filter=="Pass"),x="log2_FC",y="log10_P",color="Type",palette="jco",size=1.2,xlab="log2 (Fold change)",ylab="-log10 (P value)")
#--: pie for DEGs type
table(HCC_TN_tpm_ranktest$Type[HCC_TN_tpm_ranktest$Filter=="Pass"])->test_;
paste(test_," (",round(test_/sum(test_),4)*100,"%)",sep="")->x;
pie(test_,labels=x,col=myd_colors[1:3])
legend("top",legend=names(test_),fill=myd_colors[1:3])



######################################################################################################################
#：venn plot 
data.frame("DataSet"=rep(c("up_DEGs","down_DEGs","hypoM_genes","hyperM_genes"),c(length(up_DEGs),length(down_DEGs),length(HCC_97samples_TN_CBS_merged_mapped_hypoM_genes),length(HCC_97samples_TN_CBS_merged_mapped_hyperM_genes))),"gName"=c(up_DEGs,down_DEGs,HCC_97samples_TN_CBS_merged_mapped_hypoM_genes,HCC_97samples_TN_CBS_merged_mapped_hyperM_genes),stringsAsFactors=F)->diff_genes.df;
table(diff_genes.df$gName,diff_genes.df$DataSet)->diff_genes.df_table;
vennCounts(diff_genes.df_table)->deseq_list.vennCount;
vennDiagram(deseq_list.vennCount,names=colnames(diff_genes.df_table),cex=1.5,lwd=1.5,circle.col=brewer.pal(9,"Set1"));
#---------------------------------------------------------------------------
#：correlation analysis for CpG methylation and gene expression
calculate_CpG_expression_corr(c("chr19_40842446_40851447"),"CYP2A6",HCC_97samples_keep_probes_DATA,HCC_tpm_objRef$DATA)->CYP2A6_exp_corr_res;
#：plot 
as.numeric(unlist(lapply(CYP2A6_exp_corr_res$Probe,function(px){unlist(strsplit(px,split="_"))[2]})))->CYP2A6_exp_corr_res$Position;
1:nrow(CYP2A6_exp_corr_res)->CYP2A6_exp_corr_res$Index;
"Sig"->CYP2A6_exp_corr_res$SigType;
"Not"->CYP2A6_exp_corr_res$SigType[which(CYP2A6_exp_corr_res$Pvalue>0.05)];
#:
"Pos"->CYP2A6_exp_corr_res$Direction;
"Neg"->CYP2A6_exp_corr_res$Direction[which(CYP2A6_exp_corr_res$Corr<0)];
"Not"->CYP2A6_exp_corr_res$Direction[which(CYP2A6_exp_corr_res$Pvalue>0.05)]
ggscatter(data=CYP2A6_exp_corr_res,x="Position",y="Corr",shape="SigType",color="Direction",size=2,palette=c(get_palette("aaas",k=2)[1],"gray",get_palette("aaas",2)[2]),xlab="CYP2A6")+scale_shape_manual(values=c(16,8))+geom_hline(yintercept=c(0))+scale_x_continuous(breaks=c(40842000,40844000,40846000,40848000,40850000,40852000))->p1_scatter;
#---------------------
#: CYP2C8
HCC_97samples_TN_CBS_merged_mapped_hyperM[grep("CYP2C8",HCC_97samples_TN_CBS_merged_mapped_hyperM$gName),]
calculate_CpG_expression_corr(c("chr10_95036772_95069497"),"CYP2C8",HCC_97samples_keep_probes_DATA,HCC_tpm_objRef$DATA)->CYP2C8_exp_corr_res;
#：plot
as.numeric(unlist(lapply(CYP2C8_exp_corr_res$Probe,function(px){unlist(strsplit(px,split="_"))[2]})))->CYP2C8_exp_corr_res$Position;
1:nrow(CYP2C8_exp_corr_res)->CYP2C8_exp_corr_res$Index;
"Sig"->CYP2C8_exp_corr_res$SigType;
"Not"->CYP2C8_exp_corr_res$SigType[which(CYP2C8_exp_corr_res$Pvalue>0.05)];
#:
"Pos"->CYP2C8_exp_corr_res$Direction;
"Neg"->CYP2C8_exp_corr_res$Direction[which(CYP2C8_exp_corr_res$Corr<0)];
"Not"->CYP2C8_exp_corr_res$Direction[which(CYP2C8_exp_corr_res$Pvalue>0.05)]
ggscatter(data=CYP2C8_exp_corr_res,x="Position",y="Corr",shape="SigType",color="Direction",size=2,palette=c(get_palette("aaas",k=2)[1],"gray",get_palette("aaas",2)[2]),xlab="CYP2C8")+scale_shape_manual(values=c(16,8))+geom_hline(yintercept=c(0))->p2_scatter;
#---------------------
#: CYP2C19
HCC_97samples_TN_CBS_merged_mapped_hyperM[grep("CYP2C19",HCC_97samples_TN_CBS_merged_mapped_hyperM$gName),]
calculate_CpG_expression_corr(c("chr10_94758925_94855547"),"CYP2C19",HCC_97samples_keep_probes_DATA,HCC_tpm_objRef$DATA)->CYP2C19_exp_corr_res;
#：plot
as.numeric(unlist(lapply(CYP2C19_exp_corr_res$Probe,function(px){unlist(strsplit(px,split="_"))[2]})))->CYP2C19_exp_corr_res$Position;
1:nrow(CYP2C19_exp_corr_res)->CYP2C19_exp_corr_res$Index;
"Sig"->CYP2C19_exp_corr_res$SigType;
"Not"->CYP2C19_exp_corr_res$SigType[which(CYP2C19_exp_corr_res$Pvalue>0.05)];
#:
"Pos"->CYP2C19_exp_corr_res$Direction;
"Neg"->CYP2C19_exp_corr_res$Direction[which(CYP2C19_exp_corr_res$Corr<0)];
"Not"->CYP2C19_exp_corr_res$Direction[which(CYP2C19_exp_corr_res$Pvalue>0.05)]
ggscatter(data=CYP2C19_exp_corr_res,x="Position",y="Corr",shape="SigType",color="Direction",size=2,palette=c(get_palette("aaas",k=2)[1],"gray",get_palette("aaas",2)[2]),xlab="CYP2C19")+scale_shape_manual(values=c(16,8))+geom_hline(yintercept=c(0))+scale_x_continuous(breaks=c(94760000,94780000,94800000,94820000,94840000,94860000))->p3_scatter;
#---------------------
#: ADH1A
HCC_97samples_TN_CBS_merged_mapped_hyperM[grep("ADH1A",HCC_97samples_TN_CBS_merged_mapped_hyperM$gName),]
calculate_CpG_expression_corr(c("chr4_99275165_99291003"),"ADH1A",HCC_97samples_keep_probes_DATA,HCC_tpm_objRef$DATA)->ADH1A_exp_corr_res;
#：plot
as.numeric(unlist(lapply(ADH1A_exp_corr_res$Probe,function(px){unlist(strsplit(px,split="_"))[2]})))->ADH1A_exp_corr_res$Position;
1:nrow(ADH1A_exp_corr_res)->ADH1A_exp_corr_res$Index;
"Sig"->ADH1A_exp_corr_res$SigType;
"Not"->ADH1A_exp_corr_res$SigType[which(ADH1A_exp_corr_res$Pvalue>0.05)];
#:
"Pos"->ADH1A_exp_corr_res$Direction;
"Neg"->ADH1A_exp_corr_res$Direction[which(ADH1A_exp_corr_res$Corr<0)];
"Not"->ADH1A_exp_corr_res$Direction[which(ADH1A_exp_corr_res$Pvalue>0.05)]
ggscatter(data=ADH1A_exp_corr_res,x="Position",y="Corr",shape="SigType",color="Direction",size=2,palette=c(get_palette("aaas",k=2)[1],"gray",get_palette("aaas",2)[2]),xlab="ADH1A")+scale_shape_manual(values=c(16,8))+geom_hline(yintercept=c(0))->p4_scatter



######################################################################################################################
######################################################################################################################
#--- GO / KEGG
#-------------------------------------------------------------------------------------------------------------------
library(R.utils);
R.utils::setOption("clusterProfiler.download.method",'auto')#or options(clusterProfiler.download.method = "wininet")
library(clusterProfiler)
library(org.Hs.eg.db);
#---------------------------------------------------------------------------------------
#-------------------------GO+KEGG: downregulated DEGs  
unique(HCC_TN_tpm_ranktest_filter_down$gName)->HCC_TN_tpm_ranktest_filter_down_genes;
enrichGO(gene=HCC_TN_tpm_ranktest_filter_down_genes,OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05)->GO_KEGG.selected_genes.GO;
bitr(HCC_TN_tpm_ranktest_filter_down_genes,fromType="SYMBOL",toType=c("ENSEMBL","ENTREZID"),OrgDb=org.Hs.eg.db)->myd.gName.df;
enrichKEGG(gene=myd.gName.df[,3],organism="hsa",keyType="ncbi-geneid",pvalueCutoff=0.05)->myd.gName.df.KEGG
dotplot(GO_KEGG.selected_genes.GO, showCategory=20,label_format = 100) + ggtitle("GO terms for Down genes")->p1_GO;
dotplot(myd.gName.df.KEGG, showCategory=20,label_format = 100) + ggtitle("KEGG for Down genes")->p1_KEGG;

#-- upregulated DEGs
unique(HCC_TN_tpm_ranktest_filter_up$gName)->HCC_TN_tpm_ranktest_filter_up_genes;
enrichGO(gene=HCC_TN_tpm_ranktest_filter_up_genes,OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05)->GO_KEGG.selected_genes.GO;
bitr(HCC_TN_tpm_ranktest_filter_up_genes,fromType="SYMBOL",toType=c("ENSEMBL","ENTREZID"),OrgDb=org.Hs.eg.db)->myd.gName.df;
enrichKEGG(gene=myd.gName.df[,3],organism="hsa",keyType="ncbi-geneid",pvalueCutoff=0.05)->myd.gName.df.KEGG
dotplot(GO_KEGG.selected_genes.GO, showCategory=20,label_format = 100) + ggtitle("GO terms for Up genes")->p2_GO;
dotplot(myd.gName.df.KEGG, showCategory=20,label_format = 100) + ggtitle("KEGG for Up genes")->p2_KEGG;
#------------------------------------------------------
#-------------------------GO+KEGG: hypo+up DMGs
enrichGO(gene=HCC_97samples_TN_CBS_merged_mapped_hypoM_up_genes,OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05)->GO_KEGG.selected_genes.GO;
write.table(as.data.frame(GO_KEGG.selected_genes.GO),"GO-KEGG/hypoM_up.GOEnrichment.txt",quote=F,row.names=F,sep="\t")
bitr(HCC_97samples_TN_CBS_merged_mapped_hypoM_up_genes,fromType="SYMBOL",toType=c("ENSEMBL","ENTREZID"),OrgDb=org.Hs.eg.db)->myd.gName.df;
enrichKEGG(gene=myd.gName.df[,3],organism="hsa",keyType="ncbi-geneid",pvalueCutoff=0.05)->myd.gName.df.KEGG
write.table(as.data.frame(myd.gName.df.KEGG),"GO-KEGG/hypoM_up.KEGGEnrichment.txt",quote=F,row.names=F,sep="\t")
dotplot(GO_KEGG.selected_genes.GO, showCategory=20,label_format = 100) + ggtitle("GO terms for hypo-up genes")->p5_GO;
dotplot(myd.gName.df.KEGG, showCategory=20,label_format = 100) + ggtitle("KEGG for hyper DMGs")->p5_KEGG;
#------------------------------------------------------
#-------------------------GO+KEGG: hyper+down DMGs
enrichGO(gene=HCC_97samples_TN_CBS_merged_mapped_hyperM_down_genes,OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont="BP",pAdjustMethod="BH",pvalueCutoff=0.05)->GO_KEGG.selected_genes.GO;
write.table(as.data.frame(GO_KEGG.selected_genes.GO),"GO-KEGG/hyper_down.GOEnrichment.txt",quote=F,row.names=F,sep="\t")
bitr(HCC_97samples_TN_CBS_merged_mapped_hyperM_down_genes,fromType="SYMBOL",toType=c("ENSEMBL","ENTREZID"),OrgDb=org.Hs.eg.db)->myd.gName.df;
enrichKEGG(gene=myd.gName.df[,3],organism="hsa",keyType="ncbi-geneid",pvalueCutoff=0.05)->myd.gName.df.KEGG
write.table(as.data.frame(myd.gName.df.KEGG),"GO-KEGG/hyper_down.KEGGEnrichment.txt",quote=F,row.names=F,sep="\t")
dotplot(GO_KEGG.selected_genes.GO, showCategory=20,label_format = 100) + ggtitle("GO terms for hyper-down genes")->p6_GO;
dotplot(myd.gName.df.KEGG, showCategory=20,label_format = 100) + ggtitle("KEGG for hyper-down genes")->p6_KEGG;
#----------------dot plot 
plot_grid(p1_GO, p1_KEGG,rel_widths = c(4,4));
plot_grid(p2_GO,p2_KEGG,rel_widths = c(4,3));
plot_grid(p5_GO,p5_KEGG,rel_widths=c(5.5,4))
plot_grid(p6_GO,p6_KEGG,rel_widths=c(4,4))
