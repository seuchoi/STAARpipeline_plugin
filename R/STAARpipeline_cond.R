STAARpipeline_cond<-function(chr,gene_name,gene_name_cond,category=c("plof","plof_ds","missense","disruptive_missense","synonymous","downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"),cond_category=c("plof","plof_ds","missense","disruptive_missense","synonymous","downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"),
genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,method_cond = c("optimal","naive"),QC_label="annotation/filter",variant_type=c("SNV"),geno_missing_imputation=c("mean","minor"),Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,Use_annotation_weights=TRUE,Annotation_name,silent=FALSE){
    
res<-coding_extract(chr=chr,gene_name=gene_name_cond,genofile=genofile,obj_nullmodel=obj_nullmodel,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,silent=silent)

######
varlist0<-subset(res,type==cond_category)
cond_variants<-varlist0[,c("chr","pos","ref","alt")]


##### conditional analyses

if(category %in% c("plof","plof_ds","missense","disruptive_missense","synonymous")){

res<-Gene_Centric_Coding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,obj_nullmodel=obj_nullmodel,known_loci=cond_variants,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,method_cond = method_cond,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)



}else if(category %in% c("downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS")){

res<-Gene_Centric_Noncoding_cond(chr=chr,gene_name=gene_name,category=category,genofile=genofile,obj_nullmodel=obj_nullmodel,known_loci=cond_variants,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,method_cond = method_cond,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)

}else{}

results<-list(result_cond=res,cond_variants=cond_variants)

return(results)
}
