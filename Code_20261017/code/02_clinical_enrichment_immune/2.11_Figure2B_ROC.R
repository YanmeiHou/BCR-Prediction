# Figure 2B: cohort-specific ROC curves, preserving the original five-model comparison.
# Reuses saved probabilities; no training, sample insertion, score modification or reselection.
suppressPackageStartupMessages({library(pROC);library(ggplot2);library(patchwork)})
out <- project_path('outputs','figure2b')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source_root <- project_path('data','external','historical_models')
input_csv <- data_path('figure2b','Figure2B_predictions.csv')
models <- c('Stepglm[both]+RF','Stepglm[backward]+RF','glmBoost+RF','RF','Lasso+RF')
cohorts <- c('TCGA','DKFZ','Taylor','CancerMap')
if(!file.exists(input_csv)) {
  scores <- read.delim(file.path(source_root,'BCR_ML_4genes_v2/RS_mat.txt'),check.names=FALSE)
  prep <- readRDS(file.path(source_root,'BCR_ML_Publication/prepared_data.RDS'))
  stopifnot(!anyDuplicated(scores$ID), !anyDuplicated(rownames(prep$label)))
  ii <- match(scores$ID,rownames(prep$label));stopifnot(!anyNA(ii))
  labels <- prep$label[ii,,drop=FALSE]
  predictions <- do.call(rbind,lapply(models,function(m){
    column <- make.names(m);stopifnot(sum(names(scores)==column)==1)
    data.frame(sample_id=scores$ID,cohort=labels$Cohort,
      BCR=as.integer(as.character(labels$outcome)),model=m,prediction=scores[[column]],check.names=FALSE)
  }))
  write.csv(predictions,input_csv,row.names=FALSE)
} else predictions <- read.csv(input_csv,check.names=FALSE)
reference_file <- data_path('figure2b','original_AUC_mat.txt')
if(!file.exists(reference_file)) reference_file <- file.path(source_root,'BCR_ML_Publication/AUC_mat.txt')
reference <- read.delim(reference_file,check.names=FALSE)
rownames(reference)<-reference$ID
metrics <- list();curves <- list()
for(cn in cohorts)for(m in models){
  d <- predictions[predictions$cohort==cn & predictions$model==m,]
  stopifnot(nrow(d)>0,!anyNA(d),all(d$BCR %in% 0:1),all(is.finite(d$prediction)))
  rr <- pROC::roc(d$BCR,d$prediction,levels=c(0,1),direction='<',quiet=TRUE)
  auc_value <- as.numeric(pROC::auc(rr));ci <- as.numeric(pROC::ci.auc(rr,method='delong'))
  stored <- reference[m,cn]
  stopifnot(abs(auc_value-stored)<1e-10)
  rank_auc <- (sum(rank(d$prediction)[d$BCR==1])-sum(d$BCR==1)*(sum(d$BCR==1)+1)/2)/(sum(d$BCR==1)*sum(d$BCR==0))
  stopifnot(abs(auc_value-rank_auc)<1e-12)
  metrics[[paste(cn,m)]]<-data.frame(cohort=cn,model=m,n=nrow(d),BCR=sum(d$BCR),nonBCR=sum(d$BCR==0),
      AUC=auc_value,CI_lower=ci[1],CI_upper=ci[3],stored_AUC=stored,difference=auc_value-stored,
      evaluation=if(cn=='TCGA') 'training_apparent' else 'cross_cohort_evaluation')
  curves[[paste(cn,m)]]<-data.frame(cohort=cn,model=m,false_positive_rate=1-rr$specificities,sensitivity=rr$sensitivities)
}
metrics<-do.call(rbind,metrics);curves<-do.call(rbind,curves)
write.csv(metrics,file.path(out,'Figure2B_AUC_all_models.csv'),row.names=FALSE)
write.csv(metrics[metrics$model==models[1],],file.path(out,'Figure2B_AUC_selected_model.csv'),row.names=FALSE)
write.csv(curves,file.path(out,'Figure2B_ROC_coordinates.csv'),row.names=FALSE)
palette<-setNames(c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#FF7F00'),models)
plots<-list()
for(cn in cohorts){
  d<-curves[curves$cohort==cn,];d$model<-factor(d$model,levels=models)
  mm<-metrics[metrics$cohort==cn,];mm<-mm[match(models,mm$model),]
  labels<-setNames(sprintf('%s  AUC=%.3f',mm$model,mm$AUC),models)
  role<-if(cn=='TCGA') 'Training' else 'Validation'
  title<-sprintf('ROC - %s (%s, n=%d)',cn,role,mm$n[1])
  p<-ggplot(d,aes(false_positive_rate,sensitivity,color=model,group=model))+
    geom_abline(intercept=0,slope=1,color='#E8E8E8',linewidth=.35)+
    geom_path(linewidth=.43)+
    scale_color_manual(values=palette,breaks=models,labels=labels,drop=FALSE)+
    scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.2),labels=function(x)sprintf('%.1f',x),expand=expansion(mult=.025))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2),labels=function(x)sprintf('%.1f',x),expand=expansion(mult=.025))+
    coord_fixed()+labs(title=title,x='1 - Specificity',y='Sensitivity',color='Model')+
    theme_bw(base_size=10,base_family='Arial')+
    theme(panel.grid=element_blank(),panel.border=element_rect(color='#888888',linewidth=.35),
      axis.ticks=element_line(color='#888888',linewidth=.3),axis.text=element_text(color='black',size=9),
      axis.text.y=element_text(angle=90,hjust=.5),axis.title=element_text(color='black',size=10),
      plot.title=element_text(hjust=.5,face='bold',color='black',size=10),
      legend.position='inside',legend.position.inside=c(.985,.025),legend.justification=c(1,0),
      legend.title=element_text(size=7.5,color='black'),legend.text=element_text(size=7.2,color='black'),
      legend.background=element_rect(fill='white',color=NA),legend.key=element_rect(fill='white',color=NA),
      legend.key.height=grid::unit(.30,'cm'),legend.key.width=grid::unit(.38,'cm'),
      legend.spacing.y=grid::unit(0,'pt'),legend.margin=margin(2,2,2,2),plot.margin=margin(8,8,6,6))
  plots[[cn]]<-p
  ggsave(file.path(out,paste0('Figure2B_',cn,'.png')),p,width=5.1,height=5.1,units='in',dpi=600,bg='white')
  ggsave(file.path(out,paste0('Figure2B_',cn,'.pdf')),p,width=5.1,height=5.1,units='in',device=cairo_pdf,bg='white')
}
combined<-wrap_plots(plots,nrow=1)
ggsave(file.path(out,'Figure2B_four_cohorts.png'),combined,width=20.4,height=5.1,units='in',dpi=600,bg='white')
ggsave(file.path(out,'Figure2B_four_cohorts.pdf'),combined,width=20.4,height=5.1,units='in',device=cairo_pdf,bg='white')
ggsave(file.path(out,'Figure2B_preview.png'),combined,width=20.4,height=5.1,units='in',dpi=120,bg='white')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
print(metrics[metrics$model==models[1],c('cohort','n','AUC','CI_lower','CI_upper')],row.names=FALSE)
cat('All 20 cohort/model AUCs agree with original AUC_mat.txt.\n')
