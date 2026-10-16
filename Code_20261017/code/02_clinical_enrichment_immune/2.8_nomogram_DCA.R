# 2.8 Nomogram/DCA audit and reconstruction from original expression/clinical data.
# No duplicated-record demonstration inputs are used.
options(stringsAsFactors=FALSE,width=150)
root<-project_path('outputs','clinical')
for(p in c('Biobase','rms','pROC','ggplot2','patchwork','ragg'))stopifnot(requireNamespace(p,quietly=TRUE))
suppressPackageStartupMessages(library(ggplot2));suppressPackageStartupMessages(library(patchwork));suppressPackageStartupMessages(library(rms))
for(d in c('input','tables','figures','models','audit'))dir.create(file.path(root,d),recursive=TRUE,showWarnings=FALSE)
wr<-function(x,n)write.csv(x,file.path(root,'tables',n),row.names=FALSE,na='')
genes<-c(ENSG00000124721='DNAH8',ENSG00000131747='TOP2A',ENSG00000101938='CHRDL1',ENSG00000154330='PGM5')
gn<-unname(genes); cohorts<-c('TCGA','DKFZ','Taylor','CancerMap'); thresholds<-seq(.01,.50,.01)
source_dir<-Sys.getenv('PCA_RAW_DATA_DIR',unset=project_path('data','external','bulk'))

# 2.8.1 Extract data, retain provenance, audit exclusions and primary-tumor patient units
finput<-project_path('data','derived','clinical','original_four_gene_clinical.csv')
if(!file.exists(finput)) {
  data_list<-lapply(seq_along(cohorts),function(i) {
    f<-file.path(source_dir,paste0(c('TCGA-PRAD','DKFZ','Taylor','CancerMap')[i],'_eSet.RDS'))
    e<-readRDS(f);pd<-Biobase::pData(e);x<-Biobase::exprs(e)
    stopifnot(identical(colnames(x),rownames(pd)))
    xx<-as.data.frame(t(x[names(genes),,drop=FALSE]));names(xx)<-gn
    data.frame(cohort=cohorts[i],sample_id=rownames(pd),patient_id=pd$patient_id,sample_type=pd$sample_type,
      BCR=pd$bcr_status,Age=suppressWarnings(as.numeric(pd$age_at_diagnosis)),PSA=pd$preop_psa,
      Gleason=pd$gleason_score,time_to_bcr=pd$time_to_bcr,xx)
  })
  raw<-do.call(rbind,data_list);write.csv(raw,finput,row.names=FALSE,na='')
}else raw<-read.csv(finput)
raw$key<-paste(raw$cohort,raw$patient_id,sep='::')
wr(do.call(rbind,lapply(split(raw,raw$cohort),function(d)data.frame(cohort=d$cohort[1],
  n_total=nrow(d),n_BCR_known=sum(d$BCR %in% 0:1),Age_missing=sum(is.na(d$Age)),PSA_missing=sum(is.na(d$PSA)),Gleason_missing=sum(is.na(d$Gleason))))),'clinical_availability.csv')
raw$eligibility_reason<-ifelse(!raw$BCR %in% 0:1,'BCR unavailable',ifelse(!raw$sample_type %in% c('Primary','Tumor','Low Tumor'),'Not primary tumor','Eligible'))
raw$eligibility_reason[is.na(raw$patient_id)|raw$patient_id=='']<-'Patient ID unavailable'
d<-raw[raw$eligibility_reason=='Eligible',]
conflict<-names(which(tapply(d$BCR,d$key,function(z)length(unique(z)))>1))
raw$eligibility_reason[raw$key %in% conflict & raw$eligibility_reason=='Eligible']<-'Conflicting BCR in eligible primary-tumor records'
d<-raw[raw$eligibility_reason=='Eligible',]
# Pick a single primary tumor specimen per patient without reference to outcome or model performance.
d<-d[order(match(d$cohort,cohorts),d$patient_id,match(d$sample_type,c('Primary','Tumor','Low Tumor')),d$sample_id),]
dup<-d$sample_id[duplicated(d$key)]
raw$eligibility_reason[paste(raw$cohort,raw$sample_id) %in% paste(d$cohort[duplicated(d$key)],dup)]<-'Additional specimen from same patient'
d<-d[!duplicated(d$key),];rownames(d)<-NULL
stopifnot(!anyDuplicated(d$key),all(is.finite(as.matrix(d[,gn]))))
wr(raw,'sample_inclusion_audit.csv');wr(d,'primary_patient_data.csv')
counts<-do.call(rbind,lapply(split(d,d$cohort),function(z)data.frame(cohort=z$cohort[1],n_primary_patients=nrow(z),
  BCR=sum(z$BCR),complete_seven_predictors=sum(complete.cases(z[,c(gn,'Age','PSA','Gleason')])))))
wr(counts,'analysis_population.csv');print(counts)

# 2.8.2 Common evaluation utilities
nb<-function(y,p,t=thresholds)vapply(t,function(a)mean((p>=a)&y==1)-mean((p>=a)&y==0)*a/(1-a),numeric(1))
auc<-function(y,p)as.numeric(pROC::auc(pROC::roc(y,p,levels=c(0,1),direction='<',quiet=TRUE)))
metric<-function(y,p,model,cohort,evaluation) {
  rr<-pROC::roc(y,p,levels=c(0,1),direction='<',quiet=TRUE);ci<-as.numeric(pROC::ci.auc(rr,method='delong'))
  data.frame(model=model,cohort=cohort,evaluation=evaluation,n=length(y),events=sum(y),AUC=as.numeric(pROC::auc(rr)),
    AUC_CI95_low=ci[1],AUC_CI95_high=ci[3],Brier=mean((p-y)^2))
}
benefit<-function(y,p,model,cohort,evaluation) data.frame(model=model,cohort=cohort,evaluation=evaluation,
  threshold=thresholds,net_benefit=nb(y,p),treat_all=mean(y)-(1-mean(y))*thresholds/(1-thresholds),treat_none=0)
allmetrics<-list();alldca<-list();allpred<-list();allcoef<-list();allcal<-list();fits<-list()
store_evaluation<-function(model,cohort,evaluation,ds,p) {
  id<-paste(model,cohort,evaluation,sep='::')
  allmetrics[[id]]<<-metric(ds$BCR,p,model,cohort,evaluation)
  alldca[[id]]<<-benefit(ds$BCR,p,model,cohort,evaluation)
  allpred[[id]]<<-data.frame(model=model,cohort=cohort,evaluation=evaluation,
    ds[,c('sample_id','patient_id','BCR')],prediction=as.numeric(p))
  breaks<-unique(quantile(p,seq(0,1,.1)));groups<-cut(p,breaks,include.lowest=TRUE,labels=FALSE)
  allcal[[id]]<<-do.call(rbind,lapply(sort(unique(groups)),function(k)data.frame(model=model,cohort=cohort,evaluation=evaluation,
    bin=k,n=sum(groups==k),mean_predicted=mean(p[groups==k]),observed=mean(ds$BCR[groups==k]))))
}

# 2.8.3 Reproduce historical four-gene scripts for auditing, without reusing their claims
hist<-raw[raw$BCR %in% 0:1,];hscale<-as.data.frame(scale(hist[,gn]));hscale$BCR<-hist$BCR
hc<-glm(BCR~DNAH8+TOP2A+CHRDL1+PGM5,hscale,family=binomial())
for(cn in cohorts) {
  ix<-hist$cohort==cn;p<-predict(hc,newdata=hscale[ix,],type='response')
  store_evaluation('Historical_combined_four_gene',cn,'Apparent_pooled_training',hist[ix,],p)
}
fits$Historical_combined_four_gene<-hc

# 2.8.4 Seven-variable clinical nomogram and aligned clinical/gene-only comparators
complete<-complete.cases(d[,c(gn,'Age','PSA','Gleason')]);cc<-d[complete,]
tr<-cc[cc$cohort=='TCGA',];stopifnot(nrow(tr)>100)
center<-colMeans(tr[,gn]);spread<-sapply(tr[,gn],sd)
scale_clinical<-function(z) {z[,gn]<-as.data.frame(scale(z[,gn],center=center,scale=spread));z}
tr<-scale_clinical(tr);cc<-scale_clinical(cc)
forms<-list(GeneClinical=BCR~DNAH8+TOP2A+CHRDL1+PGM5+Age+PSA+Gleason,
  ClinicalOnly=BCR~Age+PSA+Gleason,FourGene_matched=BCR~DNAH8+TOP2A+CHRDL1+PGM5)
boot_summary<-list(); B<-500L
for(nm in names(forms)) {
  fit<-glm(forms[[nm]],tr,family=binomial());fits[[nm]]<-fit
  cf<-summary(fit)$coefficients
  allcoef[[nm]]<-data.frame(model=nm,term=rownames(cf),coefficient=cf[,1],SE=cf[,2],p=cf[,4],OR=exp(cf[,1]),
    OR_CI95_low=exp(cf[,1]-1.96*cf[,2]),OR_CI95_high=exp(cf[,1]+1.96*cf[,2]))
  for(cn in cohorts[cohorts %in% cc$cohort]) {
    ds<-cc[cc$cohort==cn,];p<-predict(fit,newdata=ds,type='response')
    store_evaluation(nm,cn,if(cn=='TCGA')'Apparent_training' else 'Fixed_TCGA_model_evaluation',ds,p)
  }
  # Bootstrap refits only the prespecified formula, not historical gene discovery.
  set.seed(20261015); optimism<-matrix(NA_real_,B,2+length(thresholds));colnames(optimism)<-c('AUC','Brier',paste0('NB_',thresholds))
  for(b in seq_len(B)) {
    ix<-sample.int(nrow(tr),nrow(tr),replace=TRUE);bt<-tr[ix,]
    bf<-glm(forms[[nm]],bt,family=binomial());pb<-predict(bf,bt,type='response');po<-predict(bf,tr,type='response')
    optimism[b,]<-c(auc(bt$BCR,pb)-auc(tr$BCR,po),mean((pb-bt$BCR)^2)-mean((po-tr$BCR)^2),nb(bt$BCR,pb)-nb(tr$BCR,po))
  }
  original<-predict(fit,tr,type='response');correction<-colMeans(optimism)
  boot_summary[[nm]]<-data.frame(model=nm,bootstrap_refits=B,apparent_AUC=auc(tr$BCR,original),
    optimism_corrected_AUC=auc(tr$BCR,original)-correction[1],apparent_Brier=mean((original-tr$BCR)^2),
    optimism_corrected_Brier=mean((original-tr$BCR)^2)-correction[2])
  dc<-benefit(tr$BCR,original,nm,'TCGA','Bootstrap_optimism_corrected')
  dc$net_benefit<-dc$net_benefit-correction[-c(1,2)]
  alldca[[paste0(nm,'_corrected')]]<-dc
  wr(data.frame(replicate=seq_len(B),optimism,check.names=FALSE),paste0('bootstrap_optimism_',nm,'.csv'))
}
wr(data.frame(gene=gn,center=center,scale=spread),'GeneClinical_training_scaling.csv')
saveRDS(list(models=fits,training_ids=tr[,c('cohort','sample_id','patient_id')],center=center,scale=spread),file.path(root,'models','nomogram_models.rds'))

# 2.8.5 Four-gene-only TCGA model for all four cohorts (separate supplementary analysis)
t4<-d[d$cohort=='TCGA',];m4<-colMeans(t4[,gn]);s4<-sapply(t4[,gn],sd)
x4<-d;x4[,gn]<-as.data.frame(scale(x4[,gn],center=m4,scale=s4))
fit4<-glm(BCR~DNAH8+TOP2A+CHRDL1+PGM5,x4[x4$cohort=='TCGA',],family=binomial())
for(cn in cohorts) {
  ds<-x4[x4$cohort==cn,]
  store_evaluation('FourGene_primary',cn,if(cn=='TCGA')'Apparent_training' else 'Fixed_TCGA_model_evaluation',ds,predict(fit4,ds,type='response'))
}
saveRDS(list(model=fit4,center=m4,scale=s4),file.path(root,'models','four_gene_primary_model.rds'))

# 2.8.6 Exact grid ranges; never bridge gaps between favorable thresholds
mt<-do.call(rbind,allmetrics);dc<-do.call(rbind,alldca);pred<-do.call(rbind,allpred);cal<-do.call(rbind,allcal)
dc$better_than_both<-dc$net_benefit>pmax(dc$treat_all,0)+1e-12
ranges<-do.call(rbind,lapply(split(dc,paste(dc$model,dc$cohort,dc$evaluation)),function(z) {
  z<-z[order(z$threshold),];run<-rle(z$better_than_both);ends<-cumsum(run$lengths);starts<-ends-run$lengths+1
  valid<-which(run$values)
  if(!length(valid))return(data.frame(model=z$model[1],cohort=z$cohort[1],evaluation=z$evaluation[1],from=NA,to=NA,n_grid_points=0))
  do.call(rbind,lapply(valid,function(j)data.frame(model=z$model[1],cohort=z$cohort[1],evaluation=z$evaluation[1],
    from=z$threshold[starts[j]],to=z$threshold[ends[j]],n_grid_points=run$lengths[j])))
}))
wr(mt,'model_performance.csv');wr(dc,'DCA_all_thresholds.csv');wr(ranges,'DCA_positive_intervals.csv')
wr(pred,'model_predictions.csv');wr(cal,'calibration_bins.csv');wr(do.call(rbind,allcoef),'model_coefficients.csv')
wr(do.call(rbind,boot_summary),'bootstrap_corrected_performance.csv')
wr(dc[round(dc$threshold*100) %in% c(5,10,20,30,40,50),],'DCA_selected_thresholds.csv')

# 2.8.7 Figures: clinical units retained; TCGA-fitted coefficients held fixed outside TCGA
theme_set(theme_classic(base_family='Arial',base_size=10)+theme(text=element_text(color='black'),
 axis.text=element_text(color='black'),plot.title=element_text(hjust=.5,face='bold'),plot.subtitle=element_text(hjust=.5,size=9),
 strip.background=element_blank(),strip.text=element_text(face='bold',color='black'),legend.title=element_blank(),legend.position='bottom'))
pal<-c(GeneClinical='#D55E00',ClinicalOnly='#0072B2',FourGene_matched='#009E73',FourGene_primary='#D55E00',Historical_combined_four_gene='#D55E00')
plot_dca<-function(sel,name,title,subtitle) {
  z<-dc[sel,];z$cohort<-factor(z$cohort,levels=cohorts)
  ref<-z[!duplicated(paste(z$cohort,z$threshold)),]
  l<-rbind(data.frame(cohort=ref$cohort,threshold=ref$threshold,NB=ref$treat_all,reference='Treat all'),
    data.frame(cohort=ref$cohort,threshold=ref$threshold,NB=0,reference='Treat none'))
  p<-ggplot(z,aes(threshold,net_benefit,color=model))+geom_line(linewidth=.8)+
    geom_line(data=l,aes(y=NB,linetype=reference),inherit.aes=TRUE,color='black',linewidth=.45)+
    scale_color_manual(values=pal,labels=c(GeneClinical='Four genes + clinical',ClinicalOnly='Clinical only',FourGene_matched='Four genes (matched)',
       FourGene_primary='Four genes',Historical_combined_four_gene='Four genes (pooled fit)'))+
    scale_linetype_manual(values=c('Treat all'=2,'Treat none'=3))+
    facet_wrap(~cohort,ncol=2)+scale_x_continuous(breaks=seq(0,.5,.1))+
    labs(title=title,subtitle=subtitle,x='Threshold probability',y='Net benefit',
      caption='Threshold grid: 0.01-0.50 in steps of 0.01. Treat-all and treat-none references are shown.')
  for(ext in c('pdf','png'))ggsave(file.path(root,'figures',paste0(name,'.',ext)),p,width=210,height=if(length(unique(z$cohort))>2)170 else 105,
    units='mm',device=if(ext=='pdf')cairo_pdf else ragg::agg_png,dpi=600)
}
plot_dca(dc$model %in% names(forms) & ((dc$cohort=='TCGA' & dc$evaluation=='Bootstrap_optimism_corrected') |
  (dc$cohort!='TCGA' & dc$evaluation=='Fixed_TCGA_model_evaluation')),'Supplement_SevenVariable_DCA',
  'Decision curve analysis of the clinical nomogram','TCGA: bootstrap-corrected | DKFZ: fixed TCGA model')
plot_dca(dc$model=='FourGene_primary','Fig2D_FourGene_DCA','Four-gene model across cohorts',
  'TCGA: apparent training performance | Other cohorts: fixed TCGA model')
plot_dca(dc$model=='Historical_combined_four_gene','Audit_Historical_Combined_DCA','Historical pooled-model DCA reproduction',
  'All four cohorts contributed to training; these are apparent cohort-specific results')
ct<-cal[cal$model=='GeneClinical',];ct$cohort<-factor(ct$cohort,levels=cohorts)
pc<-ggplot(ct,aes(mean_predicted,observed))+geom_abline(slope=1,intercept=0,linetype=2,color='grey55')+
  geom_line(color='#0072B2',linewidth=.6)+geom_point(aes(size=n),color='#0072B2')+
  facet_wrap(~cohort,ncol=2)+coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  labs(title='Clinical nomogram calibration',subtitle='TCGA apparent calibration | DKFZ fixed-model calibration',x='Mean predicted probability',y='Observed BCR proportion')+guides(size='none')
ggsave(file.path(root,'figures','Supplement_GeneClinical_Calibration.pdf'),pc,width=200,height=105,units='mm',device=cairo_pdf)
ggsave(file.path(root,'figures','Supplement_GeneClinical_Calibration.png'),pc,width=200,height=105,units='mm',dpi=600,device=ragg::agg_png)
dd<-rms::datadist(tr[,c(gn,'Age','PSA','Gleason')]);options(datadist='dd')
nomfit<-rms::lrm(forms$GeneClinical,data=tr,x=TRUE,y=TRUE)
stopifnot(max(abs(as.numeric(predict(nomfit,type='fitted'))-predict(fits$GeneClinical,type='response')))<1e-5)
nom<-rms::nomogram(nomfit,fun=plogis,fun.at=c(.01,.05,.10,.20,.40,.60,.80,.90),funlabel='Probability of recorded BCR',lp=FALSE,
  DNAH8=pretty(range(tr$DNAH8),n=5),TOP2A=pretty(range(tr$TOP2A),n=5),CHRDL1=pretty(range(tr$CHRDL1),n=5),PGM5=pretty(range(tr$PGM5),n=5),
  Age=c(40,50,60,70,80),PSA=c(0,10,20,40,60,80,100),Gleason=6:10)
draw_nom<-function() {par(mar=c(2,1,3,1),family='Arial');plot(nom,xfrac=.23,cex.var=.9,cex.axis=.8,col.grid=c('grey85','grey94'),
  total.fun=function() axis(3,at=tail(nom$PSA$points,1),labels='100',pos=3,cex.axis=.8,mgp=c(0,.4,0),tcl=-.25,padj=0));
 title(main=paste0('BCR nomogram: four genes + clinical variables (TCGA, n=',nrow(tr),')'),cex.main=1);
 mtext('Genes: TCGA-standardized expression; age: years; PSA: ng/mL; Gleason: score.',side=1,line=.5,cex=.7)}
cairo_pdf(file.path(root,'figures','Fig2C_GeneClinical_Nomogram.pdf'),width=9,height=9);draw_nom();dev.off()
ragg::agg_png(file.path(root,'figures','Fig2C_GeneClinical_Nomogram.png'),width=9,height=9,units='in',res=400);draw_nom();dev.off()
saveRDS(list(model=nomfit,nomogram=nom),file.path(root,'models','clinical_nomogram_rms.rds'))
capture.output(sessionInfo(),file=file.path(root,'sessionInfo.txt'))
cat('\nMETRICS\n');print(mt[mt$model=='GeneClinical',]);print(do.call(rbind,boot_summary))
cat('\nDCA RANGES\n');print(ranges[ranges$model=='GeneClinical',])
