if(.Platform$OS.type=='windows')invisible(Sys.setlocale('LC_CTYPE','.UTF-8'))
CODE_ROOT<-normalizePath('.',winslash='/')
source('code/00_setup/0.1_project_paths.R',encoding='UTF-8')
check_table<-function(actual,reference,tol=1e-8){
 a<-read.csv(actual,check.names=FALSE);b<-read.csv(reference,check.names=FALSE)
 stopifnot(identical(names(a),names(b)),nrow(a)==nrow(b))
 for(n in names(b)){
  if(is.numeric(b[[n]])){
   stopifnot(identical(is.na(a[[n]]),is.na(b[[n]])))
   ok<-is.finite(b[[n]])
   stopifnot(all(abs(a[[n]][ok]-b[[n]][ok]) <= tol*pmax(1,abs(b[[n]][ok]))))
  }else stopifnot(identical(as.character(a[[n]]),as.character(b[[n]])))
 }
 data.frame(actual=actual,reference=reference,rows=nrow(a),status='PASS')
}
checks<-list()
for(n in c('model_performance.csv','bootstrap_corrected_performance.csv','DCA_positive_intervals.csv')){
 checks[[paste0('clinical_',n)]]<-check_table(file.path('outputs/clinical/tables',n),file.path('data/reference/clinical',n))
}
for(n in c('performance_metrics.csv','all_495_patient_predictions.csv','OOB_confusion_matrix_050.csv')){
 checks[[paste0('tcga_',n)]]<-check_table(file.path('outputs/tcga495/tables',n),file.path('data/reference/tcga495',n))
}
for(n in c('outer_performance.csv','outer_fold_performance.csv','final_diagnostic_performance.csv','paired_AUC_differences_conditional.csv')){
 checks[[paste0('nested_',n)]]<-check_table(file.path('outputs/nested_tuning/tables',n),file.path('data/reference/nested_tuning',n))
}
for(n in c('Original_DCA_threshold_data.csv','Original_DCA_benefit_intervals.csv','Original_DCA_model_coefficients.csv')){
 checks[[paste0('dca_',n)]]<-check_table(file.path('outputs/original_pooled_dca',n),file.path('data/reference/original_pooled_dca',n))
}
checks[['figure2b']]<-check_table('outputs/figure2b/Figure2B_AUC_all_models.csv','data/reference/Figure2B_AUC_all_models.csv')
dir.create('validation',showWarnings=FALSE)
write.csv(do.call(rbind,checks),'validation/output_checks.csv',row.names=FALSE)
capture.output(sessionInfo(),file='validation/sessionInfo.txt')
pkgs<-read.csv('DEPENDENCIES.csv')$package
installed<-installed.packages()
idx<-match(pkgs,installed[,'Package'])
write.csv(data.frame(package=pkgs,installed=!is.na(idx),version=unname(installed[idx,'Version']),license=unname(installed[idx,'License']),row.names=NULL),
 'validation/package_versions.csv',row.names=FALSE,na='')
cat('PASS: ',length(checks),' output tables match saved references.\n',sep='')
