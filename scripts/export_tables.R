#!/usr/bin/env Rscript
a<-commandArgs(TRUE)
if(length(a)!=2)stop('Usage: Rscript export_tables.R RUN_DIR PACKAGE_DIR')
run<-normalizePath(a[1]);pkg<-normalizePath(a[2]);data<-file.path(run,'work/figures/data')
out<-file.path(run,'Reproducibility_Data','Analysis_Tables');dir.create(out,recursive=TRUE,showWarnings=FALSE)
read<-function(name)read.delim(file.path(data,name),check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('NA','NaN','nan',''))
source_table<-function(name)read.delim(file.path(run,'R_export/source_data',paste0(name,'.tsv')),check.names=FALSE,stringsAsFactors=FALSE)
put<-function(x,name)write.table(x,file.path(out,paste0(name,'.tsv')),sep='\t',quote=TRUE,row.names=FALSE,na='NA')
gaps<-list()
optional<-function(input,name){
 if(file.exists(file.path(data,input))){put(read(input),name);return(TRUE)}
 gaps[[length(gaps)+1]]<<-data.frame(table=name,missing_input=input);FALSE
}
meta<-read('metadata.tsv');prov<-read('metadata_provenance.tsv');mat<-read('matrix.tsv');hs<-read('haplotype_samples.tsv')
meta<-meta[order(meta$sample),];j<-match(meta$sample,prov$sample);stopifnot(!anyNA(j))
for(k in c('bioproject','biosample','study'))meta[[k]]<-prov[[k]][j]
meta$included_primary_phylogeny<-TRUE
for(g in c('mrcA','pbp1','pbp2','tp47'))meta[[paste0('included_',g,'_selection_alignment')]]<-meta$sample%in%readLines(file.path(data,paste0(g,'_ids.txt')))
for(m in c('strict','relaxed'))meta[[paste0('included_',m,'_dating')]]<-meta$sample%in%ape::read.nexus(file.path(data,paste0(m,'.tree')))$tip.label
meta$complete_mrcA_sequence<-tolower(as.character(hs$gene_complete_acgt[match(meta$sample,hs$sample)]))=='true'
meta$gene_haplotype<-hs$gene_haplotype[match(meta$sample,hs$sample)]
assignment<-source_table('Fig4C_protein_type_assignments')
meta$mrcA_protein_type<-assignment$display_type[match(meta$sample,assignment$sample)]
put(meta,'sample_metadata');put(read('source_project_crosswalk.tsv'),'source_project_crosswalk')
put(read('coverage.tsv'),'analysis_scope')
keep<-c('sample',names(mat)[sub(':.*','',names(mat))%in%c('mrcA','pbp1','pbp2','tp47','lytB')])
matrix_public<-mat[,keep,drop=FALSE]
calls<-do.call(rbind,lapply(names(matrix_public)[-1],function(k)do.call(rbind,lapply(c('All','Nichols','SS14'),function(l){
 ids<-if(l=='All')meta$sample else meta$sample[meta$lineage==l];z<-mat[match(ids,mat$sample),k]
 data.frame(locus=k,lineage=l,total=length(z),callable=sum(!is.na(z)),carriers=sum(z==1,na.rm=TRUE),noncarriers=sum(z==0,na.rm=TRUE),missing=sum(is.na(z)))
}))))
calls$carrier_percent_callable<-100*calls$carriers/calls$callable
put(calls,'variant_callability');put(matrix_public,'observed_genotypes')
focal<-read('focal.tsv')
focal_fields<-c('locus','gene','site','callable_codon_sequences','observed_codon_types','amino_acid_counts_in_selection_alignment','matrix_indicators','n_mutant','n_wildtype','n_missing','simmap_gain_mean','simmap_gain_q025','simmap_gain_q975','minimum_gains_among_optima','maximum_gains_among_optima','optimal_root_states','delta_AIC_ARD_ER','AIC_evidence','ARD_root_Pmut')
focal<-focal[,names(focal)%in%focal_fields|startsWith(names(focal),'FUBAR_'),drop=FALSE]
for(k in c('simmap_gain_mean','simmap_gain_q025','simmap_gain_q975'))names(focal)[names(focal)==k]<-paste0('archived_full_tree_ASR_',k)
displayed<-source_table('Fig3_simmap_gain_loss_intervals');displayed<-displayed[displayed$event=='Gain',]
for(k in c('mean','low','high'))focal[[paste0('Figure2E_ER_known_tips_gain_',k)]]<-displayed[[k]][match(focal$locus,displayed$locus)]
focal$Figure2E_mapping_scope<-'ER; known-state tips per variant; 100 archived stochastic maps'
put(focal,'recurrence_and_selection_summary')
lineage<-read.delim(file.path(run,'lineage_history/lineage_summary.tsv'),check.names=FALSE)
put(lineage,'A506T_lineage_parsimony')
put(read('er_ard.tsv'),'ER_ARD_model_sensitivity')
codon_selection<-read('selection.tsv')
codon_fields<-c('gene','site','callable_codon_sequences','observed_codon_types','amino_acid_counts_in_selection_alignment','matrix_indicators')
codon_selection<-codon_selection[,names(codon_selection)%in%codon_fields|startsWith(names(codon_selection),'FUBAR_'),drop=FALSE]
put(codon_selection,'all_codon_selection_results')
if(file.exists(file.path(data,'topology_stability.tsv'))){
 t<-read('topology_stability.tsv');t<-t[paste(t$gene,t$codon)%in%c('mrcA 506','mrcA 625','mrcA 708','pbp1 564','tp47 394'),]
 map<-c(gene='gene',codon='codon',primary_control_fubar_posterior='primary_fubar_posterior',
 bootstrap_trees_complete='complete_bootstrap_trees',fubar_posterior_gt_090_count='fubar_gt_09_count',
 bootstrap_fubar_posterior_min='fubar_posterior_positive_min',bootstrap_fubar_posterior_median='fubar_posterior_positive_median',bootstrap_fubar_posterior_max='fubar_posterior_positive_max')
 t<-t[,unname(map),drop=FALSE];names(t)<-names(map);stopifnot(nrow(t)==5);put(t,'focal_topology_stability')
}else gaps[[length(gaps)+1]]<-data.frame(table='focal_topology_stability',missing_input='topology_stability.tsv')
sm<-source_table('Fig6_reported_structure_metrics')
sm$SASA_unit<-'angstrom_squared';sm$SASA_sphere_points<-sm$sasa_n_points;sm$SASA_probe_radius_A<-sm$probe_radius_A
sm<-sm[,c('gene','position','wt_aa','pdb_aa','wt_match','residue_label','sasa','plddt','SASA_unit','SASA_sphere_points','SASA_probe_radius_A')]
stopifnot(nrow(sm)==3L,all(sm$gene=='mrcA'),all(sm$SASA_sphere_points==2880L))
put(sm,'reference_structural_metrics');put(read('sasa_resolution.tsv'),'MrcA_SASA_resolution')
es<-read('esm_observed.tsv');put(es[,c('gene','locus_tag','locus','replacement','position','reference_aa','alternate_aa','carriers','callable_genomes','carrier_percent_callable','ESM2_LLR')],'ESM2_observed_replacements')
put(source_table('Fig5_DNA_to_protein_types'),'protein_type_sequences')
put(source_table('Fig5_display_number_reference_crosswalk'),'protein_type_labels')
put(source_table('Fig5_type_lineages'),'protein_type_lineage_frequencies')
isolates<-read.delim(file.path(pkg,'config/published_isolates.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
j<-match(isolates$accession,assignment$sample);isolates$local_lineage<-assignment$lineage[j];isolates$local_protein_type<-assignment$protein_type[j]
isolates$local_display_type<-assignment$display_type[j];isolates$local_assignment<-assignment$assignment[j]
isolates$interpretation<-'Published focal states and local complete-protein assignments are recorded separately'
isolates<-isolates[,c('strain','accession','published_focal_states','local_lineage','local_protein_type','local_display_type','local_assignment','source_publication','interpretation')]
put(isolates,'published_isolate_crosswalk')
put(read('temporal_signal.tsv'),'temporal_randomization')
optional('temporal_extended_historical.tsv','extended_temporal_validation')
optional('TableS12A_clustered_date_tests.tsv','clustered_date_tests')
optional('TableS12B_subset_stability.tsv','subset_stability')
put(read('clock_rates.tsv'),'clock_rate_summaries')
optional('TableS13B_strict_ESS.tsv','strict_ESS')
optional('TableS13C_relaxed_ESS.tsv','relaxed_ESS')
put(read('dates.tsv'),'cross_clock_event_matching')
put(read('strict_origins.tsv'),'strict_origin_events');put(read('relaxed_origins.tsv'),'relaxed_origin_events')
files<-list.files(out,pattern='[.]tsv$',full.names=TRUE)
index<-do.call(rbind,lapply(files,function(f){x<-read.delim(f,check.names=FALSE);data.frame(part=basename(f),rows=nrow(x),columns=ncol(x))}))
index<-index[order(index$part),]
titles<-read.delim(file.path(pkg,'config/analysis_titles.tsv'),stringsAsFactors=FALSE,check.names=FALSE)
j<-match(index$part,titles$part);stopifnot(!anyNA(j))
index$title<-titles$title[j];index$title_zh<-titles$title_zh[j]
index<-index[,c('part','title','title_zh','rows','columns')];put(index,'TABLE_INDEX')
columns<-do.call(rbind,lapply(files,function(f){x<-read.delim(f,check.names=FALSE);data.frame(table=basename(f),field=names(x),type=vapply(x,function(z)class(z)[1],character(1)))}))
put(columns,'COLUMN_INDEX')
if(length(gaps))write.table(do.call(rbind,gaps),file.path(run,'OPTIONAL_TABLE_GAPS.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
jsonlite::write_json(list(analysis_files=nrow(index),table_files=nrow(index),optional_table_gaps=length(gaps)),file.path(run,'TABLE_EXPORT_STATUS.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(c('Descriptive analysis files archived from this cluster run; these are not numbered manuscript supplementary tables.',
 'TABLE_INDEX.tsv records component names, titles and dimensions; COLUMN_INDEX.tsv records fields and data types.',
 'Variant frequencies use callable genomes; protein-type frequencies use successfully typed genomes.',
 'FUBAR values are posterior probabilities that the nonsynonymous substitution rate exceeds the synonymous rate.',
 'Parsimony ranges, SIMMAP percentiles, posterior HPD intervals and equal-tail intervals retain their source definitions.',
 'The ESM-2 table contains the 40 observed replacements in the four investigated genes.',
 'The structural figure and its quantitative source data contain only MrcA residues A506, M625 and G708; all use 2880 sphere points per atom.',
 'SOURCE_MANIFEST.json and PUBLICATION_OUTPUT_MANIFEST.json link source and exported files to checksums.'),file.path(out,'README.txt'))
cat('Exported',nrow(index),'analysis archive components. Optional missing sources:',length(gaps),'\n')

source(file.path(pkg,'curate_tables.R'))
