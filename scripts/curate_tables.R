# Compact manuscript supplementary tables; full analysis tables remain in the archive.
archive_out<-out
out<-file.path(run,'Supplementary_Tables');dir.create(out,showWarnings=FALSE)
archived<-function(name)read.delim(file.path(archive_out,paste0(name,'.tsv')),check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('NA','NaN','nan',''))
metadata_spec<-file.path(run,'TABLE_S1_SCHEMA.json')
metadata_status<-system2(Sys.which('python3'),c(shQuote(file.path(pkg,'supplementary_metadata.py')),
 '--metadata',shQuote(file.path(archive_out,'sample_metadata.tsv')),
 '--matrix',shQuote(file.path(data,'matrix.tsv')),'--outdir',shQuote(out),
 '--spec',shQuote(metadata_spec)))
if(metadata_status!=0L)stop('Supplementary Table 1 amino-acid state export failed')
source(file.path(pkg,'metadata_workbook.R'))
write_metadata_workbook(file.path(out,'TableS1_sample_metadata.tsv'),metadata_spec,
 file.path(out,'TableS1_sample_metadata.xlsx'))

variants<-archived('ESM2_observed_replacements')
calls<-archived('variant_callability')
for(lineage in c('Nichols','SS14')){
 z<-calls[calls$lineage==lineage,];j<-match(variants$locus,z$locus);stopifnot(!anyNA(j))
 for(field in c('carriers','callable','carrier_percent_callable'))variants[[paste0(lineage,'_',field)]]<-z[[field]][j]
}
stopifnot(nrow(variants)==40,!anyDuplicated(variants$locus))
put(variants,'TableS2_PBP_variants_and_ESM2_scores')

topology<-archived('focal_topology_stability')
topology_fields<-c('gene','codon','primary_control_fubar_posterior','bootstrap_trees_complete','fubar_posterior_gt_090_count','bootstrap_fubar_posterior_min','bootstrap_fubar_posterior_median','bootstrap_fubar_posterior_max')
stopifnot(all(topology_fields%in%names(topology)))
put(topology[,topology_fields,drop=FALSE],'TableS3_selection_topology_sensitivity')
titles<-c('Genome metadata, amino-acid states and mrcA protein types',
 'PBP amino-acid variants and substitution-specific ESM-2 scores',
 'Sensitivity of focal FUBAR results to genome-wide phylogenetic topology')
titles_zh<-c('基因组元数据、氨基酸状态及mrcA蛋白型','PBP氨基酸变异及替换特异性ESM-2评分',
 '重点密码子FUBAR结果对全基因组系统发育拓扑的敏感性')
files<-sort(list.files(out,pattern='^TableS[1-3]_.*[.]tsv$',full.names=TRUE));stopifnot(length(files)==3)
index<-do.call(rbind,lapply(seq_along(files),function(i){
 x<-read.delim(files[i],check.names=FALSE);data.frame(table=i,part=basename(files[i]),
 spreadsheet=if(i==1L)'TableS1_sample_metadata.xlsx' else NA_character_,
 title=titles[i],title_zh=titles_zh[i],rows=nrow(x),columns=ncol(x))
}))
put(index,'TABLE_INDEX')
columns<-do.call(rbind,lapply(files,function(f){x<-read.delim(f,check.names=FALSE);data.frame(table=basename(f),field=names(x),type=vapply(x,function(z)class(z)[1],character(1)))}))
put(columns,'COLUMN_INDEX')
writeLines(c('Three manuscript supplementary tables are listed in TABLE_INDEX.tsv.',
 'Table S1 is supplied as a grouped-header Excel workbook and a flat-header TSV; its footnotes define per-codon callability.',
 'Table S1 retains 42 audited amino-acid positions across mrcA, pbp1, tp47, lytB and pbp2; untyped protein types are distinct from uncallable codons.',
 'Figure values and complete primary codon results are provided in Figure_Source_Data.',
 'Full intermediate analysis tables are retained in Reproducibility_Data/Analysis_Tables.',
 'Table S2 contains 40 observed substitutions; frequencies use variant- and lineage-specific callable denominators.',
 'The full ER-ARD model-fit and mapping results are included in Supplementary Figure 1 Source Data.',
 'Temporal diagnostics and clock parameters/ESS remain in the corresponding figure source data and reproducibility archive.'),file.path(out,'README.txt'))
jsonlite::write_json(list(table_groups=3,table_files=3,spreadsheet_files=1,optional_table_gaps=0,full_analysis_archive='Reproducibility_Data/Analysis_Tables'),file.path(run,'TABLE_EXPORT_STATUS.json'),auto_unbox=TRUE,pretty=TRUE)
cat('Manuscript supplementary tables curated: 3 tables.\n')
