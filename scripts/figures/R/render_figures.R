#!/usr/bin/env Rscript
options(warn=1)
FONT <- Sys.getenv("PBP_FIGURE_FONT", "Nimbus Sans")
a <- commandArgs(trailingOnly=TRUE)
filearg <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
PKG <- normalizePath(file.path(dirname(filearg),'..'))
OUT <- if(length(a)) a[1] else file.path(PKG,'preview')
if(dir.exists(OUT)&&length(list.files(OUT))) stop('Use an empty output directory')
dir.create(OUT,recursive=TRUE);OUT<-normalizePath(OUT)
for(d in c('figures','source_data','qa')) dir.create(file.path(OUT,d))
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(dplyr);library(tidyr);library(ggrepel)})
source(file.path(PKG,'R/style.R'));source(file.path(PKG,'qa/panel_alignment.R'));source(file.path(PKG,'qa/measure_layout.R'))
source(file.path(PKG,'R/panel_sources.R'));source(file.path(PKG,'R/original_plot_layers.R'));update_geom_defaults('text',list(family=FONT));update_geom_defaults('label',list(family=FONT));update_geom_defaults(ggrepel::GeomTextRepel,list(family=FONT));theme_set(common_theme()+theme(panel.grid=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank()));set.seed(20261023)
read_tsv <- function(n) read.delim(file.path(PKG,'data',n),check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('NA','NaN','nan',''),quote='"')
write_tsv <- function(x,n) write.table(x,file.path(OUT,'source_data',paste0(n,'.tsv')),sep='\t',quote=TRUE,row.names=FALSE,na='NA')
num<-function(x)suppressWarnings(as.numeric(x));truth<-function(x)tolower(as.character(x))=='true'
label_locus<-function(x){aa<-c(Ala='A',Arg='R',Asn='N',Asp='D',Cys='C',Gln='Q',Glu='E',Gly='G',His='H',Ile='I',Leu='L',Lys='K',Met='M',Phe='F',Pro='P',Ser='S',Thr='T',Trp='W',Tyr='Y',Val='V',Ter='*');s<-sub('.*:p[.]','',x);for(k in names(aa))s<-gsub(k,aa[k],s,fixed=TRUE);s}
FOCAL<-c('mrcA:p.Ala506Thr','mrcA:p.Ala506Val','mrcA:p.Gly708Ser','mrcA:p.Met625Val','pbp1:p.Pro564Leu','tp47:p.Ser394Arg');LAB<-label_locus(FOCAL)
meta<-read_tsv('metadata.tsv');mat<-read_tsv('matrix.tsv');stopifnot(!anyDuplicated(meta$sample),!anyDuplicated(mat$sample),setequal(meta$sample,mat$sample));meta<-meta[match(mat$sample,meta$sample),]
meta$year<-num(meta$year);meta$country<-gsub('_',' ',meta$country)
meta$lineage[is.na(meta$lineage)]<-'Other/unknown';meta$continent[is.na(meta$continent)|!meta$continent%in%names(CONTINENT_COL)]<-'Other/unknown'
for(k in names(mat)[-1]){mat[[k]]<-num(mat[[k]]);stopifnot(all(is.na(mat[[k]])|mat[[k]]%in%c(0,1)))}
stopifnot(nrow(mat)==1861)
focal<-read_tsv('focal.tsv');focal<-focal[match(FOCAL,focal$locus),];focal$label<-LAB
selection<-read_tsv('selection.tsv');stopifnot(nrow(selection)==2594,all(c('gene','site','FUBAR_Ppositive')%in%names(selection)),all(is.finite(selection$FUBAR_Ppositive)),all(selection$FUBAR_Ppositive>=0&selection$FUBAR_Ppositive<=1))
selection<-selection[,names(selection)%in%c('gene','site','callable_codon_sequences','observed_codon_types','amino_acid_counts_in_selection_alignment','matrix_indicators')|startsWith(names(selection),'FUBAR_'),drop=FALSE]
GENES<-c('mrcA','pbp1','pbp2','tp47','lytB');GLEN<-c(mrcA=884,pbp1=624,pbp2=652,tp47=434,lytB=376)
vars<-do.call(rbind,lapply(names(mat)[-1],function(k){z<-mat[[k]];data.frame(locus=k,gene=sub(':.*','',k),label=label_locus(k),n=sum(z==1,na.rm=TRUE),callable=sum(!is.na(z)),missing=sum(is.na(z)))}));vars$freq<-100*vars$n/vars$callable;vars$focal<-vars$locus%in%FOCAL;vars$site<-num(gsub('[^0-9]','',vars$label))
stopifnot(setequal(head(vars$locus[order(-vars$n)][vars$gene[order(-vars$n)]%in%GENES[1:4]],6),FOCAL))
write_tsv(vars,'variant_counts_and_scope')

gene_expr <- function(x,bold=FALSE) {
 if(is.expression(x))return(x)
 as.expression(lapply(as.character(x),function(z){
  if(is.na(z))return(NA_character_)
  parts<-strsplit(gsub('\\b(mrcA|pbp1|pbp2|tp47|lytB)\\b','|GENE:\\1|',z,perl=TRUE),'|',fixed=TRUE)[[1]]
  if(!any(grepl('^GENE:',parts)))return(z)
  args<-lapply(parts,function(v)if(grepl('^GENE:',v))call(if(bold)'bolditalic' else 'italic',sub('^GENE:','',v))else if(bold)call('bold',v)else v)
  as.call(c(list(as.name('paste')),args))
 }))
}
gene_labeller <- function(labels,multi_line=TRUE) {
 out<-lapply(labels,function(x)lapply(gene_expr(x),identity));structure(out,class='list')
}
grid_dashed<-theme(panel.grid.major.x=element_line(colour='grey82',linewidth=.25,linetype='dashed'),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank())

registry<-list();objects<-list()
savefig<-function(p,id,w=280,h=230,main=TRUE){
 idmap<-c(Alt_Fig2_integrated_recurrence='Fig2_lineages_and_recurrence',Fig4_selection='Fig3_selection',Fig5_mrcA_backgrounds='Fig4_mrcA_backgrounds',Supp2_reconstruction_sensitivity='Supp1_reconstruction_sensitivity',Supp3_ESM2='Supp3_ESM2',Supp5_temporal_regressions='Supp4_temporal_regressions',Supp6_structural_context='Supp2_structural_context',Fig6_dual_clock_history='Fig5_time_calibrated_histories')
 if(id%in%names(idmap))id<-unname(idmap[id])
 main<-grepl('^Fig[1-5]_',id)
 only<-Sys.getenv('PBP_ONLY','');if(nzchar(only)&&!id%in%strsplit(only,',',fixed=TRUE)[[1]])return(invisible(NULL))
 heights<-c(Fig1_dataset_and_variants=432,Fig2_lineages_and_recurrence=505,Fig3_selection=110,Fig4_mrcA_backgrounds=330,Fig5_time_calibrated_histories=305)
 if(main){w<-386;h<-unname(heights[id])}
 if(id=='Alt_Fig2_reconstruction_detail'){w<-386;h<-145};if(id=='Alt_Fig2_distribution_only'){w<-386;h<-360}
 tidy_labels<-function(q){
  if(inherits(q,'patchwork'))q$patches$plots<-lapply(q$patches$plots,tidy_labels)
  if(is.character(q$labels$title)&&grepl('^[A-Za-z]  +',q$labels$title)){
   q$labels$tag<-toupper(substr(q$labels$title,1,1));q$labels$title<-sub('^[A-Za-z]  +','',q$labels$title)
  }
  q$labels$subtitle<-NULL
  q$labels$caption<-NULL
  for(nm in c('title','caption','x','y'))if(is.character(q$labels[[nm]])){
   z<-q$labels[[nm]];if(nm=='title')z<-sub('^([a-z])  ',paste0(toupper(substr(z,1,1)),'  '),z)
   if(grepl('\\b(mrcA|pbp1|pbp2|tp47|lytB)\\b',z,perl=TRUE))q$labels[[nm]]<-gene_expr(z,bold=nm=='title')[[1]] else q$labels[[nm]]<-z
  }
  if(inherits(q$facet,c('FacetWrap','FacetGrid')))q$facet$params$labeller<-gene_labeller
  for(sc in q$scales$scales){
   if(is.character(sc$labels)&&any(grepl('\\b(mrcA|pbp1|pbp2|tp47|lytB)\\b',sc$labels,perl=TRUE)))sc$labels<-gene_expr(sc$labels)
   if(inherits(sc,'ScaleDiscrete')&&inherits(sc$labels,'waiver'))sc$labels<-gene_expr
  }
  if(is.null(q$scales$get_scales('x'))&&!is.null(q$mapping$x)){
   z<-tryCatch(rlang::eval_tidy(q$mapping$x,q$data),error=function(e)NULL)
   if(is.character(z)||is.factor(z))q<-q+scale_x_discrete(labels=gene_expr)
  }
  if(is.null(q$scales$get_scales('y'))&&!is.null(q$mapping$y)){
   z<-tryCatch(rlang::eval_tidy(q$mapping$y,q$data),error=function(e)NULL)
   if(is.character(z)||is.factor(z))q<-q+scale_y_discrete(labels=gene_expr)
  };q
 }
 p<-tidy_labels(p)
 p<-if(inherits(p,'patchwork'))p & theme(plot.caption=element_blank(),plot.subtitle=element_blank()) else p+theme(plot.caption=element_blank(),plot.subtitle=element_blank())
 unified<-theme(text=element_text(family=FONT,size=9),plot.title=element_text(size=10,face='bold',hjust=.5),plot.tag=element_text(size=11,face='bold',hjust=1.4,vjust=1),plot.tag.position=c(0,1),plot.margin=margin(9,7,5,14),axis.title=element_text(size=9),axis.text=element_text(size=8),legend.title=element_text(size=8),legend.text=element_text(size=8),strip.text=element_text(size=9,face='bold'),strip.background=element_rect(fill='grey85',colour='grey35',linewidth=.3))
 p<-if(inherits(p,'patchwork'))p & unified else p+unified
 message('RENDER ',id);b<-file.path(OUT,'figures',id)
 if(inherits(p,'patchwork')||inherits(p$facet,c('FacetWrap','FacetGrid')))require_patchwork_panel_alignment(if(inherits(p,'patchwork'))p else wrap_plots(p),manifest_path=file.path(OUT,'qa',paste0(id,'.layout.json')),report_path=file.path(OUT,'qa',paste0(id,'.alignment.json')),width_in=w/25.4,height_in=h/25.4,audit_script=file.path(PKG,'qa/audit_panel_alignment.py'),python=Sys.which('python3'),strict=TRUE)
 atomic_draw<-function(ext,device){target<-paste0(b,ext);tmp<-tempfile(pattern=paste0(id,'_'),tmpdir=dirname(b),fileext=ext);device(tmp);print(p);dev.off();if(!file.exists(tmp)||file.info(tmp)$size==0)stop('Empty export: ',id,ext);if(!file.rename(tmp,target))stop('Cannot finalize export: ',target)}
 atomic_draw('.pdf',function(f)grDevices::cairo_pdf(f,width=w/25.4,height=h/25.4,family=FONT))
 atomic_draw('.svg',function(f)svglite::svglite(f,width=w/25.4,height=h/25.4))
 atomic_draw('.png',function(f)ragg::agg_png(f,width=w,height=h,units='mm',res=300,background='white'))
 export_plot_layers(p,id,OUT)
 registry[[id]]<<-list(width_mm=w,height_mm=h,main=main);objects[[id]]<<-p
}
geo<-read_tsv('map_coordinates.tsv');geo$Geo_Country<-gsub('_',' ',geo$Geo_Country);geo$continent<-meta$continent[match(geo$Geo_Country,meta$country)];stopifnot(sum(geo$Count)==nrow(meta),!anyNA(geo$continent))
world<-ggplot2::map_data('world',wrap=c(-30,330));world<-world[world$region!='Antarctica',]
for(k in intersect(c('longitude','Longitude','long','lon'),names(geo)))geo[[k]]<-ifelse(geo[[k]]< -30,geo[[k]]+360,geo[[k]])
co<-aggregate(list(n=rep(1,nrow(meta))),meta[,c('country','continent')],sum);co<-co[order(match(co$continent,names(CONTINENT_COL)),-co$n),];levelsco<-rev(co$country);co$country<-factor(co$country,levels=levelsco)
cy<-aggregate(list(n=rep(1,sum(is.finite(meta$year)))),meta[is.finite(meta$year),c('country','continent','year')],sum);cy$country<-factor(cy$country,levels=levelsco)
cl<-aggregate(list(n=rep(1,nrow(meta))),meta[,c('country','lineage')],sum);cl$country<-factor(cl$country,levels=levelsco)
pv<-vars[vars$gene%in%GENES & vars$n>0,];pv<-pv[order(match(pv$gene,GENES),-pv$n,pv$site),];pv$key<-paste(pv$gene,pv$label,sep=':');pv$key<-factor(pv$key,levels=rev(pv$key));pv$gene<-factor(pv$gene,levels=GENES)
vt<-do.call(rbind,lapply(seq_len(nrow(pv)),function(i){z<-mat[[pv$locus[i]]];ok<-!is.na(z)&is.finite(meta$year);d<-data.frame(year=meta$year[ok],n=as.integer(z[ok]==1),callable=1);d<-aggregate(cbind(n,callable)~year,d,sum);d$locus<-pv$locus[i];d$key<-pv$key[i];d$gene<-pv$gene[i];d$focal<-pv$focal[i];d}));vt$key<-factor(vt$key,levels=levels(pv$key));vt$gene<-factor(vt$gene,levels=GENES)
write_tsv(cy,'Fig1_country_year');write_tsv(vt,'Fig1_variant_year_denominators');write_tsv(pv,'Fig1_displayed_variants')
oe<-new.env(parent=globalenv());oe$meta<-meta;oe$meta$Geo_Country<-meta$country;oe$meta$Continent<-meta$continent;oe$meta$TPA_Lineage<-meta$lineage
oe$country_order_df<-oe$meta %>% group_by(Continent,Geo_Country) %>% summarise(Count=n(),.groups='drop') %>% arrange(Continent,Count,Geo_Country)
oe$country_levels<-oe$country_order_df$Geo_Country
oe$Continent.cols<-data.frame(Continent=names(CONTINENT_COL),Continent.col=unname(CONTINENT_COL));oe$CONTINENT_DISPLAY_LABELS<-setNames(names(CONTINENT_COL),names(CONTINENT_COL));oe$TPA_Lineage.cols<-data.frame(Lineage=names(LINEAGE_COL),Lineage.col=unname(LINEAGE_COL))
oe$world<-world;oe$map_df<-geo;oe$map_label_df<-geo;oe$map_label_df$label_country<-geo$Geo_Country
pmap<-original_map(oe)+coord_quickmap(xlim=c(-30,330),ylim=c(-58,82),expand=FALSE)+theme(plot.title=element_text(size=10),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())
pmap$layers[[1]]$aes_params$colour<-NA
pmap$layers[[2]]$aes_params$stroke<-0
pmap$layers[[2]]$aes_params$colour<-'transparent'
pmap$layers[[2]]$aes_params$alpha<-NA_real_
pmap$layers[[2]]$aes_params$fill<-'#EC7014B8'
pmap$layers[[3]]$geom_params$min.segment.length<-Inf
pmap$layers[[3]]$geom_params$seed<-20261023
pco<-original_country_count(oe)+guides(fill='none')+theme(panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),strip.background=element_blank(),strip.text.y=element_blank())
pcy<-original_country_time(oe)+labs(size='Genomes (n)')+theme(panel.grid.minor=element_blank(),strip.background=element_blank(),strip.text.y=element_blank())
pcl<-original_country_lineage(oe)+theme(panel.grid.minor=element_blank())
tm<-vars[vars$gene%in%GENES,];tm<-tm[order(match(tm$gene,GENES),-tm$n,tm$locus),];tm$gene<-factor(tm$gene,levels=GENES);tm$axis_label<-factor(paste(tm$gene,tm$label,sep=':'),levels=rev(paste(tm$gene,tm$label,sep=':')));tm$n_carrier<-tm$n;tm$gene_category<-ifelse(tm$gene=='lytB','Peptidoglycan-pathway gene (lytB)','PBP genes');oe$top_mut<-tm
pvc<-original_variant_count(oe)+theme(panel.grid.major.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank())
pvc$layers[[2]]<-geom_text(aes(label=paste0(n_carrier,'/',callable)),hjust=-.10,size=1.9,colour='grey25')
pvc<-pvc+coord_cartesian(xlim=c(0,max(tm$n_carrier)*1.42))+labs(x='Carrier genomes (n)',caption='Labels: carriers / callable genomes')
vl<-do.call(rbind,lapply(seq_len(nrow(tm)),function(i){z<-mat[[tm$locus[i]]];ok<-!is.na(z)&is.finite(meta$year)&z==1;if(any(ok)){d<-as.data.frame(table(meta$year[ok]));names(d)<-c('year','n');d$year<-num(as.character(d$year))}else d<-data.frame(year=NA_real_,n=0L);d$axis_label<-tm$axis_label[i];d$gene<-tm$gene[i];d$focal_status<-ifelse(tm$focal[i],'Focal PBP variants','Other PBP variants');d}));vl$axis_label<-factor(vl$axis_label,levels=levels(tm$axis_label));vl$gene<-factor(vl$gene,levels=GENES);oe$pbp_long<-vl
pvt<-original_variant_time(oe)+labs(y=NULL,size='Carrier genomes (n)',caption='Counts among callable genomes with collection dates')+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.grid.minor=element_blank())
pvt_size<-pvt$scales$get_scales('size');pvt_size$name<-'Carrier genomes (n)'
pmap<-pmap+labs(title='A  Geographic distribution');pco<-pco+labs(title='B  Genome counts');pcy<-pcy+labs(title='C  Sampling distribution');pcl<-pcl+labs(title='D  Lineage composition');pvc<-pvc+labs(title='E  Variant carriage');pvt<-pvt+labs(title='F  Temporal distribution of carriers')
write_tsv(tm,'Fig1_original_style_all_target_indicators');write_tsv(vl,'Fig1_original_style_timeline')
savefig(pmap/((pco+pcy+pcl+plot_layout(widths=c(2.3,3.2,1.6),guides='collect')) & theme(legend.position='right',legend.box='vertical'))/(pvc+pvt+plot_layout(widths=c(1,1.25)))+plot_layout(heights=c(.68,1.05,1.30)),'Fig1_dataset_and_variants',386,432)
coords<-function(tree){tree<-ape::reorder.phylo(tree,'cladewise');nt<-length(tree$tip.label);nn<-nt+tree$Nnode;ord<-tree$edge[tree$edge[,2]<=nt,2];y<-rep(NA_real_,nn);y[ord]<-seq_len(nt);po<-ape::reorder.phylo(tree,'postorder');for(p in unique(po$edge[,1]))y[p]<-mean(y[po$edge[po$edge[,1]==p,2]]);list(tree=tree,x=ape::node.depth.edgelength(tree),y=y,nt=nt,nn=nn)}
tr<-ape::read.tree(file.path(PKG,'data/primary.treefile'));tr<-ape::ladderize(tr);stopifnot(setequal(tr$tip.label,mat$sample));tc<-coords(tr);tr<-tc$tree;ix<-match(tr$tip.label,mat$sample);rr<-tc$x/max(tc$x);ang<-(tc$y-1)/tc$nt*2*pi
edge_poly<-list();for(k in seq_len(nrow(tr$edge))){pa<-tr$edge[k,1];ch<-tr$edge[k,2];aa<-seq(ang[pa],ang[ch],length.out=max(2,ceiling(abs(ang[ch]-ang[pa])*40)));r<-c(rep(rr[pa],length(aa)),rr[ch]);aa<-c(aa,ang[ch]);edge_poly[[k]]<-data.frame(x=r*cos(aa),y=r*sin(aa),edge=k,child=ch)};edges<-do.call(rbind,edge_poly)
lineages<-meta$lineage[ix];sector<-list();labsector<-list()
for(l in c('Nichols','SS14')){tt<-sort(tc$y[seq_len(tc$nt)][lineages==l]);runs<-split(tt,cumsum(c(1,diff(tt)!=1)));for(j in seq_along(runs)){z<-runs[[j]];aa<-seq((min(z)-1.5)/tc$nt*2*pi,(max(z)-.5)/tc$nt*2*pi,length.out=200);sector[[length(sector)+1]]<-data.frame(x=c(0,1.02*cos(aa),0),y=c(0,1.02*sin(aa),0),lineage=l,id=paste(l,j));};big<-runs[[which.max(lengths(runs))]];mid<-(mean(range(big))-1)/tc$nt*2*pi;labsector[[l]]<-data.frame(x=1.14*cos(mid),y=1.14*sin(mid),label=l)}
sector<-do.call(rbind,sector);labsector<-do.call(rbind,labsector)
ct<-do.call(rbind,lapply(seq_along(FOCAL),function(i)do.call(rbind,lapply(c('Nichols','SS14'),function(l){z<-mat[meta$lineage==l,FOCAL[i]];data.frame(locus=FOCAL[i],label=LAB[i],lineage=l,n=sum(z==1,na.rm=TRUE),callable=sum(!is.na(z)))}))))
trees<-lapply(seq_along(FOCAL),function(i){state<-mat[ix,FOCAL[i]];nc<-ntot<-rep(0,tc$nn);nc[seq_len(tc$nt)]<-as.integer(!is.na(state)&state==1);ntot[seq_len(tc$nt)]<-1;po<-ape::reorder.phylo(tr,'postorder');for(p in unique(po$edge[,1])){kids<-po$edge[po$edge[,1]==p,2];nc[p]<-sum(nc[kids]);ntot[p]<-sum(ntot[kids])};ed<-edges;ed$colour<-ifelse(nc[ed$child]==ntot[ed$child],FOCAL_COL[LAB[i]],'#C5C5C5');tips<-data.frame(x=rr[seq_len(tc$nt)]*cos(ang[seq_len(tc$nt)]),y=rr[seq_len(tc$nt)]*sin(ang[seq_len(tc$nt)]),present=!is.na(state)&state==1)
 z<-ct[ct$locus==FOCAL[i],];subt<-paste(paste0(z$lineage,' ',z$n,'/',z$callable),collapse=' | ')
 ggplot()+geom_polygon(data=sector,aes(x,y,group=id,fill=lineage,alpha=lineage),colour=NA)+scale_alpha_manual(values=c(SS14=.02,Nichols=.03),guide='none')+scale_fill_manual(values=LINEAGE_COL,guide='none')+geom_path(data=ed,aes(x,y,group=edge,colour=colour),linewidth=.12)+geom_point(data=tips[tips$present,],aes(x,y),colour=FOCAL_COL[LAB[i]],size=.35)+scale_colour_identity()+geom_text(data=labsector,aes(x,y,label=label),size=2.5,colour='grey25')+coord_equal(xlim=c(-1.28,1.28),ylim=c(-1.26,1.26))+theme_void(base_size=9)+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())+labs(title=paste0(if(i==1)'a  ' else '',LAB[i]),subtitle=subt)+theme(plot.title=element_text(size=10,face='bold',hjust=.5),plot.subtitle=element_text(size=7,hjust=.5),plot.margin=margin(3,3,3,3))})
countryrows<-do.call(rbind,lapply(seq_along(FOCAL),function(i){z<-mat[[FOCAL[i]]];valid<-!is.na(z)&meta$lineage%in%c('SS14','Nichols')&!is.na(meta$country);d<-data.frame(country=meta$country[valid],carrier=as.integer(z[valid]==1),ss14=as.integer(meta$lineage[valid]=='SS14'),n=1);d<-aggregate(cbind(carrier,ss14,n)~country,d,sum);d<-d[d$n>=3,];d$x<-100*d$ss14/d$n;d$y<-100*d$carrier/d$n;d$label<-LAB[i];d}))
cs<-do.call(rbind,lapply(LAB,function(l){d<-countryrows[countryrows$label==l,];data.frame(label=l,n_countries=nrow(d),r=cor(d$x,d$y))}));countryrows$label<-factor(countryrows$label,levels=LAB)
write_tsv(ct,'Fig2_lineage_counts');write_tsv(countryrows,'Fig2_country_matched_denominators');write_tsv(cs,'Fig2_Pearson_r_unweighted')
oe$lineage_corr_long<-countryrows;oe$lineage_corr_long$SS14_pct<-countryrows$x;oe$lineage_corr_long$variant_pct<-countryrows$y;oe$lineage_corr_long$Geo_Country<-countryrows$country;oe$lineage_corr_long$variant<-factor(countryrows$label,levels=LAB);oe$lineage_corr_long$fig2_group<-as.character(countryrows$label)
oe$lineage_corr_stats<-data.frame(variant=factor(cs$label,levels=LAB),fig2_group=cs$label,label=sprintf('r = %+.2f',cs$r))
oe$Fig2.cols<-data.frame(fig2_group=LAB,fig2_col=unname(FOCAL_COL[LAB]))
pcorr<-original_country_correlation(oe)+labs(title='B  Geographical association',size='Callable genomes (n)',x='SS14 among callable genomes (%)',y='Carriage among callable genomes (%)')+theme(panel.grid=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
pcorr_size<-pcorr$scales$get_scales('size');pcorr_size$name<-'Callable genomes (n)'
pcorr$layers[[1]]<-geom_segment(data=data.frame(x=0,y=0,xend=100,yend=100),aes(x=x,y=y,xend=xend,yend=yend),inherit.aes=FALSE,linetype='dashed',colour='grey72',linewidth=.3)
variant_labels<-paste0(LAB,'\nr = ',sprintf('%+.2f',cs$r[match(LAB,cs$label)]))
country_labels<-do.call(rbind,lapply(seq_along(LAB),function(ki){
 k<-LAB[ki];d<-oe$lineage_corr_long[oe$lineage_corr_long$variant==k,];d<-d[order(-nchar(d$Geo_Country),-d$n,d$Geo_Country),];fit<-coef(lm(variant_pct~SS14_pct,data=d));placed<-list()
 for(j in seq_len(nrow(d))){z<-d[j,];w<-nchar(z$Geo_Country)*2.65+1;hh<-2.8
  cand<-expand.grid(x=seq(-3+w/2,103-w/2,by=1.4),y=seq(-1,101,by=1.4));cand$cost<-(cand$x-z$SS14_pct)^2+1.4*(cand$y-z$variant_pct)^2;cand<-cand[order(cand$cost),]
  chosen<-NULL
  for(ii in seq_len(nrow(cand))){xx<-cand$x[ii];yy<-cand$y[ii];xl<-xx-w/2;xr<-xx+w/2;yl<-yy-hh;yr<-yy+hh
   intersects<-function(vals)min(vals)<=yr+.5&&max(vals)>=yl-.5
   if(intersects(c(xl,xr))||intersects(fit[1]+fit[2]*c(xl,xr)))next
   rx<-1+4*sqrt(d$n/max(d$n));ry<-1+4*sqrt(d$n/max(d$n))
   if(any(abs(d$SS14_pct-xx)<w/2+rx & abs(d$variant_pct-yy)<hh+ry))next
   if(length(placed)&&any(vapply(placed,function(p)abs(xx-p$x)<(w+p$w)/2+.7&&abs(yy-p$y)<hh+p$hh+.6,logical(1))))next
   chosen<-data.frame(variant=factor(variant_labels[ki],levels=variant_labels),x=xx,y=yy,label=z$Geo_Country,w=w,hh=hh,px=z$SS14_pct,py=z$variant_pct);break
  }
  if(is.null(chosen))stop('No collision-free country label position: ',k,' ',z$Geo_Country)
  placed[[length(placed)+1]]<-chosen
 };do.call(rbind,placed)
}))
country_leaders<-list()
for(i in seq_len(nrow(country_labels))){z<-country_labels[i,];dx<-z$x-z$px;dy<-z$y-z$py;intervals<-list(c(0,1));boxes<-country_labels[country_labels$variant==z$variant,]
 for(j in seq_len(nrow(boxes))){q<-boxes[j,];lo<-0;hi<-1
  for(ax in 1:2){v0<-if(ax==1)z$px else z$py;dv<-if(ax==1)dx else dy;mn<-if(ax==1)q$x-q$w/2-.4 else q$y-q$hh-.4;mx<-if(ax==1)q$x+q$w/2+.4 else q$y+q$hh+.4
   if(abs(dv)<1e-9){if(v0<mn||v0>mx){hi<- -1;break}}else{rr<-sort(c((mn-v0)/dv,(mx-v0)/dv));lo<-max(lo,rr[1]);hi<-min(hi,rr[2])}
  }
  if(hi<=lo)next
  nextints<-list();for(v in intervals){if(hi<=v[1]||lo>=v[2]){nextints[[length(nextints)+1]]<-v}else{if(lo>v[1])nextints[[length(nextints)+1]]<-c(v[1],lo);if(hi<v[2])nextints[[length(nextints)+1]]<-c(hi,v[2])}};intervals<-nextints
 }
 for(v in intervals)if(diff(v)>1e-4)country_leaders[[length(country_leaders)+1]]<-data.frame(variant=z$variant,x=z$px+v[1]*dx,y=z$py+v[1]*dy,xend=z$px+v[2]*dx,yend=z$py+v[2]*dy)
}
country_leaders<-do.call(rbind,country_leaders)
pcorr$layers[[4]]<-geom_text(data=country_labels,aes(x,y,label=label),inherit.aes=FALSE,size=2.05,colour='grey30')
pcorr$layers[[5]]<-NULL
levels(pcorr$data$variant)<-variant_labels
pcorr$layers<-append(pcorr$layers,list(geom_segment(data=country_leaders,aes(x,y,xend=xend,yend=yend),inherit.aes=FALSE,linewidth=.15,linetype='dashed',colour='grey65')),after=1)
pcorr<-pcorr+guides(fill='none',colour='none')+facet_wrap(~variant,nrow=1)+scale_y_continuous(limits=c(-5,105),breaks=c(0,25,50,75,100))+scale_x_continuous(breaks=c(0,50,100),sec.axis=sec_axis(~100-.,breaks=c(0,50,100),name='Nichols among callable genomes (%)'))+coord_cartesian(xlim=c(-8,108))+theme(axis.text.x.top=element_text(margin=margin(b=5)),axis.title.x.top=element_text(size=9,margin=margin(b=4)),aspect.ratio=1)
stopifnot(nrow(country_labels)==nrow(countryrows))
write_tsv(country_labels,'Fig2_country_label_positions');write_tsv(country_leaders,'Fig2_country_leader_segments')
for(i in seq_along(trees)){z<-ct[ct$locus==FOCAL[i],];trees[[i]]<-trees[[i]]+labs(title=paste0(if(i==1)'A  ' else '',sub(':.*','',FOCAL[i]),' ',LAB[i]),subtitle=NULL,caption=paste0('Carriers / callable: ',paste(paste0(z$lineage,' ',z$n,'/',z$callable),collapse='; ')))+theme(plot.title=element_text(hjust=0),plot.caption=element_text(hjust=.5,size=8),plot.margin=margin(1,1,2,1))}

draws<-read_tsv('simmap_draws.tsv');er<-read_tsv('er_ard.tsv');edraw<-draws[draws$model=='ER',];gainmeans<-aggregate(cbind(gain,loss)~locus,edraw,mean);v3<-merge(vars[vars$gene%in%GENES[1:4],],gainmeans,by='locus',all.x=TRUE);v3$gene<-factor(v3$gene,levels=GENES[1:4]);v3$lane<-1
for(g in levels(v3$gene)){idx<-which(v3$gene==g);for(pos in unique(v3$site[idx])){j<-idx[v3$site[idx]==pos];v3$lane[j]<-if(length(j)==1)0 else seq(-.25,.25,length.out=length(j))}}
v3$lane[v3$locus=='mrcA:p.Ala506Thr']<-.25;v3$lane[v3$locus=='mrcA:p.Ala506Val']<--.25;v3$lane[v3$locus=='mrcA:p.Ala506Met']<-0
v3$lab<-ifelse(v3$focal,v3$label,'');lens<-data.frame(gene=factor(GENES[1:4],levels=GENES[1:4]),start=1,end=GLEN[1:4])
fd<-edraw[edraw$locus%in%FOCAL,];joint<-aggregate(list(n=rep(1,nrow(fd))),fd[,c('locus','gain','loss')],sum);den<-table(fd$locus);joint$percent<-100*joint$n/as.numeric(den[joint$locus]);joint$label<-factor(LAB[match(joint$locus,FOCAL)],levels=LAB)
profile<-v3[is.finite(v3$gain),];profile$plot_group<-ifelse(profile$focal,profile$label,'Other analysed variants')
gq<-do.call(rbind,lapply(split(edraw,edraw$locus),function(d)data.frame(locus=d$locus[1],lo=quantile(d$gain,.025),hi=quantile(d$gain,.975))))
profile<-merge(profile,gq,by='locus');profile$tag<-ifelse(profile$focal,paste(profile$gene,profile$label),'')
md<-rbind(data.frame(locus=fd$locus,event='Gain',count=fd$gain),data.frame(locus=fd$locus,event='Loss',count=fd$loss))
md<-aggregate(list(n=rep(1,nrow(md))),md,sum);md$percent<-md$n;md$variant<-factor(paste(sub(':.*','',md$locus),label_locus(md$locus)),levels=rev(paste(sub(':.*','',FOCAL),LAB)));md$label<-label_locus(md$locus);md$event<-factor(md$event,levels=c('Gain','Loss'))
write_tsv(profile,'Fig3_frequency_and_gain');write_tsv(md,'Fig3_gain_loss_mapping_frequencies');write_tsv(joint,'Fig3_joint_mapping_distribution')
rowlabels<-paste(sub(':.*','',FOCAL),LAB)
f3prev<-ct;f3prev$row<-7-match(f3prev$locus,FOCAL);f3prev$y<-f3prev$row+ifelse(f3prev$lineage=='Nichols',.16,-.16);f3prev$pct<-100*f3prev$n/f3prev$callable
f3prev$text<-sprintf('%d/%d (%.1f%%)',f3prev$n,f3prev$callable,f3prev$pct)
ps<-read_tsv('parsimony.tsv');ps<-ps[ps$tree_scope=='primary_full'&ps$root_policy=='unrestricted'&ps$locus%in%FOCAL,];stopifnot(nrow(ps)==6,!anyDuplicated(ps$locus))
ps$row<-7-match(ps$locus,FOCAL);ps$range<-ifelse(ps$minimum_gains_among_optima==ps$maximum_gains_among_optima,as.character(ps$minimum_gains_among_optima),paste0(ps$minimum_gains_among_optima,'-',ps$maximum_gains_among_optima))
f3s<-do.call(rbind,lapply(split(fd,fd$locus),function(d)do.call(rbind,lapply(c('gain','loss'),function(k)data.frame(locus=d$locus[1],event=ifelse(k=='gain','Gain','Loss'),mean=mean(d[[k]]),low=unname(quantile(d[[k]],.025)),high=unname(quantile(d[[k]],.975)),n_maps=nrow(d))))))
stopifnot(all(f3s$n_maps==100));f3s$row<-7-match(f3s$locus,FOCAL);f3s$y<-f3s$row+ifelse(f3s$event=='Gain',.16,-.16);f3s$text<-sprintf('%.2f [%.2f, %.2f]',f3s$mean,f3s$low,f3s$high)
rowscale<-function(labels=FALSE)scale_y_continuous(limits=c(.5,6.5),breaks=1:6,labels=if(labels)gene_expr(rev(rowlabels))else NULL,expand=c(0,0))
rowlines<-function()geom_hline(yintercept=seq(.5,6.5,1),colour='grey88',linewidth=.2)
p3a<-ggplot(f3prev)+rowlines()+geom_col(aes(x=pct,y=y,fill=lineage),orientation='y',width=.21,alpha=1)+geom_text(aes(x=104,y=y,label=text,colour=lineage),hjust=0,size=2.7,show.legend=FALSE)+scale_fill_manual(values=LINEAGE_COL)+scale_colour_manual(values=LINEAGE_COL,guide='none')+scale_x_continuous(limits=c(0,183),breaks=c(0,50,100),expand=c(0,0))+rowscale(TRUE)+labs(title='A  Lineage-specific carriage',x='Carriage among callable genomes (%)',y=NULL,fill='Lineage',caption='n/N: carriers / callable genomes')+theme(axis.ticks.y=element_blank())
p3p<-ggplot(ps)+rowlines()+geom_segment(aes(x=minimum_gains_among_optima,xend=maximum_gains_among_optima,y=row,yend=row),linewidth=1.7,colour='#7570B3')+geom_point(aes(minimum_gains_among_optima,row),shape=21,fill='#C6BDE0',colour='#7570B3',size=2.4)+geom_point(aes(maximum_gains_among_optima,row),shape=21,fill='#C6BDE0',colour='#7570B3',size=2.4)+geom_text(aes(x=maximum_gains_among_optima+.6,y=row,label=range),hjust=0,size=2.8)+scale_x_continuous(limits=c(-.5,19),breaks=c(0,5,10,15),expand=c(0,0))+rowscale()+labs(title='B  Parsimony gain counts',x='Number of gains',y=NULL,caption='Gain-count bounds across optimal histories')+theme(axis.ticks.y=element_blank())
p3c<-ggplot(f3s)+rowlines()+geom_segment(aes(x=low,xend=high,y=y,yend=y,colour=event),linewidth=1.1,alpha=.55)+geom_point(aes(mean,y,fill=event,shape=event),colour='grey30',stroke=.25,size=2.3)+geom_text(aes(x=17.8,y=y,label=text,colour=event),hjust=0,size=2.7,show.legend=FALSE)+scale_fill_manual(values=c(Gain='#D95F02',Loss='#1B9E77'))+scale_colour_manual(values=c(Gain='#D95F02',Loss='#1B9E77'))+scale_shape_manual(values=c(Gain=21,Loss=24))+scale_x_continuous(limits=c(-.5,25.5),breaks=c(0,5,10,15),expand=c(0,0))+rowscale()+labs(title='C  Stochastic gain and loss reconstruction',x='Number of events',y=NULL,colour='Event',fill='Event',shape='Event',caption='ER SIMMAP: mean [95% simulation interval]')+theme(axis.ticks.y=element_blank())
write_tsv(f3prev,'Fig3_lineage_prevalence');write_tsv(ps,'Fig3_parsimony_gain_bounds');write_tsv(f3s,'Fig3_simmap_gain_loss_intervals')

savefig(wrap_plots(trees,ncol=3)/pcorr/((p3a+labs(title='C  Lineage-specific carriage'))+(p3p+labs(title='D  Parsimony gain counts'))+(p3c+labs(title='E  Stochastic gain and loss reconstruction'))+plot_layout(widths=c(1.25,.85,1.35)))+plot_layout(heights=c(3.145,.85,1.2)),'Alt_Fig2_integrated_recurrence',386,455,FALSE)
selection$gene<-factor(selection$gene,levels=GENES[1:4]);selection$FUBAR_Ppositive_above_09<-selection$FUBAR_Ppositive>.9
keys<-selection[paste(selection$gene,selection$site)%in%c('mrcA 506','mrcA 625','mrcA 708','pbp1 564','tp47 394'),];keys$txt<-as.character(keys$site)
pf<-ggplot(selection,aes(site,FUBAR_Ppositive))+geom_hline(yintercept=.9,linetype=2,colour='grey55',linewidth=.25)+geom_point(aes(fill=FUBAR_Ppositive),shape=21,size=.9,colour='grey35',stroke=.16,alpha=.7)+geom_point(data=keys,shape=21,size=2.6,fill='#2166AC',colour='grey20',stroke=.4)+geom_text(data=keys,aes(x=site+ifelse(gene=='mrcA'&site==625,-30,ifelse(gene=='mrcA'&site==708,30,0)),y=1.055,label=as.character(site)),size=2.5,fontface='bold')+facet_wrap(~gene,nrow=1,scales='free_x')+scale_fill_gradient(low='grey90',high='#2166AC',limits=c(0,1),breaks=c(0,.5,1),name='Posterior\nprobability')+scale_y_continuous(limits=c(0,1.10),breaks=c(0,.5,.9))+labs(title='A  Pervasive selection in PBP-associated genes',x='Amino-acid position',y='Pr(beta > alpha)')+common_theme()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position='right')+guides(fill=guide_colourbar(barheight=grid::unit(24,'mm')))
write_tsv(selection,'Fig3_all_codons');savefig(pf,'Fig4_selection',386,110)

hap <- read_tsv('nucleotide_haplotypes.tsv');hc <- read_tsv('haplotype_lineages.tsv');hs <- read_tsv('haplotype_samples.tsv')
bases <- c('T','C','A','G');codons <- as.vector(outer(as.vector(outer(bases,bases,paste0)),bases,paste0))
gencode <- setNames(strsplit('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG','')[[1]],
 unlist(lapply(bases,function(b1) unlist(lapply(bases,function(b2) paste0(b1,b2,bases))))))
translate <- function(s) {
 stopifnot(nchar(s)%%3==0,grepl('^[ACGT]+$',s))
 s <- paste(rev(strsplit(chartr('ACGT','TGCA',s),'')[[1]]),collapse='')
 starts <- seq(1,nchar(s),by=3);aa <- paste(gencode[substring(s,starts,starts+2)],collapse='')
 stopifnot(substr(aa,nchar(aa),nchar(aa))=='*',!grepl('*',substr(aa,1,nchar(aa)-1),fixed=TRUE))
 substr(aa,1,nchar(aa)-1)
}
hap$protein <- vapply(hap$sequence_forward_reference_orientation,translate,character(1));stopifnot(all(nchar(hap$protein)==884))
seqs <- sort(unique(hap$protein));hap$protein_type <- sprintf('P%02d',match(hap$protein,seqs))
hc <- merge(hc,hap[,c('gene_haplotype','protein_type')],by='gene_haplotype',all.x=TRUE);stopifnot(!anyNA(hc$protein_type))
pc <- aggregate(n_samples~protein_type+lineage,hc,sum);tot <- aggregate(n_samples~protein_type,pc,sum)
stopifnot(sum(tot$n_samples)==sum(truth(hs$gene_complete_acgt)))
pdbref<-readLines(file.path(PKG,'data/mrcA_model.pdb'));aa3<-trimws(substr(pdbref[substr(pdbref,1,4)=='ATOM'&trimws(substr(pdbref,13,16))=='CA'],18,20))
aa1<-setNames(strsplit('ARNDCQEGHILKMFPSTWYV','')[[1]],strsplit('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL',' ')[[1]])
refprotein<-paste(aa1[aa3],collapse='');wt_id<-sprintf('P%02d',which(seqs==refprotein));stopifnot(length(wt_id)==1,nchar(refprotein)==884)
ord<-c(wt_id,setdiff(tot$protein_type[order(-tot$n_samples,tot$protein_type)],wt_id))
idmap<-data.frame(protein_type=ord,display_id=seq_along(ord),status=c('WT','Mutation'),reference='AF-O83703-F1-model_v6 sequence; exact full-length match defines reference WT')
idmap$status<-ifelse(idmap$display_id==1,'WT','Mutation')
idmap$display_label<-ifelse(idmap$display_id==1,'WT',paste0('Mutation #',idmap$display_id-1))
pm<-do.call(rbind,strsplit(seqs,''));variable<-which(apply(pm,2,function(x)length(unique(x)))>1);ref<-strsplit(refprotein,'')[[1]]
heat<-expand.grid(site=variable,protein_type=sprintf('P%02d',seq_along(seqs)),stringsAsFactors=FALSE)
heat$aa<-mapply(function(p,s)pm[match(p,sprintf('P%02d',seq_along(seqs))),s],heat$protein_type,heat$site)
heat$display<-factor(ifelse(heat$aa==ref[heat$site],'WT','Mut'),levels=c('WT','Mut'))
heat$protein_type<-factor(heat$protein_type,levels=rev(ord));pc$protein_type<-factor(pc$protein_type,levels=rev(ord));type_labels<-setNames(idmap$display_label,idmap$protein_type)
seqpal<-c(WT='#EEEEEE',Mut='#8B73AE')
ph<-ggplot(heat,aes(factor(site),protein_type,fill=display))+geom_tile(colour='grey80',linewidth=.15)+geom_text(aes(label=aa,colour=ifelse(display=='Mut','white','black')),size=2.7)+scale_colour_identity()+scale_fill_manual(values=seqpal)+scale_y_discrete(labels=type_labels)+labs(title='A  mrcA protein sequence types',x='Variable amino-acid position',y=NULL,fill=NULL,caption='Full-length callable mrcA: 1,328/1,861 genomes')+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=45,hjust=1))
pc$percent<-100*pc$n_samples/as.numeric(tapply(pc$n_samples,pc$lineage,sum)[pc$lineage]);pc$lineage<-factor(pc$lineage,levels=c('Nichols','SS14'));pc$display_label<-type_labels[as.character(pc$protein_type)]
pbubble<-ggplot(pc,aes(lineage,protein_type))+geom_hline(yintercept=seq(.5,16.5,1),linetype='dashed',colour='grey85',linewidth=.2)+geom_point(aes(size=n_samples,fill=lineage),shape=21,colour='grey30',stroke=.25,alpha=.95)+geom_text(aes(label=n_samples),hjust=-1.3,size=2.7)+scale_fill_manual(values=LINEAGE_COL,guide='none')+scale_size(range=c(1.8,8),breaks=c(1,10,100,500))+scale_x_discrete(expand=expansion(add=.9))+scale_y_discrete(limits=rev(ord))+labs(title='B  Lineage distribution',x=NULL,y=NULL,size='Samples',caption='Nichols: n = 387; SS14: n = 941')+guides(size=guide_legend(override.aes=list(fill='grey65',colour='grey40')))+theme(axis.text.y=element_blank(),panel.grid.major.x=element_line(colour='grey83',linewidth=.2,linetype='dashed'),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank())
write_tsv(hap[,c('gene_haplotype','protein_type','protein')],'Fig5_DNA_to_protein_types');write_tsv(idmap,'Fig5_display_number_reference_crosswalk');write_tsv(heat,'Fig5_amino_acid_matrix');write_tsv(pc,'Fig5_type_lineages')
nt<-tc$nt;tx<-tc$x/max(tc$x);tree_ids<-tr$tip.label;ti<-match(tree_ids,mat$sample)
stopifnot(nt==1861,!anyNA(ti),setequal(tree_ids,meta$sample))
tsg<-data.frame(x=tx[tr$edge[,1]],xend=tx[tr$edge[,2]],y=tc$y[tr$edge[,2]],yend=tc$y[tr$edge[,2]])
tvr<-do.call(rbind,lapply(unique(tr$edge[,1]),function(pa)data.frame(x=tx[pa],xend=tx[pa],y=min(tc$y[tr$edge[tr$edge[,1]==pa,2]]),yend=max(tc$y[tr$edge[tr$edge[,1]==pa,2]]))))
tsg<-rbind(tsg,tvr)
tiporder<-data.frame(sample=tree_ids,y=tc$y[seq_len(nt)],lineage=meta$lineage[ti])
single_loci<-c(FOCAL[1:2],'mrcA:p.Ala506Met',FOCAL[3:4])
stopifnot(!anyDuplicated(hs$sample),setequal(hs$sample,tree_ids),length(ord)==16)
hs$protein_type<-hap$protein_type[match(hs$gene_haplotype,hap$gene_haplotype)]
stopifnot(all(is.na(hs$protein_type)==!truth(hs$gene_complete_acgt)))
type_per_tip<-hs$protein_type[match(tree_ids,hs$sample)]
stopifnot(sum(!is.na(type_per_tip))==1328,sum(is.na(type_per_tip))==533)
for(k in ord)stopifnot(sum(type_per_tip==k,na.rm=TRUE)==tot$n_samples[match(k,tot$protein_type)])
type_x<-1.2+seq(0,by=.068,length.out=16);unknown_x<-max(type_x)+.11
track<-data.frame(sample=tree_ids,y=tiporder$y,x=1.09,state=tiporder$lineage,column='Lineage',kind='Lineage',raw_state=tiporder$lineage)
for(i in seq_along(ord)){
 st<-ifelse(is.na(type_per_tip),'Unassigned',ifelse(type_per_tip==ord[i],ifelse(i==1,'Reference type','Variant type'),'Other assigned type'))
 track<-rbind(track,data.frame(sample=tree_ids,y=tiporder$y,x=type_x[i],state=st,column=unname(type_labels[ord[i]]),kind='Protein sequence type',raw_state=type_per_tip))
}
track<-rbind(track,data.frame(sample=tree_ids,y=tiporder$y,x=unknown_x,state=ifelse(is.na(type_per_tip),'Unassigned','Other assigned type'),column='Unassigned',kind='Assignment status',raw_state=ifelse(is.na(type_per_tip),'Unassigned','Assigned')))
stopifnot(nrow(track)==nt*18,!anyNA(track$state))
track$state<-factor(track$state,levels=c('Nichols','SS14','Reference type','Variant type','Other assigned type','Unassigned'))
trackcols<-c(Nichols=unname(LINEAGE_COL['Nichols']),SS14=unname(LINEAGE_COL['SS14']),'Reference type'='#607D8B','Variant type'='#8B73AE','Other assigned type'='#FFFFFF','Unassigned'='#E6E6E6')
headers<-data.frame(x=c(1.09,type_x,unknown_x),y=nt+90,label=c('Lineage',type_labels[ord],'Unassigned'))
typecounts<-data.frame(x=c(type_x,unknown_x),y=nt+40,n=c(tot$n_samples[match(ord,tot$protein_type)],sum(is.na(type_per_tip))),type=c(type_labels[ord],'Unassigned'))
type_marks<-track[track$state%in%c('Reference type','Variant type'),]
scalevalue<-.2*max(tc$x)
ptree5<-ggplot()+geom_segment(data=tsg,aes(x=x,xend=xend,y=y,yend=yend),linewidth=.08,colour='grey55')+geom_tile(data=track,aes(x,y,fill=state),width=.047,height=1.02,colour=NA)+geom_point(data=type_marks,aes(x,y,fill=state),shape=22,size=.75,stroke=.1,colour='#58508D')+geom_text(data=headers,aes(x,y,label=label),angle=90,hjust=0,size=2.5)+geom_text(data=typecounts,aes(x,y,label=n),size=2.5,colour='grey30')+annotate('text',x=mean(type_x),y=nt+380,label='Protein sequence types',fontface='bold',size=2.8)+annotate('segment',x=.03,xend=.23,y=-45,yend=-45,linewidth=.4)+annotate('text',x=.13,y=-75,label=paste0(format(scalevalue,scientific=TRUE,digits=2),' substitutions/site'),size=2.6)+scale_fill_manual(values=trackcols,drop=FALSE)+scale_x_continuous(limits=c(-.03,unknown_x+.08),breaks=NULL,expand=c(0,0))+scale_y_continuous(limits=c(-95,nt+430),breaks=NULL,expand=c(0,0))+labs(title='C  Phylogenetic distribution of mrcA protein types',x=NULL,y=NULL,fill=NULL,caption='1,328 assigned genomes; 533 unassigned genomes retained')+theme(panel.border=element_blank(),axis.line=element_blank(),axis.ticks=element_blank(),panel.grid=element_blank(),legend.position='bottom',legend.key=element_rect(colour='grey60',linewidth=.15))+guides(fill=guide_legend(nrow=1,byrow=TRUE,override.aes=list(colour='grey60',linewidth=.15)))
write_tsv(typecounts,'Fig4C_protein_type_counts')
write_tsv(data.frame(sample=tree_ids,lineage=tiporder$lineage,protein_type=type_per_tip,display_type=unname(type_labels[type_per_tip]),assignment=ifelse(is.na(type_per_tip),'Unassigned incomplete sequence','Assigned full-length sequence')),'Fig4C_protein_type_assignments')
write_tsv(tiporder,'Fig5_primary_tree_tip_order');write_tsv(tsg,'Fig5_primary_tree_segments');write_tsv(track,'Fig5_full_tree_state_matrix')
savefig((ph+pbubble+plot_layout(widths=c(2.3,1)))/ptree5+plot_layout(heights=c(1.15,2.3)),'Fig5_mrcA_backgrounds')
pdb<-readLines(file.path(PKG,'data/mrcA_model.pdb'));ca<-pdb[substr(pdb,1,4)=='ATOM'&trimws(substr(pdb,13,16))=='CA'];ca<-data.frame(site=as.integer(substr(ca,23,26)),aa=trimws(substr(ca,18,20)),x=num(substr(ca,31,38)),y=num(substr(ca,39,46)),z=num(substr(ca,47,54)),plddt=num(substr(ca,61,66)))
stopifnot(nrow(ca)==884,ca$aa[ca$site==506]=='ALA',ca$aa[ca$site==625]=='MET',ca$aa[ca$site==708]=='GLY')
sitecols<-c('506'='#B2182B','625'='#EF8A62','708'='#EF8A62');keyca<-ca[ca$site%in%c(506,625,708),]
pconf<-ggplot(ca,aes(site,plddt))+geom_line(colour='grey60',linewidth=.35)+geom_point(data=keyca,aes(colour=factor(site)),size=2.8)+geom_text(data=keyca,aes(y=108,label=paste0(c(ALA='A',MET='M',GLY='G')[aa],site)),size=2.8)+scale_colour_manual(values=sitecols,guide='none')+scale_y_continuous(limits=c(0,115),breaks=c(0,50,100))+labs(title='C  Reference-model confidence',x='Amino-acid position',y='pLDDT',caption='Source: per-residue confidence stored in the supplied AlphaFold PDB file')
read_ppm <- function(name,crop=c(0,1,0,1),stride=1L){
 con<-file(file.path(PKG,'data',name),'rb');on.exit(close(con));stopifnot(readLines(con,n=1)=='P6')
 wh<-as.integer(strsplit(readLines(con,n=1),' ')[[1]]);stopifnot(readLines(con,n=1)=='255')
 raw<-readBin(con,'raw',n=prod(wh)*3);stopifnot(length(raw)==prod(wh)*3)
 ar<-array(raw,dim=c(3,wh[1],wh[2]));xx<-seq(max(1,floor(crop[1]*wh[1])),min(wh[1],ceiling(crop[2]*wh[1])),by=stride);yy<-seq(max(1,floor(crop[3]*wh[2])),min(wh[2],ceiling(crop[4]*wh[2])),by=stride)
 ar<-ar[,xx,yy,drop=FALSE];cols<-rgb(as.integer(ar[1,,]),as.integer(ar[2,,]),as.integer(ar[3,,]),maxColorValue=255)
 as.raster(t(matrix(cols,nrow=length(xx),ncol=length(yy))))
}
structure_panel<-function(img,title,subtitle){ggplot()+annotation_custom(grid::rasterGrob(img,interpolate=TRUE),xmin=0,xmax=1,ymin=0,ymax=1)+coord_fixed(ratio=nrow(img)/ncol(img),xlim=c(0,1),ylim=c(0,1),expand=FALSE)+theme_void(base_size=9)+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())+labs(title=title,subtitle=subtitle)+theme(plot.title=element_text(size=11,face='bold'),plot.subtitle=element_text(size=8),plot.margin=margin(6,6,6,6))}
pwhole<-structure_panel(read_ppm('MrcA_reference_structure.ppm'),'A  mrcA reference structure','Reference-model rendering; red: A506; orange: M625 and G708')+labs(caption='AlphaFold reference model · UniProt O83703')
plocal<-structure_panel(read_ppm('MrcA_local_context.ppm'),'B  Local structural context','Reference-model local view')+labs(caption='A506: red; M625 and G708: orange')+theme(plot.caption=element_text(size=8,hjust=.5))
sm<-read_tsv('structural_metrics.tsv')
sm<-sm[sm$gene=='mrcA' & sm$position%in%c(506L,625L,708L),,drop=FALSE]
sm<-sm[!duplicated(sm$position),,drop=FALSE]
sm<-sm[match(c(506L,625L,708L),sm$position),,drop=FALSE]
sr<-read_tsv('sasa_resolution.tsv');sr<-sr[sr$n_points==2880,,drop=FALSE]
stopifnot(nrow(sm)==3L,nrow(sr)==3L,setequal(sr$position,c(506L,625L,708L)),all(sm$gene=='mrcA'))
sm$sasa<-sr$SASA_A2[match(sm$position,sr$position)]
sm$plddt<-ca$plddt[match(sm$position,ca$site)]
sm$sasa_n_points<-2880L;sm$probe_radius_A<-1.4
sm$residue_label<-paste0(sm$wt_aa,sm$position)
stopifnot(identical(as.character(sm$residue_label),c('A506','M625','G708')),all(is.finite(sm$sasa)))
sm<-sm[,c('gene','position','wt_aa','pdb_aa','wt_match','residue_label','sasa','plddt','sasa_n_points','probe_radius_A')]
sm$label<-factor(sm$residue_label,levels=c('G708','M625','A506'))
pmetric<-ggplot(sm,aes(sasa,label))+geom_point(size=3,colour='#2166AC',shape=16)+geom_text(aes(label=sprintf('%.2f',sasa)),hjust=-.35,size=2.8,colour='grey20')+scale_x_continuous(expand=expansion(mult=c(.07,.30)))+labs(title='D  Residue solvent accessibility',x=expression(plain(SASA)~(ring(plain(A))^2)),y=NULL)
write_tsv(ca,'Fig6_reference_CA_coordinates');write_tsv(sm,'Fig6_reported_structure_metrics')
savefig((pwhole+plocal+plot_layout(widths=c(1.6,1)))/(pconf+pmetric+plot_layout(widths=c(1.6,1)))+plot_layout(heights=c(1.6,1)),'Supp6_structural_context',310,270,FALSE)
dates<-read_tsv('dates.tsv');keep<-dates$matching_status=='reciprocal_unique_match'&is.finite(dates$strict_posterior_simmap_median)&is.finite(dates$relaxed_posterior_simmap_median)&dates$strict_posterior_event_matching_support>=.5&dates$relaxed_posterior_event_matching_support>=.5;keep[is.na(keep)]<-FALSE;dd<-dates[keep,];events<-dd$strict_event;stopifnot(length(events)==9)
chronoplots<-list();tipyears<-list()
for(m in c('strict','relaxed')){tt<-ape::read.nexus(file.path(PKG,'data',paste0(m,'.tree')));tt<-ape::ladderize(tt);cc<-coords(tt);tt<-cc$tree;cal<-read_tsv(paste0(m,'_calibration.tsv'));xx<-cal$mcc_root_date_median[1]+cc$x;nt<-cc$nt;stopifnot(nt==533,all(tt$tip.label%in%mat$sample));tipyears[[m]]<-setNames(xx[seq_len(nt)],tt$tip.label)
 seg<-data.frame(x=xx[tt$edge[,1]],xend=xx[tt$edge[,2]],y=cc$y[tt$edge[,2]],yend=cc$y[tt$edge[,2]])
 vert<-do.call(rbind,lapply(unique(tt$edge[,1]),function(p)data.frame(x=xx[p],xend=xx[p],y=min(cc$y[tt$edge[tt$edge[,1]==p,2]]),yend=max(cc$y[tt$edge[tt$edge[,1]==p,2]]))));seg<-rbind(seg,vert)
 ti<-match(tt$tip.label,mat$sample);track<-do.call(rbind,lapply(seq_along(FOCAL),function(j){z<-mat[ti,FOCAL[j]];data.frame(x=2029+(j-1)*2.8,y=cc$y[seq_len(nt)],colour=ifelse(is.na(z),'#BBBBBB',ifelse(z==1,FOCAL_COL[LAB[j]],'#F6F6F6')),sample=tt$tip.label,locus=FOCAL[j],state=z)}))
 pp<-ggplot()+geom_segment(data=seg,aes(x=x,xend=xend,y=y,yend=yend,linetype='Phylogenetic branch'),linewidth=.14,colour='#666666')+scale_linetype_manual(values=c('Phylogenetic branch'=1),name=NULL)+geom_tile(data=track,aes(x,y,fill=colour),width=2.3,height=1.03)+scale_fill_identity(name='Variant state',breaks=c('#F6F6F6','#BBBBBB'),labels=c('Noncarrier','Uncallable'),guide='legend')+annotate('text',x=2029+(0:5)*2.8,y=nt+12,label=LAB,angle=90,hjust=0,size=2.4)+geom_vline(xintercept=1943,linetype=2,linewidth=.25,colour='grey65')+scale_x_continuous(breaks=c(1925,1950,1975,2000,2025),limits=c(1918,2046),expand=c(0,0))+scale_y_continuous(limits=c(-12,nt+63),breaks=NULL,expand=c(0,0))+labs(title=paste0(if(m=='strict')'A  Strict-clock time-scaled phylogeny' else 'B  Relaxed-clock time-scaled phylogeny'),caption='SS14: n = 533; MCC topology; CA node heights',x='Calendar year',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),panel.grid=element_blank())
 chronoplots[[m]]<-pp;write_tsv(track,paste0('Fig7_',m,'_tip_tracks'));write_tsv(seg,paste0('Fig7_',m,'_tree_segments'))
}
long<-do.call(rbind,lapply(c('strict','relaxed'),function(m)data.frame(event=dd$strict_event,source_event=dd[[paste0(m,'_event')]],clock=ifelse(m=='strict','Strict','Relaxed'),median=dd[[paste0(m,'_posterior_simmap_median')]],low=dd[[paste0(m,'_posterior_simmap_hpd_low')]],high=dd[[paste0(m,'_posterior_simmap_hpd_high')]],p_after=dd[[paste0(m,'_posterior_simmap_p_after_1943')]])))
event_locus<-sub('_origin_.*','',events);event_map<-data.frame(source_event=events,display_event=paste0(event_locus,' #',ave(seq_along(events),event_locus,FUN=seq_along)));write_tsv(event_map,'Fig6_event_number_crosswalk');long$display_event<-event_map$display_event[match(long$event,event_map$source_event)]
long$y<-match(long$event,rev(events))+ifelse(long$clock=='Strict',.18,-.18);ylabels<-rev(event_map$display_event)
orig<-read_tsv('strict_origins.tsv');orig$event_id<-paste0(orig$label,'_origin_',sprintf('%02d',orig$independent_origin_index));orig<-orig[orig$event=='gain'&orig$event_id%in%events,];stopifnot(!anyDuplicated(orig$event_id),setequal(orig$event_id,events))
bub<-do.call(rbind,lapply(seq_len(nrow(orig)),function(i){ids<-strsplit(orig$carrier_tip_set[i],';',fixed=TRUE)[[1]];yrs<-tipyears$strict[ids];stopifnot(!anyNA(yrs));tb<-table(floor(yrs+1e-6));data.frame(event=orig$event_id[i],year=as.numeric(names(tb)),n=as.numeric(tb),label=orig$label[i])}));bub$y<-match(bub$event,rev(events))
pd<-ggplot()+geom_vline(xintercept=1943,linetype=2,colour='grey65',linewidth=.3)+geom_point(data=bub,aes(year,y,size=n),colour='grey72',alpha=.65)+geom_segment(data=long,aes(x=low,xend=high,y=y,yend=y,colour=clock),linewidth=2,alpha=.5)+geom_point(data=long,aes(median,y,colour=clock,shape=clock),size=2)+scale_colour_manual(values=CLOCK_COL)+scale_shape_manual(values=c(Strict=16,Relaxed=17))+scale_size_area(max_size=5,breaks=c(1,10,30))+scale_y_continuous(breaks=seq_along(events),labels=ylabels,limits=c(.35,length(events)+.55))+labs(title='C  Reconstructed acquisition dates',subtitle='Median and 95% HPD; grey bubbles: strict-event carrier sets',x='Calendar year',y='Reconstructed gain (event ID)',caption='Event dates: median and 95% HPD',colour='Clock',shape='Clock',size='Sampled carriers (n)')
pp<-ggplot(long,aes(p_after,y,fill=clock))+geom_col(orientation='y',width=.27,alpha=.7)+geom_text(aes(x=p_after+.025,label=sprintf('%.3f',p_after),colour=clock),hjust=0,size=2.6,show.legend=FALSE)+scale_fill_manual(values=CLOCK_COL)+scale_colour_manual(values=CLOCK_COL,guide='none')+scale_y_continuous(breaks=seq_along(events),labels=ylabels,limits=c(.35,length(events)+.55))+scale_x_continuous(limits=c(0,1.20),breaks=c(0,.5,1),expand=c(0,0))+labs(title='D  Acquisition after 1943',x='Pr(gain date > 1943)',y=NULL,fill='Clock',caption=NULL)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
write_tsv(long,'Fig7_matched_dates');write_tsv(bub,'Fig7_carrier_sampling_years');dates$included<-keep;write_tsv(dates,'Fig7_matching_and_exclusions')
savefig((chronoplots$strict+chronoplots$relaxed)/(pd+pp+plot_layout(widths=c(1.8,1)))+plot_layout(heights=c(1.65,1)),'Fig6_dual_clock_history',340,320)
er$label<-paste(er$gene,label_locus(er$locus));er$label<-factor(er$label,levels=rev(unique(er$label)));erlong<-do.call(rbind,lapply(c('ER','ARD'),function(m)data.frame(label=er$label,model=m,mean=er[[paste0(m,'_simmap_gain_mean')]],lo=er[[paste0(m,'_simmap_gain_q025')]],hi=er[[paste0(m,'_simmap_gain_q975')]])));erlong$model<-factor(erlong$model,levels=c('ER','ARD'));erlong$estimate<-ifelse(is.finite(erlong$mean),sprintf('%.2f [%.2f, %.2f]',erlong$mean,erlong$lo,erlong$hi),'Fit unresolved')
limitER<-max(erlong$hi,na.rm=TRUE);erlong$textx<-limitER*1.06
pER<-ggplot(erlong,aes(y=label,colour=model))+geom_segment(data=erlong[is.finite(erlong$mean),],aes(x=lo,xend=hi,yend=label),linewidth=.6)+geom_point(aes(mean),size=1.9,na.rm=TRUE)+geom_text(aes(x=textx,label=estimate),hjust=0,size=2.35,show.legend=FALSE)+facet_wrap(~model,nrow=1)+scale_colour_manual(values=c(ER='#2166AC',ARD='#B2182B'))+scale_x_continuous(limits=c(0,limitER*1.83),breaks=seq(0,floor(limitER/5)*5,5))+labs(title='A  Gain-count model sensitivity',x='Gain events',y=NULL,colour='Model',caption='Blue: ER (equal gain/loss rates); red: ARD (different rates). Points: means; intervals: 2.5th–97.5th percentiles of stochastic mappings.\nLabels report mean [lower, upper].')+grid_dashed
savefig(pER+labs(title='ER and ARD reconstruction sensitivity'),'Supp2_reconstruction_sensitivity',386,245,FALSE);write_tsv(erlong,'Supp2_ER_ARD_numeric_labels');write_tsv(focal,'Supp2_parsimony_bounds')
es<-read_tsv('esm_site_means.tsv');ex<-read_tsv('esm_specific.tsv');names(es)[names(es)%in%c('site_mean_llr','score')]<-'mean_llr';es$site<-num(es$site);es$mean_llr<-num(es$mean_llr);es<-es[es$gene%in%GENES[1:4]&is.finite(es$site)&is.finite(es$mean_llr),];es$gene<-factor(es$gene,levels=GENES[1:4]);stopifnot(nrow(es)==2594)
obs<-vars[vars$gene%in%GENES[1:4]&vars$n>0,];marks<-aggregate(label~gene+site,obs,function(x)paste(x,collapse='/'));marks<-merge(marks,es[,c('gene','site','mean_llr')],by=c('gene','site'));marks$gene<-factor(marks$gene,levels=GENES[1:4])
marks<-marks[order(marks$gene,marks$site),];marks$labelx<-marks$site
for(g in unique(marks$gene)){ii<-which(marks$gene==g);minsep<-GLEN[as.character(g)]/60;for(j in seq_along(ii))if(j>1)marks$labelx[ii[j]]<-max(marks$site[ii[j]],marks$labelx[ii[j-1]]+minsep)}
pes<-ggplot(es,aes(site,mean_llr))+geom_hline(yintercept=0,colour='grey70',linewidth=.25)+geom_segment(aes(xend=site,y=0,yend=mean_llr),colour='grey78',linewidth=.18,alpha=.5)+geom_point(colour='grey65',size=.25,alpha=.45)+geom_segment(data=marks,aes(xend=labelx,yend=4.2),colour='grey55',linewidth=.15,linetype='dotted')+geom_point(data=marks,shape=21,fill='white',colour='grey25',size=1.8,stroke=.3)+geom_text(data=marks,aes(x=labelx,y=5,label=label),angle=90,hjust=0,size=2.2)+facet_wrap(~gene,scales='free_x',ncol=2)+scale_y_continuous(limits=c(-22,18),breaks=c(-20,-10,0),expand=expansion(mult=c(.02,.02)))+labs(title='A  ESM-2 sequence landscape',x='Amino-acid position',y='Mean LLR across scored substitutions',caption='Open circles locate observed variants on the site-mean landscape; labels at a shared codon are combined. They are not exact substitution scores.')

ex$key<-factor(paste(ex$gene,ex$label),levels=rev(paste(ex$gene,ex$label)));ex$sign<-ifelse(ex$llr<0,'Lower than reference','Higher than reference')
pex<-ggplot(ex,aes(llr,key))+geom_vline(xintercept=0,colour='grey60',linetype='dashed',linewidth=.3)+geom_segment(aes(x=0,xend=llr,yend=key),colour='grey80',linewidth=.4)+geom_point(aes(fill=sign),shape=21,colour='grey30',size=2.8,stroke=.3)+geom_text(aes(x=llr+ifelse(llr<0,-.18,.18),label=sprintf('%.3f',llr),hjust=ifelse(llr<0,1,0)),size=2.8)+scale_fill_manual(values=c('Lower than reference'='#2166AC','Higher than reference'='#B2182B'))+scale_x_continuous(limits=c(min(ex$llr)-1.4,max(ex$llr)+1.4))+labs(title='B  Substitution-specific ESM-2 scores',x='ESM-2 log-likelihood ratio',y=NULL,fill='Model compatibility',caption='Negative: lower model compatibility than reference. Positive: higher model compatibility. LLR is not a mutation rate or a measured fitness effect.')+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
savefig(pes/pex+plot_layout(heights=c(2,1)),'Supp3_ESM2',386,330,FALSE);write_tsv(es,'Supp3_ESM2_clean_landscape');write_tsv(ex,'Supp3_ESM2_exact_scores');write_tsv(marks,'Supp3_all_observed_mutation_annotations')
rtt<-read_tsv('rtt_points.tsv');summ<-read_tsv('rtt_summary.tsv');ts<-read_tsv('temporal_signal.tsv');stopifnot(!anyDuplicated(ts$n),setequal(ts$n,summ$n))
si<-match(summ$n,ts$n);summ$signal<-ifelse(ts$slope[si]>0&ts$R2[si]>ts$perm95[si],'Temporal signal supported','No supported temporal signal');summ$perm95<-ts$perm95[si]
summ$name<-ifelse(summ$n==1819,'Full cohort',ifelse(summ$n==1311,'SS14 lineage',ifelse(summ$n==508,'Nichols lineage','Selected SS14 subset')))
summ$label<-paste0(summ$name,'; n=',summ$n,'\nR-squared=',sprintf('%.4f',summ$r2),'; slope=',sprintf('%.2e',summ$slope),'\nRandomized R-squared: 95th percentile=',sprintf('%.4f',summ$perm95))
rtt$panel<-factor(summ$label[match(rtt$source_file,summ$source_file)],levels=summ$label);rtt$signal<-summ$signal[match(rtt$source_file,summ$source_file)];stopifnot(!anyNA(rtt$panel))
prtt<-ggplot(rtt,aes(year,root_to_tip))+geom_point(size=.5,alpha=.35,colour='grey35')+geom_smooth(aes(colour=signal,fill=signal),method='lm',se=TRUE,linewidth=.55,alpha=.12)+facet_wrap(~panel,ncol=2,scales='free')+scale_colour_manual(values=c('Temporal signal supported'='#2166AC','No supported temporal signal'='#B2182B'))+scale_fill_manual(values=c('Temporal signal supported'='#2166AC','No supported temporal signal'='#B2182B'))+labs(title='Temporal signal across sampling schemes',x='Collection year',y='Root-to-tip distance',colour=NULL,fill=NULL,caption='Blue requires a positive slope and observed R² above the randomized 95th percentile. Red fails at least one criterion.\nRandomization results are shown with the regressions; n denotes the actual number of genomes in each diagnostic data set, not the final dated subset.')+theme(strip.text=element_text(size=8))
savefig(prtt,'Supp5_temporal_regressions',320,330,FALSE);write_tsv(summ,'Supp5_temporal_signal_labels');write_tsv(ts,'Supplementary_table_temporal_randomization')
if(file.exists(file.path(PKG,'data/temporal_extended_historical.tsv')))write_tsv(read_tsv('temporal_extended_historical.tsv'),'Supplementary_table_historical_temporal_checks');write_tsv(read_tsv('clock_rates.tsv'),'Supplementary_table_clock_rates')
publication_order<-c('Fig1_dataset_and_variants','Fig2_lineages_and_recurrence','Fig3_selection','Fig4_mrcA_backgrounds','Fig5_time_calibrated_histories','Supp1_reconstruction_sensitivity','Supp2_structural_context','Supp3_ESM2','Supp4_temporal_regressions')
stopifnot(setequal(names(registry),publication_order))
registry<-registry[publication_order]
jsonlite::write_json(registry,file.path(OUT,'FIGURE_INDEX.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(OUT,'sessionInfo.txt'))
export_panel_sources(PKG,OUT)
writeLines('Five main and four supplementary figures rendered from frozen inputs.',file.path(OUT,'RUN_COMPLETE.txt'))
message('DONE ',OUT)
