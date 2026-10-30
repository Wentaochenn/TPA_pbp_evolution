export_plot_layers <- function(p, id, out) {
  dest <- file.path(out, 'rendered_layers', id)
  dir.create(dest, recursive=TRUE, showWarnings=FALSE)
  index <- list(); k <- 0L
  put <- function(d, name) {
    if (!is.data.frame(d) || !nrow(d)) return(invisible(NULL))
    for (n in names(d)) {
      if (is.list(d[[n]])) d[[n]] <- vapply(d[[n]], function(z) paste(z,collapse=';'),character(1))
    }
    write.table(d,file.path(dest,paste0(name,'.tsv')),sep='\t',quote=TRUE,row.names=FALSE,na='NA')
  }
  walk <- function(q) {
    if (inherits(q,'patchwork')) {
      lapply(q$patches$plots,walk)
      q$patches <- NULL; class(q) <- setdiff(class(q),'patchwork')
    }
    if (!inherits(q,'ggplot')) return(invisible(NULL))
    k <<- k+1L; prefix <- sprintf('plot_%02d',k)
    title <- paste(deparse(q$labels$title),collapse=' ')
    index[[length(index)+1L]] <<- data.frame(plot=prefix,title=title,tag=paste(q$labels$tag,collapse=''))
    put(q$data,paste0(prefix,'_input'))
    for(i in seq_along(q$layers)) put(q$layers[[i]]$data,paste0(prefix,'_layer_',i,'_input'))
    built <- ggplot_build(q)
    for(i in seq_along(built$data)) put(built$data[[i]],paste0(prefix,'_layer_',i,'_rendered'))
  }
  walk(p)
  write.table(do.call(rbind,index),file.path(dest,'LAYER_INDEX.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
}

export_panel_sources <- function(pkg,out) {
  root <- file.path(out,'panel_data'); dir.create(root,recursive=TRUE)
  entries <- list()
  emit <- function(fig,panel,files,meaning) {
    dest <- file.path(root,fig,panel);dir.create(dest,recursive=TRUE)
    for (f in files) {
      input <- if(startsWith(f,'data/'))file.path(pkg,f)else file.path(out,'source_data',paste0(f,'.tsv'))
      name <- sub('^(Fig[0-9]+[A-Z]?|Supp[0-9]+|Supplementary_table)_','',basename(input))
      name <- sub('^TableS[0-9]+[A-Z]?_','',name)
      replacements <- c(original_style_all_target_indicators.tsv='variant_carriage_counts.tsv',
                        original_style_timeline.tsv='carrier_counts_by_year.tsv')
      if(name %in% names(replacements)) name <- unname(replacements[name])
      stopifnot(file.exists(input),!file.exists(file.path(dest,name)))
      stopifnot(file.copy(input,file.path(dest,name),overwrite=FALSE))
      entries[[length(entries)+1L]] <<- data.frame(figure=fig,panel=panel,file=file.path(fig,panel,name),description=meaning)
    }
  }
  emit('Figure1','A',c('data/map_coordinates.tsv','data/metadata.tsv'),'Country counts and map locations; map geometry is supplied with the R maps package and exported in rendered_layers.')
  emit('Figure1','B',c('data/metadata.tsv'),'Country genome counts; aggregation code in render_figures.R.')
  emit('Figure1','C',c('Fig1_country_year'),'Country-by-year genome counts.')
  emit('Figure1','D',c('data/metadata.tsv'),'Country-by-lineage counts and proportions.')
  emit('Figure1','E',c('Fig1_original_style_all_target_indicators'),'Carrier, callable and missing counts per displayed variant.')
  emit('Figure1','F',c('Fig1_original_style_timeline','Fig1_variant_year_denominators'),'Year-specific carrier counts and callable denominators. Categories without dated carriers have missing years and no plotted point.')
  emit('Figure2','A',c('data/primary.treefile','data/matrix.tsv','data/metadata.tsv','Fig2_lineage_counts'),'Full phylogeny, observed tip genotypes and lineage-specific denominators.')
  emit('Figure2','B',c('Fig2_country_matched_denominators','Fig2_Pearson_r_unweighted','Fig2_country_label_positions','Fig2_country_leader_segments'),'Callable-matched country frequencies, Pearson correlations and final label coordinates.')
  emit('Figure2','C',c('Fig3_lineage_prevalence'),'Lineage-specific carrier and callable counts.')
  emit('Figure2','D',c('Fig3_parsimony_gain_bounds','data/lineage_summary.tsv'),'Full-tree parsimony gain bounds and lineage-specific A506T histories.')
  emit('Figure2','E',c('Fig3_simmap_gain_loss_intervals','data/simmap_draws.tsv'),'ER stochastic-map gains and losses and empirical quantiles.')
  emit('Figure3','A',c('Fig3_all_codons'),'FUBAR posterior probabilities for all 2,594 codons across four genes; focal codons are labelled.')
  emit('Figure4','A',c('Fig5_DNA_to_protein_types','Fig5_display_number_reference_crosswalk','Fig5_amino_acid_matrix'),'Full-length protein sequence definitions, reference identity and displayed residues.')
  emit('Figure4','B',c('Fig5_type_lineages'),'Protein type counts and percentages within the complete-sequence lineage subsets.')
  emit('Figure4','C',c('data/primary.treefile','Fig5_primary_tree_tip_order','Fig5_primary_tree_segments','Fig5_full_tree_state_matrix','Fig4C_protein_type_assignments','Fig4C_protein_type_counts'),'All 1861 tree tips; assigned and unassigned protein types remain distinct.')
  for(m in c('strict','relaxed')) emit('Figure5',if(m=='strict')'A'else'B',c(paste0('data/',m,'.tree'),paste0('data/',m,'_calibration.tsv'),paste0('Fig7_',m,'_tree_segments'),paste0('Fig7_',m,'_tip_tracks'),'data/clock_rates.tsv',if(m=='strict')'data/TableS13B_strict_ESS.tsv'else'data/TableS13C_relaxed_ESS.tsv'),'Time-scaled CA-MCC phylogeny, observed terminal genotypes, posterior rate summaries and scalar ESS.')
  emit('Figure5','C',c('Fig7_matched_dates','Fig7_carrier_sampling_years','Fig7_matching_and_exclusions','Fig6_event_number_crosswalk'),'Matched event-date medians and HPDs; sampled carriers are observations, not inferred origin dates.')
  emit('Figure5','D',c('Fig7_matched_dates','Fig7_matching_and_exclusions'),'Posterior probabilities of matched event dates after 1943, by clock model.')
  emit('FigureS1','A',c('Supp2_ER_ARD_numeric_labels','data/er_ard.tsv'),'Completed ER/ARD comparisons for 23 eligible variants, including tp47 R107G.')
  emit('FigureS2','A',c('data/MrcA_reference_structure.png','data/mrcA_model.pdb','data/structure_camera_views.json','data/structure_render_provenance.json','data/structure_marked_residues.tsv'),'Reference PDB, CPU-rendered overview, fixed camera parameters and marked-residue coordinates.')
  emit('FigureS2','B',c('data/MrcA_local_context.png','data/mrcA_model.pdb','data/structure_camera_views.json','data/structure_render_provenance.json','data/structure_marked_residues.tsv'),'Reference PDB, independently rendered local view, fixed camera parameters and marked-residue coordinates.')
  emit('FigureS2','C',c('Fig6_reference_CA_coordinates'),'Reference residue pLDDT from the archived coordinate file.')
  emit('FigureS2','D',c('Fig6_reported_structure_metrics','data/mrcA_model.pdb'),'MrcA reference residues A506, M625 and G708; uniform 2880-point SASA, 1.4-A probe.')
  emit('FigureS3','A',c('Supp3_ESM2_clean_landscape','Supp3_all_observed_mutation_annotations'),'Site means across 19 substitutions; annotations identify observed variants.')
  emit('FigureS3','B',c('Supp3_ESM2_exact_scores'),'Archived substitution-specific LLRs; these are not site means.')
  emit('FigureS4','A',c('data/rtt_points.tsv','Supp5_temporal_signal_labels','Supplementary_table_temporal_randomization','data/TableS12A_clustered_date_tests.tsv','data/TableS12B_subset_stability.tsv'),'Root-to-tip data and archived permutation diagnostics, faceted by sampling scheme.')
  write.table(do.call(rbind,entries),file.path(root,'PANEL_DATA_INDEX.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
}
