#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop('Usage: Rscript run_figures.R PREPARED_WORK_DIR OUTPUT_DIR')
root <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
required <- c('ggplot2','patchwork','dplyr','tidyr','ggrepel','ape','jsonlite','maps','svglite','ragg','systemfonts')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop('Missing R packages: ', paste(missing, collapse = ', '))
if (as.character(packageVersion('ggplot2')) != '3.4.4' ||
    as.character(packageVersion('patchwork')) != '1.2.0') {
    stop('This plotting code retains the validated ggplot2 3.4.4 / patchwork 1.2.0 interface. Use the dedicated figure environment; do not change the analytical tp_wgs environment.')
}
font <- Sys.getenv('PBP_FIGURE_FONT', 'Nimbus Sans')
available <- unique(systemfonts::system_fonts()$family)
if (!tolower(font) %in% tolower(available)) {
    stop('Font not installed: ', font, '. Install it, or explicitly set PBP_FIGURE_FONT to the agreed alternative.')
}
Sys.setenv(PBP_FIGURE_FONT = font)
Sys.unsetenv('DISPLAY')
script <- file.path(root, 'figures', 'R', 'render_figures.R')
if (!file.exists(file.path(root, 'INPUT_MANIFEST.json'))) stop('Prepared cluster input manifest is missing')
code <- system2(file.path(R.home('bin'), 'Rscript'), c(shQuote(script), shQuote(out)))
if (code != 0L) stop('Main and supplementary figure rendering failed, exit code ', code)
if (!file.exists(file.path(out,'panel_data','PANEL_DATA_INDEX.tsv'))) stop('Panel source-data index was not exported')
index <- jsonlite::read_json(file.path(out,'FIGURE_INDEX.json'))
expected <- c('Fig1_dataset_and_variants','Fig2_lineages_and_recurrence','Fig3_selection',
              'Fig4_mrcA_backgrounds','Fig5_time_calibrated_histories',
              'Supp1_reconstruction_sensitivity','Supp2_structural_context',
              'Supp3_ESM2','Supp4_temporal_regressions')
if (!setequal(names(index), expected)) stop('Unexpected figure allocation')
for (name in names(index)) {
    for (ext in c('pdf','svg','png')) {
        path <- file.path(out,'figures',paste0(name,'.',ext))
        if (!file.exists(path) || file.info(path)$size == 0) stop('Missing figure: ',path)
    }
}
cat('Exported 5 main figures, 4 supplementary figures and panel source data.\n')
