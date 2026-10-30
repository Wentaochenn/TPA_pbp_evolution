# Render the publication metadata table from its canonical TSV and header specification.
write_metadata_workbook <- function(tsv, spec, destination) {
  if (!requireNamespace('openxlsx', quietly = TRUE)) {
    stop('Install openxlsx in the publication R environment before exporting Table S1.')
  }
  specification <- jsonlite::fromJSON(spec)
  x <- read.delim(tsv, check.names = FALSE, stringsAsFactors = FALSE,
                  colClasses = 'character', na.strings = 'NA')
  fields <- specification$columns$field
  stopifnot(identical(names(x), fields), nrow(x) == specification$row_count,
            ncol(x) == specification$column_count)
  years <- suppressWarnings(as.numeric(x$year))
  stopifnot(all(is.na(x$year) | !is.na(years)))
  x$year <- years
  sheet <- 'Metadata'
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheet, gridLines = FALSE)
  openxlsx::modifyBaseFont(wb, fontName = 'Arial', fontSize = 11,
                           fontColour = '#000000')
  openxlsx::writeData(wb, sheet, x, startRow = 2, colNames = TRUE,
                      rowNames = FALSE, withFilter = TRUE, keepNA = FALSE)
  columns <- seq_len(ncol(x))
  data_rows <- seq.int(3L, nrow(x) + 2L)
  base <- openxlsx::createStyle(fontName = 'Arial', fontSize = 11,
                               fontColour = '#000000', valign = 'center')
  header <- openxlsx::createStyle(fontName = 'Arial', fontSize = 11,
                                 fontColour = '#000000', textDecoration = 'bold',
                                 fgFill = '#E7E6E6', halign = 'center',
                                 valign = 'center', wrapText = TRUE,
                                 border = 'bottom', borderColour = '#B7B7B7')
  gene_header <- openxlsx::createStyle(fontName = 'Arial', fontSize = 11,
                                      fontColour = '#000000',
                                      textDecoration = c('bold', 'italic'),
                                      fgFill = '#E7E6E6', halign = 'center',
                                      valign = 'center', border = 'bottom',
                                      borderColour = '#B7B7B7')
  openxlsx::addStyle(wb, sheet, base, rows = seq_len(nrow(x) + 2L),
                     cols = columns, gridExpand = TRUE)
  openxlsx::addStyle(wb, sheet, header, rows = 1:2, cols = columns,
                     gridExpand = TRUE, stack = TRUE)
  metadata_columns <- which(specification$columns$kind != 'amino_acid')
  for (j in metadata_columns) {
    openxlsx::writeData(wb, sheet, '', startCol = j, startRow = 2,
                        colNames = FALSE)
    openxlsx::writeData(wb, sheet, specification$columns$label[j],
                        startCol = j, startRow = 1, colNames = FALSE)
    openxlsx::mergeCells(wb, sheet, cols = j, rows = 1:2)
  }
  site_columns <- which(specification$columns$kind == 'amino_acid')
  for (gene in unique(specification$columns$group[site_columns])) {
    gene_columns <- which(specification$columns$group == gene)
    stopifnot(identical(gene_columns, seq.int(min(gene_columns), max(gene_columns))))
    openxlsx::writeData(wb, sheet, gene, startCol = min(gene_columns),
                        startRow = 1, colNames = FALSE)
    openxlsx::mergeCells(wb, sheet, cols = gene_columns, rows = 1)
    openxlsx::addStyle(wb, sheet, gene_header, rows = 1,
                       cols = gene_columns, gridExpand = TRUE, stack = TRUE)
    openxlsx::writeData(wb, sheet,
                        matrix(specification$columns$label[gene_columns], nrow = 1),
                        startCol = min(gene_columns), startRow = 2, colNames = FALSE)
  }
  reference_style <- openxlsx::createStyle(fgFill = '#F2F2F2', halign = 'center')
  alternative_style <- openxlsx::createStyle(fgFill = '#FFF2CC', halign = 'center')
  uncallable_style <- openxlsx::createStyle(fgFill = '#FFFFFF', halign = 'center')
  summary_rows <- nrow(x) + c(4L, 5L)
  summary_values <- matrix(NA_real_, nrow = 2, ncol = length(site_columns))
  for (k in seq_along(site_columns)) {
    j <- site_columns[k]
    field <- fields[j]
    values <- x[[field]]
    reference <- substr(specification$columns$label[j], 1, 1)
    stopifnot(!anyNA(values))
    reference_rows <- data_rows[values == reference]
    alternate_rows <- data_rows[values != reference & values != 'uncallable']
    missing_rows <- data_rows[values == 'uncallable']
    if (length(reference_rows)) openxlsx::addStyle(wb, sheet, reference_style,
        rows = reference_rows, cols = j, gridExpand = TRUE, stack = TRUE)
    if (length(alternate_rows)) openxlsx::addStyle(wb, sheet, alternative_style,
        rows = alternate_rows, cols = j, gridExpand = TRUE, stack = TRUE)
    if (length(missing_rows)) openxlsx::addStyle(wb, sheet, uncallable_style,
        rows = missing_rows, cols = j, gridExpand = TRUE, stack = TRUE)
    summary_values[, k] <- c(length(alternate_rows), nrow(x) - length(missing_rows))
    expected <- specification$summaries[specification$summaries$field == field, ]
    stopifnot(nrow(expected) == 1L,
              summary_values[1, k] == expected$carriers,
              summary_values[2, k] == expected$callable)
  }
  openxlsx::writeData(wb, sheet, summary_values, startCol = min(site_columns),
                      startRow = summary_rows[1], colNames = FALSE)
  labels <- c('Non-reference carriers (n)', 'Callable genomes (n)')
  for (k in seq_along(summary_rows)) {
    openxlsx::writeData(wb, sheet, labels[k], startCol = 1,
                        startRow = summary_rows[k], colNames = FALSE)
    openxlsx::mergeCells(wb, sheet, cols = seq_len(min(site_columns) - 1L),
                         rows = summary_rows[k])
  }
  openxlsx::addStyle(wb, sheet, header, rows = summary_rows, cols = columns,
                     gridExpand = TRUE, stack = TRUE)
  note_style <- openxlsx::createStyle(fontName = 'Arial', fontSize = 11,
                                     fontColour = '#000000', halign = 'left',
                                     valign = 'top', wrapText = TRUE)
  note_rows <- max(summary_rows) + 2L + seq_along(specification$footnotes) - 1L
  for (k in seq_along(specification$footnotes)) {
    note <- paste0(k, '. ', specification$footnotes[k])
    openxlsx::writeData(wb, sheet, note, startCol = 1,
                        startRow = note_rows[k], colNames = FALSE)
    openxlsx::mergeCells(wb, sheet, cols = 1:8, rows = note_rows[k])
    openxlsx::addStyle(wb, sheet, note_style, rows = note_rows[k],
                       cols = 1:8, gridExpand = TRUE, stack = TRUE)
    openxlsx::setRowHeights(wb, sheet, rows = note_rows[k],
                            heights = max(30, 15 * ceiling(nchar(note) / 140)))
  }
  openxlsx::setColWidths(wb, sheet, cols = metadata_columns, widths = 24)
  openxlsx::setColWidths(wb, sheet, cols = site_columns, widths = 13)
  openxlsx::setRowHeights(wb, sheet, rows = 1:2, heights = 24)
  openxlsx::setRowHeights(wb, sheet, rows = data_rows, heights = 18)
  openxlsx::setRowHeights(wb, sheet, rows = summary_rows, heights = 24)
  openxlsx::freezePane(wb, sheet, firstActiveRow = 3, firstActiveCol = 2)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, destination, overwrite = TRUE)
  invisible(destination)
}
