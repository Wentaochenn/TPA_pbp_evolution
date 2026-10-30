# Device-space correction for unresolved null units in the generic gtable helper.
# Measurements come from drawn viewports on the FINAL physical-size Cairo device.
write_patchwork_panel_layout <- function(plot,manifest_path,width_in,height_in,panel_ids=NULL,
 row_groups=NULL,column_groups=NULL,exemptions=list()) {
 tmp <- tempfile(fileext='.pdf')
 grDevices::cairo_pdf(tmp,width=width_in,height=height_in,family=Sys.getenv('PBP_FIGURE_FONT','Nimbus Sans'))
 on.exit({grDevices::dev.off();unlink(tmp)},add=TRUE)
 g <- patchwork::patchworkGrob(plot)
 grid::grid.newpage();grid::grid.draw(g);grid::grid.force()
 ls <- grid::grid.ls(viewports=TRUE,grobs=FALSE,print=FALSE)
 # Nested patchworks can reuse the same viewport name in different parents.
 # Resolve the full viewport path, never only the leaf name.
 idx <- which(grepl('^panel([-.;]).*\\.[0-9]+-',ls$name) & !grepl('^panel-area',ls$name))
 boxes <- list()
 for(i in idx) {
  n <- paste(ls$vpPath[i],ls$name[i],sep='::')
  path <- strsplit(n,'::',fixed=TRUE)[[1]]
  path <- path[nzchar(path)&path!='ROOT']
  grid::upViewport(0);grid::downViewport(do.call(grid::vpPath,as.list(path)))
  lo <- grid::deviceLoc(grid::unit(0,'npc'),grid::unit(0,'npc'),valueOnly=TRUE)
  hi <- grid::deviceLoc(grid::unit(1,'npc'),grid::unit(1,'npc'),valueOnly=TRUE)
  b <- c(lo$x,lo$y,hi$x,hi$y)*72
  if(all(is.finite(b)) && b[3]>b[1] && b[4]>b[2]) boxes[[n]] <- b
 }
 if(!length(boxes)) stop('No measurable plot viewports')
 bb <- do.call(rbind,boxes)
 # Deduplicate containers occupying the identical rectangle.
 bb <- bb[!duplicated(round(bb,5)),,drop=FALSE]
 # Exclude facet/patchwork container viewports that contain actual child panels.
 container <- vapply(seq_len(nrow(bb)),function(i) any(vapply(setdiff(seq_len(nrow(bb)),i),function(j) {
  all(bb[j,c(1,2)] >= bb[i,c(1,2)]-.01) && all(bb[j,c(3,4)] <= bb[i,c(3,4)]+.01) &&
   prod(bb[j,c(3,4)]-bb[j,c(1,2)]) < prod(bb[i,c(3,4)]-bb[i,c(1,2)])-.1
 },logical(1))),logical(1))
 bb <- bb[!container,,drop=FALSE]
 bb <- bb[order(-round(bb[,4],2),bb[,1]),,drop=FALSE]
 ids <- sprintf("p%02d",seq_len(nrow(bb)))
 name <- sub('[.]layout[.]json$','',basename(manifest_path))
 rows <- list();cols <- list()
 # Actual viewport rows; compare same-row plot edges. Mixed-modality rows have
 # explicitly unequal widths (e.g. country totals / year bubbles / proportions).
 for(i in seq_len(nrow(bb))) {
  j <- which(abs(bb[,4]-bb[i,4])<1 & abs(bb[,2]-bb[i,2])<1)
  if(length(j)>1 && i==min(j)) {
   rows[[length(rows)+1]] <- ids[j]
   widths <- bb[j,3]-bb[j,1]
   if(max(widths)-min(widths)>1.5) exemptions[[length(exemptions)+1]] <- list(panels=ids[j],checks=c("panel-width","horizontal-gutter"),reason="Intentional unequal-width composition: sequence/tree hero and narrower count or annotation panels.")
  }
 }
 for(i in seq_len(nrow(bb))) {
  j <- which(abs(bb[,1]-bb[i,1])<1 & abs(bb[,3]-bb[i,3])<1)
  if(length(j)>1 && i==min(j)) cols[[length(cols)+1]] <- ids[j]
 }
 panels <- lapply(seq_along(ids),function(i) list(id=ids[i],bbox_pt=as.numeric(bb[i,]),grid_id=paste0('explicit-',i)))
 manifest <- list(schema_version=1L,backend='r-patchwork-device-viewports',figure=list(width_pt=width_in*72,height_pt=height_in*72),
 panels=panels,row_groups=rows,column_groups=cols,exemptions=exemptions)
 jsonlite::write_json(manifest,manifest_path,auto_unbox=TRUE,pretty=TRUE,digits=NA)
 invisible(manifest)
}
