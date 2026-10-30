#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ape))
a <- commandArgs(TRUE)
stopifnot(length(a)%in%c(3L,4L))
tree <- read.nexus(a[1]); if(inherits(tree,"multiPhylo")) stop("Expected one MCC tree")
if(length(a)==4L) {
    target<-read.nexus(a[4])
    if(!isTRUE(all.equal(tree,target,use.edge.length=FALSE)))stop("CA output changed the target topology")
}
if(length(tree$tip.label)!=as.integer(Sys.getenv("PBP_EXPECT_TIPS","533")))stop("Unexpected MCC sample count")
dates <- read.delim(a[2],check.names=FALSE,stringsAsFactors=FALSE)
pick <- function(options) {j<-match(tolower(options),tolower(names(dates)));j<-j[!is.na(j)];if(!length(j))stop("Date column missing");names(dates)[j[1]]}
id <- pick(c("taxon","sample","sample_id","strain","name"))
dc <- pick(c("date","year","decimal_date","decimal_year"))
dd <- setNames(suppressWarnings(as.numeric(dates[[dc]])),dates[[id]])[tree$tip.label]
if(anyDuplicated(tree$tip.label)||any(!is.finite(dd)))stop("MCC tips not uniquely dated")
if(any(!is.finite(tree$edge.length))||any(tree$edge.length < -1e-10))stop("MCC has invalid/negative branches")
depth<-node.depth.edgelength(tree)[seq_along(dd)]
roots<-dd-depth
if(diff(range(roots))>1e-5)stop("MCC branch lengths disagree with tip dates: root-date range=",diff(range(roots)))
write.table(data.frame(file=a[1],n_tips=length(dd),root_date=median(roots),
  root_date_range=diff(range(roots)),zero_edges=sum(tree$edge.length<=0),
  min_edge=min(tree$edge.length),height_method="common_ancestor"),a[3],sep="\t",row.names=FALSE,quote=FALSE)
cat("PASS: valid chronogram and exact tip-date calibration:",a[1],"\n")
