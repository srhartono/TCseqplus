#!/usr/bin/env Rscript

if (commandArgs()[1] != "RStudio") {

  ARGS <- c(
    "tlxfile", "character", "",
    "output","character",""
  )

  OPTS <- c(
    "as.filter","logical",TRUE,"set to FALSE to use filters as selectors instead",
    "remove.adapter","logical",TRUE,"remove adapter junctions (not when first junction - baitonly)",
    "f.unaligned","character","","",
    "f.baitonly","character","","",
    "f.uncut","character","","",
    "f.misprimed","character","","",
    "f.freqcut","character","","",
    "f.largegap","character","","",
    "f.mapqual","character","","",
    "f.breaksite","character","","",
    "f.sequential","character","","",
    "f.repeatseq","character","","",
    "f.duplicate","character","",""
  )

  source_local <- function(fname){
    argv <- commandArgs(trailingOnly = FALSE)
    base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
    source(paste(base_dir, fname, sep="/"))
  }

  source_local("Rsub.R")
  parseArgs("TranslocFilter.R", ARGS, OPTS)

} else {
  source("~/AltLab/transloc_pipeline/R/Rsub.R")
  source("~/AltLab/transloc_pipeline/R/TranslocHelper.R")


  tlxfile <- "~/AltLab/Data/Alt055/results-new3/RF204_Alt055/RF204_Alt055.tlx"
  output <- "~/AltLab/Data/Alt055/results-new3/RF204_Alt055/RF204_Alt055_filtered.txt"

  as.filter <- TRUE
  remove.adapter <- TRUE
  f.unaligned <- "1"
  f.baitonly <- "1"
  f.uncut <- "1"
  f.misprimed <- "L10"
  f.freqcut <- "1"
  f.largegap <- "G30"
  f.mapqual <- "1"
  f.breaksite <- "1"
  f.sequential <- "1"
  f.repeatseq <- "1"
  f.duplicate <- "1"
}

suppressPackageStartupMessages(library(readr, quietly=TRUE))
suppressPackageStartupMessages(library(dplyr, quietly=TRUE))

stats.file <- paste(gsub("^(.+)\\.[a-zA-Z0-9]+$","\\1",output,perl=T),"_stats.txt",sep="")
#stats.file <- sub(paste(".",file_ext(output),sep=""),"_stats.txt",output)


filter.names <- c("unaligned","baitonly","uncut","misprimed","freqcut","largegap","mapqual","breaksite","sequential","repeatseq","duplicate")
filter.values <- c()

for (filter.name in filter.names) {

  tmp.value <- get(paste("f.",filter.name,sep=""))

  if (grepl("^[0-9]+$",tmp.value)) {

#	curr.tmp.value = FALSE
#	if (tmp.value == 1) {curr.tmp.value = TRUE}
#    filter.values[filter.name] <- paste("==",curr.tmp.value,sep="")
#    filter.values[filter.name] <- paste("==",tmp.value,sep="")
    filter.values[filter.name] <- tmp.value

  } else if (grepl("^[GL]?[E]?[0-9]+$",tmp.value)) {

    tmp.value <- sub("[Gg][Ee]",">=",tmp.value)
    tmp.value <- sub("[Ll][Ee]","<=",tmp.value)
    tmp.value <- sub("[Gg]",">",tmp.value)
    tmp.value <- sub("[Ll]","<",tmp.value)
    tmp.value <- sub("[Ee]","==",tmp.value)

    filter.values[filter.name] <- tmp.value

  } else if (tmp.value == "") {
    filter.values[filter.name] <- tmp.value
  } else {
    stop(paste("Error:",filter.name,"filter entered in wrong format"))
  }
}

# tlx <- fread(tlxfile,sep="\t",header=T,select=c("Qname","JuncID","Rname",filter.names))

filter.cols <- rep("i", length(filter.names))
names(filter.cols) <- filter.names
tlx <- read_tsv(tlxfile,
                col_types=do.call(cols_only, c(list("Qname"="c",
                                                    "JuncID"="i",
                                                    "Rname"="c"),
                                               filter.cols)))

#print(tlx)
if (remove.adapter) {
#  tlx <- filter(tlx,!(JuncID > 1 & Rname == "Adapter"))
  tlx <- filter(tlx,!(tlx$JuncID > 1 & tlx$Rname == "Adapter"))
}

stats.names <- c("total",filter.names,"result")
filter.stats <- data.frame(reads=rep(0,length(stats.names)),junctions=rep(0,length(stats.names)),row.names=stats.names)

reads.total <- n_distinct(tlx$Qname)
#junctions.total <- nrow(filter(tlx, tlx$Rname != "" & tlx$Rname != "Adapter"))
junctions.total <- nrow(filter(tlx, Rname != "" & Rname != "Adapter"))

filter.stats["total",] <- c(reads.total,junctions.total)

tlx.filt.list <- list()
myreads.current = reads.total
myind = 0
for (filter.name in filter.names) {
	myind = myind + 1;
	myindprint = paste(myind,".",sep="")
  filter.value <- filter.values[filter.name]
  if (filter.value == "") {
		print(paste(myindprint,"filter.name =",filter.name,"; filter.value =",filter.value,"(SKIPPED) ; reads.total =",reads.total,"; reads.current =",myreads.current,"; reads.filtered.out =",0))#mytotalfilterline))
		next
	}
  filter.text <- paste(filter.name,filter.value,sep="")


#  tlx.filt.list[[filter.name]] <- filter(tlx, filter.text) %>% select(tlx$Qname,tlx$JuncID,tlx$Rname)
	mytotalfilterline = dim(tlx)
#	print(paste("Current filter text = ",head(filter.text)))
#	print(paste("head tlx = ",head(tlx)))
	currtlx = filter(tlx,tlx[,colnames(tlx) == filter.name] == filter.value)
  tlx.filt.list[[filter.name]] <- currtlx %>% select(Qname,JuncID,Rname)
#  tlx.filt.list[[filter.name]] <- filter(tlx, filter.text) %>% select(Qname,JuncID,Rname)

#  junctions.count <- nrow(filter(tlx.filt.list[[filter.name]],Rname != "" & Rname != "Adapter"))
  junctions.count <- nrow(filter(tlx.filt.list[[filter.name]],Rname != "" & Rname != "Adapter"))
  reads.count <- as.integer(summarize(tlx.filt.list[[filter.name]],n_distinct(Qname)))
	myreads.current = myreads.current - reads.count
	print(paste(myindprint,"filter.name =",filter.name,"; filter.value =",filter.value,"; reads.total =",reads.total,"; reads.current =",myreads.current,"; reads.filtered.out =",reads.count))#mytotalfilterline))
  filter.stats[filter.name,] <- c(reads.count,junctions.count)

}
tlx.filt <- bind_rows(tlx.filt.list)

if (as.filter) {
  tlx <- anti_join(tlx,tlx.filt,by=c("Qname","JuncID"))

} else {
  tlx <- semi_join(tlx,tlx.filt,by=c("Qname","JuncID"))
}

reads.count <- as.integer(summarize(tlx,n_distinct(Qname)))
junctions.count <- nrow(filter(tlx,Rname != "" & Rname != "Adapter"))
filter.stats["result",] <- c(reads.count,junctions.count)

filter.stats <- filter.stats %>% mutate(filter = rownames(.)) %>%
    select(filter, reads, junctions)

tlx <- select(tlx,Qname,JuncID)

print(paste("Output =",output))
print(paste("Filter stats =",stats.file))
write_tsv(tlx, output, col_names=F)
write_tsv(filter.stats, stats.file)
