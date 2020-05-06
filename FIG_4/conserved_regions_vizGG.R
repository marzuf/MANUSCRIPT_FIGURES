# IGV style

# Rscript conserved_regions_vizGG.R
# Rscript conserved_regions_vizGG.R norm_vs_tumor
# Rscript conserved_regions_vizGG.R subtypes
# Rscript conserved_regions_vizGG.R wt_vs_mut

# Rscript conserved_regions_viz.R <cmpType>

cat("> START ", "conserved_regions_viz.R", "\n")

startTime <- Sys.time()

plotType <- "svg"



source("../settings.R")

require(ggrepel)
require(ggsci)
tad_col <- pal_d3()(3)[1]
gene_col <- pal_d3()(3)[2]
syn_col <- pal_d3()(3)[3]

consAll_col <- "red"
consAbove_col <- "orange"
consBelow_col <- "black"


myWidth <- myWidth*1.2

outFolder <- file.path("CONSERVED_REGIONS_VIZGG")
dir.create(outFolder, recursive = TRUE)


args <- commandArgs(trailingOnly = TRUE)

if(length(args) == 0) {
  cmpType <- ""
  filePrefix <- ""
  cmpTit <- paste0("all")
} else if(length(args) == 1) {
  cmpType <- args[1]  
  filePrefix <- paste0(cmpType, "_")
  cmpTit <- cmpType
}else {
  stop("---error\n")
}
signif_column <- "adjPvalComb"
signifThresh <- 0.01
signifcol <- paste0(signif_column, "_", signifThresh)
minOverlapBpRatio <- 0.8
minIntersectGenes <- 3
gene_matching_fuse_threshold <- 0.8


setDir <- "/media/electron"
setDir <- ""
entrezDT_file <- paste0(setDir, "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt")
gff_dt <- read.delim(entrezDT_file, header = TRUE, stringsAsFactors = FALSE)
gff_dt$entrezID <- as.character(gff_dt$entrezID)

cond_fc_dt <- get(load(file.path(runFolder, "CREATE_COND_MEANFC", "all_dt.Rdata")))

result_dt <- get(load(file.path(runFolder, "CREATE_FINAL_TABLE/all_result_dt.Rdata")))
nDS <- length(unique(file.path(result_dt$hicds, result_dt$exprds)))

inFolder <- file.path(runFolder, "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2")
outFile <- file.path(inFolder, paste0(filePrefix, "conserved_regions_with_genes_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
stopifnot(file.exists(outFile))
conserved_dt <- get(load(outFile))
conserved_dt$conserved_region <- as.character(conserved_dt$conserved_region)

outFile <- file.path(inFolder, paste0(filePrefix, "conserved_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
stopifnot(file.exists(outFile))
conserved_list <- get(load(outFile))
nConserved <- lengths(conserved_list)

stopifnot(length(nConserved) == nrow(conserved_dt))

maxConserved <- names(which.max(nConserved))

max_dt <- conserved_dt[conserved_dt$conserved_region == maxConserved,,drop=F]
stopifnot(nrow(max_dt) == 1)

all_max_regions <- unlist(strsplit(max_dt$corresp_tads, split=","))
stopifnot(length(all_max_regions) == nConserved[maxConserved])

all_max_entrez <- unlist(strsplit(max_dt$intersect_genes_entrez, split=","))
stopifnot(all_max_entrez %in% gff_dt$entrezID)

nDScons <- length(all_max_regions)

colConsThresh <- ceiling(nDScons/2)

# retrieve all the genes in the corresponding regions

all_symbs <- as.character(unlist(sapply(all_max_regions, function(ds) {
  reg_symb <- unlist(strsplit(x=result_dt$region_genes[result_dt$hicds == dirname(dirname(ds)) & result_dt$exprds == basename(dirname(ds)) & result_dt$region == basename(ds)], split=","))
  stopifnot(length(reg_symb) > 0 )
  reg_symb
})))
symb_count <- setNames(as.numeric(table(all_symbs)), names(table(all_symbs)))

all_genes_starts_ends <- sapply(unique(all_symbs), function(x) {
  nSymb <- as.numeric(symb_count[paste0(x)])
  c(start = gff_dt$start[gff_dt$symbol == x],
    end = gff_dt$end[gff_dt$symbol == x],
    symbol = x,
    count = nSymb,
    col = ifelse(nSymb == nDScons, consAll_col, ifelse(nSymb >= colConsThresh, consAbove_col, consBelow_col))
  )
})
gene_plot_dt <- data.frame(t(all_genes_starts_ends))
rownames(gene_plot_dt) <- NULL
gene_plot_dt$start <- as.numeric(as.character(gene_plot_dt$start))
gene_plot_dt$end <- as.numeric(as.character(gene_plot_dt$end))
stopifnot(is.numeric(gene_plot_dt$start))
stopifnot(is.numeric(gene_plot_dt$end))
stopifnot(gene_plot_dt$chromo == gene_plot_dt$chromo[1])


# all_genes_starts_ends <- all_genes_starts_ends[,order(plot_dt$start)]


all_regions_starts_ends <- sapply(all_max_regions, function(x) {
  g2t_dt <- read.delim(file.path(runFolder, dirname(dirname(x)), "genes2tad", "all_assigned_regions.txt"), header=FALSE, stringsAsFactors=FALSE, col.names=c("chromo", "region", "start", "end"))
  stopifnot(sum(g2t_dt$region == basename(x)) == 1)
  g2t_dt$start[g2t_dt$region == basename(x)]
  c(start=g2t_dt$start[g2t_dt$region == basename(x)], end=g2t_dt$end[g2t_dt$region == basename(x)], chromo = g2t_dt$chromo[g2t_dt$region == basename(x)])
})

all_exprds <- basename(dirname(colnames(all_regions_starts_ends)))
all_regions_starts_ends <- all_regions_starts_ends[,rev(colnames(all_regions_starts_ends)[order(all_exprds)])]

region_plot_dt <- data.frame(t(all_regions_starts_ends))
region_plot_dt$region_id <- rownames(region_plot_dt)
rownames(region_plot_dt) <- NULL
region_plot_dt$hicds <- dirname(dirname(region_plot_dt$region_id))
region_plot_dt$exprds <- basename(dirname(region_plot_dt$region_id))
region_plot_dt$region <- basename(region_plot_dt$region_id)
region_plot_dt$start <- as.numeric(as.character(region_plot_dt$start))
region_plot_dt$end <- as.numeric(as.character(region_plot_dt$end))
stopifnot(is.numeric(region_plot_dt$start))
stopifnot(is.numeric(region_plot_dt$end))
stopifnot(region_plot_dt$chromo == region_plot_dt$chromo[1])

region_plot_dt$hicds_lab <- hicds_names[paste0(region_plot_dt$hicds)]
region_plot_dt$exprds_lab <- exprds_names[paste0(region_plot_dt$exprds)]
stopifnot(!is.na(region_plot_dt$hicds_lab))
stopifnot(!is.na(region_plot_dt$exprds_lab))

region_plot_dt$raw_labels <- paste0(as.character(region_plot_dt$hicds_lab), " - ", as.character(region_plot_dt$exprds_lab))

# save(all_regions_starts_ends,file= "all_regions_starts_ends.Rdata", version=2)
# save(all_genes_starts_ends,file= "all_genes_starts_ends.Rdata", version=2)

# load("CONSERVED_REGIONS_VIZGG/region_plot_dt.Rdata")
# load("CONSERVED_REGIONS_VIZGG/gene_plot_dt.Rdata")

cat(nrow(region_plot_dt), "\n")

region_plot_dt <- merge(region_plot_dt, cond_fc_dt[, c("hicds", "exprds","region", "meanFC") ], by=c("hicds", "exprds", "region"), all.x=T, all.y=F)


cat(nrow(region_plot_dt), "\n")

chromo <- unique(as.character(all_regions_starts_ends["chromo",]))
stopifnot(length(chromo) == 1)

tad_gene_space <- 0.5
gene_line_offset <- 0.1

region_plot_dt <- region_plot_dt[order(region_plot_dt$hicds_lab, region_plot_dt$exprds_lab, decreasing = TRUE),]
region_plot_dt$ds_rank <- 1:nrow(region_plot_dt)

gene_plot_dt <- gene_plot_dt[order(gene_plot_dt$start, gene_plot_dt$end),]
gene_plot_dt$gene_rank <- max(region_plot_dt$ds_rank) + tad_gene_space + 1:nrow(gene_plot_dt)

gene_plot_dt$gene_pos <- 0.5*(gene_plot_dt$start+gene_plot_dt$end)

TADlinecol <- "darkblue"
TADlt <- 1
TADlw <- 2
geneLt <- 1
geneLw <- 1
geneDelLt <- 2
geneDelLw <- 0.5

region_p <- ggplot() + 
  labs(x="", y="")+
  # lines for the TADs
  geom_segment( aes(x = region_plot_dt$start, y = region_plot_dt$ds_rank, 
                    xend=region_plot_dt$end, yend = region_plot_dt$ds_rank), 
                colour = TADlinecol, linetype = TADlt, size = TADlw,
                inherit.aes = F) + 
  # lines for the genes
  geom_segment( aes(x = gene_plot_dt$start, y = gene_plot_dt$gene_rank,
                    xend=gene_plot_dt$end, yend = gene_plot_dt$gene_rank),
                colour = gene_plot_dt$col, linetype = geneLt, size = geneLw,
                inherit.aes = F) +
  # vertical lines for the gene delimiters
  geom_segment( aes(x = c(gene_plot_dt$start, gene_plot_dt$end), 
                    xend= c(gene_plot_dt$start, gene_plot_dt$end),
                    # y = rep(min(region_plot_dt$ds_rank)-gene_line_offset, 2),
                    y = rep(min(region_plot_dt$ds_rank)-gene_line_offset, 1),
                    # yend=rep(gene_plot_dt$gene_rank, 2)), 
                yend=rep(gene_plot_dt$gene_rank, 2)), 
                colour = rep(gene_plot_dt$col,2), linetype = geneDelLt, size = geneDelLw,
                inherit.aes = F) +
  theme_void()

region_p2 <-  region_p + 
geom_text_repel(
  aes(x = gene_plot_dt$gene_pos, y =  gene_plot_dt$gene_rank, label = gene_plot_dt$symbol), inherit.aes = F,
  font_face="italic",
  # nudge_x      = 0.05,
  # direction    = "y",
  nudge_y      = 2,
  direction    = "y",
  hjust        = 0.5,
  segment.size = 1
) 

region_plot_dt$raw_labels <- paste0(as.character(region_plot_dt$hicds_lab), " - ", as.character(region_plot_dt$exprds_lab))
region_plot_dt$ds_lab <- sapply(1:nrow(region_plot_dt), function(i) {
  label_part1 <- gsub("(.+) (.+) vs. (.+)", "\\1", region_plot_dt$raw_labels[i])
  label_part2 <- gsub("(.+) (.+) vs. (.+)", "\\2", region_plot_dt$raw_labels[i])
  label_part3 <- gsub("(.+) (.+) vs. (.+)", "\\3", region_plot_dt$raw_labels[i])
  if(region_plot_dt$meanFC[i] > 0){
    mylab <- gsub(" ","~", paste0(label_part1, "~", label_part2, "~vs.~bold(", label_part3, ")"))
  }else {
    mylab <- gsub(" ","~", paste0(label_part1, "~bold(", label_part2, ")~vs.~", label_part3))
  }
  
  mylab <- gsub("_", "[{\"-\"}]*", mylab)
  
  mylab
})

# region_plot_dt$ds_lab[27]=region_plot_dt$ds_lab[10]
# region_plot_dt$ds_lab[27]="786[{\"-\"}]*O~-~KICH~bold(normal)~vs.~tumor"
# region_plot_dt$ds_lab[1:26]="."

axisOffset <- 2*min(region_plot_dt$end-region_plot_dt$start)#90000

region_p3 <- region_p2 + xlim(min(c(region_plot_dt$start, gene_plot_dt$start)- axisOffset), NA) +
        geom_text_repel(
          aes(x = region_plot_dt$start, y =  region_plot_dt$ds_rank, label = region_plot_dt$ds_lab), inherit.aes = F,
          nudge_x       = 3.5 - region_plot_dt$start,
          direction     = "y",
          hjust         = 1, parse = T, size=3,
          force_pull   = 0
        ) +
  
  theme(plot.margin =margin(t = 0, r = 0, b = 0, l = 10, unit = "pt") )



outFile <- file.path(outFolder, paste0(maxConserved, "_viz.", plotType))
ggsave(region_p3, filename=outFile, height=myHeightGG, width=myWidthGG *1.2)
cat(paste0("... written: ", outFile, "\n"))


stop("-ok")


for(i in 1:nrow(region_plot_dt)) {
  
  label_part1 <- gsub("(.+) (.+) vs. (.+)", "\\1", region_plot_dt$raw_labels[i])
  label_part2 <- gsub("(.+) (.+) vs. (.+)", "\\2", region_plot_dt$raw_labels[i])
  label_part3 <- gsub("(.+) (.+) vs. (.+)", "\\3", region_plot_dt$raw_labels[i])
  
  if(region_plot_dt$meanFC[i] > 0){
    mylab <-bquote(.(label_part1)~.(label_part2)~' vs. '~bold(.(label_part3)))  
  }else {
    mylab <- bquote(.(label_part1)~bold(.(label_part2))~' vs. '~.(label_part3))
  }
  # cat(mylab,"\n")
  
  text(
    x = region_plot_dt$start[i],
    y = dsPos[i] + textOffset,
    # labels = colnames(all_regions_starts_ends),
    # labels = paste0(region_plot_dt$hicds_lab, " - ", region_plot_dt$exprds_lab) , #dirname(colnames(all_regions_starts_ends)),
    labels = mylab,
    # cex = 0.5,
    cex = 0.7,
    pos=2,
    col = tad_col
  ) 
  
  
}








































dsSpace <- 0.5
geneSpace <- 0.8

ySynt <- 0.1
syntLwd <- 2

tadOffset <- 50000

textOffset <- 0.05

xStart <- min(region_plot_dt$start) - tadOffset
xEnd <- max(region_plot_dt$end) + tadOffset

yOffset <- 0.3   
synmatch_start <- NULL
synmatch_end <- NULL

dsPos <- seq(from=0, by=dsSpace, length.out=ncol(all_regions_starts_ends))
genePos <- seq(from = max(dsPos) + dsSpace, by=geneSpace, length.out=ncol(all_genes_starts_ends))
# genePos <- seq(from = max(dsPos) + dsSpace, by=0, length.out=ncol(all_genes_starts_ends))  # all genes on same y position
yStart <- min(c(dsPos, genePos)) - yOffset
yEnd <- max(c(dsPos, genePos)) + yOffset

subTit <- paste0("(# DS conserv. = ", nDScons, ")")

outFile <- file.path(outFolder, paste0(maxConserved, "_viz.", plotType))
do.call(plotType, list(outFile, height=myHeight, width=myWidth*2))
dev.control(displaylist="enable")
initMar <- par()$mar
par(mar=initMar+c(0,10,0,0))
par(family="Heshley")
par(xpd=TRUE)
plot(NULL,
     main = paste0(maxConserved),
     xlim = c(xStart, xEnd),
     # ylim = c(yStart, yEnd),
     ylim = c(0, yEnd),
     xlab = "",
     # xlab = paste0(gsub("chr", "chromosome ", chromo)),
     ylab = "",
     axes = FALSE,
     cex.main = plotCex
)
mtext(side = 3, text = subTit)
mtext(side = 1, text= paste0(gsub("chr", "chromosome ", chromo)), font = 2, cex = 1, line=2)
axis(1,
     at = unique(sort(c(region_plot_dt$start, region_plot_dt$end))), 
     cex = 0.6)
# draw the tads
segments(
  x0 = region_plot_dt$start,
  y0 = dsPos,
  x1 = region_plot_dt$end,
  y1 = dsPos,
  col=tad_col
)
# text(
#   x = region_plot_dt$start,
#   y = dsPos + textOffset,
#   # labels = colnames(all_regions_starts_ends),
#   # labels = paste0(region_plot_dt$hicds_lab, " - ", region_plot_dt$exprds_lab) , #dirname(colnames(all_regions_starts_ends)),
#   labels = region_plot_dt$labels,
#   # cex = 0.5,
#   cex = 0.7,
#   pos=2,
#   col = tad_col
# )

cat(nrow(region_plot_dt), "\n")

for(i in 1:nrow(region_plot_dt)) {
  
  label_part1 <- gsub("(.+) (.+) vs. (.+)", "\\1", region_plot_dt$raw_labels[i])
  label_part2 <- gsub("(.+) (.+) vs. (.+)", "\\2", region_plot_dt$raw_labels[i])
  label_part3 <- gsub("(.+) (.+) vs. (.+)", "\\3", region_plot_dt$raw_labels[i])
  
  if(region_plot_dt$meanFC[i] > 0){
    mylab <-bquote(.(label_part1)~.(label_part2)~' vs. '~bold(.(label_part3)))  
  }else {
    mylab <- bquote(.(label_part1)~bold(.(label_part2))~' vs. '~.(label_part3))
  }
  # cat(mylab,"\n")
  
  text(
    x = region_plot_dt$start[i],
    y = dsPos[i] + textOffset,
    # labels = colnames(all_regions_starts_ends),
    # labels = paste0(region_plot_dt$hicds_lab, " - ", region_plot_dt$exprds_lab) , #dirname(colnames(all_regions_starts_ends)),
    labels = mylab,
    # cex = 0.5,
    cex = 0.7,
    pos=2,
    col = tad_col
  ) 
  
  
}





segments(
  x0 = gene_plot_dt$start,
  y0 = genePos,
  x1 = gene_plot_dt$end,
  y1 = genePos,
  col=all_genes_starts_ends["col",]
)
text(
  x = 0.5*(gene_plot_dt$start + gene_plot_dt$end),
  y = genePos - 5*textOffset,
  # y = genePos + 0,
  labels = as.character(all_genes_starts_ends["symbol",]),
  cex = 1,
  pos=3,
  col=as.character(all_genes_starts_ends["col",])
)

segments(x0=c(gene_plot_dt$start, gene_plot_dt$end),
         y0=0-yOffset,
         x1=c(gene_plot_dt$start, gene_plot_dt$end),
         y1 = rep(genePos,2),
         lty=2, 
         col = as.character(all_genes_starts_ends["col",])
)


legend(
  "bottom",
  horiz = TRUE,
  legend = c(
    paste0("# cons. = ", nDScons), 
    paste0("# cons. >= ", colConsThresh),
    paste0("# cons. < ", colConsThresh)
  ),
  text.col = c(consAll_col, consAbove_col, consBelow_col), 
  inset = c(-0.0, -0.2),
  bty = "n",
  xpd = TRUE
)

vizplot <- recordPlot()


invisible(dev.off())
