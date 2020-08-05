# IGV style
#Rscript conserved_regions_vizGG_v2_withImmune.R all conserved_region_130 
# Rscript conserved_regions_vizGG_v2.R all <conserved_region_130>  <conserved_region_42> 
# Rscript conserved_regions_vizGG_v2.R norm_vs_tumor <conserved_region_22> <conserved_region_9
# Rscript conserved_regions_vizGG_v2.R subtypes
# Rscript conserved_regions_vizGG_v2.R wt_vs_mut

# Rscript conserved_regions_viz2_withImmune.R cmpType <regions_to_plot>

cat("> START ", "conserved_regions_vizGG_v2_withImmune.R", "\n")

startTime <- Sys.time()

plotType <- "svg"

source("../settings.R")

require(ggrepel)
require(ggsci)
require(foreach)
require(GenomicRanges)
require(doMC)
registerDoMC(40)
tad_col <- pal_d3()(3)[1]
gene_col <- pal_d3()(3)[2]
syn_col <- pal_d3()(3)[3]

# consAll_col <- "red"
consAll_col <- "darkgreen"
consAbove_col <- "orange"
consBelow_col <- "black"


myWidth <- myWidth*1.2

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) > 0)

if(args[1] == "all"){
  cmpType <- ""
  filePrefix <- ""
  cmpTit <- paste0("all")
} else {
  cmpType <- args[1]  
  filePrefix <- paste0(cmpType, "_")
  cmpTit <- cmpType
  stopifnot(cmpType %in% c("subtypes", "wt_vs_mut", "norm_vs_tumor"))
}

outFolder <- file.path("CONSERVED_REGIONS_VIZGG_V2_WITHIMMUNE", cmpType)
dir.create(outFolder, recursive = TRUE)

if(length(args) > 1){
  conserved_regions_to_plot <- args[2:length(args)]  
} else {
  conserved_regions_to_plot <- NA
}


immune_ds <- c("AWS_GM12878",  "AWS_erythroid_progenitor","HG00514_B_Lymphocyte", 
               "HG00512_B_Lymphocyte","HG00733_B_Lymphocyte",  "HG00732_B_Lymphocyte","HG00731_B_Lymphocyte", 
               "GM19240_B_Lymphocyte","GM19238_B_Lymphocyte", "AWS_T_cell",  "GSM3112403_MDM_6h_mock_r1")
immune_ds <- gsub("AWS_", "", immune_ds)
immune_ds <- paste0(immune_ds, "_40kb")


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

tmp_dt <- result_dt
tmp_dt$cmpType <- all_cmps[tmp_dt$exprds]
if(cmpTit != "all") {
  tmp_dt <- tmp_dt[tmp_dt$cmpType == cmpType,]
  stopifnot(nrow(tmp_dt) > 0)
}
nCmps <- length(unique(file.path(tmp_dt$hicds, tmp_dt$exprds)))

inFolder <- file.path(runFolder, "TAD_MATCHING_SIGNIF_ACROSS_HICDS_ALLMATCH_v2", cmpType)
outFile <- file.path(inFolder, paste0(filePrefix, 
                                      "conserved_regions_with_genes_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
stopifnot(file.exists(outFile))
conserved_dt <- get(load(outFile))
conserved_dt$conserved_region <- as.character(conserved_dt$conserved_region)

outFile <- file.path(inFolder, 
                     paste0(filePrefix, "conserved_signif_tads", signif_column, signifThresh, "_minBpRatio", minOverlapBpRatio, "_minInterGenes", minIntersectGenes, ".Rdata"))
stopifnot(file.exists(outFile))
conserved_list <- get(load(outFile))
nConserved <- lengths(conserved_list)

stopifnot(length(nConserved) == nrow(conserved_dt))

if(is.na(conserved_regions_to_plot)) {
  conserved_regions_to_plot <- names(nConserved)[nConserved == max(nConserved)]
} 


for(maxConserved in conserved_regions_to_plot) {
  
  
  max_dt <- conserved_dt[conserved_dt$conserved_region == maxConserved,,drop=F]
  stopifnot(nrow(max_dt) == 1)
  
  all_max_regions <- unlist(strsplit(max_dt$corresp_tads, split=","))
  stopifnot(length(all_max_regions) == nConserved[maxConserved])
  
  all_max_entrez <- unlist(strsplit(max_dt$intersect_genes_entrez, split=","))
  stopifnot(all_max_entrez %in% gff_dt$entrezID)
  
  nDScons <- length(all_max_regions)
  
  
  
  colConsThresh <- ceiling(nDScons/2)
  
  # geneColSetting <- setNames( c(consAll_col, consAbove_col, consBelow_col), 
  #   c(paste0("# cons. = ", nDScons), 
  #   paste0("# cons. >= ", colConsThresh),
  #   paste0("# cons. < ", colConsThresh)
  #   )
  # )
  geneColSetting <- setNames( c(consAll_col, consAbove_col, consBelow_col), 
                              c(paste0("= ", nDScons, "/", nDScons), 
                                paste0(">= ", colConsThresh, "/", nDScons), 
                                paste0("< ", colConsThresh, "/", nDScons)
                              )
  )
  
  # geneColSetting <- c("black"# cons." = "black", "red" = "red")  # name is legend name
  
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
      # col = ifelse(nSymb == nDScons, paste0("# cons. = ", nDScons), 
      #              ifelse(nSymb >= colConsThresh, paste0("# cons. >= ", colConsThresh),   paste0("# cons. < ", colConsThresh)))
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
  
  region_plot_dt <- merge(region_plot_dt, cond_fc_dt[, c("hicds", "exprds","region", "meanFC") ], 
                          by=c("hicds", "exprds", "region"), all.x=T, all.y=F)
  
  
  
  gene_reg_dt <- merge(result_dt[,c("hicds", "exprds", "region", "region_genes", "meanLogFC", "adjPvalComb"),], 
                       region_plot_dt[,c("hicds", "exprds", "region")], by=c("hicds", "exprds", "region"),
                       all.x=F, all.y=T)
  stopifnot(!is.na(gene_reg_dt))
  gene_reg_dt <- gene_reg_dt[order(gene_reg_dt$adjPvalComb),]
  outFile <- file.path(outFolder, paste0(filePrefix, maxConserved, "_allDatasets_withSymbs.txt"))
  write.table(gene_reg_dt, file = outFile, col.names = T, row.names=F, quote=F, sep="\t")
  
  
  cat(nrow(region_plot_dt), "\n")
  
  chromo <- unique(as.character(all_regions_starts_ends["chromo",]))
  stopifnot(length(chromo) == 1)
  
  
  myTit <- paste0(chromo, ": ", maxConserved , " - ", cmpTit , " datasets (n=", nCmps, ")")
  subTit <- paste0("conserved in ", nDScons, " datasets")
  
  
  tad_gene_space <- 0.5
  gene_line_offset <- 0.1
  
  # region_plot_dt <- region_plot_dt[order(region_plot_dt$hicds_lab, region_plot_dt$exprds_lab, decreasing = TRUE),]
  # region_plot_dt$ds_rank <- 1:nrow(region_plot_dt)
  # ### CHANGE HERE TO HAVE RANK BY CMP TYPE and space between them
  region_plot_dt$tad_size <- region_plot_dt$end- region_plot_dt$start
  region_plot_dt$cmpType <- all_cmps[paste0(region_plot_dt$exprds)]
  stopifnot(!is.na(region_plot_dt$cmpType))
  region_plot_dt$cmpType <- factor(region_plot_dt$cmpType, levels=sort(unique(as.character(all_cmps))))
  stopifnot(!is.na(region_plot_dt$cmpType))
  region_plot_dt$cmpType_num <- as.numeric(region_plot_dt$cmpType)-1
  # region_plot_dt <- region_plot_dt[order(-region_plot_dt$cmpType_num, region_plot_dt$hicds_lab, region_plot_dt$exprds_lab, decreasing = TRUE),]
  region_plot_dt <- region_plot_dt[order(-region_plot_dt$cmpType_num, region_plot_dt$tad_size, decreasing=TRUE),]
  # region_plot_dt$ds_rank <- 1:nrow(region_plot_dt) ## CHANGED HERE
  region_plot_dt$ds_rank <- 1:nrow(region_plot_dt) + region_plot_dt$cmpType_num
  
  # for label colors
  region_plot_dt$ds_col <- all_cols[as.character(region_plot_dt$cmpType)]
  stopifnot(!is.na(region_plot_dt$ds_col))
  
  
  
  
  
  
  
  
  
  
  
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
    
    mylab <- gsub("_", "[{\"-\"}]*", mylab)  # underscore are not recognize -> replace with dash lowerscript
    
    # if the first caracter are numbers -> need to separate ! -> should match 22Rv1 but not 786[]
    # mylab <- gsub("^(\\d+)", "\\1*", mylab)
    mylab <- gsub("^(\\d+)([a-zA-Z])", "\\1*\\2", mylab)
    
    mylab
  })
  
  #### ADD HERE THE IMMUNE DATA FOR PLOTTING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ds=immune_ds[1]
  immuneFolder <- file.path("../../IMMUNE_DATA_v2_Yuanlong")
  immuneCol <- "blueviolet"
  
  plstart <- min(region_plot_dt$start)
  plend <- max(region_plot_dt$end)
  stopifnot(plend > plstart)
  plchr <- unique(region_plot_dt$chromo)
  stopifnot(length(plchr) == 1)
  immune_region_plot_dt <- foreach(ds = immune_ds, .combine='rbind') %dopar% {
    
    region_dt <- read.delim(file.path(immuneFolder, ds, "genes2tad", "all_assigned_regions.txt"),
                            stringsAsFactors = FALSE, header=F,
                            col.names = c("chromo", "region", "start", "end"))
    region_dt <- region_dt[grepl("_TAD", region_dt$region),]
    sub_dt <- region_dt[region_dt$start <= plend & region_dt$end >= plstart & region_dt$chromo == plchr,]
    # should have following columns:
    # [1] "hicds"       "exprds"      "region"      "start"       "end"        
    # [6] "chromo"      "region_id"   "hicds_lab"   "exprds_lab"  "raw_labels" 
    # [11] "meanFC"      "tad_size"    "cmpType"     "cmpType_num" "ds_rank"    
    # [16] "ds_col"      "ds_lab"     
    sub_dt$hicds <- ds
    sub_dt$exprds <- ""
    sub_dt$region_id <- file.path(sub_dt$hicds, sub_dt$exprds, sub_dt$region)
    sub_dt$hicds_lab <- gsub("_40kb", "", sub_dt$hicds)
    sub_dt$exprds_lab <- sub_dt$exprds
    sub_dt$raw_labels <- paste0(sub_dt$hicds, sub_dt$exprds)
    sub_dt$meanFC <- NA
    sub_dt$tad_size <- sub_dt$end - sub_dt$start + 1
    sub_dt$cmpType <- "immuneCells"
    sub_dt$cmpType_num <- max(region_plot_dt$cmpType_num) + 1
    sub_dt$ds_rank <- max(region_plot_dt$ds_rank) + which(immune_ds == ds)
    sub_dt$ds_col <- immuneCol
    sub_dt$ds_lab <- paste(sub_dt$hicds_lab, sub_dt$exprds_lab)
    
    ### TAKE ONLY THE MOST OVERLAPPING
    # sub_dt <- sub_dt[sub_dt$tad_size == max(sub_dt$tad_size),]
    
    ## Using pairs to find intersection of overlapping ranges
    query <- IRanges(start=plstart, end=plend) 
    subject <- IRanges(start = sub_dt$start, end=sub_dt$end)
    hits <- findOverlaps(query, subject)
    p <- Pairs(query, subject, hits=hits)
    v0 <- which.max(   width(pintersect(p)))
    
    v1 <- which.max( sapply(1:nrow(sub_dt), function(x) {
      ## Using pairs to find intersection of overlapping ranges
      query <- IRanges(start=plstart, end=plend) 
      subject <- IRanges(start = sub_dt[x,"start"], end=sub_dt[x,"end"])
      hits <- findOverlaps(query, subject)
      p <- Pairs(query, subject, hits=hits)
      width(pintersect(p))
    }))
    stopifnot(v0 == v1)
    
    sub_dt[v0,]
  }
  
  save(immune_region_plot_dt, file="immune_region_plot_dt.Rdata", version=2)
  save(region_plot_dt, file="region_plot_dt.Rdata", version=2)
  
  
  stopifnot(setequal(colnames(immune_region_plot_dt), colnames(region_plot_dt)))
  
  immune_region_plot_dt <- immune_region_plot_dt[,colnames(region_plot_dt)]
  
  region_plot_dt <- rbind(region_plot_dt, immune_region_plot_dt)
  
  
  
  
  
  
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
  
  gene_plot_dt$col <- factor(gene_plot_dt$col, levels = as.character(geneColSetting))
  
  xscale <- seq(from=min(region_plot_dt$start) , to=max(region_plot_dt$end) , length.out=10)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  region_p <- ggplot() + 
    
    ggtitle(myTit, subtitle=subTit)+
    
    labs(x="", y="")+
    # lines for the TADs
    geom_segment( aes(x = region_plot_dt$start, y = region_plot_dt$ds_rank, 
                      xend=region_plot_dt$end, yend = region_plot_dt$ds_rank), 
                  colour = TADlinecol, linetype = TADlt, size = TADlw,
                  inherit.aes = F) + 
    # lines for the genes
    geom_segment( aes(x = gene_plot_dt$start, y = gene_plot_dt$gene_rank,
                      xend=gene_plot_dt$end, yend = gene_plot_dt$gene_rank,
                      
                      # colour = geneColSetting[paste0(gene_plot_dt$col)]
                      colour = gene_plot_dt$col
                      
                      
                      ),
                  # colour = gene_plot_dt$col, 
                  linetype = geneLt, size = geneLw,
                  
                  show.legend=TRUE,
                  
                  inherit.aes = F) +
    
    scale_color_manual(values= setNames(as.character(geneColSetting),as.character(geneColSetting)), 
                       labels = setNames(names(geneColSetting), as.character(geneColSetting))) +
    labs(color="# conserved")+
  
    
    # vertical lines for the gene delimiters
    geom_segment( aes(x = c(gene_plot_dt$start, gene_plot_dt$end), 
                      xend= c(gene_plot_dt$start, gene_plot_dt$end),
                      # y = rep(min(region_plot_dt$ds_rank)-gene_line_offset, 2),
                      y = rep(min(region_plot_dt$ds_rank)-gene_line_offset, 1),
                      # yend=rep(gene_plot_dt$gene_rank, 2)), 
                  yend=rep(gene_plot_dt$gene_rank, 2)), 
                  colour = rep(as.character(gene_plot_dt$col),2), linetype = geneDelLt, size = geneDelLw,
                  inherit.aes = F) +
    
    
    theme_void()
  
  save(gene_plot_dt, file="gene_plot_dt.Rdata", version=2)
  
  region_p2 <-  region_p + 
  geom_text_repel(
    aes(x = gene_plot_dt$gene_pos, y =  gene_plot_dt$gene_rank, label = gene_plot_dt$symbol,
        color=gene_plot_dt$col), inherit.aes = F, show.legend=F,
    font_face="italic",
    # nudge_x      = 0.05,
    # direction    = "y",
    nudge_y      = 2,
    direction    = "y",
    hjust        = 0.5,
    segment.size = 1
  )

  
  
  # region_plot_dt$ds_lab[27]=region_plot_dt$ds_lab[10]
  # region_plot_dt$ds_lab[27]="786[{\"-\"}]*O~-~KICH~bold(normal)~vs.~tumor"
  # region_plot_dt$ds_lab[1:26]="."
  
  # region_plot_dt$ds_lab[6]="22Rv1~-~PRAD~bold(normal)~vs.~tumor"
  # should be 
  # region_plot_dt$ds_lab[6]="22*Rv1~-~PRAD~bold(normal)~vs.~tumor"
  
  axisOffset <- min(region_plot_dt$end-region_plot_dt$start)#90000
  
  region_plot_dt$ds_lab[6]
  


#  geom_label_repel(
#    arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "first"),
#    force = 10,
#    xlim  = x_limits  ### => all labels start at 3 if x_limits = c(3, NA)
#  ) +


  region_p3 <- region_p2 + xlim(min(c(region_plot_dt$start, gene_plot_dt$start)- axisOffset), NA) +
          geom_text_repel(
            aes(x = region_plot_dt$start, y =  region_plot_dt$ds_rank, label = region_plot_dt$ds_lab), inherit.aes = F,
            nudge_x       = 3.5 - region_plot_dt$start,
            direction     = "y",
            color = region_plot_dt$ds_col, # ADDED 
            hjust         = 1, parse = T, size=3,
            force_pull   = 0
          ) +
    
    theme(plot.margin =margin(t = 25, r = 25, b = 10, l = 10, unit = "pt"),
          text = element_text(family=fontFamily),
          plot.title = element_text(hjust=0.5, size = 16, face="bold"),
          plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
          legend.text = element_text(size=10),
          legend.title =  element_text(size=12, face ="bold")
          ) +
    # scale_x_continuous(labels =xscale, breaks=xscale , name="") +
    # theme(axis.line.x = element_()) +
    expand_limits(x = c(NA, max(xscale)+10000))+
  geom_segment(aes(x=min(xscale), xend=max(xscale),y=0, yend=0), lineend="round", colour="darkgrey", size=2)+  
  # annotate("text", x=c(range(xscale)), y = -0.8, vjust=1, 
  #          size=5,
  #          label=range(xscale) , colour="darkgrey", fontface="italic") +
    annotate("text", x=c(min(xscale), 0.5*(min(xscale)+max(xscale)), max(xscale)), y = -0.8, vjust=1, 
             size=5,
             label=c(min(xscale), chromo, max(xscale)) , colour="darkgrey", fontface="italic") 
    
  # region_p3+scale_color_manual(values= c(consAll_col, consAbove_col, consBelow_col),  guide = 'legend')
    
  save(gene_plot_dt, file="gene_plot_dt.Rdata", version=2)
  save(region_plot_dt, file="region_plot_dt.Rdata", version=2)

genes_plot_dt <- gene_plot_dt
saveFile <- file.path(outFolder, paste0("fig4C_", maxConserved, "_genes_plot_dt.Rdata"))
save(genes_plot_dt, file=saveFile, version=2)
cat(paste0("... written:" , saveFile, "\n"))

tads_plot_dt <- region_plot_dt
saveFile <- file.path(outFolder, paste0("fig4C_", maxConserved, "_tads_plot_dt.Rdata"))
save(tads_plot_dt, file=saveFile, version=2)
cat(paste0("... written:" , saveFile, "\n"))


  
  save(region_p3, file="region_p3.Rdata", version=2)
  
  outFile <- file.path(outFolder, paste0(filePrefix, maxConserved, "_viz.", plotType))
  ggsave(region_p3, filename=outFile, height=myHeightGG*1.2, width=myWidthGG *2)
  cat(paste0("... written: ", outFile, "\n"))
  
}
