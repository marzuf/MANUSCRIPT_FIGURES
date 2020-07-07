


load("../data/conserved_region_130_tads_plot_dt.RData")
load("../data/conserved_region_130_genes_plot_dt.RData")

function(genes_plot_dt, tads_plot_dt, myTit="Conservation across datasets"){}

genes_plot_dt$start <- as.numeric(as.character(genes_plot_dt$start))
stopifnot(!is.na(genes_plot_dt$start))
genes_plot_dt$end <- as.numeric(as.character(genes_plot_dt$end))
stopifnot(!is.na(genes_plot_dt$end))
genes_plot_dt$chromo <- as.character(genes_plot_dt$chromo)
genes_plot_dt$count <- as.numeric(as.character(genes_plot_dt$count))
stopifnot(!is.na(genes_plot_dt$count))

tads_plot_dt$chromo <- as.character(tads_plot_dt$chromo)
tads_plot_dt$start <- as.numeric(as.character(tads_plot_dt$start))
stopifnot(!is.na(tads_plot_dt$start))
tads_plot_dt$end <- as.numeric(as.character(tads_plot_dt$end))
stopifnot(!is.na(tads_plot_dt$end))

stopifnot(length(unique(genes_plot_dt$chromo)) == 1)
stopifnot(length(unique(tads_plot_dt$chromo)) == 1)

chromo <- unique(genes_plot_dt$chromo)
  
stopifnot(!duplicated(file.path(tads_plot_dt$dataset, tads_plot_dt$cond1, tads_plot_dt$cond2)))

nDScons <- length(file.path(tads_plot_dt$dataset, tads_plot_dt$cond1, tads_plot_dt$cond2))

stopifnot(max(genes_plot_dt$count) <= nDScons)

if(!suppressPackageStartupMessages(require("ggplot2"))) stop("-- ggplot2 package required\n")  
if(!suppressPackageStartupMessages(require("ggsci"))) stop("-- ggsci package required\n")  
if(!suppressPackageStartupMessages(require("ggrepel"))) stop("-- ggrepel package required\n")  


tad_col <- pal_d3()(3)[1]
gene_col <- pal_d3()(3)[2]
syn_col <- pal_d3()(3)[3]
colConsThresh <- ceiling(nDScons/2)
geneColSetting <- setNames( c(consAll_col, consAbove_col, consBelow_col), 
                            c(paste0("= ", nDScons, "/", nDScons), 
                              paste0(">= ", colConsThresh, "/", nDScons), 
                              paste0("< ", colConsThresh, "/", nDScons)
                            )
)

# consAll_col <- "red"
consAll_col <- "darkgreen"
consAbove_col <- "orange"
consBelow_col <- "black"
tad_gene_space <- 0.5
gene_line_offset <- 0.1
TADlinecol <- "darkblue"
TADlt <- 1
TADlw <- 2
geneLt <- 1
geneLw <- 1
geneDelLt <- 2
geneDelLw <- 0.5

  
  subTit <- paste0("conserved in ", nDScons, " datasets")
  
  all_cols=NULL

  

    # ### CHANGE HERE TO HAVE RANK BY CMP TYPE and space between them
  tads_plot_dt$tad_size <- tads_plot_dt$end- tads_plot_dt$start

  stopifnot(!is.na(tads_plot_dt$cmpType))
  tads_plot_dt$cmpType <- factor(tads_plot_dt$cmpType, levels=sort(unique(as.character(tads_plot_dt$cmpType))))
  stopifnot(!is.na(tads_plot_dt$cmpType))
  tads_plot_dt$cmpType_num <- as.numeric(tads_plot_dt$cmpType)-1
  # tads_plot_dt <- tads_plot_dt[order(-tads_plot_dt$cmpType_num, tads_plot_dt$hicds_lab, tads_plot_dt$exprds_lab, decreasing = TRUE),]
  tads_plot_dt <- tads_plot_dt[order(-tads_plot_dt$cmpType_num, tads_plot_dt$tad_size, decreasing=TRUE),]
  # tads_plot_dt$ds_rank <- 1:nrow(tads_plot_dt) ## CHANGED HERE
  tads_plot_dt$ds_rank <- 1:nrow(tads_plot_dt) + tads_plot_dt$cmpType_num
  
  # for label colors
  if(is.null(all_cols)) {
    tads_plot_dt$ds_col <- "black"
  } else {
    tads_plot_dt$ds_col <- all_cols[as.character(tads_plot_dt$cmpType)]
  }

  stopifnot(!is.na(tads_plot_dt$ds_col))
  
  genes_plot_dt <- genes_plot_dt[order(genes_plot_dt$start, genes_plot_dt$end),]
  genes_plot_dt$gene_rank <- max(tads_plot_dt$ds_rank) + tad_gene_space + 1:nrow(genes_plot_dt)
  genes_plot_dt$gene_pos <- 0.5*(genes_plot_dt$start+genes_plot_dt$end)
  
  genes_plot_dt$col <- ifelse(genes_plot_dt$count == nDScons, consAll_col, ifelse(genes_plot_dt$count >= colConsThresh, consAbove_col, consBelow_col))
  # col = ifelse(nSymb == nDScons, paste0("# cons. = ", nDScons), 
  #              ifelse(nSymb >= colConsThresh, paste0("# cons. >= ", colConsThresh),   paste0("# cons. < ", colConsThresh)))
  )
  
  genes_plot_dt$col <- factor(genes_plot_dt$col, levels = as.character(geneColSetting))
  stopifnot(!is.na(genes_plot_dt$col))
  
  xscale <- seq(from=min(tads_plot_dt$start) , to=max(tads_plot_dt$end) , length.out=10)
  
  region_p <- ggplot() + 
    
    ggtitle(myTit, subtitle=subTit)+
    
    labs(x="", y="")+
    # lines for the TADs
    geom_segment( aes(x = tads_plot_dt$start, y = tads_plot_dt$ds_rank, 
                      xend=tads_plot_dt$end, yend = tads_plot_dt$ds_rank), 
                  colour = TADlinecol, linetype = TADlt, size = TADlw,
                  inherit.aes = F) + 
    # lines for the genes
    geom_segment( aes(x = genes_plot_dt$start, y = genes_plot_dt$gene_rank,
                      xend=genes_plot_dt$end, yend = genes_plot_dt$gene_rank,
                      
                      # colour = geneColSetting[paste0(genes_plot_dt$col)]
                      colour = genes_plot_dt$col
                      
                      
                      ),
                  # colour = genes_plot_dt$col, 
                  linetype = geneLt, size = geneLw,
                  
                  show.legend=TRUE,
                  
                  inherit.aes = F) +
    
    scale_color_manual(values= setNames(as.character(geneColSetting),as.character(geneColSetting)), 
                       labels = setNames(names(geneColSetting), as.character(geneColSetting))) +
    labs(color="# conserved")+
  
    
    # vertical lines for the gene delimiters
    geom_segment( aes(x = c(genes_plot_dt$start, genes_plot_dt$end), 
                      xend= c(genes_plot_dt$start, genes_plot_dt$end),
                      # y = rep(min(tads_plot_dt$ds_rank)-gene_line_offset, 2),
                      y = rep(min(tads_plot_dt$ds_rank)-gene_line_offset, 1),
                      # yend=rep(genes_plot_dt$gene_rank, 2)), 
                  yend=rep(genes_plot_dt$gene_rank, 2)), 
                  colour = rep(as.character(genes_plot_dt$col),2), linetype = geneDelLt, size = geneDelLw,
                  inherit.aes = F) +
    
    
    theme_void()
  
  save(genes_plot_dt, file="genes_plot_dt.Rdata", version=2)
  
  region_p2 <-  region_p + 
  geom_text_repel(
    aes(x = genes_plot_dt$gene_pos, y =  genes_plot_dt$gene_rank, label = genes_plot_dt$symbol,
        color=genes_plot_dt$col), inherit.aes = F, show.legend=F,
    font_face="italic",
    # nudge_x      = 0.05,
    # direction    = "y",
    nudge_y      = 2,
    direction    = "y",
    hjust        = 0.5,
    segment.size = 1
  )
  

  tads_plot_dt$raw_labels <- paste0(as.character(tads_plot_dt$dataset), " ", as.character(tads_plot_dt$cond1),  " - ",  as.character(tads_plot_dt$cond2))
  tads_plot_dt$ds_lab <- sapply(1:nrow(tads_plot_dt), function(i) {
    label_part1 <- tads_plot_dt$dataset[i]
    label_part2 <- tads_plot_dt$cond1[i]
    label_part3 <- tads_plot_dt$cond2[i]
    if(tads_plot_dt$upCond[i] == tads_plot_dt$cond1[i]){
      mylab <- gsub(" ","~", paste0(label_part1, "~", label_part2, "~bold(", label_part3, ")"))
    }else {
      mylab <- gsub(" ","~", paste0(label_part1, "~bold(", label_part2, ")~", label_part3))
    }
    
    mylab <- gsub("_", "[{\"-\"}]*", mylab)  # underscore are not recognize -> replace with dash lowerscript
    
    mylab <- gsub("^(\\d+)([a-zA-Z])", "\\1*\\2", mylab)
    
    mylab
  })
  
  
  axisOffset <- min(tads_plot_dt$end-tads_plot_dt$start)#90000
  

  
  

  region_p3 <- region_p2 + xlim(min(c(tads_plot_dt$start, genes_plot_dt$start)- axisOffset), NA) +
          geom_text_repel(
            aes(x = tads_plot_dt$start, y =  tads_plot_dt$ds_rank, label = tads_plot_dt$ds_lab), inherit.aes = F,
            nudge_x       = 3.5 - tads_plot_dt$start,
            direction     = "y",
            color = tads_plot_dt$ds_col, # ADDED 
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
    expand_limits(x = c(NA, max(xscale)+10000))+
  geom_segment(aes(x=min(xscale), xend=max(xscale),y=0, yend=0), lineend="round", colour="darkgrey", size=2)+  

        annotate("text", x=c(min(xscale), 0.5*(min(xscale)+max(xscale)), max(xscale)), y = -0.8, vjust=1, 
             size=5,
             label=c(min(xscale), chromo, max(xscale)) , colour="darkgrey", fontface="italic") 
    
  
  