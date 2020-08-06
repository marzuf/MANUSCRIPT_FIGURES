
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#' Plot TADs of a conserved region
#'
#' Function to plot the TADs and genes of a conserved region.
#'
#' @param genes_dt The text to be printed. Should contain at least the following columns: symbol/chromo/start/end/<count>
#' @param tads_dt Dataframe with information regarding conserved TADs to plot. Should contain at least the following columns: dataset/dsCat/cond1/cond2/chromo/start/end/<upCond>. Datasets will be grouped according to dsCat.
#' @param dsCat_cols Should be a named vector of colors. The names should match the levels of tads_dt$dsCat.
#' @param ... Other parameters for fine-tuning the plot (might be to look a bit in the code to adapt to your data).
#' @return The plot (ggplot).
#' @export
#' 
plot_conservedRegions <- function(genes_dt, tads_dt,
                                  TADlinecol1="darkblue",  TADlinecol2="darkturquoise",
                                  consAll_col="darkgreen",
                                  consAbove_col="orange",
                                  consBelow_col="black",
                                  tad_gene_space= 0.5,
                                  gene_line_offset= 0.1,
                                  TADlt=1,
                                  TADlw=2,
                                  geneLt=1,
                                  geneLw=1,
                                  geneDelLt=2, # vertical lines
                                  geneDelLw=0.5, # vertical lines
#                                  fontFamily="Hershey",
                                  dsCat_cols=NULL,
                                  colConsThresh=NULL,
                                  subTit=NULL,
                                  myTit=NULL
                                  ){


  tad_bicols <- c(TADlinecol1, TADlinecol2)

  
  if(!suppressPackageStartupMessages(require("ggplot2"))) stop("-- ggplot2 package required\n")  
  if(!suppressPackageStartupMessages(require("ggrepel"))) stop("-- ggrepel package required\n")  
  
  if(!is.null(dsCat_cols)) {
    stopifnot(as.character(tads_dt$dsCat) %in% names(dsCat_cols))
  }
  
  # stopifnot(c("symbol", "count", "chromo", "start", "end") %in% colnames(genes_dt))
  stopifnot(c("symbol", "chromo", "start", "end") %in% colnames(genes_dt))
  # stopifnot(c("dataset","dsCat", "upCond", "cond1", "cond2", "chromo", "start", "end") %in% colnames(tads_dt))
  stopifnot(c("dataset","dsCat", "cond1", "cond2", "chromo", "start", "end") %in% colnames(tads_dt))
  
  genes_dt$start <- as.numeric(as.character(genes_dt$start))
  stopifnot(!is.na(genes_dt$start))
  genes_dt$end <- as.numeric(as.character(genes_dt$end))
  stopifnot(!is.na(genes_dt$end))
  genes_dt$chromo <- as.character(genes_dt$chromo)
  if("count" %in% colnames(genes_dt)) {
    genes_dt$count <- as.numeric(as.character(genes_dt$count))
    stopifnot(!is.na(genes_dt$count))
  }
  tads_dt$chromo <- as.character(tads_dt$chromo)
  tads_dt$start <- as.numeric(as.character(tads_dt$start))
  stopifnot(!is.na(tads_dt$start))
  tads_dt$end <- as.numeric(as.character(tads_dt$end))
  stopifnot(!is.na(tads_dt$end))
  
  chromo <- unique(genes_dt$chromo)
  nDScons <- length(file.path(tads_dt$dataset, tads_dt$cond1, tads_dt$cond2))
  
  stopifnot(!duplicated(file.path(tads_dt$dataset, tads_dt$cond1, tads_dt$cond2)))
  stopifnot(length(unique(genes_dt$chromo)) == 1)
  stopifnot(length(unique(tads_dt$chromo)) == 1)
  if("count" %in% colnames(genes_dt)) stopifnot(max(genes_dt$count) <= nDScons)
  
  
  if(is.null(colConsThresh)) {
    colConsThresh <- ceiling(nDScons/2)  
  } else {
    stopifnot(!is.numeric(colConsThresh))
  }
  
  geneColSetting <- setNames( c(consAll_col, consAbove_col, consBelow_col), 
                              c(paste0("= ", nDScons, "/", nDScons), 
                                paste0(">= ", colConsThresh, "/", nDScons), 
                                paste0("< ", colConsThresh, "/", nDScons)
                              ))
  
  if(is.null(subTit)) subTit <- paste0("conserved in ", nDScons, " datasets")
  if(is.null(myTit)) myTit <- paste0("Conservation across datasets")
  

  tads_dt$tad_size <- tads_dt$end- tads_dt$start
  stopifnot(!is.na(tads_dt$dsCat))
  tads_dt$dsCat <- factor(tads_dt$dsCat, levels=sort(unique(as.character(tads_dt$dsCat))))
  stopifnot(!is.na(tads_dt$dsCat))
  tads_dt$cmpType_num <- as.numeric(tads_dt$dsCat)-1
  tads_dt <- tads_dt[order(-tads_dt$cmpType_num, tads_dt$tad_size, decreasing=TRUE),]
  tads_dt$ds_rank <- 1:nrow(tads_dt) + tads_dt$cmpType_num
  
  # for label colors
  if(is.null(dsCat_cols)) {
    tads_dt$ds_col <- "black"
  } else {
    tads_dt$ds_col <- dsCat_cols[as.character(tads_dt$dsCat)]
  }
  stopifnot(!is.na(tads_dt$ds_col))
  
  genes_dt <- genes_dt[order(genes_dt$start, genes_dt$end),]
  genes_dt$gene_rank <- max(tads_dt$ds_rank) + tad_gene_space + 1:nrow(genes_dt)
  genes_dt$gene_pos <- 0.5*(genes_dt$start+genes_dt$end)
  if("count" %in% colnames(genes_dt)) {
    genes_dt$col <- ifelse(genes_dt$count == nDScons, consAll_col, ifelse(genes_dt$count >= colConsThresh, consAbove_col, consBelow_col))
    genes_dt$col <- factor(genes_dt$col, levels = as.character(geneColSetting))
  } else {
    genes_dt$col <- consAll_col
  }
  stopifnot(!is.na(genes_dt$col))
  
  
  tads_dt$raw_labels <- paste0(as.character(tads_dt$dataset), " ", as.character(tads_dt$cond1),  " - ",  as.character(tads_dt$cond2))
  tads_dt$ds_lab <- sapply(1:nrow(tads_dt), function(i) {
    label_part1 <- tads_dt$dataset[i]
    label_part2 <- tads_dt$cond1[i]
    label_part3 <- tads_dt$cond2[i]
    if("upCond" %in% colnames(tads_dt)) {
      if(tads_dt$upCond[i] == tads_dt$cond1[i]){
        mylab <- gsub(" ","~", paste0(label_part1, "~", label_part2, "~bold(", label_part3, ")"))
      }else {
        mylab <- gsub(" ","~", paste0(label_part1, "~bold(", label_part2, ")~", label_part3))
      }  
    } else {
      mylab <- gsub(" ","~", paste0(label_part1, "~", label_part2, "~", label_part3))
    }
    mylab <- gsub("_", "[{\"-\"}]*", mylab)  # underscore are not recognize -> replace with dash lowerscript
    mylab <- gsub("^(\\d+)([a-zA-Z])", "\\1*\\2", mylab)
    mylab
  })


####>>>>> UPDATE HERE IN CASE THERE ARE SEVERAL TADs BY DATASET
  all_tad_cols <- unlist(sapply(as.numeric(table(tads_dt$ds_rank)), function(x) rep(tad_bicols, ceiling(x/2) )[1:x]))
  stopifnot(length(all_tad_cols) == nrow(tads_dt))
  tads_dt$tadLineCol <- all_tad_cols
####>>>>>


  
  axisOffset <- min(tads_dt$end-tads_dt$start)#90000
  xscale <- seq(from=min(tads_dt$start) , to=max(tads_dt$end) , length.out=10)
  
  region_p <- ggplot() + 
    ggtitle(myTit, subtitle=subTit)+
    labs(x="", y="")+
    # lines for the TADs
    geom_segment( aes(x = tads_dt$start, y = tads_dt$ds_rank, 
                      xend=tads_dt$end, yend = tads_dt$ds_rank), 
####>>>>>   UPDATE HERE IN CASE THERE ARE SEVERAL TADs BY DATASET               colour = TADlinecol, 
                  colour = tads_dt$tadLineCol, 
				linetype = TADlt, size = TADlw,
                  inherit.aes = F) + 
    # lines for the genes
    geom_segment( aes(x = genes_dt$start, y = genes_dt$gene_rank,
                      xend=genes_dt$end, yend = genes_dt$gene_rank,
                      colour = genes_dt$col),
                  linetype = geneLt, size = geneLw, show.legend=TRUE,inherit.aes = F) 
    
    
    if("count" %in% colnames(genes_dt))    {     region_p <- region_p + 
    
    scale_color_manual(values= setNames(as.character(geneColSetting),as.character(geneColSetting)), 
                       labels = setNames(names(geneColSetting), as.character(geneColSetting))) +
      labs(color="# conserved")
      
    } else {
      region_p <- region_p + guides(color=FALSE)
    }


####>>>>> UPDATE HERE IN CASE THERE ARE SEVERAL TADs BY DATASET
  ds_tads_dt <- tads_dt[,c("start", "ds_rank", "ds_lab", "ds_col")]
  ds_tads_dt <- aggregate(start~., data=ds_tads_dt, FUN=min)
  stopifnot(!duplicated(ds_tads_dt$ds_lab))




    
    region_p <- region_p + 
    
    # vertical lines for the gene delimiters
    geom_segment( aes(x = c(genes_dt$start, genes_dt$end), 
                      xend= c(genes_dt$start, genes_dt$end),
                      # y = rep(min(tads_dt$ds_rank)-gene_line_offset, 2),
                      #y = rep(min(tads_dt$ds_rank)-gene_line_offset, 1), 
                      y = rep(min(ds_tads_dt$ds_rank)-gene_line_offset, 1), ####>>>>> UPDATE HERE IN CASE THERE ARE SEVERAL TADs BY DATASET
                      # yend=rep(genes_dt$gene_rank, 2)), 
                      yend=rep(genes_dt$gene_rank, 2)), 
                  colour = rep(as.character(genes_dt$col),2), linetype = geneDelLt, size = geneDelLw,
                  inherit.aes = F) +
    theme_void() + 
    geom_text_repel(
      aes(x = genes_dt$gene_pos, y =  genes_dt$gene_rank, label = genes_dt$symbol,
          color=genes_dt$col), inherit.aes = F, show.legend=F,
      #font_face="italic",
      # nudge_x      = 0.05,
      # direction    = "y",
      nudge_y      = 2,
      direction    = "y",
      hjust        = 0.5,
      segment.size = 1
    )+ 
#    xlim(min(c(tads_dt$start, genes_dt$start)- axisOffset), NA) +
#    geom_text_repel(
#      aes(x = tads_dt$start, y =  tads_dt$ds_rank, label = tads_dt$ds_lab), inherit.aes = FALSE,
#      nudge_x       = 3.5 - tads_dt$start,
#      direction     = "y",
#      color = tads_dt$ds_col, # ADDED 
#      hjust         = 1, parse = T, size=3#,
#      # force_pull   = 0
#    ) +
 ####>>>>> UPDATE HERE IN CASE THERE ARE SEVERAL TADs BY DATASET
    xlim(min(c(ds_tads_dt$start, genes_dt$start)- axisOffset), NA) +
    geom_text_repel(
      aes(x = ds_tads_dt$start, y =  ds_tads_dt$ds_rank, label = ds_tads_dt$ds_lab), inherit.aes = FALSE,
      nudge_x       = 3.5 - ds_tads_dt$start,
      direction     = "y",
      color = ds_tads_dt$ds_col, # ADDED 
      hjust         = 1, parse = T, size=3#,
      # force_pull   = 0
    ) +
    theme(plot.margin =margin(t = 25, r = 25, b = 10, l = 10, unit = "pt"),
#          text = element_text(family=fontFamily),
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
  region_p
  return(region_p)
}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#' Volcano plot TAD meanLogFC and meanCorr 
#'
#' Function that returns a volcano-like plot of TAD meanLogFC and intra-TAD meanCorr.
#'
#' @param meanCorr Vector of intra-TAD mean correlations (vector names should be TAD IDs).
#' @param meanFC Vector of TAD meanLogFC (vector names should be TAD IDs).
#' @param comb_pval Vector of TAD p-values (vector names should be TAD IDs).
#' @param tads_to_annot (optional) TAD that should be annotated with a label.
#' @param padjusted If the p-values are already adjusted (BH method).
#' @param ... Other parameters for fine-tuning the plot.
#' @return The plot (ggplot).
#' @export


plot_volcanoTADsCorrFC <- function(meanCorr, meanFC, comb_pval, 
                                   padjusted=FALSE,
                                   tads_to_annot=NULL,
                                   x_lab="TAD mean LogFC",
                                   y_lab="TAD mean intraCorr",
                                   plotTit=NULL,
                                   subTit=NULL,
                                   fcUpDownThresh=0,
                                   signifThresh=0.05,
                                   verySignifThresh=0.01,
                                   strongUp="#EE0011FF",
                                   lightUp="#EE001132",
                                   strongDown="#0C5BB0FF",
                                   lightDown="#0C5BB032"
                                   ) {
  
  if(!suppressPackageStartupMessages(require("ggplot2"))) stop("-- ggplot2 package required\n")  
  if(!suppressPackageStartupMessages(require("ggrepel"))) stop("-- ggrepel package required\n")  


  if(!padjusted) {
    cat("... adjust p-values with BH method\n")
    adj_comb_pval <- p.adjust(comb_pval, method="BH")
  } else {
    adj_comb_pval <- comb_pval
  }
  
  interReg <- Reduce(intersect, list(names(meanCorr) , names(meanFC), names(adj_comb_pval)))
  keepCorr <- meanCorr[interReg]
  keepFC <- meanFC[interReg]
  keepComb <- adj_comb_pval[interReg]
  
  cat(paste0("... kept meanCorr values:\t", length(keepCorr), "/", length(meanCorr), "\n"))
  cat(paste0("... kept meanFC values:\t", length(keepFC), "/", length(meanFC), "\n"))
  cat(paste0("... kept comb. p-val values:\t", length(keepComb), "/", length(adj_comb_pval), "\n"))
  
  
  ex_DT <- data.frame(
    region = interReg,
    meanLogFC = keepFC[interReg],
    meanCorr = keepCorr[interReg],
    adjPvalComb = keepComb[interReg],
    stringsAsFactors = FALSE
  )
  stopifnot(!is.na(ex_DT))
  ex_DT$adjPvalComb_log10 <- -log10(ex_DT$adjPvalComb)
  
  
  nSignif <- sum(ex_DT$adjPvalComb<=signifThresh)
  nVerySignif <- sum(ex_DT$adjPvalComb<=verySignifThresh)
  
  if(is.null(plotTit)) plotTit <- paste0("TAD meanLogFC and meanCorr")
  if(is.null(subTit)) subTit <- paste0("# TADs=", nrow(ex_DT), "; # signif.<=", signifThresh, "=",nSignif,"; # signif.<=", verySignifThresh, "=",nVerySignif)
  
  ex_DT$adjPvalComb_log10_col <- ex_DT$adjPvalComb_log10
  ex_DT$adjPvalComb_log10_col[ex_DT$adjPvalComb_log10_col <= -log10(signifThresh)] <- -log10(signifThresh)
  
  ex_DT$adjPvalComb_crop <- ex_DT$adjPvalComb 
  ex_DT$adjPvalComb_crop[ex_DT$adjPvalComb_crop >= signifThresh] <- signifThresh
  
  ex_DT$adjPvalCombSignif <- ex_DT$adjPvalComb <= signifThresh
  
  
  ex_DT$dotcol_labs1 <- ifelse(ex_DT$meanLogFC > fcUpDownThresh, paste0("logFC>", fcUpDownThresh), paste0("logFC<=", fcUpDownThresh))
  ex_DT$dotcol_labs2 <- ifelse(ex_DT$adjPvalComb <= verySignifThresh, paste0("adjP<=", verySignifThresh),
                               ifelse(ex_DT$adjPvalComb <= signifThresh, paste0("adjP<=", signifThresh), NA))
  
  ex_DT$dotcol_lab <- interaction(ex_DT$dotcol_labs1, ex_DT$dotcol_labs2)
  ex_DT$dotcol_lab <- as.character(ex_DT$dotcol_lab)
  stopifnot(which(is.na(ex_DT$dotcol_labs2)) == which(ex_DT$adjPvalComb > signifThresh))
  ex_DT$dotcol_lab[is.na(ex_DT$dotcol_labs2)] <- paste0("adjP>", signifThresh)
  stopifnot(!is.na(ex_DT$dotcol_lab))
  
  mycols <- setNames(c(strongUp, lightUp, strongDown, lightDown, "darkgrey"), 
                     c(paste0("logFC>", fcUpDownThresh, ".adjP<=", verySignifThresh),
                       paste0("logFC>", fcUpDownThresh, ".adjP<=", signifThresh),
                       paste0("logFC<=", fcUpDownThresh, ".adjP<=", verySignifThresh),
                       paste0("logFC<=",fcUpDownThresh, ".adjP<=", signifThresh), 
                       paste0("adjP>", signifThresh)))
  
  mycols_labs <- setNames(c(paste0("mean LogFC > ", fcUpDownThresh, " &\nadj. p-val <=", verySignifThresh),
                            paste0("mean LogFC > ", fcUpDownThresh, " &\nadj. p-val <=", signifThresh),
                            paste0("mean LogFC <= ", fcUpDownThresh, " &\nadj. p-val <=", verySignifThresh),
                            paste0("mean LogFC <= ", fcUpDownThresh, " &\nadj. p-val <=", signifThresh), 
                            paste0("adj. p-val > ", signifThresh)),
                          c(
                            paste0("logFC>",fcUpDownThresh, ".adjP<=", verySignifThresh),
                            paste0("logFC>", fcUpDownThresh, ".adjP<=", signifThresh),
                            paste0("logFC<=", fcUpDownThresh, ".adjP<=", verySignifThresh),
                            paste0("logFC<=", fcUpDownThresh, ".adjP<=", signifThresh), 
                            paste0("adjP>", signifThresh)))
  
  fc_corr_p <- ggplot(ex_DT, aes(x=meanLogFC, y = meanCorr, fill = dotcol_lab, color=dotcol_lab )) +
    ggtitle(plotTit, subtitle = subTit)+
    scale_fill_manual(values=mycols, labels = mycols_labs)+  
    scale_color_manual(values=mycols,
                       labels=setNames(names(mycols), mycols), guide=F)+
    guides(color=F)+
    
    geom_point(data=ex_DT[ex_DT$dotcol_lab==paste0("adjP>", signifThresh),],
               shape=21,
               aes(size = adjPvalComb_log10_col, color = dotcol_lab), stroke=1) +
    geom_point(data=ex_DT[ex_DT$dotcol_lab!=paste0("adjP>", signifThresh),],
               shape=21,
               aes(size = adjPvalComb_log10_col, color = dotcol_lab), stroke=1)+
    
    labs( size= "TAD adj. p-val", fill="")+
    geom_label_repel(data = ex_DT[ex_DT$region %in% tads_to_annot,],
                     aes(x= meanLogFC, 
                         y=meanCorr, label=region),
                     inherit.aes = F,
                     min.segment.length = unit(0, 'lines'),
                     force = 10,
                     segment.size = 0.5,
                     nudge_y = 0.1)+
    scale_x_continuous(name=paste0(x_lab),breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name=paste0(y_lab),breaks = scales::pretty_breaks(n = 10))+
    scale_size(range=c(0.5,8),#expand=c(2,0),
               breaks=c(-log10(0.05), -log10(0.01), -log10(0.005), -log10(0.001)),
               labels=c(paste0(">= 0.05"),"0.01"," 0.005","0.001"),
               guide="legend") +
    guides(size = guide_legend(order=1))+
    guides(fill = guide_legend(override.aes = list(color="white", size=4),order=2))+
    theme(
#      text = element_text(family=fontFamily),
      panel.background = element_rect(fill = "transparent"),
      panel.grid.major.x =  element_blank(),
      panel.grid.minor.x =  element_blank(),
      axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
      axis.title.y = element_text(size=14, hjust=0.5, vjust=0.5),
      plot.title = element_text(hjust=0.5, size = 16, face="bold"),
      plot.subtitle = element_text(hjust=0.5, size = 14, face="italic"),
      axis.line = element_line(),
      axis.text.y = element_text(size=12, hjust=0.5, vjust=0.5, color="black"),
      axis.text.x = element_text(size=12, hjust=0.5, vjust=0.5, color="black"),
      panel.grid = element_blank(),
      legend.key = element_rect(fill = NA),
      legend.title =  element_text(size=12,face="bold"),
      panel.grid.major.y =  element_blank(),
      panel.grid.minor.y =  element_blank(),
      legend.text = element_text(margin = margin(b = 0.1, unit = 'in'), size=10)
    ) +
    geom_hline(yintercept = 0, linetype=2, color="darkgrey")+
    geom_vline(xintercept = 0, linetype=2, color="darkgrey")
  fc_corr_p
  return(fc_corr_p)
}



