# Rscript cmp_nSignif_sampleSize.R GSE105381_HepG2_40kb TCGAlihc_wt_mutCTNNB1
# Rscript cmp_nSignif_sampleSize.R ENCSR312KHQ_SK-MEL-5_40kb TCGAskcm_lowInf_highInf
# Rscript cmp_nSignif_sampleSize.R ENCSR079VIJ_G401_40kb TCGAkich_norm_kich


setDir <- "/media/electron"
setDir <- ""

require(foreach)
require(doMC)
registerDoMC(40)
require(ggplot2)
require(ggpubr)
require(reshape2)

plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 9


tadSignifThresh <- 0.01
geneSignifThresh <- 0.01

init_hicds <- "GSE105381_HepG2_40kb"
init_exprds <- "TCGAlihc_wt_mutCTNNB1"

source("../settings.R")

myWidthGG <- myWidthGG*1.2
myWidthGG <- 7
myHeightGG <- 5

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
init_hicds <- args[1]
init_exprds <- args[2]

stopifnot(init_hicds %in% names(hicds_names))
stopifnot(init_exprds %in% names(exprds_names))

hicds_lab <- hicds_names[paste0(init_hicds)]
exprds_lab <- exprds_names[paste0(init_exprds)]

outFolder <- file.path("CMP_NSIGNIF_SAMPLESIZE", init_hicds, init_exprds)
dir.create(outFolder, recursive = TRUE)

stopifnot(init_hicds %in% all_obs_hicds)
randomsub_hicds <- all_hicds[grepl("RANDOMSUB", all_hicds) & grepl(gsub(paste0("_40kb"), "", init_hicds), all_hicds)]

my_hicds <- c(init_hicds, randomsub_hicds)
my_exprds <- get_exprds(my_hicds)

script1_name <- "1_runGeneDE"
script11_name <- "11sameNbr_runEmpPvalCombined"

nsub=""
nsub=20


custom_plot <- function(p) {
  newp <- p+
    scale_fill_manual(values = c(col1, col2)) +
    scale_color_manual(values = c(col1, col2)) +
    theme(
      text = element_text(family=fontFamily),
      plot.title = element_text(hjust = 0.5, face = "bold", size=16),
      plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 14),
      panel.grid = element_blank(),
      panel.grid.major.y = element_line(colour = "grey"),
      panel.grid.minor.y = element_line(colour = "grey"),
      axis.line.x= element_line(size = .2, color = "black"),
      axis.line.y = element_line(size = .2, color = "black"),
      axis.text.y = element_text(color="black", hjust=1,vjust = 0.5, size=12),
      axis.text.x =element_text(color="black", hjust=0.5,vjust = 0.5, size=12),
      # axis.ticks.x = element_blank(),
      axis.title.y = element_text(color="black", size=14),
      axis.title.x = element_text(color="black", size=14),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      legend.key.size = unit(1.0, 'cm'),
      legend.background =  element_rect(),
      legend.text = element_text(size=12),
      legend.key = element_blank(),
      legend.title = element_text(face="bold", size=12, hjust=0.5)
    )
  return(newp)
}


all_tadSignif_dt <- foreach(hicds = my_hicds, .combine='rbind') %dopar% {
      exprds = all_exprds[[paste0(hicds)]][1]
      hicds_dt <- foreach(exprds = my_exprds[[paste0(hicds)]], .combine='rbind') %do% {
        # 
        # 
        # 
        # if(nsub == "") {
        #   hicds <- paste0(init_hicds,"_", nsub, "40kb")
        # } else{
        #   hicds <- paste0(init_hicds,"_RANDOMSUB", nsub, "_40kb")
        # }
        
        if(init_exprds != exprds) return(NULL)
        
        settingF <- file.path(runFolder, "PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R"))
        stopifnot(file.exists(settingF))
        source(settingF)  
        s1 <- get(load(file.path(setDir, sample1_file)))
        s2 <- get(load(file.path(setDir, sample2_file)))
        # stopifnot(length(s1) == length(s2))
        # if(nsub != "") stopifnot(length(s1) == nsub)
        
        pval_file <- file.path(pipFolder,  hicds, exprds, script11_name, "emp_pval_combined.Rdata")
#        if(!file.exists(pval_file)) return(NULL)
		stopifnot(file.exists(pval_file))
        pval <- get(load(pval_file))
        adj_pval <- p.adjust(pval, method="BH")
        
        data.frame(
          hicds = hicds,
          exprds = exprds,
          nSamp1 = length(s1),
          nSamp2 = length(s2),
          region = names(adj_pval),
          adjCombPval = as.numeric(adj_pval),
          stringsAsFactors = FALSE
        )
        
      }
      hicds_dt
      
}
outFile <- file.path(outFolder, "all_tadSignif_dt.Rdata")
save(all_tadSignif_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
# all_tadSignif_dt <- get(load(file.path(outFolder, "all_tadSignif_dt.Rdata")))


all_tadSignif_dt$nSamp <- all_tadSignif_dt$nSamp1 +  all_tadSignif_dt$nSamp2
all_tadSignif_dt <- all_tadSignif_dt[order(all_tadSignif_dt$nSamp),]
all_tadSignif_dt$sampLab <- paste0(all_tadSignif_dt$nSamp1, "+", all_tadSignif_dt$nSamp2 )

labLevels <- unique(as.character(all_tadSignif_dt$sampLab))


tadSignif_agg_dt <- aggregate(adjCombPval~hicds+exprds+sampLab,data=all_tadSignif_dt,FUN=function(x) mean(x <= tadSignifThresh))
colnames(tadSignif_agg_dt)[colnames(tadSignif_agg_dt) == "adjCombPval"] <- "ratioSignifTADs"

outFile <- file.path(outFolder, "tadSignif_agg_dt.Rdata")
save(tadSignif_agg_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))


all_geneSignif_dt <- foreach(hicds = my_hicds, .combine='rbind') %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  hicds_dt <- foreach(exprds = my_exprds[[paste0(hicds)]], .combine='rbind') %do% {
    #
    if(init_exprds != exprds) return(NULL)
    settingF <- file.path(runFolder, "PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R"))
    stopifnot(file.exists(settingF))
    source(settingF)  
    s1 <- get(load(file.path(setDir, sample1_file)))
    s2 <- get(load(file.path(setDir, sample2_file)))

    pval_file <- file.path(pipFolder,  hicds, exprds, script1_name, "DE_topTable.Rdata")
#    if(!file.exists(pval_file)) return(NULL)
    stopifnot(file.exists(pval_file))
    de_dt <- get(load(pval_file))
    
    data.frame(
      hicds = hicds,
      exprds = exprds,
      nSamp1 = length(s1),
      nSamp2 = length(s2),
      gene = as.character(de_dt$genes) ,
      adjPval = de_dt$adj.P.Val,
      stringsAsFactors = FALSE
    )
  }
  hicds_dt
}
outFile <- file.path(outFolder, "all_geneSignif_dt.Rdata")
save(all_geneSignif_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
# all_geneSignif_dt <- get(load(file.path(outFolder, "all_geneSignif_dt.Rdata")))

all_geneSignif_dt$nSamp <- all_geneSignif_dt$nSamp1 +  all_geneSignif_dt$nSamp2
all_geneSignif_dt$nSamp <- all_geneSignif_dt[order(all_geneSignif_dt$nSamp),]
all_geneSignif_dt$sampLab <- paste0(all_geneSignif_dt$nSamp1, "+", all_geneSignif_dt$nSamp2 )


geneSignif_agg_dt <- aggregate(adjPval~hicds+exprds+sampLab,data=all_geneSignif_dt,FUN=function(x) mean(x <= geneSignifThresh))
colnames(geneSignif_agg_dt)[colnames(geneSignif_agg_dt) == "adjPval"] <- "ratioSignifGenes"

outFile <- file.path(outFolder, "geneSignif_agg_dt.Rdata")
save(geneSignif_agg_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

allSignif_agg_dt <- merge(geneSignif_agg_dt, tadSignif_agg_dt, by=c("hicds", "exprds", "sampLab"), all=TRUE)
stopifnot(!is.na(allSignif_agg_dt))


plot_dt <- melt(allSignif_agg_dt, id.vars = c("hicds", "exprds", "sampLab"))

plot_dt$variable <- as.character(plot_dt$variable)

plot_dt$variable[plot_dt$variable == "ratioSignifTADs"] <- paste0("TADs\n(adj. p-val<=", tadSignifThresh, ")")
plot_dt$variable[plot_dt$variable == "ratioSignifGenes"] <- paste0("genes\n(adj. p-val<=", geneSignifThresh, ")")

plot_dt$sampLab <- factor(plot_dt$sampLab, levels=labLevels)

subTit <- "(variable sample size)"
subTit <- paste0(hicds_lab, " - ", exprds_lab, " - variable sample size")
subTit <- paste0(hicds_lab, " - ", exprds_lab)

col1 <- "darkblue"
col2 <- "darkorange"

outFile <- file.path(outFolder, "plot_dt.Rdata")
save(plot_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))

###############################################################
############################################################### # detected features
###############################################################

init_p <- custom_plot(ggplot(data = plot_dt,aes_string(x="sampLab",y="value", fill = "variable", color="variable"))+
                        ggtitle("Ratio signif. and sample size", subtitle = subTit)+
  scale_x_discrete(name= "# of samples")+
  scale_y_continuous(name="Ratio signif. features", breaks=scales::pretty_breaks(n = 10))+
  labs(fill  = paste0("Ratio signif.") , color  = paste0("Ratio signif.") ) )

signif_p <- init_p +  
  geom_boxplot(notch = FALSE, outlier.shape=NA, fill=NA)+
  geom_point(aes(color=variable), position=position_jitterdodge(),  alpha=0.5)    +
  # guides(color = FALSE)
  guides(fill = FALSE)



outFile <- file.path(outFolder, paste0(init_hicds, "_", init_exprds, "_signifSubSamples_boxplot.", plotType))
ggsave(plot = signif_p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


signif_p_strip <- init_p + 
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 1.2
  ) +
  stat_summary(
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.4,
    position = position_dodge(0.8)
  ) +
  guides(fill = FALSE)



outFile <- file.path(outFolder, paste0(init_hicds, "_", init_exprds, "_signifSubSamples_stripchart.", plotType))
ggsave(plot = signif_p_strip, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


###############################################################
############################################################### # consistency
###############################################################


all_geneSignif <- foreach(hicds = my_hicds) %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  hicds_data <- foreach(exprds = my_exprds[[paste0(hicds)]]) %do% {
    
    if(init_exprds != exprds) return(NULL)
    settingF <- file.path(runFolder, "PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R"))
    stopifnot(file.exists(settingF))
    source(settingF)  
    s1 <- get(load(file.path(setDir, sample1_file)))
    s2 <- get(load(file.path(setDir, sample2_file)))
    
    pval_file <- file.path(pipFolder,  hicds, exprds, script1_name, "DE_topTable.Rdata")
    if(!file.exists(pval_file)) return(NULL)
    de_dt <- get(load(pval_file))
    as.character(de_dt$genes)[de_dt$adj.P.Val <= geneSignifThresh]
  }
  names(hicds_data) <- my_exprds[[paste0(hicds)]]
  hicds_data_filt <- Filter(f=function(x) !is.null(x), hicds_data)
  hicds_data_filt
}
names(all_geneSignif) <- my_hicds
outFile <- file.path(outFolder, "all_geneSignif.Rdata")
save(all_geneSignif, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
# all_geneSignif <- get(load(file.path(outFolder, "all_geneSignif.Rdata")))

all_tadSignif <- foreach(hicds = my_hicds) %dopar% {
  exprds = all_exprds[[paste0(hicds)]][1]
  hicds_data <- foreach(exprds = my_exprds[[paste0(hicds)]]) %do% {

    if(init_exprds != exprds) return(NULL)
    settingF <- file.path(runFolder, "PIPELINE", "INPUT_FILES", hicds, paste0("run_settings_", exprds, ".R"))
    stopifnot(file.exists(settingF))
    source(settingF)  
    s1 <- get(load(file.path(setDir, sample1_file)))
    s2 <- get(load(file.path(setDir, sample2_file)))
    # stopifnot(length(s1) == length(s2))
    # if(nsub != "") stopifnot(length(s1) == nsub)
    
    pval_file <- file.path(pipFolder,  hicds, exprds, script11_name, "emp_pval_combined.Rdata")
    if(!file.exists(pval_file)) return(NULL)
    pval <- get(load(pval_file))
    adj_pval <- p.adjust(pval, method="BH")
    
    names(adj_pval)[adj_pval <= tadSignifThresh]
  }
  names(hicds_data) <- my_exprds[[paste0(hicds)]]
  hicds_data_filt <- Filter(f=function(x) !is.null(x), hicds_data)
  hicds_data_filt
}
names(all_tadSignif) <- my_hicds
outFile <- file.path(outFolder, "all_tadSignif.Rdata")
save(all_tadSignif, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
# all_tadSignif <- get(load(file.path(outFolder, "all_tadSignif.Rdata")))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
#### PLOT1:  consist. to ref
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

ref_geneSignif <- all_geneSignif[[paste0(init_hicds)]][[paste0(init_exprds)]]
ref_tadSignif <- all_tadSignif[[paste0(init_hicds)]][[paste0(init_exprds)]]

refConst_geneSignif_dt <- foreach(i = 1:nrow(allSignif_agg_dt), .combine='rbind') %dopar% {
  
  i_hicds <- allSignif_agg_dt$hicds[i]
  i_exprds <- allSignif_agg_dt$exprds[i]
  i_geneSignif <- all_geneSignif[[paste0(i_hicds)]][[paste0(i_exprds)]]
  i_tadSignif <- all_tadSignif[[paste0(i_hicds)]][[paste0(i_exprds)]]
  
  geneRatioDetected <- mean(ref_geneSignif %in% i_geneSignif)
  tadRatioDetected <- mean(ref_tadSignif %in% i_tadSignif)
  
  data.frame(
    hicds=i_hicds,
    exprds=i_exprds,
    ratioGeneRefSignif = geneRatioDetected,
    ratioTadRefSignif = tadRatioDetected,
    stringsAsFactors = FALSE
  )
  
}
outFile <- file.path(outFolder, "refConst_geneSignif_dt.Rdata")
save(refConst_geneSignif_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
# refConst_geneSignif_dt <- get(load(file.path(outFolder, "refConst_geneSignif_dt.Rdata")))


foo_dt  <- plot_dt
foo_dt$dataset <- file.path(foo_dt$hicds, foo_dt$exprds)
getSampLab <- setNames(as.character(foo_dt$sampLab), as.character(foo_dt$dataset))

plot_dt2 <- melt(refConst_geneSignif_dt, id=c("hicds", "exprds") )
# plot_dt2$dataset <- file.path(plot_dt2$hicds, plot_dt2$exprds)
# stopifnot(table(plot_dt2$dataset) == 2)
# plot_dt2$dataset <- NULL
# plot_dt2$variable <- as.character(plot_dt2$variable)
# plot_dt2 <- merge(plot_dt[,c("hicds", "exprds", "sampLab")], plot_dt2, by=c("hicds", "exprds"), all.x=F, all.y=T)
# stopifnot(!is.na(plot_dt2))
# plot_dt2$dataset <- NULL
plot_dt2$dataset <- file.path(plot_dt2$hicds, plot_dt2$exprds)
stopifnot(table(plot_dt2$dataset) == 2)
plot_dt2$sampLab <- getSampLab[paste0(plot_dt2$dataset)]
plot_dt2$sampLab <- factor(as.character(plot_dt2$sampLab), levels=labLevels)
stopifnot(table(plot_dt2$dataset) == 2)

plot_dt2$variable <- as.character(plot_dt2$variable)
plot_dt2$variable[plot_dt2$variable == "ratioGeneRefSignif"] <- paste0("genes\n(adj. p-val<=", geneSignifThresh, ")")
plot_dt2$variable[plot_dt2$variable == "ratioTadRefSignif"] <- paste0("TADs\n(adj. p-val<=", tadSignifThresh, ")")

init_p2 <- custom_plot(ggplot(data = plot_dt2,aes_string(x="sampLab",y="value", fill = "variable", color="variable"))+
                         ggtitle("Ratio of ref. signif. and sample size", subtitle = subTit)+
  scale_x_discrete(name= "# of samples")+
  scale_y_continuous(name="Ratio of ref. signif. features", breaks=scales::pretty_breaks(n = 10))+
  labs(fill  = paste0("Ratio from\nref. signif.") , color  = paste0("Ratio from\nref. signif.") ))

ratioRef_p <- init_p2 +  
  geom_boxplot(notch = FALSE, outlier.shape=NA, fill=NA)+
  geom_point(aes(color=variable), position=position_jitterdodge(),  alpha=0.5)    +
  # guides(color = FALSE)
guides(fill = FALSE)



outFile <- file.path(outFolder, paste0(init_hicds, "_", init_exprds, "_signifSubSamplesRatioRef_boxplot.", plotType))
ggsave(plot = ratioRef_p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


ratioRef_p_strip <- init_p2 + 
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 1.2
  ) +
  stat_summary(
    fun.data="mean_sdl",  fun.args = list(mult=1), 
    geom = "pointrange",  size = 0.4,
    position = position_dodge(0.8)
  ) +
  guides(fill = FALSE)



outFile <- file.path(outFolder, paste0(init_hicds, "_", init_exprds, "_signifSubSamplesRatioRef_stripchart.", plotType))
ggsave(plot = ratioRef_p_strip, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### PLOT2: constit. among each level
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
tmp_dt <- allSignif_agg_dt
sampLevels <- unique(as.character(allSignif_agg_dt$sampLab))

stopifnot(allSignif_agg_dt$exprds == init_exprds)
nsamp = sampLevels[2]
acrossConst_geneSignif_dt <- foreach(nsamp = sampLevels, .combine='rbind') %dopar% {
  
  matching_hicds <- allSignif_agg_dt$hicds[as.character(allSignif_agg_dt$sampLab) == nsamp]

  # if(length(matching_hicds) == 1) {
  #   matching_geneSignif <- lapply(all_geneSignif[matching_hicds], function(x)x[[1]])
  # } else {
  #   
  # } 
  matching_geneSignif <- unlist(all_geneSignif[matching_hicds],recursive=FALSE)
  # matching_geneSignif <- lapply(all_geneSignif[matching_hicds], function(x)x[[1]])
  intersect_geneMatch <- Reduce(intersect, matching_geneSignif)
  union_geneMatch <- Reduce(union, matching_geneSignif)
  geneConsistRatio <- length(intersect_geneMatch)/length(union_geneMatch)
  
  matching_tadSignif <- unlist(all_tadSignif[matching_hicds], recursive=FALSE)
    # matching_tadSignif <- lapply(all_tadSignif[matching_hicds], function(x)x[[1]])
  intersect_tadMatch <- Reduce(intersect, matching_tadSignif)
  union_tadMatch <- Reduce(union, matching_tadSignif)
  tadConsistRatio <- length(intersect_tadMatch)/length(union_tadMatch)
  
  data.frame(
    hicds = init_hicds,
    exprds=init_exprds,
    geneConsistRatio=geneConsistRatio,
    tadConsistRatio=tadConsistRatio,
    sampLab = nsamp,
    stringsAsFactors = FALSE
  )
  
}
outFile <- file.path(outFolder, "acrossConst_geneSignif_dt.Rdata")
save(acrossConst_geneSignif_dt, file=outFile, version=2)
cat(paste0("... written: ", outFile, "\n"))
# acrossConst_geneSignif_dt <- get(load(file.path(outFolder, "acrossConst_geneSignif_dt.Rdata")))


plot_dt3 <- melt(acrossConst_geneSignif_dt, id=c("hicds", "exprds", "sampLab") )
stopifnot(!is.na(plot_dt3))

plot_dt3$variable <- as.character(plot_dt3$variable)
plot_dt3$variable[plot_dt3$variable == "geneConsistRatio"] <- paste0("genes\n(adj. p-val<=", geneSignifThresh, ")")
plot_dt3$variable[plot_dt3$variable == "tadConsistRatio"] <- paste0("TADs\n(adj. p-val<=", tadSignifThresh, ")")

plot_dt3$sampLab <- factor(plot_dt3$sampLab, levels=labLevels)


init_p3 <- custom_plot(ggplot(data = plot_dt3,aes_string(x="sampLab",y="value", fill = "variable", color="variable"))+
                         ggtitle("Ratio consist. signif. and sample size", subtitle = subTit)+
  scale_x_discrete(name= "# of samples")+
  scale_y_continuous(name="Ratio of consist. signif. features", breaks=scales::pretty_breaks(n = 10))+
  labs(fill  = paste0("Ratio consist. signif.") , color  = paste0("Ratio consist. signif.") ) )

ratioMatch_p <- init_p3 +  
  geom_bar(stat="identity", position="dodge")+
  guides(color = FALSE)


outFile <- file.path(outFolder, paste0(init_hicds, "_", init_exprds, "_signifSubSamplesRatioMatch_boxplot.", plotType))
ggsave(plot = ratioMatch_p, filename = outFile, height=myHeightGG, width = myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


# ratioMatch_p_strip <- init_p3 + 
#   geom_jitter(
#     position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
#     size = 1.2
#   ) +
#   stat_summary(
#     fun.data="mean_sdl",  fun.args = list(mult=1), 
#     geom = "pointrange",  size = 0.4,
#     position = position_dodge(0.8)
#   ) +
#   guides(fill = FALSE)
# 
# 
# 
# outFile <- file.path(outFolder, paste0(init_hicds, "_", init_exprds, "_signifSubSamplesRatioMatch_stripchart.", plotType))
# ggsave(plot = ratioMatch_p_strip, filename = outFile, height=myHeightGG, width = myWidthGG)
# cat(paste0("... written: ", outFile, "\n"))


