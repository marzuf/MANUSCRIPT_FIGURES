# Rscript fcc_density_example.R LG1_40kb/TCGAluad_mutKRAS_mutEGFR ENCSR312KHQ_SK-MEL-5_40kb/TCGAskcm_lowInf_highInf

options(scipen = 100)

SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

script_name <- "fcc_density_example.R"
cat("> START ", script_name, "\n")
startTime <- Sys.time()

script0_name <- "0_prepGeneData"

require(foreach)
require(doMC)
require(ggpubr)
require(ggplot2)
registerDoMC(40)

plotType <- "svg"
myHeight <- 400
myWidth <- 400
myHeightGG <- 7
myWidthGG <- 7

plotCex <- 1.4


all_ds <- commandArgs(trailingOnly = T)
# hicds <- args[seq(from=1, to=length(args), by=2)]
# exprds <- args[seq(from=2, to=length(args), by=2)]
stopifnot(length(all_ds) >= 1)

source("../../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../settings.R")

pipFolder<- runFolder
stopifnot(dir.exists(pipFolder))

source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")

pipOutFolder <- file.path(pipFolder, "PIPELINE", "OUTPUT_FOLDER")
stopifnot(dir.exists(pipOutFolder))

myHeightGG <- 7
myWidthGG <- 10

outFolder <- "FCC_DENSITY_EXAMPLE"
dir.create(outFolder, recursive = TRUE)

all_dt <- foreach(ds  = all_ds, .combine='rbind') %dopar% {
  hicds <- dirname(ds)
  exprds <- basename(ds)
  fcc_file <- file.path(pipOutFolder, hicds, exprds, "8cOnlyFCC_runAllDown", "all_obs_prodSignedRatio.Rdata")
  if(!file.exists(fcc_file)) {
    fcc_file <- file.path(pipOutFolder, hicds, exprds, "8cOnlyFCConlyObs_runAllDown", "all_obs_prodSignedRatio.Rdata")  
  }
  stopifnot(file.exists(fcc_file))
  
  data.frame(
    hicds=hicds, 
    exprds=exprds,
    fcc = as.numeric(get(load(fcc_file))),
    stringsAsFactors = FALSE
  )
}

save(all_dt, file=file.path(outFolder, "all_dt.Rdata"))

all_dt$dataset <- paste0(hicds_names[as.character(all_dt$hicds)], " \n ", exprds_names[as.character(all_dt$exprds)])
stopifnot(!is.na(all_dt$dataset))

if(length(all_ds) == 2) {
  require(paletteer)
  # mycols <- as.character(paletteer::paletteer_d("nord::lake_superior")[3:4])
  mycols <- c("#C87D4BFF" ,"#4B647DFF")
} else{
  stop("need to define color scale manually\n")
}

plotTit <- paste0("FCC distribution")
subTit <- ""

dist_p <- ggplot(all_dt, aes(x=fcc, fill=dataset, col=dataset)) + 
  ggtitle(plotTit, subtitle = subTit)+
  labs(fill="Dataset", x = paste0("TAD FCC"), y ="Density") + 
  guides(color=FALSE)+
  geom_histogram(aes(y=..density..), alpha=0.5,
                 position="identity")+
  geom_density(alpha=.2) +
  # geom_histogram(aes(y=..density..))+
  # geom_density(alpha=.3)+
  # eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
  # eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
  scale_fill_manual(values=mycols)+
  scale_color_manual(values=mycols)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8))+
  my_box_theme+
  theme(
    axis.line = element_line(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16),
    legend.title =element_text(size=14),
    legend.text = element_text(margin = margin(b = 0.1, unit = 'in'), size=12)
  )

outFile <- file.path(outFolder, paste0("FCC_hist_and_density_", paste0(gsub("/", "_", all_ds), collapse="_"), ".", plotType))
ggsave(dist_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



dist_p <- ggplot(all_dt, aes(x=fcc, fill=dataset, col=dataset)) + 
  ggtitle(plotTit, subtitle = subTit)+
  labs(fill="Dataset", x = paste0("TAD FCC"), y ="Density") + 
  guides(color=FALSE)+
  geom_histogram(aes(y=..density..), alpha=0.5,
                 position="identity")+
  # geom_density(alpha=.2) +
  # geom_histogram(aes(y=..density..))+
  # geom_density(alpha=.3)+
  # eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
  # eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
  scale_fill_manual(values=mycols)+
  scale_color_manual(values=mycols)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8))+
  my_box_theme+
  theme(
    axis.line = element_line(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16),
    legend.title =element_text(size=14),
    legend.text = element_text(margin = margin(b = 0.1, unit = 'in'), size=12)
  )

outFile <- file.path(outFolder, paste0("FCC_histOnly_", paste0(gsub("/", "_", all_ds), collapse="_"), ".", plotType))
ggsave(dist_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



dist_p <- ggplot(all_dt, aes(x=fcc, fill=dataset, col=dataset)) + 
  ggtitle(plotTit, subtitle = subTit)+
  labs(fill="Dataset", x = paste0("TAD FCC"), y ="Density") + 
  guides(color=FALSE)+
  # geom_histogram(aes(y=..density..), alpha=0.5,
  #                position="identity")+
  geom_density(alpha=.2) +
  # geom_histogram(aes(y=..density..))+
  # geom_density(alpha=.3)+
  # eval(parse(text=paste0("scale_fill_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
  # eval(parse(text=paste0("scale_color_", ggsci_pal, "(", ggsci_subpal, ")"))) + 
  scale_fill_manual(values=mycols)+
  scale_color_manual(values=mycols)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8))+
  my_box_theme+
  theme(
    axis.line = element_line(),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16),
    legend.title =element_text(size=14),
    legend.text = element_text(margin = margin(b = 0.1, unit = 'in'), size=12)
  )

outFile <- file.path(outFolder, paste0("FCC_densityOnly_", paste0(gsub("/", "_", all_ds), collapse="_"), ".", plotType))
ggsave(dist_p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))





