cmp_names <- c(

    "wt_vs_mut" = "wt. vs. mut.",
    "norm_vs_tumor" = "normal vs. tumor",
    "subtypes" = "subtypes"

) 

hicds_names <- c( 
"Barutcu_MCF-10A_40kb" = "Barutcu_MCF-10A",            
  "Barutcu_MCF-7_40kb" = "Barutcu_MCF-7",               
 "ENCSR079VIJ_G401_40kb" = "G401",
 "ENCSR312KHQ_SK-MEL-5_40kb"= "SK-MEL-5",
 "ENCSR346DCU_LNCaP_40kb" = "LNCaP",           
 "ENCSR401TBQ_Caki2_40kb" = "Caki2",
  "ENCSR444WCZ_A549_40kb" = "A549", 
 "ENCSR489OCU_NCI-H460_40kb" = "NCI-H460", 
 "ENCSR504OTV_transverse_colon_40kb" = "Transv. colon",
 "ENCSR549MGQ_T47D_40kb" = "T47D",
 "ENCSR862OGI_RPMI-7951_40kb" = "RPMI-7951",  
 "GSE105194_cerebellum_40kb" = "Cerebellum",
"GSE105194_spinal_cord_40kb" = "Spinal cord",
 "GSE105318_DLD1_40kb" = "DLD1",
 "GSE105381_HepG2_40kb" ="HepG2",    
 "GSE109229_BT474_40kb" = "BT474",   
 "GSE109229_SKBR3_40kb" ="SKBR3",  
 "GSE118514_22Rv1_40kb" ="22Rv1",   
 "GSE118514_RWPE1_40kb"  = "RWPE1",          
 "GSE118588_Panc_beta_40kb" ="Panc. \u03B2-cells" ,
 "GSE99051_786_O_40kb" = "786_O",              
 "HMEC_40kb" = "HMEC",      
 "K562_40kb" ="K562",              
 "LG1_40kb" ="Lung (1)",               
 "LG2_40kb" = "Lung (2)",              
 "LI_40kb" = "Liver",        
"PA2_40kb" = "Pancreas (1)",                 
 "PA3_40kb" = "Pancreas (2)",                          
 "Panc1_rep12_40kb" = "Pancreas (3)",
 "Rao_HCT-116_2017_40kb"   = "HCT-116"
)         

exprds_names <- c(
  "TCGAbrca_lum_bas" = "BRCA luminal vs. basal",             
  "TCGAkich_norm_kich" = "KICH normal vs. tumor",           
 "TCGAskcm_lowInf_highInf" = "SKCM low inf. vs. high inf.", 
 "TCGAskcm_wt_mutBRAF" = "SKCM BRAFwt vs. BRAFmut",          
 "TCGAskcm_wt_mutCTNNB1" = "SKCM CTNNB1wt vs. CTNNB1mut",          
 "TCGAprad_norm_prad" = "PRAD normal vs. tumor",           
 "TCGAluad_mutKRAS_mutEGFR" = "LUAD KRASmut vs. EGFRmut",
 "TCGAluad_nonsmoker_smoker" = "LUAD nonsmoker vs. smoker",
 "TCGAluad_norm_luad" = "LUAD normal vs. tumor",
 "TCGAluad_wt_mutKRAS" = "LUAD KRASwt vs. KRASmut",          
 "TCGAlusc_norm_lusc" = "LUSC normal vs. tumor",
 "TCGAcoad_msi_mss" = "COAD MSI vs. MSS",  
 "TCGAgbm_classical_mesenchymal" = "GBM classical vs. mesenchymal",
 "TCGAgbm_classical_neural"      = "GBM classical vs. neural",
 "TCGAgbm_classical_proneural"  = "GBM classical vs. proneural",
 "TCGAlgg_IDHwt_IDHmutnc" = "LGG IDHwt vs. IDHmut",
 "TCGAlihc_norm_lihc" = "LIHC normal vs. tumor",
 "TCGAlihc_wt_mutCTNNB1" = "LIHC CTNNB1wt vs. CTNNB1mut",
 "TCGApaad_wt_mutKRAS" = "PAAD KRASwt vs. KRASmut",
 "TCGAlaml_wt_mutFLT3" = "LAML FLT3wt vs. FLT3mut"
)

cond1_names <- sapply(exprds_names, function(x) gsub(".+ (.+) vs. (.+)", "\\1", x))
names(cond1_names) <- sapply(names(exprds_names), function(x) gsub("TCGA.+_(.+)_(.+)", "\\1", x))
cond2_names <- sapply(exprds_names, function(x) gsub(".+ (.+) vs. (.+)", "\\2", x))
names(cond2_names) <- sapply(names(exprds_names), function(x) gsub("TCGA.+_(.+)_(.+)", "\\2", x))



