## DEG Analysis Code Reformatting -Ongoing Original Working code still there 889 to 1367

#### COMPARISON BETWEEN GROUPS [VOLCANO PLOT]----
# Inputs: Sample_Data, Counts_Data, Feature_Data, Normalized Counts Data
#Sample_Data should have had the SlideName alrerady modified by now, and ROI Column added,
#Ensure that it is indeed added.
# Sample_Data <- Sample_Data %>% mutate(SlideName = gsub("BBP", "BFP", SlideName))
# Sample_Data[, "ROI"] <-  rownames(Sample_Data)
#Normalized counts data should be imported from the previously prepared data
Sample_Data <- Sample_Data #Make sure Sample_Data has ROI column extracted from the row names previously
Normalized_Counts_Data <- Normalized_Counts_Data_TMM

# Ensure the input files exist
if (exists ("Sample_Data") & exists ("Counts_Data") & exists ("Feature_Data") &
    exists ("Normalized_Counts_Data_TMM")) {
  print("All Required Files Exist")
} else {
  print("All Input Files Don't Exist, Load Them")
}


COMPARE_TABLE <- as.data.frame(cbind(
  COMPARE_GROUP_NAME=c("Samples", "Tumor", "Tumor"),
          GROUP1_NAME=c("LUCAT1_KD","TUMOR_WHOLE", "TUMOR_WHOLE"),
          GROUP2_NAME=c("CONTROL","TUMOR_EDGE", "TUMOR_INSIDE")
  ))

comp =3 #[INPUT_NEEDED]
for (comp in 3:3) { #[INPUT_NEEDED]
  print (comp)
  COMPARE_GROUP_NAME = COMPARE_TABLE[comp,"COMPARE_GROUP_NAME"] 
  GROUP1_NAME = COMPARE_TABLE[comp,"GROUP1_NAME"] 
  GROUP2_NAME <-  COMPARE_TABLE[comp,"GROUP2_NAME"] 
  
  print(paste(COMPARE_GROUP_NAME,":",GROUP1_NAME,"_Vs_",GROUP2_NAME))
  #Sample Data
  Experiment_Sample_Data <- Sample_Data
  Experiment_Sample_Data[, "ROI"] <-  rownames(Experiment_Sample_Data)
  Comparision_Sample_Data <- Experiment_Sample_Data[,c("ROI", "ScanLabel","Samples","Tumor")]
  Comparision_Sample_Data[,"Group"] <- Experiment_Sample_Data[,COMPARE_GROUP_NAME]
  print(table(Comparision_Sample_Data[,"Group"]))
  
  Comparision_Sample_Data <- Comparision_Sample_Data[which(Comparision_Sample_Data[,"Group"]==GROUP1_NAME | Comparision_Sample_Data[,"Group"]==GROUP2_NAME),]
  rownames(Comparision_Sample_Data) <- Comparision_Sample_Data[,"ROI"]
  
  print(table(Comparision_Sample_Data[,"Group"]))
  
  GROUP1_SAMPLES <- as.character(Comparision_Sample_Data[which(Comparision_Sample_Data[,"Group"]==GROUP1_NAME),"ROI"])
  GROUP2_SAMPLES <- as.character(Comparision_Sample_Data[which(Comparision_Sample_Data[,"Group"]==GROUP2_NAME),"ROI"])
  
  #Counts Data
  Experiment_Counts_Data <- Counts_Data
  Comparision_Counts_Data <- Experiment_Counts_Data[,c(GROUP1_SAMPLES,GROUP2_SAMPLES)]
  print(dim(Comparision_Counts_Data))
  
  #Normalized Counts Data
  Experiment_Normalised_Counts <- Normalized_Counts_Data_TMM
  Comparision_Normalised_Counts_Data <- Experiment_Normalised_Counts[,c(GROUP1_SAMPLES,GROUP2_SAMPLES)]
  print(dim(Comparision_Normalised_Counts_Data))
  
  #Feature Data
  Comparision_Feature_Data <- Feature_Data[,c("GeneID","TargetGene")]
  colnames(Comparision_Feature_Data) <- c("GeneID","Gene")
  rownames(Comparision_Feature_Data) <- Comparision_Feature_Data[,"Gene"]
  #######################################################################
  ######################################################################
  EXP_GeoMX_DSP_Data_FOR_COMPARISION <- NULL
  EXP_GeoMX_DSP_Data_FOR_COMPARISION <- SpatialExperiment(assays = list(counts = Comparision_Counts_Data,logcounts=Comparision_Normalised_Counts_Data))
  colData(EXP_GeoMX_DSP_Data_FOR_COMPARISION) <- DataFrame(Comparision_Sample_Data)
  rowData(EXP_GeoMX_DSP_Data_FOR_COMPARISION) <- DataFrame(Comparision_Feature_Data)
  ########################################################################
  COMPARISION_NAME <- paste(COMPARE_GROUP_NAME,"_",GROUP1_NAME,"_Vs_",GROUP2_NAME,sep="")
  
  DEG_parent_dir <- paste(output_dir, "/", "DEG", sep="")
  DEG_test_dir <- paste(DEG_parent_dir, "/", COMPARISION_NAME, sep="")
  
  DEG_TEST <- "EDGER"
  # Using Limma-Voom Pipeline to do a DEG Analysis, it requires a DGEList Framework, SE2DGEList function converts the same. 
  ##############################################################################################################
  dge <- SE2DGEList(EXP_GeoMX_DSP_Data_FOR_COMPARISION)
  ##############################################################################################################
  # Create a Design Matrix for GROUP1 vs GROUP2 Analysis
  ####################################
  design <- model.matrix(~0 + Group , data = colData(EXP_GeoMX_DSP_Data_FOR_COMPARISION))
  colnames(design) <- gsub("^Group","",colnames(design))
  colnames(design) <- gsub(" ","_",colnames(design))
  write.table(design,paste(output_dir, "/", GROUP1_NAME,"_VS_",GROUP2_NAME,"_DESIGN_MATRIX.txt",sep=""),sep="\t",row.names = T,col.names = T,quote = T)
  print("EDGER")
  print(dim(design))
  print(colSums(design))
  
  # Using the Limma framework create a Contrast levels between GROUP1 and GROUP2. 
  #makeContrasts (Limma) = Construct the contrast matrix corresponding to specified contrasts of a set of parameters.
  ####################################
  contr.matrix <- makeContrasts(GROUP1_Vs_GROUP2 = get(GROUP1_NAME) - get(GROUP2_NAME),levels = colnames(design))
  write.table(contr.matrix,paste(output_dir, "/", GROUP1_NAME,"_VS_",GROUP2_NAME,"_CONTRAST_MATRIX.txt",sep=""),sep="\t",row.names = T,col.names = T,quote = T)
  
  #Based on the suggestions, we remove the filter out genes with low coverage in the dataset to allow a more accurate mean-variance relationship 
  #and reduce the number of statistical tests. Here we use the filterByExpr function from the edgeR package to filter genes based on the model matrix, 
  #keeping as many genes as possible with reasonable counts.
  ########################################################
  keep <- filterByExpr(dge, design)
  #######################################################3
  LowCoverage_Genes <- as.data.frame(as.character(rownames(dge)[!keep]))
  colnames(LowCoverage_Genes) <- "Gene"
  write.table(LowCoverage_Genes,paste(output_dir, "/", GROUP1_NAME,"_VS_",GROUP2_NAME,"_LowCoverage_Genes.txt",sep=""),sep="\t",row.names = F,col.names = T,quote = T)
  ######################################################
  ############################################################
  # Filter the Low expressing genes for the analysis
  ###########################################################
  #dge_all <- dge[keep, ]
  dge_all <- dge
  
  ##### BCV check
  # Biological CV (BCV) is the coefficient of variation with which the (unknown) true abundance of the gene varies 
  # between replicate RNA samples. For more detail about dispersion and BCV calculation
  ###########################################################################
  dge_all <- estimateDisp(dge_all, design = design, robust = TRUE)
  ##############################################################
  bcv_df <- data.frame(
    'BCV' = sqrt(dge_all$tagwise.dispersion),
    'AveLogCPM' = dge_all$AveLogCPM,
    'gene_id' = rownames(dge_all)
  )
  bcv_df[,"HighBCV"] <- ifelse(bcv_df[,"BCV"]>0.8,"HIGH-BCV","LOW-BCV")
  bcv_df <- bcv_df[order(bcv_df[,"BCV"],decreasing = T),]
  write.table(bcv_df,paste(output_dir, "/", GROUP1_NAME,"_VS_",GROUP2_NAME,"_Genes_Biological_Variation_Coefficient.txt",sep=""),sep="\t",row.names = F,col.names = T,quote = T)
  ######################################################
  highbcv <- bcv_df$BCV > 0.8
  highbcv_df <- bcv_df[highbcv, ]
  ###############################################
  # BCV Plot
  # Change it to ggplot
  ###############################################
  graphics.off()
  pdf(paste(output_dir, "/", GROUP1_NAME,"_VS_",GROUP2_NAME,"_Genes_Biological_Variation_Coefficient.pdf",sep=""),width = 10,height = 10)
  plotBCV(dge_all, legend.position = "topleft", ylim = c(0, 1.3))
  points(highbcv_df$AveLogCPM, highbcv_df$BCV, col = "#FF3158")
  text(highbcv_df$AveLogCPM, highbcv_df$BCV, labels = highbcv_df$gene_id, pos = 4,cex = 0.5)
  graphics.off()
  
  ### Differential Expression Analysis using   limma-voom pipeline
  ############################################################################
  graphics.off()
  pdf(paste(output_dir, "/", GROUP1_NAME,"_VS_",GROUP2_NAME,"_Mean_Variance_Trend.pdf",sep=""),width = 10,height = 10)
  limma_voom_v <- voom(dge_all, design, plot = TRUE) 
  graphics.off()
  
  # Using EdgeR
  ####################################################################################
  ########################
  #glmFit fits genewise negative binomial glms, all with the same design matrix but possibly different dispersions, offsets and weights. 
  #When the design matrix defines a one-way layout, or can be re-parametrized to a one-way layout, the glms are fitting very quickly using 
  ########################
  edgeR_gfit <- glmFit(dge_all, design = design)
  #####################
  EdgeR_Normalised_Counts_Data <- as.data.frame(edgeR_gfit[[2]])
  ###########################################################################
  GROUP1_VS_GROUP2_DE_EDGER_MEANS <- as.data.frame(cbind(rownames(EdgeR_Normalised_Counts_Data),
                                     GROUP1_Mean  = rowMeans(log(EdgeR_Normalised_Counts_Data[,GROUP1_SAMPLES],2),na.rm = T),
                                     GROUP2_Mean  = rowMeans(log(EdgeR_Normalised_Counts_Data[,GROUP2_SAMPLES],2),na.rm = T),
                                     MEANVAL = rowMeans(log(EdgeR_Normalised_Counts_Data[,c(GROUP2_SAMPLES,GROUP1_SAMPLES)],2),na.rm = T),
                                     LOG2FC = rowMeans(log(EdgeR_Normalised_Counts_Data[,GROUP1_SAMPLES],2),na.rm = T) - rowMeans(log(EdgeR_Normalised_Counts_Data[,GROUP2_SAMPLES],2),na.rm = T)
                                                         
  )) 
  
  colnames(GROUP1_VS_GROUP2_DE_EDGER_MEANS)[1] <- "Gene"
  
  #####################
  #glmLRT conducts likelihood ratio tests for one or more coefficients in the linear model. 
  #If coef is used, the null hypothesis is that all the coefficients indicated by coef are equal to zero. 
  #If contrast is non-null, then the null hypothesis is that the specified contrasts of the coefficients are equal to zero. 
  #For example, a contrast of c(0,1,-1), assuming there are three coefficients, would test the hypothesis that the second and third coefficients are equal.
  ####################
  edgeR_glrt <- glmLRT(edgeR_gfit, design, contrast = contr.matrix)
  edgeR_glrt_DEG <- as.data.frame(edgeR_glrt)
  edgeR_glrt_DEG[,"adj.P.Val"] <- p.adjust(as.numeric(as.character(edgeR_glrt_DEG[,"PValue"])),method = "BH",n = as.numeric(length(edgeR_glrt_DEG[,1])))
  ########################################
  #logFC	:log2-fold change of expression between conditions being tested.
  #logCPM	 :average log2-counts per million, the average taken over all libraries in y.
  #LR	:likelihood ratio statistics.
  #PValue	:p-values.
  #########################################
  GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS <- edgeR_glrt_DEG
  #######################################################################
  GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS <- GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,c("GeneID","Gene","logCPM","logFC","PValue","adj.P.Val")]
  colnames(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS) <- c("GeneID","Gene","EDGER_MEANVAL","EDGER_LOG2FC","EDGER_PVAL","EDGER_PVAL_ADJUST")
  # GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS <- plyr:::join(GROUP1_VS_GROUP2_DE_EDGER_MEANS,GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS,by="Gene",type="full",match="all")
  GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS <- GROUP1_VS_GROUP2_DE_EDGER_MEANS %>%
    full_join(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS, by = "Gene")
  ##################################################
  ##################################################
  GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"EDGER_PVAL_ADJUST"] <- as.numeric(as.character(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"EDGER_PVAL_ADJUST"]))
  GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"EDGER_PVAL"] <- as.numeric(as.character(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"EDGER_PVAL"]))
  GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"EDGER_LOG2FC"] <- as.numeric(as.character(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"EDGER_LOG2FC"]))
  # ################################################
  GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"PVALADJ_SIG_GENES"] <- ifelse(((GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"EDGER_PVAL_ADJUST"]<0.05) & (abs(as.numeric(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"EDGER_LOG2FC"]))>=log(2,2))),as.numeric(sign(as.numeric(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"EDGER_LOG2FC"]))),NA)
  GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"PVAL_SIG_GENES"] <- ifelse(((GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"EDGER_PVAL"]<0.05) & (abs(as.numeric(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"EDGER_LOG2FC"]))>=log(2,2))),as.numeric(sign(as.numeric(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"EDGER_LOG2FC"]))),NA)
  GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS <- GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[order(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"EDGER_LOG2FC"],decreasing = T),]
  GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS <- GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[order(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"PVAL_SIG_GENES"],decreasing = T),]
  GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS <- GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[order(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"PVALADJ_SIG_GENES"],decreasing = T),]  
  GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS <- GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,c("Gene","GeneID",setdiff(colnames(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS),c("Gene","GeneID")))]
  ##########################################
  colnames(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS) <- sapply(colnames(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS),function(x) gsub("GROUP1",GROUP1_NAME,as.character(x)))
  colnames(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS) <- sapply(colnames(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS),function(x) gsub("GROUP2",GROUP2_NAME,as.character(x)))
  ###########################################
  write.table(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS,paste(output_dir, "/", GROUP1_NAME,"_VS_",GROUP2_NAME,"_",DEG_TEST,"_RESULTS.txt",sep=""),row.names = F,col.names = T,quote = T,sep="\t")
  ########################################################################
  print(table(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"PVALADJ_SIG_GENES"]))
  print(table(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"PVAL_SIG_GENES"]))
  
  ##########################################################################
  
  
  FC_MIN <- floor(min(c(unlist(as.numeric(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"EDGER_LOG2FC"])))))#,
  # unlist(as.numeric(GROUP1_VS_GROUP2_DE_LIMMA_ANALYSIS[,"LIMMA_LOG2FC"])),
  # unlist(as.numeric(GROUP1_VS_GROUP2_DE_TTEST_ANALYSIS[,"TTEST_LOG2FC"])))))
  
  FC_MAX <- ceiling(max(c(unlist(as.numeric(GROUP1_VS_GROUP2_DE_EDGER_ANALYSIS[,"EDGER_LOG2FC"])))))#,
  # unlist(as.numeric(GROUP1_VS_GROUP2_DE_LIMMA_ANALYSIS[,"LIMMA_LOG2FC"])),
  # unlist(as.numeric(GROUP1_VS_GROUP2_DE_TTEST_ANALYSIS[,"TTEST_LOG2FC"])))))
  
  
  L2FC_MIN <- -1*(max(abs(c(FC_MIN,FC_MAX))))
  L2FC_MAX <- 1*(max(abs(c(FC_MIN,FC_MAX))))
  
  
  
  
  DE_Analysis = paste(COMPARE_GROUP_NAME,":",GROUP1_NAME,"_Vs_",GROUP2_NAME,sep="")
  
  
  #### Nested Loop Will Start
  DEG_TESTS <- c("EDGER")#,"LIMMA","TTEST")
  Df_Positive_List <- NULL
  Df_Negative_List <- NULL
  test = 1
  COMP_SIG_GENES <- NULL
  
  for(test in 1:length(DEG_TESTS)){
    
    TEST = DEG_TESTS[test]
    print(TEST)
    # TEST_DIR <- paste(output_dir, "/", TEST,"/",sep="")
    TEST_DATA <- read.delim(paste(output_dir, "/", GROUP1_NAME,"_VS_",GROUP2_NAME,"_",TEST,"_RESULTS.txt",sep=""),sep="\t",header = T,check.names = F,stringsAsFactors = F)
    ############################################################
    Df_Positive_List[[TEST]] <- as.character(TEST_DATA[which(TEST_DATA[,"PVAL_SIG_GENES"]==1),"Gene"])
    Df_Negative_List[[TEST]] <- as.character(TEST_DATA[which(TEST_DATA[,"PVAL_SIG_GENES"]== -1),"Gene"])
    #############################################################
    PLOT_DATA <- NULL
    
    LOG2FC_COL <- paste(TEST,"_LOG2FC",sep="")
    MEANEXP_COL <- paste(TEST,"_MEANVAL",sep="")
    PVAL_COL <- paste(TEST,"_PVAL",sep="")
    PVALADJ_COL <- paste(TEST,"_PVAL_ADJUST",sep="")
    
    PLOT_DATA <- NULL
    PLOT_DATA <- TEST_DATA[,c("Gene",LOG2FC_COL,PVAL_COL,PVALADJ_COL)]
    colnames(PLOT_DATA) <- c("GENE","LOG2FC","PVAL","PVALADJ")
    PLOT_DATA$THRESHOLD <- 0
    PLOT_DATA$THRESHOLD <- ifelse((PLOT_DATA[,"PVAL"]<0.05),1,PLOT_DATA[,"THRESHOLD"])
    ###############################
    PLOT_DATA$THRESHOLD <- ifelse((PLOT_DATA[,"PVAL"]<0.05 & abs(PLOT_DATA[,"LOG2FC"])>=as.numeric(log(2,2))),2,PLOT_DATA[,"THRESHOLD"])
    PLOT_DATA$THRESHOLD <- ifelse((PLOT_DATA[,"PVALADJ"]<0.05 & abs(PLOT_DATA[,"LOG2FC"])>=as.numeric(log(2,2))),3,PLOT_DATA[,"THRESHOLD"])
    ################################
    
    
    #Milan Dont transform before filtering; do this later
    
    # PLOT_DATA[,"PVAL"] <- -log10(PLOT_DATA[,"PVAL"])
    ################################################
    PVAL_th <- -log10(0.05)
    PVALADJ_th <- -log10(0.01)
    ################################################
    Max_Th <- L2FC_MAX
    Min_Th <- L2FC_MIN
    #################################################
    PLOT_DATA <- PLOT_DATA[order(PLOT_DATA[,"LOG2FC"],decreasing = T),]
    rownames(PLOT_DATA) <- NULL
    PLOT_DATA[,"Rank"] <- as.numeric(rownames(PLOT_DATA))
    ###~~~~~~~~~~~~
    #SSSS Highlight significant genes if they have a THRESHOLD value of 2 or higher and are among the top 25 ranked genes based on their LOG2FC values.
    PLOT_DATA[,"SIG_GENES"] <- NA
    PLOT_DATA[,"SIG_GENES"] <- ifelse(((abs(PLOT_DATA[,"THRESHOLD"])>=2) & (PLOT_DATA[,"Rank"]<=25)),PLOT_DATA[,"GENE"],NA)
    PLOT_DATA <- PLOT_DATA[order(PLOT_DATA[,"LOG2FC"],decreasing = F),]
    rownames(PLOT_DATA) <- NULL
    PLOT_DATA[,"Rank"] <- as.numeric(rownames(PLOT_DATA))
    PLOT_DATA[,"SIG_GENES"] <- ifelse(((abs(PLOT_DATA[,"THRESHOLD"])>=2) & (PLOT_DATA[,"Rank"]<=25)),PLOT_DATA[,"GENE"],PLOT_DATA[,"SIG_GENES"])
    ###~~~~~~~~~~~~
    #Keep ~~~ OR ^^^ Segment
    ###^^^^^^^^^^^^^^^
    # Milan Set significance thresholds #[INPUT_NEEDED]
    hypoxia_genes <- c("HIF1A", "VEGFA", "PGK1", "LDHA", "ADM", "CA9", "EPO", "NDRG1", "BNIP3", "SLC2A1", "ENO1", "COL1A1", "TPI1", "ALDOA", "PGAM1", "HBEGF", "ANXA1", "P4HA1", "LOX", "PFN2")
    wnt_genes_mouse_test <- c("Fzd10","Wnt3","Wnt6","Wnt5b","Wnt9b","Wnt11","Rspo3","Dkk4", "Draxin", "Ngf","Snai2","Sox2","Sox17","Adamts5","Adam11")
    log2fc_threshold <- log2(2)  # 2-fold change
    pval_threshold <- 0.05 #-log10(0.05)
    
    # Milan Identify significant genes from the curated list
    # PLOT_DATA[,"SIG_GENES"] <- NA
    # PLOT_DATA[,"SIG_GENES"] <- ifelse(
    #   PLOT_DATA[,"GENE"] %in% hypoxia_genes & #[INPUT_NEEDED] immune_pathways_genes_unique_c2, wnt_genes_mouse_test or wnt_genes_mouse, immune_pathways_genes_unique, naChannel_pathways_genes_unique
    #     abs(PLOT_DATA[,"LOG2FC"]) >= log2fc_threshold &
    #     PLOT_DATA[,"PVAL"] < pval_threshold, # Note: '>' instead of '<' because of -log10 transformation
    #   PLOT_DATA[,"GENE"],
    #   NA
    # )
    ###^^^^^^^^^^^^^^^
    
    
    # Transform p-values to -log10 scale AFTER filtering
    PLOT_DATA[,"PVAL"] <- -log10(PLOT_DATA[,"PVAL"])
    # Create a new column for coloring points
    # PLOT_DATA[,"COLOR"] <- ifelse(PLOT_DATA[,"GENE"] %in% wnt_genes_mouse & 
    #                                 abs(PLOT_DATA[,"LOG2FC"]) >= log2fc_threshold &
    #                                 PLOT_DATA[,"PVAL"] > pval_threshold,
    #                               "Highlighted", "Not Highlighted")
    
    ## Optional: Limit number of labeled genes to prevent overcrowding (Top 100 significant genes) #[INPUT_NEEDED] or comment out below section
    ###****###
    # sig_genes <- PLOT_DATA[!is.na(PLOT_DATA[,"SIG_GENES"]),]
    # sig_genes <- sig_genes[order(sig_genes[,"PVAL"], decreasing = TRUE),]  # Note: decreasing = TRUE because of -log10 transformation
    # top_n <- min(100, nrow(sig_genes))
    # top_genes <- sig_genes[1:top_n, "GENE"]
    # PLOT_DATA[,"SIG_GENES"] <- ifelse(PLOT_DATA[,"GENE"] %in% top_genes, PLOT_DATA[,"GENE"], NA)
    ###^^^^^^^^^^^^^^^^^^^^^^^^^^
    
    ###################################################
    x_limits <- c(NA, -2) #Volcanoplot label restricting at a location.
    TEST_SIG_GENES <- as.character(PLOT_DATA[complete.cases(PLOT_DATA[,"SIG_GENES"]),"SIG_GENES"])
    COMP_SIG_GENES <- c(COMP_SIG_GENES,unlist(TEST_SIG_GENES))
    Volcano_Plot <- NULL
    Volcano_Plot <- ggplot(PLOT_DATA,aes(x = LOG2FC, y = PVAL,label=SIG_GENES, size = factor(THRESHOLD)))+
      #geom_label_repel(aes(label=SIG_GENES),box.padding= 0.3,point.padding = 0.5,segment.color="grey80",na.rm =T,colour = "black",fill="white",size = 3,force = 4,direction = "both",max.overlaps = Inf,nudge_x = 0,nudge_y = 0)+
      geom_point(aes(x = LOG2FC, y = PVAL,size = factor(THRESHOLD),colour = factor(THRESHOLD),alpha=factor(THRESHOLD))) +
      geom_point(data=PLOT_DATA[!is.na(PLOT_DATA[,"SIG_GENES"]),],aes(size = factor(THRESHOLD)),alpha = 0.6,colour="black",shape=1,pch = 21,stroke = 1,show.legend = F)+
      
      # geom_label_repel(aes(label=SIG_GENES),box.padding= 0.3,point.padding = 0.5,segment.color="grey80",na.rm =T,colour = "black",fill="white",size = 1.5,force = 4,direction = "both",max.overlaps = Inf,nudge_x = 0,nudge_y = 0)+
      geom_vline(xintercept = 0,color = "black", linetype='dashed',color = "#4268F4")+
      geom_vline(xintercept = as.numeric(log(2,2)),color = "black", linetype='dashed',color = "#4268F4",alpha = 0.5)+
      geom_vline(xintercept = -as.numeric(log(2,2)),color = "black", linetype='dashed',color = "#4268F4",alpha = 0.5)+
      geom_vline(xintercept = as.numeric(log(2,2)),color = "black", linetype='dashed',color = "#4268F4",alpha = 0.5)+
      geom_vline(xintercept = -as.numeric(log(2,2)),color = "black", linetype='dashed',color = "#4268F4",alpha = 0.5)+
      geom_hline(yintercept = PVAL_th, linetype='dashed',color = "#4268F4")+
      # geom_hline(yintercept = PVALADJ_th, linetype='dashed',color = "#4268F4")+
      geom_hline(yintercept = 0,color = "black")+
      # geom_text(x=L2FC_MIN+1, y=ceiling(max(PLOT_DATA[,"PVAL"],na.rm = T)), label=gsub("_","\n",paste(GROUP2_NAME,"_(N:",length(GROUP2_SAMPLES),")",sep="")),show.legend = F,size=5,color = "#4268F4") +
      # geom_text(x=L2FC_MAX-1, y=ceiling(max(PLOT_DATA[,"PVAL"],na.rm = T)), label=gsub("_","\n",paste(GROUP1_NAME,"_(N:",length(GROUP1_SAMPLES),")",sep="")),show.legend = F,size=5,color = "#303030") +
      #[INPUT_NEEDED] for BFP, since the name is BBP replace appropriate group name GROUP2_NAME to "BFP"
      geom_text(x=L2FC_MIN+1, y=ceiling(max(PLOT_DATA[,"PVAL"],na.rm = T)), label=GROUP2_NAME,show.legend = F,size=20, color = "#4268F4", check_overlap = TRUE) +
      geom_text(x=L2FC_MAX-1, y=ceiling(max(PLOT_DATA[,"PVAL"],na.rm = T)), label=GROUP1_NAME,show.legend = F,size=20,color = "#4268F4", check_overlap = TRUE) +
      
      xlim(c(L2FC_MIN,L2FC_MAX))+
      scale_color_manual(values = c("grey50","#FF3158","#4268F4","#42B858"),breaks = c(0,1,2,3,4)) +
      # scale_color_manual(values = c("Highlighted" = "red", "Not Highlighted" = "grey")) +
      
      # geom_text_repel(aes(label = SIG_GENES),
      #                 box.padding = 0.5,
      #                 point.padding = 0.3,
      #                 segment.color = "grey50",
      #                 show.legend = FALSE,
      #                 na.rm = TRUE) +
      # All labels should be to the left of 3.
      
      geom_text_repel(
        aes(label = SIG_GENES), 
        # xlim  = x_limits, #To limit positioning of labels
        size=10.0,
        box.padding = 0.5, 
        point.padding = 0.3,
        segment.color = "grey50",
        show.legend = FALSE,
        na.rm = TRUE,
        force = 10,
        max.overlaps = Inf,
        min.segment.length = 0,
        max.time = 15,
        max.iter = 100000
      ) +
      
      scale_alpha_manual(values = c(0.5,1,1,1),breaks = c(0,1,2,3,4)) +
      scale_size_manual(values = c(0,1,2,3,4),breaks = c(0,1,2,3,4)) +
      #geom_text_repel(aes(x = G1_BE_VS_SQ_Like_Nuclei_DIFF, y = PVAL,label = SIG_GENES)) +
      # ggtitle(paste(DE_Analysis,"\n",TEST,"\nLog 2 FoldChange",sep = ""))+
      ggtitle(paste("Significant Differentially Expressed Genes",sep = ""))+ #[INPUT_NEEDED]
      
      # xlab(paste(DE_Analysis,"\n Log2 Fold Change",sep="")) +
      xlab(paste("Log2 Fold Change",sep="")) +
      ylab(paste("-Log10(","P-Value)",sep="")) +
      ylim(c(0,(ceiling(max(PLOT_DATA[,"PVAL"],na.rm = T))+1)))+
      coord_cartesian(xlim = c(Min_Th,Max_Th),ylim = c(0,(ceiling(max(PLOT_DATA[,"PVAL"],na.rm = T))+0.5)))+
      theme_bw() + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
      theme(plot.title = element_text(hjust = 0.5))+
      theme(text = element_text(face = "bold",size = 10),
            axis.text.x=element_text(angle=0, hjust=1, vjust=0,size=20), #From10
            axis.title = element_text(size=24,face="bold"), #From12
            axis.text.y = element_text(size = 20, angle = 0, hjust = 0, vjust = 0, face = "bold"), #From16
            axis.title.y.right = element_text(size = 12),
            legend.text=element_text(size=12),
            legend.title=element_text(size=12),
            axis.line = element_line(size=2),
            legend.position = "bottom")+
      theme(plot.title = element_text(size=28))+ #From 8
      guides(fill=guide_legend(nrow=2, byrow=TRUE))+
      theme(legend.position = "none")
    
    #### Mean Plot
    MEAN_COLS <- c(paste(GROUP1_NAME,"_Mean",sep=""),paste(GROUP2_NAME,"_Mean",sep=""),"MEANVAL",MEANEXP_COL)
    #######################################################################
    MEANEXP_MIN <- floor(min(TEST_DATA[,MEAN_COLS])) 
    MEANEXP_MAX <- ceiling(max(TEST_DATA[,MEAN_COLS]))
    
    mean_col = 1
    Mean_Plot_List <- NULL
    for(mean_col in 1:length(MEAN_COLS)){
      
      MEAN_COL = MEAN_COLS[mean_col] 
      print(MEAN_COL)
      
      PLOT_DATA <- TEST_DATA[,c("Gene",LOG2FC_COL,MEAN_COL,PVAL_COL,PVALADJ_COL)]
      colnames(PLOT_DATA) <- c("GENE","LOG2FC","Average_Expression","PVAL","PVALADJ")
      #PLOT_DATA[,"PVAL"] <- ifelse(PLOT_DATA[,"PVAL"]<0.000000001,0.000000001,PLOT_DATA[,"PVAL"])
      PLOT_DATA$THRESHOLD <- ifelse((PLOT_DATA[,"PVAL"]<0.05),1,0)
      PLOT_DATA$THRESHOLD <- ifelse((PLOT_DATA[,"PVAL"]<0.05 & abs(PLOT_DATA[,"LOG2FC"])>=as.numeric(log(2,2))),2,PLOT_DATA[,"THRESHOLD"])
      PLOT_DATA$THRESHOLD <- ifelse((PLOT_DATA[,"PVALADJ"]<0.05 & abs(PLOT_DATA[,"LOG2FC"])>=as.numeric(log(2,2))),3,PLOT_DATA[,"THRESHOLD"])
      PLOT_DATA[,"PVAL"] <- -log10(PLOT_DATA[,"PVAL"])
      
      PVAL_th <- -log10(0.05)
      Max_Th <- ceiling(as.numeric(max(PLOT_DATA[,"LOG2FC"])))
      Min_Th <- floor(as.numeric(min(PLOT_DATA[,"LOG2FC"])))
      PLOT_DATA <- PLOT_DATA[order(PLOT_DATA[,"LOG2FC"],decreasing = T),]
      rownames(PLOT_DATA) <- NULL
      PLOT_DATA[,"Rank"] <- as.numeric(rownames(PLOT_DATA))
      PLOT_DATA[,"SIG_GENES"] <- NA
      PLOT_DATA[,"SIG_GENES"] <- ifelse(((PLOT_DATA[,"THRESHOLD"]>=2) & (PLOT_DATA[,"Rank"]<=25)),PLOT_DATA[,"GENE"],NA)
      PLOT_DATA <- PLOT_DATA[order(PLOT_DATA[,"LOG2FC"],decreasing = F),]
      rownames(PLOT_DATA) <- NULL
      PLOT_DATA[,"Rank"] <- as.numeric(rownames(PLOT_DATA))
      PLOT_DATA[,"SIG_GENES"] <- ifelse(((PLOT_DATA[,"THRESHOLD"]>=2) & (PLOT_DATA[,"Rank"]<=25)),PLOT_DATA[,"GENE"],PLOT_DATA[,"SIG_GENES"])
      
      Mean_Plot <- NULL
      Mean_Plot <- ggplot(PLOT_DATA,aes(x = Average_Expression, y = LOG2FC,label=SIG_GENES, size = factor(THRESHOLD)))+
        geom_label_repel(aes(label=SIG_GENES),box.padding= 0.3,point.padding = 0.5,segment.color="grey80",na.rm =T,colour = "black",fill="white",size = 1.5,force = 4,direction = "both",max.overlaps = Inf,nudge_x = 0,nudge_y = 0)+
        geom_point(aes(x = Average_Expression, y = LOG2FC,size = factor(THRESHOLD),colour = factor(THRESHOLD),alpha=factor(THRESHOLD))) +
        geom_point(data=PLOT_DATA[!is.na(PLOT_DATA[,"SIG_GENES"]),],aes(size = factor(THRESHOLD)),alpha = 1,colour="black",shape=1,pch = 21,stroke = 1,show.legend = F)+
        geom_vline(xintercept = 0,color = "black", linetype='dashed',color = "#4268F4")+
        # geom_vline(xintercept = as.numeric(log(2,2)),color = "black", linetype='dashed',color = "#4268F4",alpha = 0.5)+
        # geom_vline(xintercept = -as.numeric(log(2,2)),color = "black", linetype='dashed',color = "#4268F4",alpha = 0.5)+
        # geom_vline(xintercept = as.numeric(log(2,2)),color = "black", linetype='dashed',color = "#4268F4",alpha = 0.5)+
        # geom_vline(xintercept = -as.numeric(log(2,2)),color = "black", linetype='dashed',color = "#4268F4",alpha = 0.5)+
        # geom_hline(yintercept = PVAL_th, linetype='dashed',color = "#4268F4")+
        geom_hline(yintercept = log(2,2), linetype='dashed',color = "#4268F4")+
        geom_hline(yintercept = -log(2,2), linetype='dashed',color = "#4268F4")+
        geom_hline(yintercept = 0,color = "black")+
        geom_text(y=L2FC_MIN+1, x=MEANEXP_MAX, label=gsub("_","\n",paste(GROUP2_NAME,"_(N:",length(GROUP2_SAMPLES),")",sep="")),show.legend = F,size=5,color = "#4268F4") +
        geom_text(y=L2FC_MAX-1, x=MEANEXP_MAX, label=gsub("_","\n",paste(GROUP1_NAME,"_(N:",length(GROUP1_SAMPLES),")",sep="")),show.legend = F,size=5,color = "#4268F4") +
        # 
        scale_color_manual(values = c("grey50","#FF3158","#4268F4","#42B858"),breaks = c(0,1,2,3,4)) +
        scale_alpha_manual(values = c(0.5,1,1,1),breaks = c(0,1,2,3,4)) +
        scale_size_manual(values = c(0,1,2,3,4),breaks = c(0,1,2,3,4)) +
        xlim(c(MEANEXP_MIN,MEANEXP_MAX))+
        ylim(c(L2FC_MIN,L2FC_MAX))+
        #geom_text_repel(aes(x = G1_BE_VS_SQ_Like_Nuclei_DIFF, y = PVAL,label = SIG_GENES)) +
        ggtitle(paste(DE_Analysis,"\n",TEST,"\n",MEAN_COL," Vs Log2Fold Change",sep = ""))+
        ylab(paste(DE_Analysis,"\n Log2 Fold Change",sep="")) +
        xlab(paste(MEAN_COL,sep="")) +
        #ylim(c(0,(ceiling(max(PLOT_DATA[,"PVAL"],na.rm = T))+5)))+
        #coord_cartesian(xlim = c(Min_Th,Max_Th),ylim = c(0,(ceiling(max(PLOT_DATA[,"PVAL"],na.rm = T))+5)))+
        theme_bw() + 
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
        theme(plot.title = element_text(hjust = 0.5))+
        theme(text = element_text(face = "bold",size = 10),
              axis.text.x=element_text(angle=0, hjust=1, vjust=0,size=10),
              axis.title = element_text(size=12,face="bold"),
              axis.text.y = element_text(size = 16, angle = 0, hjust = 0, vjust = 0, face = "bold"),
              axis.title.y.right = element_text(size = 12),
              legend.text=element_text(size=12),
              legend.title=element_text(size=12),
              axis.line = element_line(size=2),
              legend.position = "bottom")+
        theme(plot.title = element_text(size=8))+
        guides(fill=guide_legend(nrow=2, byrow=TRUE))+
        theme(legend.position = "none")
      
      Mean_Plot_List[[mean_col]] <- Mean_Plot
      
    }
    
    # Combining the plots using patchwork
    # combined_plot <- Volcano_Plot + 
    #   (Mean_Plot_List[[1]] + Mean_Plot_List[[2]] + 
    #      Mean_Plot_List[[3]] + Mean_Plot_List[[4]]) + 
    #   plot_layout(ncol = 2, widths = c(0.5, 0.5), heights = c(1, 1))
    #   
    # combined_plot <- (
    #   (Volcano_Plot / (Mean_Plot_List[[1]] + Mean_Plot_List[[2]]))/
    #     (Mean_Plot_List[[3]] + Mean_Plot_List[[4]])
    #   ) +
    #     plot_layout(ncol = 1, guides = 'collect')
    
    combined_plot <- Volcano_Plot
    graphics.off()
    pdf(paste(output_dir, "/", GROUP1_NAME,"_VS_",GROUP2_NAME,"_DEG_",TEST,"_DEG_ANALYSIS_PVAL.pdf",sep=""),width = 30,height = 20) #From 30x15
    plot(combined_plot)
    graphics.off()
    
  }
}
