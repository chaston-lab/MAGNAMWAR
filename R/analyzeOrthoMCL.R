#' Main OrthoMCL Analysis
#' 
#' Main function for analyzing the statistical association of PDG (phylogenetic distribution group) presence with phenotype data
#' @param mcl_file Path to OrthoMCL output file containing COG (cluster of orthologous groups) groupings
#' @param var_file Path to tab separated phenotype data
#' @param species_name Column name in var_file containing 4-letter species designations
#' @param model Linear Model with gene presence as fixed effect (lm), Wilcox Test with gene presence as fixed effect (wx), Linear Mixed Effect models with gene presence as fixed effect and additional variables specified as: one random effect (lmeR1); two independent random effects (lmeR2ind); two random effects with rndm2 nested in rndm1 (lmeR2nest); or two independent random effects with one additional fixed effect (lmeF2)
#' @param resp Column name in var_file containing response variable
#' @param fix2 Column name in var_file containing second fixed effect
#' @param rndm1 Column name in var_file containing first random variable
#' @param rndm2 Column name in var_file containing second random variable 
#' @param format Select format of OrthoMCL output file (default: "ortho" for version 2.0 orthologGroups, option: "groups" for version 1.0 output)
#' @return A matrix with the following columns: PDG, p-values, Bonferroni corrected p-values, mean phenotype of PDG-containing taxa, mean pheotype of PDG-lacking taxa, taxa that contain PDG, taxa that lack PDG
#' @references Some sore of reference
#' @examples 
#' #Linear Model
#' mcl_mtrx <- analyze_OrthoMCL("sample_data/orthologGroups", "sample_data/pheno_data.tsv", "Treatment", "lm", resp="RespVar")
#' 
#' #Wilcox Test
#' mcl_mtrx <- analyze_OrthoMCL("sample_data/orthologGroups", "sample_data/pheno_data.tsv", "Treatment", "wx", resp="RespVar")
#' 
#' #Linear Mixed Effect with one random effect
#' mcl_mtrx <- analyze_OrthoMCL("sample_data/orthologGroups", "sample_data/pheno_data.tsv", "Treatment", "lmeR1", resp="RespVar", rndm1="Experiment")
#' 
#' #Linear Mixed Effect with two independent random effects
#' mcl_mtrx <- analyze_OrthoMCL("sample_data/orthologGroups", "sample_data/pheno_data.tsv", "Treatment", "lmeR2ind", resp="RespVar", rndm1="Experiment", rndm2="Vial")
#' 
#' #Linear Mixed Effect with rndm2 nested in rndm1
#' mcl_mtrx <- analyze_OrthoMCL("sample_data/orthologGroups", "sample_data/pheno_data.tsv", "Treatment", "lmeR2nest", resp="RespVar", rndm1="Experiment", rndm2="Vial")
#' 
#' #Linear Mixed Effect with two independent random effects and one additional fixed effect
#' mcl_mtrx <- analyze_OrthoMCL("sample_data/orthologGroups", "sample_data/pheno_data.tsv", "Treatment", "lmeF2", resp="RespVar", fix2="Treatment", rndm1="Experiment", rndm2="Vial")
#' 
#' @export
analyze_OrthoMCL <- function(mcl_file, var_file, species_name, model, resp=NULL, fix2=NULL, rndm1=NULL, rndm2=NULL, format="ortho") {
  
  cat(" _____ _____ _____ _____ _____ _____ _ _ _ _____ _____ \n")
  cat("|     |  _  |   __|   | |  _  |     | | | |  _  | __  |\n")
  cat("| | | |     |  |  | | | |     | | | | | | |     |    -|\n")
  cat("|_|_|_|__|__|_____|_|___|__|__|_|_|_|_____|__|__|__|__|\n")
  
  if(format == "ortho") {
    haplo_data <- parse_orthologGroups(mcl_file)
  } else if(format == "groups") {
    haplo_data <- parse_orthologGroups(mcl_file, format="groups")
  } else {
   cat("format declaration invalid")
  }
  
  haplo_mtrx <- haplo_data$matrix
  haplo_names <- haplo_data$names
  
  if(model == "lmeF2")
    mtrx <- analyze.ffrr(haplo_mtrx, haplo_names, var_file, species_name, resp, fix2, rndm1, rndm2)
  else if(model == "lmeR2nest")
    mtrx <- analyze.frr.div(haplo_mtrx, haplo_names, var_file, species_name, resp, rndm1, rndm2)
  else if(model == "lmeR2ind")
    mtrx <- analyze.frr.plus(haplo_mtrx, haplo_names, var_file, species_name, resp, rndm1, rndm2)
  else if(model == "lmeR1")
    mtrx <- analyze.fr(haplo_mtrx, haplo_names, var_file, species_name, resp, rndm1)
  else if(model == "lm")
    mtrx <- analyze.f(haplo_mtrx, haplo_names, var_file, species_name, resp)
  else if(model == "wx")
    mtrx <- analyze.wilcox(haplo_mtrx, haplo_names, var_file, species_name, resp)
  else 
    cat("Error: Could not find a correct match for your model declaration\n")
  
  return(mtrx)
}

analyze.ffrr <- function(haplo_mtrx, haplo_names, var_file, species_name, resp_var, fix_var2, rndm1_eff, rndm2_eff) {
  
  library(lme4)
  library(multcomp)
  
  cat("Importing Data\n")
  
  tx <- read.table(var_file, header=T)
  
  cat("Merging Files\n")
  haplo_tx <- merge(tx, haplo_mtrx, by.y="V1", by.x=species_name, all=F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(haplo_mtrx)[2] - 1 #number of phylogenetic distribution groups
  
  count <- 0
  for(i in haplo_names) {
    for(j in haplo_names[i])
      count <- count + 1
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow=count, ncol=9)
  
  sub <- matrix(, nrow=0, ncol=0)
  sub$resp_var <- haplo_tx[,which(colnames(haplo_tx)==resp_var)]
  sub$fix2 <- haplo_tx[,which(colnames(haplo_tx)==fix_var2)]
  sub$rndm1 <- haplo_tx[,which(colnames(haplo_tx)==rndm1_eff)]
  sub$rndm2 <- haplo_tx[,which(colnames(haplo_tx)==rndm2_eff)]
  
  sub$rndm1 <- factor(sub$rndm1)
  sub$rndm2 <- factor(sub$rndm2)
  
  cat("\tPercent Complete:\n\t")
  out_cnt <- 1
  for(i in 2:(num_pdg + 1)) { 
    
    #getting column names for expiramental fixed effect
    name <- paste("V", i, sep="")
    sub$fix <- haplo_tx[, which(colnames(haplo_tx)==name)]
    sub$species <- haplo_tx[, which(colnames(haplo_tx)==species_name)]
    #View(sub)
    
    #fitting the linear mixed model
    lmm <- try(lmer(resp_var~fix + fix2 + (1|rndm1) + (1|rndm2), data=sub),T)
    l1.1 <- try(glht(lmm, mcp(fix="Tukey")),T)
    l1.2 <- try(summary(l1.1),T)
    p_val1 <- try(l1.2$test$pvalues[1],T)
    
    l2.1 <- try(glht(lmm, mcp(fix2="Tukey")),T)
    l2.2 <- try(summary(l2.1),T)
    p_val2 <- "variable not a factor"
    try(p_val2 <- l2.2$test$pvalues[1],T)
    
    p_val2_corrected <- "ibid"
    try(p_val2_corrected <- p_val2*num_pdg,T)
    
    #calculating meta-data
    mean_calc <- aggregate(resp_var~fix, sub, mean) #matrix with mean_contain and mean_missing
    mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
    mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing=F),]
    mean_contain <- mean_calc2[2,2]
    mean_missing <- mean_calc2[1,2]
    
    taxa <- aggregate(as.numeric(as.character(fix))~species, sub, mean)
    colnames(taxa) <- c("taxa_name", "presence") 
    
    taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence==1))
    taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse="|")
    
    taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence==0))
    taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse="|")
    
    #adding data to output matrix
    for( j in haplo_names[[i-1]]) {
      
      try(output[out_cnt,] <- c(j, p_val1, (p_val1*num_pdg), p_val2, p_val2_corrected, mean_contain, mean_missing, taxa_contain, taxa_missing),T)
      out_cnt <- out_cnt + 1
    }
    
    #printing progress
    if(((i-1) %% 100) == 0) 
      cat(paste(round(((i-1)/num_pdg*100), digits=2), "%__", sep=""))
  }
  cat("Finished!\n")
  
  cat("Cleaning Final Output")
  output_clean <- output[rowSums(is.na(output))!=9, ] 
  colnames(output_clean) <- c("pdg", "p-val1", "corrected_p-val1", "p-val2", "corrected_p-val2", "mean_pdgContain", "mean_pdgLack", "taxa_contain", "taxa_miss")
  
  return(output_clean)
}

analyze.frr.plus <- function(haplo_mtrx, haplo_names, var_file, species_name, resp_var, rndm1_eff, rndm2_eff) {
  
  library(lme4)
  library(multcomp)
  
  cat("Importing Data\n")
  
  tx <- read.table(var_file, header=T)
  
  cat("Merging Files\n")
  haplo_tx <- merge(tx, haplo_mtrx, by.y="V1", by.x=species_name, all=F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(haplo_mtrx)[2] - 1 #number of phylogenetic distribution groups
  
  count <- 0
  for(i in haplo_names) {
    for(j in haplo_names[i])
      count <- count + 1
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow=count, ncol=7)
  
  sub <- matrix(, nrow=0, ncol=0)
  sub$resp_var <- haplo_tx[,which(colnames(haplo_tx)==resp_var)]
  sub$rndm1 <- haplo_tx[,which(colnames(haplo_tx)==rndm1_eff)]
  sub$rndm2 <- haplo_tx[,which(colnames(haplo_tx)==rndm2_eff)]
  
  sub$rndm1 <- factor(sub$rndm1)
  sub$rndm2 <- factor(sub$rndm2)
  
  cat("\tPercent Complete:\n\t")
  out_cnt <- 1
  for(i in 2:(num_pdg + 1)) { 
    
    #getting column names for expiramental fixed effect
    name <- paste("V", i, sep="")
    sub$fix <- haplo_tx[, which(colnames(haplo_tx)==name)]
    sub$species <- haplo_tx[, which(colnames(haplo_tx)==species_name)]
    #View(sub)
    
    #fitting the linear mixed model
    lmm <- try(lmer(resp_var~fix + (1|rndm1) + (1|rndm2), data=sub),T)
    l1 <- try(glht(lmm, mcp(fix="Tukey")),T)
    l2 <- try(summary(l1),T)
    p_val <- try(l2$test$pvalues[1],T)
    
    #calculating meta-data
    mean_calc <- aggregate(resp_var~fix, sub, mean) #matrix with mean_contain and mean_missing
    mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
    mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing=F),]
    mean_contain <- mean_calc2[2,2]
    mean_missing <- mean_calc2[1,2]
    
    taxa <- aggregate(as.numeric(as.character(fix))~species, sub, mean)
    colnames(taxa) <- c("taxa_name", "presence") 
    
    taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence==1))
    taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse="|")
    
    taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence==0))
    taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse="|")
    
    #adding data to output matrix
    for( j in haplo_names[[i-1]]) {
      
      try(output[out_cnt,] <- c(j, p_val, (p_val*num_pdg), mean_contain, mean_missing, taxa_contain, taxa_missing),T)
      out_cnt <- out_cnt + 1
    }
    
    #printing progress
    if(((i-1) %% 100) == 0) 
      cat(paste(round(((i-1)/num_pdg*100), digits=2), "%__", sep=""))
  }
  cat("Finished!\n")
  
  cat("Cleaning Final Output")
  output_clean <- output[rowSums(is.na(output))!=7, ] 
  colnames(output_clean) <- c("pdg", "p-val", "corrected_p-val", "mean_pdgContain", "mean_pdgLack", "taxa_contain", "taxa_miss")
  
  return(output_clean)
}

analyze.frr.div <- function(haplo_mtrx, haplo_names, var_file, species_name, resp_var, rndm1_eff, rndm2_eff) {
  
  library(lme4)
  library(multcomp)
  
  cat("Importing Data\n")
  
  tx <- read.table(var_file, header=T)
  
  cat("Merging Files\n")
  haplo_tx <- merge(tx, haplo_mtrx, by.y="V1", by.x=species_name, all=F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(haplo_mtrx)[2] - 1 #number of phylogenetic distribution groups
  
  count <- 0
  for(i in haplo_names) {
    for(j in haplo_names[i])
      count <- count + 1
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow=count, ncol=7)
  
  sub <- matrix(, nrow=0, ncol=0)
  sub$resp_var <- haplo_tx[,which(colnames(haplo_tx)==resp_var)]
  sub$rndm1 <- haplo_tx[,which(colnames(haplo_tx)==rndm1_eff)]
  sub$rndm2 <- haplo_tx[,which(colnames(haplo_tx)==rndm2_eff)]
  
  sub$rndm1 <- factor(sub$rndm1)
  sub$rndm2 <- factor(sub$rndm2)
  
  cat("\tPercent Complete:\n\t")
  out_cnt <- 1
  for(i in 2:(num_pdg + 1)) { 
    
    #getting column names for expiramental fixed effect
    name <- paste("V", i, sep="")
    sub$fix <- haplo_tx[, which(colnames(haplo_tx)==name)]
    sub$species <- haplo_tx[, which(colnames(haplo_tx)==species_name)]
    #View(sub)
    
    #fitting the linear mixed model
    lmm <- try(lmer(resp_var~fix + (1|rndm1/rndm2), data=sub),T)
    l1 <- try(glht(lmm, mcp(fix="Tukey")),T)
    l2 <- try(summary(l1),T)
    p_val <- try(l2$test$pvalues[1],T)
    
    #calculating meta-data
    mean_calc <- aggregate(resp_var~fix, sub, mean) #matrix with mean_contain and mean_missing
    mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
    mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing=F),]
    mean_contain <- mean_calc2[2,2]
    mean_missing <- mean_calc2[1,2]
    
    taxa <- aggregate(as.numeric(as.character(fix))~species, sub, mean)
    colnames(taxa) <- c("taxa_name", "presence") 
    
    taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence==1))
    taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse="|")
    
    taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence==0))
    taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse="|")
    
    #adding data to output matrix
    for( j in haplo_names[[i-1]]) {
      
      try(output[out_cnt,] <- c(j, p_val, (p_val*num_pdg), mean_contain, mean_missing, taxa_contain, taxa_missing),T)
      out_cnt <- out_cnt + 1
    }
    
    #printing progress
    if(((i-1) %% 100) == 0) 
      cat(paste(round(((i-1)/num_pdg*100), digits=2), "%__", sep=""))
  }
  cat("Finished!\n")
  
  cat("Cleaning Final Output")
  output_clean <- output[rowSums(is.na(output))!=7, ] 
  colnames(output_clean) <- c("pdg", "p-val", "corrected_p-val", "mean_pdgContain", "mean_pdgLack", "taxa_contain", "taxa_miss")
  
  return(output_clean)
}

analyze.fr <- function(haplo_mtrx, haplo_names, var_file, species_name, resp_var, rndm_eff) {
  
  library(lme4)
  library(multcomp)
  
  cat("Importing Data\n")
  
  tx <- read.table(var_file, header=T)
  
  cat("Merging Files\n")
  haplo_tx <- merge(tx, haplo_mtrx, by.y="V1", by.x=species_name, all=F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(haplo_mtrx)[2] - 1 #number of phylogenetic distribution groups
  
  count <- 0
  for(i in haplo_names) {
    for(j in haplo_names[i])
      count <- count + 1
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow=count, ncol=7)
  
  sub <- matrix(, nrow=0, ncol=0)
  sub$resp_var <- haplo_tx[,which(colnames(haplo_tx)==resp_var)]
  sub$rndm <- haplo_tx[,which(colnames(haplo_tx)==rndm_eff)]
  
  sub$rndm <- factor(sub$rndm)
  
  cat("\tPercent Complete:\n\t")
  out_cnt <- 1
  for(i in 2:(num_pdg + 1)) { 
    
    #getting column names for expiramental fixed effect
    name <- paste("V", i, sep="")
    sub$fix <- haplo_tx[, which(colnames(haplo_tx)==name)]
    sub$species <- haplo_tx[, which(colnames(haplo_tx)==species_name)]
    #View(sub)
    
    #fitting the linear mixed model
    lmm <- try(lmer(resp_var~fix + (1|rndm), data=sub),T)
    l1 <- try(glht(lmm, mcp(fix="Tukey")),T)
    l2 <- try(summary(l1),T)
    p_val <- try(l2$test$pvalues[1],T)
    
    #calculating meta-data
    mean_calc <- aggregate(resp_var~fix, sub, mean) #matrix with mean_contain and mean_missing
    mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
    mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing=F),]
    mean_contain <- mean_calc2[2,2]
    mean_missing <- mean_calc2[1,2]
    
    taxa <- aggregate(as.numeric(as.character(fix))~species, sub, mean)
    colnames(taxa) <- c("taxa_name", "presence") 
    
    taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence==1))
    taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse="|")
    
    taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence==0))
    taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse="|")
    
    #adding data to output matrix
    for( j in haplo_names[[i-1]]) {
      
      try(output[out_cnt,] <- c(j, p_val, (p_val*num_pdg), mean_contain, mean_missing, taxa_contain, taxa_missing),T)
      out_cnt <- out_cnt + 1
    }
    
    #printing progress
    if(((i-1) %% 100) == 0) 
      cat(paste(round(((i-1)/num_pdg*100), digits=2), "%__", sep=""))
  }
  cat("Finished!\n")
  
  cat("Cleaning Final Output")
  output_clean <- output[rowSums(is.na(output))!=7, ] 
  colnames(output_clean) <- c("pdg", "p-val", "corrected_p-val", "mean_pdgContain", "mean_pdgLack", "taxa_contain", "taxa_miss")
  
  return(output_clean)
}

analyze.f <- function(haplo_mtrx, haplo_names, var_file, species_name, resp_var) {
  
  library(multcomp)
  
  cat("Importing Data\n")
  
  tx <- read.table(var_file, header=T)
  
  cat("Merging Files\n")
  haplo_tx <- merge(tx, haplo_mtrx, by.y="V1", by.x=species_name, all=F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(haplo_mtrx)[2] - 1 #number of phylogenetic distribution groups
  
  count <- 0
  for(i in haplo_names) {
    for(j in haplo_names[i])
      count <- count + 1
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow=count, ncol=7)
  
  sub <- matrix(, nrow=0, ncol=0)
  sub$resp_var <- haplo_tx[,which(colnames(haplo_tx)==resp_var)]
  
  cat("\tPercent Complete:\n\t")
  out_cnt <- 1
  for(i in 2:(num_pdg + 1)) { 
    
    #getting column names for expiramental fixed effect
    name <- paste("V", i, sep="")
    sub$fix <- haplo_tx[, which(colnames(haplo_tx)==name)]
    sub$species <- haplo_tx[, which(colnames(haplo_tx)==species_name)]
    #View(sub)
    
    #fitting the linear model
    lm <- try(lm(resp_var~fix, data=sub),T)
    l1 <- try(glht(lm, mcp(fix="Tukey")),T)
    l2 <- try(summary(l1),T)
    p_val <- try(l2$test$pvalues[1],T)
    
    #calculating meta-data
    mean_calc <- aggregate(resp_var~fix, sub, mean) #matrix with mean_contain and mean_missing
    mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
    mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing=F),]
    mean_contain <- mean_calc2[2,2]
    mean_missing <- mean_calc2[1,2]
    
    taxa <- aggregate(as.numeric(as.character(fix))~species, sub, mean)
    colnames(taxa) <- c("taxa_name", "presence") 
    
    taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence==1))
    taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse="|")
    
    taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence==0))
    taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse="|")
    
    #adding data to output matrix
    for( j in haplo_names[[i-1]]) {
      
      try(output[out_cnt,] <- c(j, p_val, (p_val*num_pdg), mean_contain, mean_missing, taxa_contain, taxa_missing),T)
      out_cnt <- out_cnt + 1
    }
    
    #printing progress
    if(((i-1) %% 100) == 0) 
      cat(paste(round(((i-1)/num_pdg*100), digits=2), "%__", sep=""))
  }
  cat("Finished!\n")
  
  cat("Cleaning Final Output")
  output_clean <- output[rowSums(is.na(output))!=7, ] 
  colnames(output_clean) <- c("pdg", "p-val", "corrected_p-val", "mean_pdgContain", "mean_pdgLack", "taxa_contain", "taxa_miss")
  
  return(output_clean)
}

analyze.wilcox <- function(haplo_mtrx, haplo_names, var_file, species_name, resp_var) {
  
  cat("Importing Data\n")
  
  tx <- read.table(var_file, header=T)
  
  cat("Merging Files\n")
  haplo_tx <- merge(tx, haplo_mtrx, by.y="V1", by.x=species_name, all=F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(haplo_mtrx)[2] - 1 #number of phylogenetic distribution groups
  
  count <- 0
  for(i in haplo_names) {
    for(j in haplo_names[i])
      count <- count + 1
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow=count, ncol=7)
  
  sub <- matrix(, nrow=0, ncol=0)
  sub$resp_var <- haplo_tx[,which(colnames(haplo_tx)==resp_var)]
  
  cat("\tPercent Complete:\n\t")
  out_cnt <- 1
  for(i in 2:(num_pdg + 1)) { 
    
    #getting column names for expiramental fixed effect
    name <- paste("V", i, sep="")
    sub$fix <- haplo_tx[, which(colnames(haplo_tx)==name)]
    sub$species <- haplo_tx[, which(colnames(haplo_tx)==species_name)]
    #View(sub)
    
    #fitting the linear model
    wix <- try(wilcox.test(resp_var~fix, data=sub),T)
    p_val <- try(wix$p.value,T)
    
    #calculating meta-data
    mean_calc <- aggregate(resp_var~fix, sub, mean) #matrix with mean_contain and mean_missing
    mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
    mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing=F),]
    mean_contain <- mean_calc2[2,2]
    mean_missing <- mean_calc2[1,2]
    
    taxa <- aggregate(as.numeric(as.character(fix))~species, sub, mean)
    colnames(taxa) <- c("taxa_name", "presence") 
    
    taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence==1))
    taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse="|")
    
    taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence==0))
    taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse="|")
    
    #adding data to output matrix
    for( j in haplo_names[[i-1]]) {
      
      try(output[out_cnt,] <- c(j, p_val, (p_val*num_pdg), mean_contain, mean_missing, taxa_contain, taxa_missing),T)
      out_cnt <- out_cnt + 1
    }
    
    #printing progress
    if(((i-1) %% 100) == 0) 
      cat(paste(round(((i-1)/num_pdg*100), digits=2), "%__", sep=""))
  }
  cat("Finished!\n")
  
  cat("Cleaning Final Output")
  output_clean <- output[rowSums(is.na(output))!=7, ] 
  colnames(output_clean) <- c("pdg", "p-val", "corrected_p-val", "mean_pdgContain", "mean_pdgLack", "taxa_contain", "taxa_miss")
  
  return(output_clean)
}














