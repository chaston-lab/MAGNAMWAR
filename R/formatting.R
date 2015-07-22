#' Format all raw GenBank fastas to single OrthoMCL compatible fasta file
#' 
#' Creates the composite fasta file for use in running OrthoMCL and/or submitting to www.orthomcl.org
#' @param fa_dir Path to the directory where all raw GenBank files are stored. Note, all file names must be changed to a 4-letter code representing each species and have ".fasta" file descriptor
#' @param genbnk_id The index of the sequence ID in the GenBank pipe-separated annotation line (default: 4)
#' @return The path to the final OrthoMCL compatible fasta file
#' @examples
#' formatted_file <- formatMCLfastas("sample_data/")
#' @export
format_MCLfastas <- function(fa_dir, genbnk_id=4) {
  
  library("seqinr")
  
  filename <- "MCLformatted_all.fasta"
  outfile <- paste(c(fa_dir, filename), collapse="")
  
  files <- dir(fa_dir, pattern = ".fasta")
  if(filename %in% files) {
    file.remove(outfile)
  }
  
  for(i in 1:length(files)) {
    
    if(files[i] != filename & files[i] != "repseq.fasta") {
      
      abs_path <- paste(c(fa_dir, files[i]), collapse="")
      cat(abs_path)
      cat("\n")
      
      id <- strsplit(files[i], split="\\.")
      
      reps_fa <- read.fasta(file=abs_path, as.string=T, forceDNAtolower=F, seqonly=F, strip.desc=T)
      
      for(j in 1:length(reps_fa)) {
        
        info <- strsplit(getAnnot(reps_fa[[j]]), split="\\|")
        
        mcl_info <- paste(c(id[[1]][1], info[[1]][genbnk_id]), collapse="|")
        
        write.fasta(reps_fa[[j]][1], mcl_info, outfile, open="a")
      }
    }
  }
  
  return(outfile)
}

#' Join Representative Sequences
#' 
#' Joins the OrthoMCL output matrix to representative sequences
#' @param reps_file A fasta file containing one representative sequence for every COG in mcl_mtrx
#' @param mcl_mtrx OrthoMCL output matrix from analyze_OrthoMCL()
#' @return Returns the original OrthoMCL output matrix with additional columns: representative sequence taxon, representative sequence id, representative sequence annotation, representative sequence 
#' @examples
#' joined_mtrx <- join_repset(sample_data/repseq.fasta", mcl_mtrx)
#' #mcl_mtrx previously derived from analyze_OrthoMCL()
#' @export
join_repset <- function(reps_file, mcl_mtrx) {
  
  library("seqinr")
  
  reps_fa <- read.fasta(file=reps_file, as.string=T, forceDNAtolower=F, seqonly=F, strip.desc=T)
  fa_mtrx <- matrix(, nrow=length(reps_fa), ncol=5)
  colnames(fa_mtrx) <- c("pdg", "rep_taxon", "rep_id", "rep_annot", "rep_seq")
  
  for(i in 1:length(reps_fa)) {
    info <- strsplit(getAnnot(reps_fa[[i]]), split = "\t")
    fa_mtrx[i,] <- c(info[[1]][1], info[[1]][2], info[[1]][3], info[[1]][4], reps_fa[[i]][1])
  }
  
  mcl_reps <- merge(mcl_mtrx, fa_mtrx, by="pdg", all=F)
  
  return(mcl_reps)
}

#' Pick Representative Sequences
#' 
#' Randomly picks a representative sequence from GenBank fasta files for every OrthoMCL COG
#' @param mcl_file OrthoMCL output file
#' @param fa_dir Path to the directory where all raw GenBank files are stored. Note, all file names must be changed to a 4-letter code representing each species and have ".fasta" file descriptor
#' @param del_mid Print intermediate fully concatinated fasta (default: FALSE)
#' @param format Select format of OrthoMCL output file (default: "ortho" for version 2.0 orthologGroups, option: "groups" for version 1.0 output)
#' @return The path to the representative sequence file in fasta format
#' @examples
#' repseqfile <- pick_repseq("sample_data/orthologGroups", "sample_data/")
#' @export
pick_repseq <- function(mcl_file, fa_dir, del_mid=T, format="ortho") {
  
  filename <- "repseq.fasta"
  repseqfile <- paste(c(fa_dir, filename), collapse="")
  
  files <- dir(fa_dir, pattern = ".fasta")
  if(filename %in% files) {
    file.remove(repseqfile)
  }
  
  cat("reading in files:\n")
  if(format == "ortho") {
    mcl_data <- parse_MCLdata(mcl_file)
  } else if(format == "groups"){
    mcl_data <- as.matrix(read.delim(mcl_file, header=F, sep=" ", row.names=1))
    rownames(mcl_data) <- sub(":", "", rownames(mcl_data))
    mcl_data <- as.matrix(t(mcl_data))
  } else {
    cat("format declaration invalid")
  }
  #View(mcl_data)
  
  annot_data <- format_MCLreps(fa_dir, 4, 5, del_mid)
  
  cat("\npicking and writing representative sequence for PDG:\n")
  for(i in 1:dim(mcl_data)[2]) {
    
    total_seq <- sum(mcl_data[, i] != "")
    rndm_seq <- sample(1:total_seq, 1)
    rep_info <- strsplit(mcl_data[rndm_seq, i], split="\\|")
    
    for(j in 1:length(annot_data)) {
      
      annot_info <- strsplit(getAnnot(annot_data[[j]]), split="\\|")
      
      if(rep_info[[1]][2] == annot_info[[1]][2]) { 
        
        annotations <- paste(c(colnames(mcl_data)[i], rep_info[[1]][1], rep_info[[1]][2], " N/A"), collapse="\t")
        try(annotations <- paste(c(colnames(mcl_data)[i], rep_info[[1]][1], rep_info[[1]][2], annot_info[[1]][3]), collapse="\t"),T)
        
        cat(colnames(mcl_data)[i])
        cat("\n")
        
        write.fasta(annot_data[[j]][1], annotations, repseqfile, open="a")
        break
      }
    }
  } 
  
  return(repseqfile)
}

format_MCLreps <- function(fa_dir, genbnk_id, annot_id, del_mid) {
  
  library(seqinr)
  
  files <- dir(fa_dir)
  if("temp.fa" %in% files) {
    file.remove("temp.fa")
  }
  files <- dir(fa_dir, pattern = ".fasta")
  
  outfile <- paste(c(fa_dir, "temp.fa"), collapse="")
  for(i in 1:length(files)) {
    
    if(files[i] != "MCLformatted_all.fasta" & files[i] != "repseq.fasta") {
      
      abs_path <- paste(c(fa_dir, files[i]), collapse="")
      cat(abs_path)
      cat("\n")
      
      id <- strsplit(files[i], split="\\.")
      
      reps_fa <- read.fasta(file=abs_path, as.string=T, forceDNAtolower=F, seqonly=F, strip.desc=T)
      
      for(j in 1:length(reps_fa)) {
        
        info <- strsplit(getAnnot(reps_fa[[j]]), split="\\|")
        
        mcl_info <- paste(c(id[[1]][1], info[[1]][genbnk_id], info[[1]][annot_id]), collapse="|")
        
        write.fasta(reps_fa[[j]][1], mcl_info, outfile, open="a")
      }
    }
  }
  
  all_fastas <- read.fasta(file=outfile, as.string=T, forceDNAtolower=F, seqonly=F, strip.desc=T)
  
  if(del_mid)
    file.remove(outfile)
  
  return(all_fastas)
}