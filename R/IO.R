#' Print analyzed matrix
#' 
#' Writes a tab separated version of the analyzed OrthoMCL data with or without the joined representative sequences
#' @param mtrx Matrix derived from analyze_OrthoMCL
#' @param filename File name to save final output
#' @return The path to the written file
#' @examples
#' path <- write_mcl(mcl_mtrx, "matrix.tsv")
#' #mcl_mtrx previously derived from analyze_OrthoMCL() or join_repset()
#' @export
write_mcl <- function(mtrx, filename) {
  write.tsv(mtrx, filename)
  
  return(filename)
}

parse_haplo_old <-function(file) {
  
  table <- as.matrix(read.table(file, header=F, sep="|")) #now it will separate into two listsls
  
  haplo <- table[,1]
  haplo_list <- strsplit(haplo, split = " ")
  haplo_matrix <- matrix(unlist(haplo_list), ncol=length(haplo_list[[1]]), byrow=T)
  haplo_matrix <- t(haplo_matrix)
  
  return(haplo_matrix)
}

parse_names_old <-function(file) {
  
  table <- as.matrix(read.table(file, header=T, sep="|"))
  
  names <- table[,2]
  names_list <- strsplit(names, split = ",")
  
  return(names_list)
}

parse_orthologGroups <- function(file, format="ortho") {
  print(file)
  
  cat("reading in OrthoMCL data...")
  if(format == "ortho") {
    
    mcl_data <- parse_MCLdata(file)
    
  } else if(format == "groups") {
    
    mcl_data <- as.matrix(read.delim(file, header=F, sep=" ", row.names=1))
    rownames(mcl_data) <- sub(":", "", rownames(mcl_data))
    mcl_data <- as.matrix(t(mcl_data))
    
  } else {
    cat("format declaration invalid")
  }
  cat("done\n")
  
  #Create static 4-letter code list for grouping presence-absence
  cat("creating reference taxon codes...")
  code_list <- c()
  for(i in 1:dim(mcl_data)[2]) {
    for(j in 1:length(mcl_data[,i])){
      code <- strsplit(mcl_data[j,i], split="\\|")[[1]][1]
      if(!(code %in% code_list) & !is.na(code)) {
        code_list[length(code_list)+1] <- code 
      }
    }
  }
  cat("done\n")
  
  #Define group definitions from code_list and group similar groupings together
  cat("grouping protein sequences (this may take a while)...")
  haplo_test <- matrix(, ncol=0, nrow=length(code_list))
  haplo_names <- list()
  for(i in 1:dim(mcl_data)[2]) {
    
    code_vect <- strsplit(mcl_data[,i], split="\\|")
    code_mtrx <- matrix(unlist(code_vect), ncol=2, byrow=T)
    all_codes <- unlist(code_mtrx[,1])
    
    haploTF <- code_list%in%all_codes
    
    contains <- 0
    for(j in 1:dim(haplo_test)[2]) {
      if(dim(haplo_test)[2] == 0)
        break
      if(sum(haploTF == haplo_test[, j]) == length(code_list)) {
        contains <- j
        break
      }
    }
    if(contains > 0) {
      haplo_names[[contains]] <- append(haplo_names[[contains]], colnames(mcl_data)[i])
    } else {
      haplo_test <- cbind(haplo_test, haploTF)
      haplo_names[length(haplo_names)+1] <- list(colnames(mcl_data)[i])
    }
  }
  cat("done\n")
  
  haplo_test <- ifelse((haplo_test == T), 1, 0)
  
  #populate final ouput matrix
  haplo_mtrx <- matrix(, ncol=dim(haplo_test)[2]+1, nrow=length(code_list))
  haplo_mtrx[,1] <- code_list
  for(i in 1:dim(haplo_test)[2]) {
    haplo_mtrx[,i+1] <- haplo_test[,i]
  }
 
  return(list("matrix" = haplo_mtrx, "names" = haplo_names))
}

parse_MCLdata <- function(file) {
  
  table <- as.matrix(read.table(file, header=F, sep="\t"))
  table <- table[order(table[,2]),]
  
  group <- ".:NA:."
  group_vals <- c()
  mcl_data <- matrix(, ncol=0, nrow=dim(table)[1])
  maxsize <- dim(table)[1]
  
  for(i in 1:dim(table)[1]) {
    
    if(group != table[i,2]) { 
      #save old grouping
      if(group != ".:NA:.") {
        length(group_vals) <- dim(table)[1] - 1
        mcl_data <- cbind(mcl_data, group_vals)
        colnames(mcl_data)[dim(mcl_data)[2]] <- group
      }
      
      #start new grouping
      group <- table[i,2]
      group_vals <- c(table[i,1])
    }
    else {
      #continue grouping
      group_vals[length(group_vals)+1] <- table[i,1]
    }
  }
  length(group_vals) <- dim(table)[1] - 1
  mcl_data <- cbind(mcl_data, group_vals)
  colnames(mcl_data)[dim(mcl_data)[2]] <- group
  
  #set NA values to empty string
  sums <- colSums(!is.na(mcl_data))
  mcl_data <- head(mcl_data, n=max(sums))
  mcl_data[is.na(mcl_data)] <- ""
  
  #cleaning steps
  delcols <- c()
  for(i in 1:dim(mcl_data)[2]) {
    rows <- sum(mcl_data[, i] != "")
    
    #remove no group values
    if(colnames(mcl_data)[i] == "NO_GROUP") {
      delcols[length(delcols)+1] <- i
    }
    
    #remove singleton groups
    if(rows <= 1) {
      delcols[length(delcols)+1] <- i
    }
    
    #remove species-only groups
    code <- strsplit(mcl_data[1,i], split="\\|")[[1]][1]
    addcol <- T
    for(j in 1:rows) {
      test <- strsplit(mcl_data[j,i], split="\\|")[[1]][1]
      if(code != test) {
        addcol <- F
        break
      }
    }
    if(addcol & !(i %in% delcols)) {
      delcols[length(delcols)+1] <- i
    }
  }
  mcl_data <- mcl_data[,-delcols]

  return(mcl_data)
}

merge_tx <-function(file, haplo_matrix) {
  
  tx <- read.table(file, header=T)
  
  haplo_tx <- merge(tx, haplo_matrix, by.y="V1", by.x="Treatment", all=F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx$tgl_fly>0))   
  
  return(haplo_tx)
}

create_haplo_matrix <-function(haplo_file, var_file) {
  
  haplo_matrix <- parse_haplo(haplo_file)
  haplo_tx <- merge_tx(var_file, haplo_matrix)
  
  return(haplo_tx)
}