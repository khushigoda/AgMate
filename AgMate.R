###########################################################################
###########################################################################
###                                                                     ###
### AgMate Code                                                         ###
### Written By: Khushi Goda, PhD Candidate, NCSU                        ###
### Tree Improvement Program`                                           ###
###########################################################################
###########################################################################

###########################################################################
### Library Packages  ###
###########################################################################

library(shinydashboard)
library(pedigree)
library(AGHmatrix)
library(plyr)
library(fastmatch)
library(ggplot2)
library(DEoptim)
library(splitstackshape)
library(shiny)
library(shinyFiles)
library(shinyjs)
library(shinyvalidate)
library(shinyalert)
library(tidyverse)

###########################################################################
### Global Variables  ###
###########################################################################

sorted_output_filename = "Pedigree.csv"
sorted_output_filename2 = "AMatrix.csv"
output_filename = "Mating_List.csv"
output_filename_optimized = "Mating_List_Optimized.csv"
output_logfile2 = "Log_Part2.txt"
output_logfile1 = "Log_Part1.txt"
next_round_file_name <- "Input_Next_Round.csv"
graph1_filename <- "Density.png"
graph2_filename <- "Index.png"
graph_stats_data <- data.frame(matrix(ncol = 2, nrow = 2))

glo_parents_list <- ""
glo_AMatrix <- ""
glo_ebvweight <- 0
glo_progenyinbreedingweight <- 0
glo_parentalcoaweight <- 0
glo_numberbabies <- 0
glo_amatrixcutoff <- 0
glo_iterations <- 0
glo_plot1 <- 0
glo_plot2 <- 0

###########################################################################
### Part 1 Functions  ###
###########################################################################

###########################################################################
### Function to check validity of files uploaded for Part 1 input       ###
###########################################################################

checkFiles_part1 <- function(selection_list, reference_list, x) {
  
  error1 <- "Candidate file needs to have the following columns at least: 'CandidateID', 'Parent1', 'Parent2', and 'Index'. Please upload an appropriate candidate file."
  error2 <- "Reference file needs to have the following columns at least: 'ID', 'Parent1', and 'Parent2'. Please upload an appropriate reference file."
  success <- "1"
  
  #case x = 1 is when only selection file is given and reference file = selection file
  if (x == 1) {
    if (ncol(selection_list) < 4) {
      return(error1)
    }
    else {
      if (names(selection_list)[1] == "CandidateID" & names(selection_list)[2] == "Parent1" & names(selection_list)[3] == "Parent2" & names(selection_list)[4] == "Index")
        return(success)
      else
        return(error1)
    }
  }
  #case x = 2 is when both selection and reference files are given
  if (x==2) {
    if (ncol(reference_list) < 3) {
      return(error2)
    }
    else if (!(ncol(selection_list) == 1 | ncol(selection_list) >= 4)) {
      return(error1)
    }
    else {
      if (names(selection_list)[1] == "CandidateID" & names(selection_list)[2] == "Parent1" & names(selection_list)[3] == "Parent2" & names(selection_list)[4] == "Index")
        if (names(reference_list)[1] == "ID" & names(reference_list)[2] == "Parent1" & names(reference_list)[3] == "Parent2")
          return(success)
        else
          return(error2)
      else
        return(error1)
      
    }
  }
}

###########################################################################
### Function to generate A Matrix for the given candidate file          ###
###########################################################################

generate_AMatrix <- function(selection_list,reference_list,progress, output) {
  if (ncol(selection_list) == 1) {
    names(selection_list)[1] <- "CandidateID"
    namevector <- c("Parent1", "Parent2", "Index")
    selection_list[, namevector] <- NA
    for (i in nrow(selection_list)) {
      j <- which(reference_list[,1] == selection_list[i,1])
      if (is.integer(j) && length(j) > 0) {
        selection_list[i,2] = as.character(reference_list[j,2])
        selection_list[i,3] = as.character(reference_list[j,3])
      }
    }
  }
  else {
    names(selection_list)[1] <- "CandidateID"
    names(selection_list)[2] <- "Parent1"
    names(selection_list)[3] <- "Parent2"
    names(selection_list)[4] <- "Index"
  }
  
  names(reference_list)[1] <- "CandidateID"
  names(reference_list)[2] <- "Parent1"
  names(reference_list)[3] <- "Parent2"
  #names(reference_list)[4] <- "Index"
  
  selection_list$Parent1[which(selection_list$Parent1==0)] <- NA
  selection_list$Parent2[which(selection_list$Parent2==0)] <- NA
  reference_list$Parent1[which(reference_list$Parent1==0)] <- NA   
  reference_list$Parent2[which(reference_list$Parent2==0)] <- NA     
  numberIndividualsOriginal <- nrow(selection_list)
  
  while (1)
  {
    #count number of individuals in the file before, and compare with after at the end
    numberIndividualsBefore <- nrow(selection_list)
    #adding missing founders/parents not in the selection col1 to the sorted input file col1
    #ensure no duplication in adding the founders/parents to col1 
    selection_list_complete <- unique(add.Inds(selection_list)[ , 1:3])
    reference_list_complete <- unique(add.Inds(reference_list)[ , 1:3])
    
    #Individuals checked with entire database/reference file to ensure missing information(in this case parents parents (grandparents)) included. 
    for(i in 1:nrow(selection_list_complete))
    {
      progress$inc(1/nrow(selection_list_complete), detail=paste("Auto Completion of Pedigree"))
      if( (is.na(selection_list_complete[i,2])) & (is.na(selection_list_complete[i,3]))) {
        j <- which(reference_list[,1] == selection_list_complete[i,1])
        if (is.integer(j) && length(j) > 0) {
          selection_list_complete[i,2] = as.character(reference_list[j,2])
          selection_list_complete[i,3] = as.character(reference_list[j,3])
        }
      }
    }
    #after adding founders/parents and filling in their parents(grandparents) sorted again
    #deciding the order and sorting sccording to the order
    selection_list_complete <- selection_list_complete[order(orderPed(selection_list_complete)),]
    reference_list_complete <- reference_list_complete[order(orderPed(reference_list_complete)),]
    if (numberIndividualsBefore == nrow(selection_list_complete)) {
      break
    }
    selection_list <- selection_list_complete
    reference_list <- reference_list_complete
  }  
  #end of loop
  
  #computing the numerator A matrix using AGHmatrix package applying Henderson
  selection_list_complete[is.na(selection_list_complete)] <-0      #NA for unknown parents changed back to 0
  #reference_list_complete[is.na(reference_list_complete)] <-0      #NA for unknown parents changed back to 0
  
  progress$set(message = "Running Part 1", value = 0)
  progress$inc(0.5, detail=paste("Generating A Matrix"))
  
  #get numerator A matrix
  selection_list_A <- Amatrix(data = selection_list_complete, ploidy = 2)
  #reference_list_A <- Amatrix(data = reference_list_complete, ploidy = 2) 
  
  progress$set(message = "Running Part 1", value = 0)
  progress$inc(0.5, detail=paste("Writing Output to Files"))
  
  write.csv(selection_list_A, row.names=FALSE, file = sorted_output_filename2 )
  
  progress$set(message = "Running Part 1", value = 0)
  progress$inc(0.5, detail=paste("Calculating Statistics"))
  
  #computing the population inbreeding coefficient and the group coancestry.
  str5 <- calculate_coefficients(selection_list_A, "")
  
  selection_list_complete$inbreeding_coeff <- 0
  selection_list_complete$coancestry_coeff <- 0
  
  for (i in 1:nrow(selection_list_complete)) 
  {
    selection_list_complete[i,4] <- selection_list_A[i,i]-1
    par1 <- selection_list_complete[i,2]
    par2 <- selection_list_complete[i,3]
    if (par1 != 0 && par2 != 0) {
      index1 <- match(par1, colnames(selection_list_A))
      index2 <- match(par2, colnames(selection_list_A))
      selection_list_complete[i,5] <- selection_list_A[index1, index2]  / 2
    }
  }
  
  #only print selected columns in pedigree file
  write.csv(selection_list_complete[,c("CandidateID", "Parent1", "Parent2")], row.names=FALSE, file = sorted_output_filename)
  
  numberIndividualsAfter <- nrow(selection_list_complete)

  strtitle <- paste("Outputs saved in the folder")
  strtitlestrong <- strong(strtitle)
  
  str0 <- paste("Complete Pedigree is saved as ", sorted_output_filename, sep="")
  str001 <- paste("Numerator Relationship Matrix (", sep="")
  str002 <- strong(paste("A"),sep="")
  str003 <- paste(" Matrix) file is saved as ", sorted_output_filename2, sep="")
  str000strong <- paste(str001, str002, str003)
  str000 <- paste("Numerator Relationship Matrix (A Matrix) file is saved as ", sorted_output_filename2, sep="")
  str00 <- paste("Log File is saved as ", output_logfile1, sep="")
  
  str1 <- paste("Candidate file statistics")
  str1strong <- strong(str1)
  str2 <- paste("Original Number of Individuals: ", numberIndividualsOriginal, sep="")
  str3 <- paste("Final Number of Individuals: ", numberIndividualsAfter, sep="")
  str4 <- paste("Number of Individuals Added: ", numberIndividualsAfter - numberIndividualsOriginal, sep="")
  estr <- paste("")
  nstr <- paste("\n")
  
  str5split <- strsplit(str5, split='<br/>')
  
  str5merged <- paste(str5split[[1]][1], nstr, str5split[[1]][2], nstr, str5split[[1]][3], nstr, str5split[[1]][5], nstr, str5split[[1]][6], nstr, str5split[[1]][7])
  
  logfile_text <- paste(nstr, str1, nstr, str2, nstr, str3, nstr, str4, nstr, str5merged, nstr, strtitle, nstr, str0, nstr, str000, nstr, str00, nstr)
  
  logfile <- file(output_logfile1)
  write(as.character(logfile_text), logfile, append=TRUE)
  close(logfile)
  
  HTML(paste(estr, str1strong, estr, str2, str3, str4, estr, str5, strtitlestrong, estr, str1, estr, str0, str000strong, str00, estr, sep='<br/>'))
}

###########################################################################
### Function to calculcate inbreeding and coacnestry coefficients       ###
###########################################################################

calculate_coefficients <- function (list_A, message_string) {
  inbreeding_coeff <- sum_coeff <- count_coeff <- 0
  inbreeding_max <- coancestry_max <- 0
  inbreeding_min <- coancestry_min <- 100
  
  for (i in 1:nrow(list_A))
  {
    inbreeding_coeff <- inbreeding_coeff + list_A[i,i]
    if (list_A[i,i] > inbreeding_max) {
      inbreeding_max <- list_A[i,i]
    }
    if (list_A[i,i] < inbreeding_min) {
      inbreeding_min <- list_A[i,i]
    }
    #lower triangle only
    for(j in 1:i) 
    {
      sum_coeff <- sum_coeff + list_A[i,j]
      count_coeff <- count_coeff + 1
      
      if ((list_A[i,j]/2) >= coancestry_max && i!=j)  {
        coancestry_max <- list_A[i,j] / 2
      }  
      if (list_A[i,j] < coancestry_min && i!=j) {
        coancestry_min <- list_A[i,j] / 2
      }  
    }
  }
  #removing values of principal diagonal
  sum_coeff <- sum_coeff - inbreeding_coeff
  count_coeff <- count_coeff - nrow(list_A)
  
  coa_matrix <- list_A[lower.tri(list_A)]
  #print(mean(coa_matrix)/2)
  #print(sd(coa_matrix)/2)
  #print((sd(coa_matrix)/2) / sqrt(nrow(list_A)))
  
  #calculating inbreeding coefficient and coancestry coefficient
  inbreeding_coeff <- as.double(inbreeding_coeff / nrow(list_A))-1
  coancestry_coeff <- (sum_coeff / count_coeff) / 2
  
  inbreeding_coeff <- signif(inbreeding_coeff,4)
  coancestry_coeff <- signif(coancestry_coeff,4)
  
  #inbreeding and coancestry of both selection and reference files, average max, min (should be 0)
  str11 <- paste("Average Inbreeding Coeffecient: ", inbreeding_coeff, sep="")
  str21 <- paste("Max Inbreeding Coeffecient: ", inbreeding_max-1, sep="")
  str31 <- paste("Min Inbreeding Coeffecient: ", inbreeding_min-1, sep="")
  
  str41 <- paste("Average Coancestry Coeffecient: ", coancestry_coeff, sep="")
  str51 <- paste("Max Coancestry Coeffecient: ", coancestry_max, sep="")
  str61 <- paste("Min Coancestry Coeffecient: ", coancestry_min, sep="")
  estr <- paste("")

  if (message_string == "Selection File")
  {
    graph_stats_data[1,1] <<- coancestry_coeff
    graph_stats_data[1,2] <<- inbreeding_coeff
    
  }  
  
  else if (message_string == "Mating List File")
  {
    graph_stats_data[2,1] <<- coancestry_coeff
    graph_stats_data[2,2] <<- inbreeding_coeff
    
  }  
  
  else
    return (HTML(paste(str11,str21,str31,estr,str41,str51,str61,estr, sep='<br/>')))
}

###########################################################################
### Function to check validity of files uploaded for Part 2 input       ###
###########################################################################

checkFiles_part2 <- function (filedata) 
{
  error1 <- "Candidate file needs to have the following columns at least: 'CandidateID', 'Parent1', 'Parent2', 'Index', 'MaxUse'. Please upload an appropriate candidate file"
  success <- "1"
  
  if (names(filedata)[1] == "CandidateID" & names(filedata)[2] == "Parent1" & names(filedata)[3] == "Parent2" & names(filedata)[4] == "Index" & names(filedata)[5] == "MaxUse")
    return(success)
  else
    return(HTML(error1))
}

run_omm_code <- function ()
{
  lower <- c(-10,-10)
  upper <- -lower
  
  iterations <- glo_iterations
  
  DEoptim(optimal_mating_solution,lower,upper,DEoptim.control(itermax = iterations))
  
}

###########################################################################
### Function to run the optimal mating solution (DE Algorithm)          ###
###########################################################################

optimal_mating_solution <- function (x)
{
  parents_list <- glo_parents_list
  AMatrix <- glo_AMatrix
  ebvweight <- glo_ebvweight
  progenyinbreedingweight <- glo_progenyinbreedingweight
  parentalcoaweight <- glo_parentalcoaweight
  numberbabies <- glo_numberbabies
  amatrixcutoff <- glo_amatrixcutoff
  
  progress <- shiny::Progress$new()
  # Make sure it closes when we exit this reactive, even if there's an error
  on.exit(progress$close())
  
  progress$set(message = "Running Optimal Mating Algorithm", value = 0)
  
  parents_list$MaxUse <- ceiling(parents_list$MaxUse)

  #editing parents_list such that we repeat entries per ID depending on the number of maxuse
  parents_list <- expandRows(parents_list, "MaxUse")
  parents_list$MaxUse <- 1
  
  file_babies <- sum(parents_list$MaxUse)
  number_babies <- numberbabies
  
  if (file_babies < number_babies)
  {
    print("Number of progeny requested in the program is more than the total of MaxUse for Parents")
    return(1)
    
  }
  diff <- babies <- 0
  
  A_list <- parents_list
  
  A_list$Sorting <- 0
  sort_col <- match("Sorting",names(A_list))
  
  for (i in 1:(number_babies*2))
  {
    A_list[i,sort_col] <- runif(1, 0, 99)  
  }
  
  A_list <- A_list[order(-A_list$Sorting),] 
  
  B_list <- subset(A_list, select = -c(Sorting))
  
  B_list <- B_list[order(-B_list$Index),]
  
  #creating data frame for mating list
  #initializing mating list matrix to 0
  mating_list <- data.frame(matrix(nrow=nrow(A_list), ncol=nrow(B_list)))
  mating_list[is.na(mating_list)] <- 0
  mating_output <- data.frame(CandidateID= character(), Parent1= character(), Parent2= character(), Index=numeric(), Inbreeding_Coefficient=numeric(), Coancestry_Coefficient=numeric(), stringsAsFactors = FALSE)
  mating_output_id <- 1
  
  A_list$used <- B_list$used <- 0
  
  #prepare mating list
  #1. figure out which parent2 to mate by taking the next available parent2 based on maxuse
  #2. figure out which parent1 to mate by checking 3 conditions. 1. parent1 used < maxuse 2. not full sibs 3. unrelated parent1 and parent2
  #3. even though this loop is for all babies that we want, it will only run till the value of min(sum of parent1_maxuse, sum of parent2_maxuse)
  #4. B_list maxuse will be reached. A_list maxuse for some individuals will not be reached in this loop
  #5. for the rest of the babies, next for loop
  
  #initiate parent2_id to 0
  #find parent2_id by taking the next available parent2 to mate
  #if a parent2 has not been used as much as the maxuse of that parent2, then that is the parent2 we will mate
  for (i in 1:nrow(A_list))
  {
    parent2_id <- 0
    
    tempi <- i
    while (parent2_id == 0)
    {
      if (A_list$used[i] < 1)
      {
        parent2_id <- i
      }  
      i <- i + 1
    }
    i<- tempi
    
    #if babies produced is more than what is required, then stop
    if (babies >= number_babies)
      break
    else
      progress$inc(1/nrow(number_babies), detail=paste("Running Optimal Mating Algorithm"))
    
    #mate that parent2_id with the next suitable parent1 
    for (j in 1:parent2_id)
    {
      #find if suitable or not. 3 conditions
      #1. if the parent1 has not been used upto the maxuse. 
      #2. if the parent1 has not mated with this parent2_id already. 
      #3. parent1 the parent2 are unrelated
    
      if (B_list$used[j] >= 1)
        next
      
      if (unrelated(as.character(B_list[j,1]), as.character(A_list[parent2_id,1]), AMatrix, amatrixcutoff) == 0)
        next
      if (check_mated(mating_output, as.character(A_list[parent2_id,1]), as.character(B_list[j,1])) == 0)
      {
        #make the parent1 and parent2 mate
        mating_list[parent2_id,j] = 1
        
        #increase the used component of the parent1 and parent2 both by 1
        B_list$used[j] = B_list$used[j] + 1
        A_list$used[parent2_id] = A_list$used[parent2_id] + 1
        
        #added to make monoecious to count it double
        A_list$used[j] = A_list$used[j] + 1
        B_list$used[parent2_id] = B_list$used[parent2_id] + 1
        
        col <- which(colnames(AMatrix)==B_list$CandidateID[j])
        row <- which(colnames(AMatrix)==A_list$CandidateID[parent2_id])
        
        #increase the count of babies we have had so far by 1
        babies <- babies + 1
        
        mating_output[mating_output_id, 1] <- paste("PB_", mating_output_id, sep='')
        mating_output[mating_output_id, 2] <- as.character(A_list[parent2_id,1])
        mating_output[mating_output_id, 3] <- as.character(B_list[j,1])
        mating_output[mating_output_id, 4] <- (A_list$Index[parent2_id] + B_list$Index[j]) / 2
        mating_output[mating_output_id, 5] <- AMatrix[row,col] / 2
        
        par0 <- mating_output[mating_output_id, 1]
        par1 <- mating_output[mating_output_id, 2]
        par2 <- mating_output[mating_output_id, 3]
        index1 <- match(par1, colnames(AMatrix))
        index2 <- match(par2, colnames(AMatrix))
        mating_output[mating_output_id, 6] <- AMatrix[index1, index2] / 2
        
        mating_output_id <- mating_output_id + 1
        
        break
      }
    }
  }

  return_list <- calculate_msi(A_list, B_list, mating_list, AMatrix, ebvweight, progenyinbreedingweight, parentalcoaweight)
  
  write.csv(mating_output, row.names=FALSE, file = output_filename)
  
  str0 <- paste("Completed creating Optimized Mating list. File saved as ", output_filename_optimized, sep="")
  estr <- paste("")
  
  return(-return_list$MSI)
}

###########################################################################
### Function to find if both parents are related or not                 ###
### Args - first parent, second parent, A matrix of progeny, cutoff     ###
### Returns 0 if related and 1 if not related                           ###
###########################################################################

unrelated <- function (first_parent, second_parent, AMatrix, amatrixcutoff)
{
  a_matrix_cutoff <- amatrixcutoff
  col <- which(colnames(AMatrix)==first_parent)
  row <- which(colnames(AMatrix)==second_parent)
  #check if that value is 0 (unrelated) or not (related)
  if (AMatrix[row,col]/2 > a_matrix_cutoff)
    return(0)
  else
    return(1)
}

###########################################################################
### Function to calculate MSI of the mating list                        ###
###########################################################################

calculate_msi <- function (A_list, B_list, mating_list, AMatrix, ebvweight, progenyinbreedingweight, parentalcoaweight)
{
  msi_total <- 0
  
  #w1
  #ebv_weight <- 1
  ebv_weight <- ebvweight
  #w2
  #progeny_inbreeding_weight <- -0.1
  progeny_inbreeding_weight <- progenyinbreedingweight
  #w3
  #parental_coa_weight <- -1
  parental_coa_weight <- parentalcoaweight
  
  #msi_part2 = (weight on progeny inbreeding) * (predicted progency mean inbreeding coefficient)
  msi_part2 <- 0
  numberOfMatings <- 0
  
  mating_list$sumRow <- rowSums(mating_list)
  
  for (i in 1:nrow(mating_list))
  {
    if (mating_list$sumRow[i] != 1)
      next
    for (j in 1:ncol(mating_list))
    {
      if(mating_list[i,j] == 1)
      {
        A_list_name <- A_list$CandidateID[i]
        B_list_name <- B_list$CandidateID[j]
        AMatrix_A_index <- which(names(AMatrix) == A_list_name)
        AMatrix_B_index <- which(names(AMatrix) == B_list_name)
        msi_part2 <- msi_part2 + (AMatrix[AMatrix_A_index, AMatrix_B_index] / 2)
        numberOfMatings <- numberOfMatings + 1
        break
      }
    }
  }
  
  msi_part2 <- msi_part2 / numberOfMatings
  
  #msi_part1 = (Weight on bgv) * (contribution of candidates vector prime) * (vector of EBV)
  msi_part1 <- 0
  for (i in 1:nrow(A_list))
    msi_part1 <- msi_part1 + (A_list$Index[i] * A_list$used[i])
  for (i in 1:nrow(B_list))
    msi_part1 <- msi_part1 + (B_list$Index[i] * B_list$used[i])
  
  #take average of EBV
  msi_part1 <- msi_part1 / (2*numberOfMatings)
  
  #msi_part3 = (weight on patrental COA / 2) * (contribution of candidates vector prime) * (A matrix of mating list) * (contribution of candidates vector)
  msi_part3 <- 0
  x_contribution <- data.frame(matrix(nrow=ncol(AMatrix),ncol=1))
  x_contribution_prime <- data.frame(matrix(ncol=ncol(AMatrix),nrow=1))
  
  for (i in 1:ncol(AMatrix))
  {
    A_colname <- colnames(AMatrix[i])
    A_used <- A_list$used[which(A_list$CandidateID == A_colname)]
    B_used <- B_list$used[which(B_list$CandidateID == A_colname)]
    if( length(A_used) > 0 )
    {
      x_contribution[i,1] <-  x_contribution_prime[1,i] <- 1
    }
    else if (length(B_used) > 0)
    {
      x_contribution[i,1] <- x_contribution_prime[1,i] <- 1
    }
    else 
    {
      x_contribution[i,1] <- x_contribution_prime[1,i] <- 0
    }
  }
  
  x_contribution <- do.call(cbind, x_contribution)
  x_contribution_prime <- do.call(cbind, x_contribution_prime)
  AMatrix <- do.call(cbind, AMatrix)
  resultant_matrix <- ((x_contribution_prime %*% AMatrix)) # %*% x_contribution)
  resultant_matrix <- resultant_matrix %*% x_contribution
  
  msi_part3 <- as.double(resultant_matrix[1,1]) * 0.5
  
  #average of msi_part3
  msi_part3 <- msi_part3 / (4*numberOfMatings*numberOfMatings)
  
  #we are multiplying with weights because in the total equation we need weights factored in
  msi_part1 <- msi_part1 * ebv_weight
  msi_part2 <- msi_part2 * progeny_inbreeding_weight
  msi_part3 <- msi_part3 * parental_coa_weight
  
  msi_total <- msi_part1 + msi_part2 + msi_part3
  return_list <- list("MSI" = msi_total, "EBV" = msi_part1, "ProgenyInbreeding" = msi_part2, "ParentalCOA" = msi_part3)
  
  return(return_list)
}

###########################################################################
### Function to check if 2 candidates have already mated                ###
###########################################################################

check_mated <- function (mating_output, x_name, y_name)
{
  return_value <- 0
  if (nrow(mating_output) > 0)
  {
    for (i in 1:nrow(mating_output))
    {
      if (is.na(mating_output[i,2]) && is.na(mating_output[i,3]))
        next
      else if (mating_output[i,3] == y_name && mating_output[i,2] == x_name)
        return(1)
      else if (mating_output[i,3] == x_name && mating_output[i,2] == y_name)
        return(1)
    }
  }
  return(return_value)
}

###########################################################################
### Function to print outputs of Part 2 of the code                     ###
###########################################################################

print_output_part2 <- function (output_filename, parents_list, duration)
{
  selection_list <- read.csv(output_filename, check.names = FALSE)
  reference_list <- parents_list
  
  selection_list$Parent1[which(selection_list$Parent1==0)] <- NA
  selection_list$Parent2[which(selection_list$Parent2==0)] <- NA
  reference_list$Parent1[which(reference_list$Parent1==0)] <- NA   
  reference_list$Parent2[which(reference_list$Parent2==0)] <- NA     
  numberIndividualsOriginal <- nrow(selection_list)

  while (1)
  {
    #count number of individuals in the file before, and compare with after at the end
    numberIndividualsBefore <- nrow(selection_list)
    #adding missing founders/parents not in the selection col1 to the sorted input file col1
    #ensure no duplication in adding the founders/parents to col1 
    selection_list_complete <- unique(add.Inds(selection_list)[ , 1:3])
    reference_list_complete <- unique(add.Inds(reference_list)[ , 1:3])
    
    #Individuals checked with entire database/reference file to ensure missing information(in this case parents parents (grandparents)) included. 
    for(i in 1:nrow(selection_list_complete))
    {
      if( (is.na(selection_list_complete[i,2])) & (is.na(selection_list_complete[i,3]))) {
        j <- which(reference_list[,1] == selection_list_complete[i,1])
        if (is.integer(j) && length(j) > 0) {
          selection_list_complete[i,2] = as.character(reference_list[j,2])
          selection_list_complete[i,3] = as.character(reference_list[j,3])
        }
      }
    }
    #after adding founders/parents and filling in their parents(grandparents) sorted again
    #deciding the order and sorting sccording to the order
    selection_list_complete <- selection_list_complete[order(orderPed(selection_list_complete)),]
    reference_list_complete <- reference_list_complete[order(orderPed(reference_list_complete)),]
    if (numberIndividualsBefore == nrow(selection_list_complete)) {
      break
    }
    selection_list <- selection_list_complete
    reference_list <- reference_list_complete
  }  
  #end of loop
  
  #computing the numerator A matrix using AGHmatrix package applying Henderson
  selection_list_complete[is.na(selection_list_complete)] <-0      #NA for unknown parents changed back to 0
  
  #get numerator A matrix
  selection_list_A <- Amatrix(data = selection_list_complete, ploidy = 2)
  
  #computing the population inbreeding coefficient and the group coancestry.
  str50 <- calculate_coefficients(selection_list_A, "")
  
  #avg_index <- round(mean(output_data$Index),digits=2)
  #avg_inbreeding <- signif(mean(output_data$Inbreeding_Coefficient),4)
  #avg_coa <- signif(mean(output_data$Coancestry_Coefficient),4)
  
  crosses <- length(selection_list$Index)
  
  str11 <- paste("Mating List statistics")
  str11strong <- strong(str11)
  str1 <- paste("Number of Crosses = ", crosses, sep="")
  str16 <- paste("Number of Iterations = ", glo_iterations, sep="")
  str12 <- paste("Index Weight = ", glo_ebvweight, sep="")
  str13 <- paste("Inbreeding Weight = ", glo_progenyinbreedingweight, sep="")
  str14 <- paste("Coancestry Weight = ", glo_parentalcoaweight, sep="")
  str15 <- paste("Coancestry Constraint = ", glo_amatrixcutoff, sep="")
 
  #str2 <- paste("Average Progeny Index = ",avg_index, sep="")
  #str3 <- paste("Average Progeny Inbreeding = ",avg_inbreeding, sep="")
  #str4 <- paste("Average Progeny Coancestry = ",avg_coa, sep="")
  
  strtitle <- paste("Outputs saved in the folder")
  strtitlestrong <- strong(strtitle)
  str0 <- paste("Mating List is saved as ", output_filename, sep="")
  str00 <- paste("Log File is saved as ", output_logfile2, sep="")
  estr <- paste("")
  nstr <- paste("\n")
  
  if (duration < 60)
      str5 <- paste("Runtime ", signif(duration,2), " secs", sep="")
  else
  {
    duration <- duration/60
    duration2 <- duration%%60
    duration <- duration + duration2/100
    str5 <- paste("Runtime ", signif(duration,2), " mins", sep="")
    
  }
  
  str5strong <- strong(str5)
  
  logfile_text <- paste(nstr, str11, nstr, str1, nstr, str16, nstr, str12, nstr, str13, nstr, str14, nstr, str15, nstr, str50, nstr, strtitle, nstr, str0, nstr, str5, nstr)
  
  logfile <- file(output_logfile2)
  write(as.character(logfile_text), logfile, append=TRUE)
  
  close(logfile)
  
  HTML(paste(estr, str11strong, estr, str1, str16, estr, str12, str13, str14, str15, estr, str50, estr, strtitlestrong, estr, str0, str00, estr, str5strong, estr, sep='<br/>'))
}

###########################################################################
### UI definition for the ShinyApp                                      ###
###########################################################################

#shinyApp UI in Dashboard style
ui <- dashboardPage(
  dashboardHeader(title = "AgMate"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Part 1 - Generating A Matrix", tabName = "part1", icon = icon("dashboard")),
      menuItem("Part 2 - Mating Algorithm", tabName = "part2", icon = icon("th")),
      menuItem("Part 3 - Visualizing Results", tabName = "part3", icon = icon("bar-chart-o"))
    )
  ),
  dashboardBody(
    tags$style(
      HTML(".shiny-notification {
              height: 100px;
              width: 400px;
              position:fixed;
              top: calc(50% - 50px);;
              right: calc(50% - 400px);;
            }
           "
      )
    ),
    tabItems(
      # First tab content
      tabItem(tabName = "part1",
        fluidRow(
          box (
            h3("Part 1"),
            hr(),
            h4("Select Candidate File"),
            h5("Description: A CSV file with list of candidates to be selected and mated."),
            actionButton("part1_sample1", "Sample"),
            fileInput("file1", "Choose CSV File",multiple = FALSE,
                    accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
            h4("Select Reference File (Optional to upload)"),
            h5("Description: A CSV file with all the individuals needed to complete the pedigree."),
            actionButton("part1_sample2", "Sample"),
            fileInput("file2", "Choose CSV File",multiple = FALSE,
                    accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
            hr(),
            actionButton("part1_button", "Run"),
            width = 6
          ),
          box (
            h3("Output from Part 1"),
            htmlOutput("contents1"),
            width = 6
          )
        ),
        
        fluidRow(
          box(
            actionButton("close1", "Exit", align="right",onclick = "setTimeout(function(){window.close();},500);"),
            width = 6
          ),
          box (
            downloadButton("downloadData_part1_1", "Download Stats"),
            downloadButton("downloadData_part1_2", "Download Pedigree File"),
            downloadButton("downloadData_part1_3", "Download A Matrix File"),
            width = 6
          )
        ),
        fluidRow(
          box(
            title = "AgMate : An Optimal Mating Software by Khushi Goda",
            width = 12 
          )
        )
      ),
      
     # Second tab content
      tabItem(tabName = "part2",
        fluidRow(
          box (
            box (
              h3("Define parameters"),
              hr(),
              numericInput("numberbabies", 
                           h4("Number of Crosses"), 
                           value = 50),
              
              numericInput("amatrixcutoff",
                           h4("Coancestry Constraint"),
                           value = 0.124,
                           min=0,
                           max=2,
                           step=0.001),
              
              sliderInput("iterations",
                          h4("Number of Iterations"),
                          min = 0, max = 100, value = 1),
            ),
            
            box (
              h3("Weights for optimization"),
              hr(),
              
              numericInput("ebvweight", 
                         h4("Index Weight"), 
                         value = 1,
                         min=-20,
                         max=20,
                         step=0.1),
            
              numericInput("progenyinbreedingweight", 
                         h4("Inbreeding Weight"), 
                         value = -0.1,
                         min=-20,
                         max=20,
                         step=0.1),
            
              numericInput("parentalcoaweight", 
                         h4("Coancestry Weight"), 
                         value = -1,
                         min=-20,
                         max=20,
                         step=0.1),
              
            ),
            width = 6
          ),
          
          box (
            h3("Part 2"),
            hr(),
            h4("Select Candidate File"),
            h5("Description: A CSV File with the list of candidates to mate and their maximum use."),
            actionButton("part2_sample1", "Sample"),
            fileInput("file3", "Choose CSV File",multiple = FALSE,
                      accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
            h4("Select Numerator Relationship Matrix"),
            h5("Description: Numerator Relationship Matrix (A Matrix) file for the whole pedigree from Part 1."),
            fileInput("file4", "Choose CSV File",multiple = FALSE,
                      accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
            hr(),
            actionButton("part2_button", "Run"),
            width = 3
          ),
          
          box (
            h3("Output from Part 2"),
            htmlOutput("contents2"),
            width = 3
          )
        ),
            
        fluidRow(
          box(
            actionButton("close2", "Exit", align="right",onclick = "setTimeout(function(){window.close();},500);"),
            width = 6
          ),
          box (
            downloadButton("downloadData_part2_1", "Download Stats"),
            downloadButton("downloadData_part2_2", "Download Mating List"),
            width = 6
          )
        ),
        fluidRow(
          box(
            title = "AgMate : An Optimal Mating Software by Khushi Goda",
            width = 12 
          )
        )
      ),
     
     # Third tab content
     tabItem(tabName = "part3",
             fluidRow(
               box (
                 h3("Part 3"),
                 hr(),
                 h4("Select Candidate File"),
                 #h5("Description: A CSV file with list of candidates to be selected and mated."),
                 fileInput("file5", "Choose CSV File",multiple = FALSE,
                           accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                 h4("Select Mating List"),
                 #h5("Description: A CSV file with all the individuals needed to complete the pedigree."),
                 fileInput("file6", "Choose CSV File",multiple = FALSE,
                           accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                 hr(),
                 actionButton("part3_button", "Plot"),
                 width = 4
               ),
               box (
                 h3("Distribution of Index"),
                 plotOutput("contents3_1"),
                 hr(),
                 downloadButton("download_graph1", "Download Graph"),
                 width = 4
               ),
               box (
                 h3("Mean Coancestry and Inbreeding"),
                 plotOutput("contents3_2"),
                 hr(),
                 downloadButton("download_graph2", "Download Graph"),
                 width = 4
               )
             ),
             
             fluidRow(
               box(
                 actionButton("close3", "Exit", align="right",onclick = "setTimeout(function(){window.close();},500);"),
                 width = 12
               )
             ),
             fluidRow(
               box(
                 title = "AgMate : An Optimal Mating Software by Khushi Goda",
                 width = 12 
               )
             )
     )
    )
  )
)

###########################################################################
### Server definition for the ShinyApp                                  ###
###########################################################################

server <- function(input, output) {
  
  iv <- InputValidator$new()
  iv$add_rule("ebvweight", sv_between(-20,20))
  iv$add_rule("parentalcoaweight", sv_between(-20,20))
  iv$add_rule("progenyinbreedingweight", sv_between(-20,20))
  iv$add_rule("amatrixcutoff", sv_between(0,2))
  iv$add_rule("iterations", sv_between(1,100))
  iv$enable()
  
  observe({
    if (input$close1 > 0) stopApp()    # stop shiny
  })
  
  observe({
    if (input$close2 > 0) stopApp()    # stop shiny
  })
  
  observe({
    if (input$close3 > 0) stopApp()    # stop shiny
  })
  
  observeEvent(input$part1_sample1, {
    showModal(modalDialog(
      tags$img(src='sample1part1.png'),
      size = c("l"),
      easyClose = TRUE
    ))
  })
  
  observeEvent(input$part1_sample2, {
    showModal(modalDialog(
      tags$img(src='sample2part1.png'),
      size = c("l"),
      easyClose = TRUE
    ))
  })
  
  observeEvent(input$part2_sample1, {
    showModal(modalDialog(
      tags$img(src='sample1part2.png'),
      size = c("l"),
      easyClose = TRUE
    ))
  })
  
  part1_reactive <- eventReactive(input$part1_button,{
    if (is.null(input$file1$datapath)) {
      return(HTML("Candidate File is Required. Please Re-run the App."))
    }
    else if (is.null(input$file2$datapath)) {
      selection_list <- read.csv(input$file1$datapath, check.names = FALSE)
      reference_list <- selection_list
      x <- 1
      checkFilesReturnMessage <- checkFiles_part1(selection_list, reference_list,x)    
    }
    else {
      selection_list <- read.csv(input$file1$datapath, check.names = FALSE)    
      reference_list <- read.csv(input$file2$datapath, check.names = FALSE)   
      x <- 2
      checkFilesReturnMessage <- checkFiles_part1(selection_list, reference_list,x)    
    }
    
    if (!(checkFilesReturnMessage == "1")) {
      return(HTML(checkFilesReturnMessage))
    }
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())
    
    progress$set(message = "Running Part 1", value = 0)
    generate_AMatrix(selection_list,reference_list,progress, output)
  
  })
  
  output$contents1 <- renderUI({
    if (input$part1_button)
      part1_reactive()
  })
  
  output$downloadData_part1_1 <- downloadHandler(
    filename = function() { 
      paste(output_logfile1, sep="")
    },
    content = function(file) {
      file.copy(output_logfile1,file)
    }
  )
  
  output$downloadData_part1_2 <- downloadHandler(
    filename = function() { 
      paste(sorted_output_filename, sep="")
    },
    content = function(file) {
      file.copy(sorted_output_filename,file)
    }
  )
  
  output$downloadData_part1_3 <- downloadHandler(
    filename = function() { 
      paste(sorted_output_filename2, sep="")
    },
    content = function(file) {
      file.copy(sorted_output_filename2,file)
    }
  )
  
  part2_reactive <- eventReactive(input$part2_button,{
    
    if (is.null(input$file3$datapath)) {
      return(HTML("Candidate File is Required. Please Re-run the App."))
    }
    else if (is.null(input$file4$datapath)) {
      return(HTML("AMatrix File is Required. Please Re-run the App."))
    }
    else {
      parents_list <- read.csv(input$file3$datapath, check.names = FALSE)
      AMatrix_file <- read.csv(input$file4$datapath, check.names = FALSE)    
    }
    
    validate_files_message <- checkFiles_part2(parents_list)
    if (validate_files_message != "1")
      return(validate_files_message)
    
    #global variables
    glo_parents_list <<- read.csv(input$file3$datapath, check.names = FALSE)
    glo_AMatrix <<- read.csv(input$file4$datapath, check.names = FALSE)    
    glo_ebvweight <<- input$ebvweight
    glo_progenyinbreedingweight <<- input$progenyinbreedingweight
    glo_parentalcoaweight <<- input$parentalcoaweight
    glo_numberbabies <<- input$numberbabies
    glo_amatrixcutoff <<- input$amatrixcutoff
    glo_iterations <<- input$iterations
    
    start_time <- Sys.time()
    run_omm_code()
    
    #optimize_EBV(output_filename, output_filename_optimized)
    #generate_input_file_next_round(input$file3$datapath, output_filename_optimized)
    
    end_time <- Sys.time()
    duration <- difftime(end_time, start_time, units=c("secs"))
    print_output_part2(output_filename, glo_parents_list, as.numeric(duration))
    
  })
  
  output$contents2 <- renderUI({
    if (input$part2_button)
      part2_reactive()
  })
  
  output$downloadData_part2_1 <- downloadHandler(
    filename = function() { 
      paste(output_logfile2, sep="")
    },
    content = function(file) {
      file.copy(output_logfile2,file)
    }
  )
  
  output$downloadData_part2_2 <- downloadHandler(
    filename = function() { 
      paste(output_filename, sep="")
    },
    content = function(file) {
      file.copy(output_filename,file)
    }
  )
  
  ###########################################################################
  ### Code for generating Index Chart                                     ###
  ###########################################################################
  
  part3_reactive_1 <- eventReactive(input$part3_button,{
      
    if (is.null(input$file5$datapath)) {
      return(HTML("Candidate File is Required. Please Re-run the App."))
    }
    else if (is.null(input$file6$datapath)) {
      return(HTML("Mating List is Required. Please Re-run the App."))  
    }
    else {
      selection_list <- read.csv(input$file5$datapath, check.names = FALSE)    
      mating_list <- read.csv(input$file6$datapath, check.names = FALSE)   
    }
    
    # read candidate list data
    can <- selection_list
    # create a column called List
    can$List <- "Candidate list"
    
    # read mating list data
    mat <- mating_list
    # create a column called List
    mat$List <- "Mating list"
    
    # data for index
    dat.index <- rbind(can[,-5], mat[,-c(5:6)])
    
    # get summary stats
    sum.index <- dat.index %>% group_by(List) %>% summarise(
      Mean = mean(Index)
    ) 
    
    # grouped density plot
    glo_plot1 <<- ggplot(dat.index, aes(x=Index, color=List, fill=List)) + theme_classic()+
      geom_density(aes(y=..count..), alpha=.3)+
      geom_vline(data=sum.index, aes(xintercept=Mean,linetype=List, color=List),size=.8)+
      scale_color_manual(values = c("blue4","red4"))+
      scale_fill_manual(values = c("white","lightyellow"))+
      scale_y_continuous(expand = c(0,0))+
      labs(x="Index value",y="Count of trees")+
      theme(legend.position = "top",
            legend.title = element_blank(),
            axis.title = element_text(size=16, color="black"),
            axis.text = element_text(size=12, colour = 'black'),
            panel.grid = element_blank(),
            axis.title.x = element_text(vjust = -1))
    
    ggsave(graph1_filename, glo_plot1)
    
    glo_plot1
      
  })
  
  output$contents3_1 <- renderPlot({
    if (input$part3_button)
      part3_reactive_1()
  })
    
  ###########################################################################
  ### Code for generating Inbreeding and Coancestry Chart                 ###
  ###########################################################################
    
  part3_reactive_2 <- eventReactive(input$part3_button,{
      
     if (is.null(input$file5$datapath)) {
       return(HTML("Candidate File is Required. Please Re-run the App."))
     }
     else if (is.null(input$file6$datapath)) {
       return(HTML("Mating List is Required. Please Re-run the App."))  
     }
     else {
       selection_list <- read.csv(input$file5$datapath, check.names = FALSE)    
       mating_list <- read.csv(input$file6$datapath, check.names = FALSE)   
     }
     
     #candidate file stats
     selection_list <- selection_list
     reference_list <- selection_list
    
     selection_list$Parent1[which(selection_list$Parent1==0)] <- NA
     selection_list$Parent2[which(selection_list$Parent2==0)] <- NA
     reference_list$Parent1[which(reference_list$Parent1==0)] <- NA   
     reference_list$Parent2[which(reference_list$Parent2==0)] <- NA     
     numberIndividualsOriginal <- nrow(selection_list)
    
     while (1)
     {
       #count number of individuals in the file before, and compare with after at the end
       numberIndividualsBefore <- nrow(selection_list)
       #adding missing founders/parents not in the selection col1 to the sorted input file col1
       #ensure no duplication in adding the founders/parents to col1 
       selection_list_complete <- unique(add.Inds(selection_list)[ , 1:3])
       reference_list_complete <- unique(add.Inds(reference_list)[ , 1:3])
       
       #Individuals checked with entire database/reference file to ensure missing information(in this case parents parents (grandparents)) included. 
       for(i in 1:nrow(selection_list_complete))
       {
         if( (is.na(selection_list_complete[i,2])) & (is.na(selection_list_complete[i,3]))) {
           j <- which(reference_list[,1] == selection_list_complete[i,1])
           if (is.integer(j) && length(j) > 0) {
             selection_list_complete[i,2] = as.character(reference_list[j,2])
             selection_list_complete[i,3] = as.character(reference_list[j,3])
           }
         }
       }
       #after adding founders/parents and filling in their parents(grandparents) sorted again
       #deciding the order and sorting sccording to the order
       selection_list_complete <- selection_list_complete[order(orderPed(selection_list_complete)),]
       reference_list_complete <- reference_list_complete[order(orderPed(reference_list_complete)),]
       if (numberIndividualsBefore == nrow(selection_list_complete)) {
         break
       }
       selection_list <- selection_list_complete
       reference_list <- reference_list_complete
     }  
     #end of loop
     
     #computing the numerator A matrix using AGHmatrix package applying Henderson
     selection_list_complete[is.na(selection_list_complete)] <-0      #NA for unknown parents changed back to 0
     #reference_list_complete[is.na(reference_list_complete)] <-0      #NA for unknown parents changed back to 0
     
     #get numerator A matrix
     selection_list_A <- Amatrix(data = selection_list_complete, ploidy = 2)
     #reference_list_A <- Amatrix(data = reference_list_complete, ploidy = 2) 
     
     #computing the population inbreeding coefficient and the group coancestry.
     calculate_coefficients(selection_list_A, "Selection File")
     
     #mating file stats
     selection_list <- mating_list
     reference_list <- selection_list
     
     selection_list$Parent1[which(selection_list$Parent1==0)] <- NA
     selection_list$Parent2[which(selection_list$Parent2==0)] <- NA
     reference_list$Parent1[which(reference_list$Parent1==0)] <- NA   
     reference_list$Parent2[which(reference_list$Parent2==0)] <- NA     
     numberIndividualsOriginal <- nrow(selection_list)
     
     while (1)
     {
       #count number of individuals in the file before, and compare with after at the end
       numberIndividualsBefore <- nrow(selection_list)
       #adding missing founders/parents not in the selection col1 to the sorted input file col1
       #ensure no duplication in adding the founders/parents to col1 
       selection_list_complete <- unique(add.Inds(selection_list)[ , 1:3])
       reference_list_complete <- unique(add.Inds(reference_list)[ , 1:3])
       
       #Individuals checked with entire database/reference file to ensure missing information(in this case parents parents (grandparents)) included. 
       for(i in 1:nrow(selection_list_complete))
       {
         if( (is.na(selection_list_complete[i,2])) & (is.na(selection_list_complete[i,3]))) {
           j <- which(reference_list[,1] == selection_list_complete[i,1])
           if (is.integer(j) && length(j) > 0) {
             selection_list_complete[i,2] = as.character(reference_list[j,2])
             selection_list_complete[i,3] = as.character(reference_list[j,3])
           }
         }
       }
       #after adding founders/parents and filling in their parents(grandparents) sorted again
       #deciding the order and sorting sccording to the order
       selection_list_complete <- selection_list_complete[order(orderPed(selection_list_complete)),]
       reference_list_complete <- reference_list_complete[order(orderPed(reference_list_complete)),]
       if (numberIndividualsBefore == nrow(selection_list_complete)) {
         break
       }
       selection_list <- selection_list_complete
       reference_list <- reference_list_complete
     }  
     #end of loop
     
     #computing the numerator A matrix using AGHmatrix package applying Henderson
     selection_list_complete[is.na(selection_list_complete)] <-0      #NA for unknown parents changed back to 0
     #reference_list_complete[is.na(reference_list_complete)] <-0      #NA for unknown parents changed back to 0
     
     #get numerator A matrix
     selection_list_A <- Amatrix(data = selection_list_complete, ploidy = 2)
     #reference_list_A <- Amatrix(data = reference_list_complete, ploidy = 2) 
     
     #computing the population inbreeding coefficient and the group coancestry.
     calculate_coefficients(selection_list_A, "Mating List File")
     
     # dataframe for inbreeding and coancestry
     inco <- data.frame(List = c("Candidate list","Mating list"),
                        Parameter = c("Coancestry","Coancestry","Inbreeding","Inbreeding"),
                        Coefficient = c(graph_stats_data[1,1],graph_stats_data[2,1],graph_stats_data[1,2],graph_stats_data[2,2]))
     
     glo_plot2 <<- ggplot(inco, aes(x=Parameter, y=Coefficient))+ theme_classic()+
       geom_bar(aes(fill=List),position = 'dodge', stat = "identity", alpha=.6,
                color='NA', size=.3, width=.7)+
       scale_fill_manual(values = c("orange2","dodgerblue"))+
       labs(x="",y="Coefficient")+
       scale_y_continuous(expand = c(0,0))+
       theme(legend.position = 'top',
             legend.title = element_blank(),
             legend.direction = "horizontal",
             legend.background = element_rect(fill = 'NA', color='NA', size=.2),
             panel.grid = element_blank(),
             axis.text = element_text(size=11, color = "black"),
             axis.title.x = element_blank(),
             axis.title.y = element_text(size=14, colour = 'black', vjust=2.2))
     
     ggsave(graph2_filename, glo_plot2)
     
     glo_plot2
      
  })
  
  output$contents3_2 <- renderPlot({
    if (input$part3_button)
      part3_reactive_2()
  })
  
  output$download_graph1 <- downloadHandler(
    filename = function() { 
      paste(graph1_filename, sep="")
    },
    content = function(file) {
      file.copy(graph1_filename,file)
    }
  )
  
  output$download_graph2 <- downloadHandler(
    filename = function() { 
      paste(graph2_filename, sep="")
    },
    content = function(file) {
      file.copy(graph2_filename,file)
    }
  )
  
}

shinyApp(ui, server)