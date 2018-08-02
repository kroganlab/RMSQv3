#! /usr/bin/env Rscript

suppressMessages(library(grDevices))
suppressMessages(library(stats))
suppressMessages(library(graphics))
suppressMessages(library(reshape2))
suppressMessages(library(limma))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(data.table))
suppressMessages(library(stringr))

theme_set(theme_bw(base_size = 15, base_family="Helvetica"))

####################################
## SMALL FUNCTIONS

#' @title Long to Wide format selecting the `Sequence` column of the evidence file
#' @description Facilitates applying the dcast function, i.e., takes long-format data and casts it into wide-format data.
#' @param d_long the data.frame in long format
#' @keywords data.frame, dcast
#' castMaxQToWide()
#' @export
castMaxQToWide <- function(d_long){
  data_w = data.table::dcast( Proteins + Sequence + Charge ~ RawFile + IsotopeLabelType, data=d_long, value.var='Intensity', fun.aggregate=sum, fill = NA)
  return(data_w)
}

#' @title Long to Wide format selecting the `Modified.sequence` column of the evidence file
#' @description Facilitates applying the dcast function, i.e., takes long-format data and casts it into wide-format data.
#' @param d_long the data.frame in long format
#' @keywords data.frame, dcast
#' castMaxQToWidePTM()
#' @export
castMaxQToWidePTM <- function(d_long){
  data_w = data.table::dcast( Proteins + Modified.sequence + Charge ~ RawFile + IsotopeLabelType, data=d_long, value.var='Intensity', fun.aggregate=sum, fill=NA)
  setnames(data_w,2,'Sequence')
  return(data_w)
}

#' @title Change a specific column name in a given data.frame
#' @description Making easier j
#' @param dataset the data.frame with the column name you want to change
#' @param oldname the old column name
#' @param neename the new name for that column
#' @keywords rename, data.frame, columns
#' changeColumnName()
#' @export
changeColumnName <- function(dataset, oldname, newname){
  if( !(oldname %in% colnames(dataset)) ){
    stop("The Column name provided <",oldname,"> was not found in the data.table provided")
  }
  names(dataset)[grep(paste0('^',oldname,'$'), names(dataset))] <- newname
  return(dataset)
}

#' @title Remove contaminants and empty proteins
#' @description Remove contaminants and erronously identified 'reverse' sequences as by MaxQuant
#' @param d_long the data.frame in long format
#' @keywords cleanup, contaminants
#' filterMaxqData()
#' @export
filterMaxqData <- function(data){
  # Remove contaminants and reversed sequences (labeled by MaxQuant)
  data_selected = data[grep("CON__|REV__",data$Proteins, invert=T),]
  # Remove empty proteins names
  blank.idx <- which(data_selected$Proteins == "")
  if(length(blank.idx)>0)  data_selected = data_selected[-blank.idx,]
  return(data_selected)
}

# Not used. Is useful?
is.uniprotAc <- function(identifier){
  grepl('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}',identifier)
}

mergeMaxQDataWithKeys <- function(data, keys, by=c('RawFile')){
  # Check if the number of RawFiles is the same.
  unique_data <- unique(data$RawFile)
  unique_keys <- unique(keys$RawFile)
  
  if (length(unique_keys) != length(unique_data)){
    keys_not_found = setdiff(unique_keys, unique_data)
    data_not_found = setdiff(unique_data, unique_keys)
    cat(sprintf("keys found: %s \t keys not in data file:\n%s\n", length(unique_keys)-length(keys_not_found), paste(keys_not_found,collapse='\t')))
    cat(sprintf("data found: %s \t data not in keys file:\n%s\n", length(unique_data)-length(data_not_found), paste(data_not_found, collapse='\t')))
  }else{
    cat("\n\t>-----+ Check point: the number of RawFiles in both keys and evidences file is identical\n\n")
  }
  
  ## select only required attributes from MQ format
  data = merge(data, keys, by=by)
  return(data)
}

#' @title Convert the SILAC evidence file to MSstats format
#' @description Converting the evidence file from a SILAC search to a format
#' compatible with MSstats. It basically modifies the Raw.files adding the
#' Heavy and Light label
#' @param filename The text filepath to the evidence file
#' @param output The text filepath of the output name
#' @keywords 
#' MQutil.SILACToLong()
#' @export
MQutil.SILACToLong = function(filename, output){
  library(data.table)
  library(reshape2)
  file = Sys.glob(filename)
  cat(sprintf('>> PROCESSING SILAC EVIDENCE FILE\n'))
  tmp = fread(file, integer64 = 'double')
  
  # reshape the data and split the heavy and light data
  tmp_long = reshape2::melt(tmp, measure.vars = c('Intensity L','Intensity H'))
  tmp_long[,Intensity:=NULL]
  setnames(tmp_long,'value','Intensity')
  setnames(tmp_long,'variable','IsotopeLabelType')
  setnames(tmp_long,'Raw file','Raw.file')
  levels(tmp_long$IsotopeLabelType) = c('L','H')
  tmp_long[!is.na(tmp_long$Intensity) && tmp_long$Intensity<1,]$Intensity=NA
  write.table(tmp_long, file=output, sep='\t', quote=F, row.names=F, col.names=T)
  cat("----- + File ",output, " has been created\n")
  return(tmp_long)
}

plotHeat <- function(mss_F, out_file, labelOrder=NULL, names='Protein', cluster_cols=F, display='log2FC'){
  heat_data = data.frame(mss_F, names=names)
  #heat_data = mss_F[,c('uniprot_id','Label','log2FC')]
  
  ## create matrix from log2FC or p-value as user defined
  if(display=='log2FC'){
    # Issues with extreme_val later if we have Inf/-Inf values.
    if( sum(is.infinite(heat_data$log2FC)) > 0 ){
      idx <- is.infinite(heat_data$log2FC)
      heat_data$log2FC[ idx ] <- NA
    }
    heat_data_w = dcast(names ~ Label, data=heat_data, value.var='log2FC') 
  }else if(display=='pvalue'){
    heat_data$adj.pvalue = -log10(heat_data$adj.pvalue+10^-16)  
    heat_data_w = dcast(names ~ Label, data=heat_data, value.var='adj.pvalue')  
  }
  
  ## try
  #gene_names = uniprot_to_gene_replace(uniprot_ac=heat_data_w$Protein)
  rownames(heat_data_w) = heat_data_w$names
  heat_data_w = heat_data_w[,-1]
  heat_data_w[is.na(heat_data_w)]=0
  max_val = ceiling(max(heat_data_w))
  min_val = floor(min(heat_data_w))
  extreme_val = max(max_val, abs(min_val))
  if(extreme_val %% 2 != 0) extreme_val=extreme_val+1
  bin_size=2
  signed_bins = (extreme_val/bin_size)
  colors_neg = rev(colorRampPalette(brewer.pal("Blues",n=extreme_val/bin_size))(signed_bins))
  colors_pos = colorRampPalette(brewer.pal("Reds",n=extreme_val/bin_size))(signed_bins)
  colors_tot = c(colors_neg, colors_pos)
  
  if(is.null(labelOrder)){
    cat("\t Saving heatmap\n")
    pheatmap(heat_data_w, scale="none", cellheight=10, cellwidth=10, filename =out_file, color=colors_tot, breaks=seq(from=-extreme_val, to=extreme_val, by=bin_size), cluster_cols=cluster_cols, fontfamily="mono")  
  }else{
    heat_data_w = heat_data_w[,labelOrder]
    pheatmap(heat_data_w, scale="none", cellheight=10, cellwidth=10, filename=out_file, color=colors_tot, breaks=seq(from=-extreme_val, to=extreme_val, by=bin_size), cluster_cols=cluster_cols, fontfamily="mono")
  }
  
  heat_data_w
}

prettyPrintHeatmapLabels <- function(uniprot_acs, uniprot_ids, gene_names){
  #uniprot_ids_trunc = gsub('([A-Z,0-9]+)_([A-Z,0-9]+)','\\1',uniprot_ids)
  #longest_id = max(nchar(uniprot_ids_trunc))
  #tmp_frame = data.frame(t=uniprot_ids_trunc, s=longest_id-nchar(uniprot_ids_trunc)+1, g=gene_names,a=uniprot_acs, stringsAsFactors=F)
  #tmp_frame[is.na(tmp_frame$t),]$t=tmp_frame[is.na(tmp_frame$t),]$a
  #result = apply(tmp_frame, 1, function(x)paste0(x[1],paste(rep(' ',x[2]),collapse=''),x[3]))
  #result = apply(tmp_frame, 1, function(x)paste0(x[4],' ',x[3],collapse=''))
  result = paste(uniprot_acs,uniprot_ids,gene_names,sep=' ')
  return(result)
}

#' @title Read in Evidence File
#' @description Read in a MaxQuant searched Evidence file using data.table. This function propperly classes each column and so fread doesn't have to guess.
#' @param evidence_file The filepath to the MaxQuant searched data (evidence) file (txt tab delimited file).
#' @keywords MaxQuant, evidence
#' read_evidence_file()
#' @export
read_evidence_file <- function(evidence_file){
  cat("Reading in evidence file...\n")
  # read in the first line to get the header names
  cols <- readLines(evidence_file, 1)
  cols <- data.frame( V1 = unlist(strsplit(cols, '\t')), stringsAsFactors = F)
  cols$idx <- 1:dim(cols)[1]
  
  # get data frame of pre-recorded column names and their respective classes
  col.classes <- as.data.frame( matrix(c("Sequence","character","Length","integer","Modifications","character","Modified sequence","character","Oxidation (M) Probabilities","character","Oxidation (M) Score Diffs","character","Acetyl (Protein N-term)","integer","Oxidation (M)","integer","Missed cleavages","integer","Proteins","character","Leading proteins","character","Leading razor protein","character","Gene names","character","Protein names","character","Type","character","Raw file","character","Experiment","character","MS/MS m/z","numeric","Charge","integer","m/z","numeric","Mass","numeric","Resolution","numeric","Uncalibrated - Calibrated m/z [ppm]","numeric","Uncalibrated - Calibrated m/z [Da]","numeric","Mass Error [ppm]","numeric","Mass Error [Da]","numeric","Uncalibrated Mass Error [ppm]","numeric","Uncalibrated Mass Error [Da]","numeric","Max intensity m/z 0","numeric","Retention time","numeric","Retention length","numeric","Calibrated retention time","numeric","Calibrated retention time start","numeric","Calibrated retention time finish","numeric","Retention time calibration","numeric","Match time difference","numeric","Match m/z difference","numeric","Match q-value","numeric","Match score","numeric","Number of data points","integer","Number of scans","integer","Number of isotopic peaks","integer","PIF","numeric","Fraction of total spectrum","numeric","Base peak fraction","numeric","PEP","numeric","MS/MS Count","integer","MS/MS Scan Number","integer","Score","numeric","Delta score","numeric","Combinatorics","integer","Intensity","numeric","Reverse","character","Potential contaminant","character","id","integer","Protein group IDs","character","Peptide ID","integer","Mod. peptide ID","integer","MS/MS IDs","character","Best MS/MS","integer","AIF MS/MS IDs","logical","Oxidation (M) site IDs","character", "Leading Proteins", "character", "Contaminant", "character"), ncol=2, byrow=T), stringsAsFactors = F)
  # merge the classes to the columns
  cols.matched = merge(cols, col.classes, by="V1", all.x=T)
  # re-order things to match the initial order
  cols.matched <- cols.matched[order(cols.matched$idx),]
  
  # Stop if there is an issue
  if(length(which(is.na(cols.matched$V2)))>0){
    stop(paste0("OH NO!! YOUR EVIDENCE FILE CONTAINS A COLUMN THAT I DON'T RECOGNIZE :( PLEASE TELL THE 'col.classes' IN THE read_evidence_file' FUNCTION AND ADD IN THIS NEW COLUMN(S) CALLED \n\t", paste(cols.matched$V1[which(is.na(cols.matched$V2))], collapse="\n\t"), "\n" ) )
  }
  
  # read in the evidence file with their classes
  x <- fread(evidence_file, integer64 = 'double', colClasses = cols.matched$V2)
  return(x)
}

removeMaxQProteinGroups <- function(data){
  data_selected = data[grep(";",data$Proteins, invert=T),]
  return(data_selected)
}

sampleCorrelationHeatmap <- function (data_w, keys, config) {
  mat = log2(data_w[,4:ncol(data_w),with=F])
  mat[is.na(mat)]=0
  mat_cor = cor(mat, method = 'pearson', use = 'everything')
  ordered_keys = keys[with(keys, order(RawFile)),] ## we want to make informarive row names so order by RawFile because that's how data_w is ordered
  mat_names = paste(ordered_keys$Condition, ordered_keys$BioReplicate, ordered_keys$Run)
  colnames(mat_cor) = mat_names
  rownames(mat_cor) = mat_names
  colors_pos = colorRampPalette(brewer.pal("Blues",n=5))(10)
  colors_neg = rev(colorRampPalette(brewer.pal("Reds",n=5))(10))
  colors_tot = c(colors_neg, colors_pos)
  pheatmap(mat = mat_cor, cellwidth = 10, cellheight = 10, scale = 'none', filename = gsub('.txt','-heatmap.pdf',config$files$output), color = colors_tot, breaks = seq(from=-1,to = 1, by=.1), fontfamily="mono")
  
}

samplePeptideBarplot <- function(data_f, config){
  # set up data into ggplot compatible format
  data_f = data.table(data_f, labels=paste(data_f$RawFile, data_f$Condition, data_f$BioReplicate))
  data_f = data_f[with(data_f, order(labels,decreasing = T)),]
  
  # plot the peptide counts for all the samples TOGETHER
  p = ggplot(data = data_f, aes(x=labels))
  p = p + geom_bar() + theme(axis.text.x = element_text(angle = 90, hjust = 1, family = 'mono')) + ggtitle('Unique peptides per run\n after filtering') + coord_flip()
  ggsave(filename = gsub('.txt','-peptidecounts.pdf',config$files$output), plot=p, width = 8, height = 10)
  
  w = 10
  h = ceiling( (7/5+2) * ceiling(length(unique(data_f$Condition))/5) )
  # plot the peptide counts for all the samples PER BAIT
  p = ggplot(data = data_f, aes(x=as.factor(BioReplicate)))
  p = p + geom_bar() + theme(axis.text.x = element_text(angle = 90, hjust = 1, family = 'mono')) + ggtitle('Unique peptides per run\n after filtering') + facet_wrap(~Condition, scales='free', ncol=5)  + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename = gsub('.txt','-peptidecounts-perBait.pdf',config$files$output), plot=p, width = w, height = h)
  
}

significantHits <- function(mss_results, labels='*', LFC=c(-2,2), FDR=0.05){
  ## get subset based on labels
  selected_results = mss_results[grep(labels,mss_results$Label), ]
  cat(sprintf('\tAVAILABLE LABELS FOR HEATMAP:\t%s\n',paste(unique(mss_results$Label), collapse=',')))
  cat(sprintf('\tSELECTED LABELS FOR HEATMAP:\t%s\n',paste(unique(selected_results$Label), collapse=',')))
  significant_proteins = selected_results[(!is.na(selected_results$log2FC) & selected_results$adj.pvalue <= FDR & (selected_results$log2FC >= LFC[2] | selected_results$log2FC <= LFC[1])) , 'Protein']
  significant_results = selected_results[selected_results$Protein %in% significant_proteins, ]
  return(significant_results)
}

simplifyAggregate <- function(str, sep=',', numeric=F){
  str_vec = unlist(str_split(str, pattern = sep))
  if(numeric){
    str_vec = sort(unique(as.numeric(str_vec)))
  }else{
    str_vec = sort(unique(str_vec))  
  }
  
  str_new = str_join(as.character(str_vec), collapse = ';')
  return(str_new)
}

simplifyOutput <- function(input){
  input$Protein = apply(input, 1, function(x) simplifyAggregate(unname(x['Protein']), sep = ';'))
  if(any(grepl('mod_sites',colnames(input)))){
    input$mod_sites = apply(input, 1, function(x) simplifyAggregate(unname(x['mod_sites'])))  
  }
  if(any(grepl('uniprot_ac',colnames(input)))){
    input$uniprot_ac = apply(input, 1, function(x) simplifyAggregate(unname(x['uniprot_ac'])))
  }
  if(any(grepl('entrezgene',colnames(input)))){
    input$entrezgene = apply(input, 1, function(x) simplifyAggregate(unname(x['entrezgene']), numeric = T))
  }
  if(any(grepl('uniprot_genename',colnames(input)))){
    input$uniprot_genename = apply(input, 1, function(x) simplifyAggregate(unname(x['uniprot_genename'])))
  }
  if(any(grepl('description',colnames(input)))){
    input$description = apply(input, 1, function(x) simplifyAggregate(unname(x['description'])))
  }
  input[,sample_1:=gsub('([A-Z,0-9,a-z,_,\\s]+)\\-{1}([A-Z,0-9,a-z,_,\\s]+)','\\1',input$Label)]
  input[,sample_2:=gsub('([A-Z,0-9,a-z,_,\\s]+)\\-{1}([A-Z,0-9,a-z,_,\\s]+)','\\2',input$Label)]
  input[,Label:=NULL]
  if(any(grepl('group_id',colnames(input)))){
    input[,group_id:=NULL]
  }
  return(input)
}


volcanoPlot <- function(mss_results_sel, lfc_upper, lfc_lower, FDR, file_name='', PDF=T, decimal_threshold=16){
  
  # handle cases where log2FC is Inf. There are no pvalues or other information for these cases :(
  # Issues with extreme_val later if we have Inf/-Inf values.
  if( sum(is.infinite(mss_results_sel$log2FC)) > 0 ){
    idx <- is.infinite(mss_results_sel$log2FC)
    mss_results_sel$log2FC[ idx ] <- NA
  }
  
  min_x = -ceiling(max(abs(mss_results_sel$log2FC), na.rm=T))
  max_x = ceiling(max(abs(mss_results_sel$log2FC), na.rm=T))
  # Deal with special cases in the data where we have pvalues = Inf,NA,0
  if( sum(is.na(mss_results_sel$adj.pvalue))>0 ) mss_results_sel <- mss_results_sel[!is.na(mss_results_sel$adj.pvalue),]
  if(nrow(mss_results_sel[mss_results_sel$adj.pvalue == 0 | mss_results_sel$adj.pvalue == -Inf,]) > 0) mss_results_sel[!is.na(mss_results_sel$adj.pvalue) & (mss_results_sel$adj.pvalue == 0 | mss_results_sel$adj.pvalue == -Inf),]$adj.pvalue = 10^-decimal_threshold
  max_y = ceiling(-log10(min(mss_results_sel[mss_results_sel$adj.pvalue > 0,]$adj.pvalue, na.rm=T))) + 1
  
  l = length(unique(mss_results_sel$Label))
  w_base = 7
  h_base = 7
  
  if(l<=2){
    w=w_base*l 
  }else{
    w=w_base*2
  }
  h = h_base*ceiling(l/2)
  
  if(PDF) pdf(file_name, width=w, height=h)
  p = ggplot(mss_results_sel, aes(x=log2FC,y=-log10(adj.pvalue)))
  print(p + geom_point(colour='grey') + 
          geom_point(data = mss_results_sel[mss_results_sel$adj.pvalue <= FDR & mss_results_sel$log2FC>=lfc_upper,], aes(x=log2FC,y=-log10(adj.pvalue)), colour='red', size=2) +
          geom_point(data = mss_results_sel[mss_results_sel$adj.pvalue <= FDR & mss_results_sel$log2FC<=lfc_lower,], aes(x=log2FC,y=-log10(adj.pvalue)), colour='blue', size=2) +
          geom_vline(xintercept=c(lfc_lower,lfc_upper), lty='dashed') + 
          geom_hline(yintercept=-log10(FDR), lty='dashed') + 
          xlim(min_x,max_x) + 
          ylim(0,max_y) + 
          facet_wrap(facets = ~Label, ncol = 2, scales = 'fixed')) 
  if(PDF) dev.off()
}


#' @title Generate the contrast matrix required by MSstats from a txt file
#' @description It simplifies the process of creating the contrast file
#' @param contrast_file The text filepath of contrasts
#' @keywords 
#' writeContrast()
#' @export
writeContrast <- function(contrast_file) {
  input_contrasts <- readLines(contrast_file, warn=F)
  
  # check if contrast_file is old-format (i.e the contrast_file is a matrix)
  headers <- unlist(strsplit(input_contrasts[1], split = "\t"))
  if (length(headers)>1) {
    newinput_contrasts <- c()
    for (i in 2:length(input_contrasts)){
      newinput_contrasts <- c(newinput_contrasts, unlist(strsplit(input_contrasts[i], split = "\t"))[1])
    }
    input_contrasts <- newinput_contrasts
  }
  
  # validate the input
  input_contrasts <- trimws(input_contrasts)
  valid <- TRUE
  accepted_chars <- c(LETTERS, letters, 0:9, '-','_')
  for (x in input_contrasts) {
    characs <- unlist(strsplit(x, split=''))
    not_allowed_count <- length(which(!(characs %in% accepted_chars)))
    if (not_allowed_count > 0) {
      valid <- FALSE
      cat(paste(x, "is not a valid input"))
      break
    }
    
    dash_count <- length(which(characs == '-'))
    if (dash_count != 1) {
      valid <- FALSE
      cat(paste(x, "needs to contain exactly 1 '-'"))
      break    
    }
  }
  
  if (valid) {
    mat <- t(as.data.frame(strsplit(input_contrasts, split = '-')))
    rownames(mat) <- NULL
    conds <- sort(unique(c(mat[,1], mat[,2])))
    contrast_matrix <- matrix(0, nrow = nrow(mat), ncol = length(conds))
    colnames(contrast_matrix) <- conds
    rownames(contrast_matrix) <- input_contrasts
    
    for (i in 1:nrow(mat)) {
      cond1 <- mat[i,1]
      cond2 <- mat[i,2]
      contrast_matrix[i, cond1] <- 1
      contrast_matrix[i, cond2] <- -1
    }
    return (contrast_matrix)
  }
  return (NA)
}