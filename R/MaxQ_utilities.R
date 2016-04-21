#! /usr/bin/Rscript --vanilla

suppressMessages(library(data.table))
suppressWarnings(library(seqinr))
suppressMessages(library(stringr))
suppressMessages(require(bit64))
suppressMessages(library(getopt))
suppressMessages(library(reshape2))
suppressMessages(library(biomaRt))
suppressMessages(library(limma))

###############################
## FILE AND LIB LOADING #######

#########################
## CONFIG LOADING #######

ALLOWED_COMMANDS = c('concat','convert-silac','keys','convert-sites','annotate','results-wide','mapback-sites','heatmap','simplify','saint-format','data-plots','spectral-counts', 'mist')

spec = matrix(c(
  'verbose', 'v', 2, "integer", "",
  'help'   , 'h', 0, "logical", "available arguments (this screen)",
  'command'  , 'c', 1, "character", sprintf("command to run. Currently supported commands: %s",paste(ALLOWED_COMMANDS,collapse=',')),
  'files'  , 'f', 1, "character", "files to feed to command. accepts regexp but needs to be quoted",
  'output'  , 'o', 1, "character", "Output file",
  'keys','k', 1, "character", "keys file",
  'proteome'  , 'p', 1, "character", "Reference Proteome FASTA file",
  'mapping'  , 'm', 1, "character", "mapping file produced by convert-sites",
  'species'  , 's', 1, "character", "species to use look for protein masses from uniprot from. Can list multiple species, separate with a '-'. Current species= HUMAN, MOUSE. Used with 'mist' option",
  'uniprot_dir','d', 1, "character", "Location of the uniprot files are located. Used with 'mist' option.",
  'mod_type','t', 1, "character", "Modification type: ub|ph.",
  'lfc_lower','l', 1, "double", "Lower Log2FC to include for heatmap plotting",
  'lfc_upper','u', 1, "double", "Upper Log2FC to include for heatmap plotting",
  'q_value','q', 1, "double", "q-value (FDR) to include for heatmap plotting",
  'labels','', 1, "character", "labels to include for heatmap plotting",
  'identifier_column','i', 1, "character", "Identifier column"),
  byrow=T, ncol=5)

opt = getopt(spec = spec, opt = commandArgs(TRUE), command = get_Rscript_filename(), usage = FALSE, debug = FALSE)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

# warnings for when arguments are forgotten
if( is.null(opt$command) ) { cat("NO COMMAND (-c) ENTERED!!\n") }
if( is.null(opt$files) ) { cat("NO FILE (-f) SEPCIFIED!!\n") }
if( is.null(opt$output) ) { cat("NO OUTPUT FILE (-o) SEPCIFIED!!\n") }
# set default values for arguments
if( is.null(opt$mod_type) ) { opt$mod_type = 'ub' }
if( is.null(opt$lfc_lower) ) { opt$lfc_lower = '-2' }
if( is.null(opt$lfc_upper) ) { opt$lfc_upper = '2' }
if( is.null(opt$q_value) ) { opt$q_value = '0.05' }
if( is.null(opt$labels) ) { opt$labels = '*' }

loadLibs = function(){
  cat(">> LOADING EXTERNAL FILES AND LIBRARIES\n")
  # set source directory
  args <- commandArgs(trailingOnly = F) 
  scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
  
  ## load all external files
  if(length(scriptPath)==0){
    source("R/MSstats_functions.R")
  }else{
    source(paste(scriptPath,"/MSstats_functions.R",sep=""))  
  }
}

MQutil.SILACToLong = function(filename, output){
  library(data.table)
  library(reshape2)
  file = Sys.glob(filename)
  cat(sprintf('\tPROCESSING:\n\t%s\n',paste(file,collapse='\n\t')))
  tmp = fread(file, integer64 = 'double')
  tmp_long = reshape2::melt(tmp, measure.vars = c('Intensity L','Intensity H'))
  tmp_long[,Intensity:=NULL]
  setnames(tmp_long,'value','Intensity')
  setnames(tmp_long,'variable','IsotopeLabelType')
  setnames(tmp_long,'Raw file','Raw.file')
  levels(tmp_long$IsotopeLabelType) = c('L','H')
  tmp_long[!is.na(tmp_long$Intensity) && tmp_long$Intensity<1,]$Intensity=NA
  write.table(tmp_long, file=output, sep='\t', quote=F, row.names=F, col.names=T)
}

MQutil.concat = function(filenames, output){
  library(data.table)
  files = Sys.glob(filenames)
  cat(sprintf('\tPROCESSING:\n\t%s\n',paste(files,collapse='\n\t')))
  
  res = NULL
  unique_files = c()
  
  for(f in 1:length(files)){
    file = files[f]
    tmp = fread(file, stringsAsFactors=F, colClasses = c(Intensity='character'), integer64 = 'double')
    tmp$Intensity = as.numeric(tmp$Intensity)
    #tmp$Intensity L = as.numeric(tmp[,'Intensity L',with=F])
    
    unique_files_current = unique(tmp[['Raw file']])
    if(!is.null(intersect(unique_files_current,unique_files)) && length(intersect(unique_files_current,unique_files))>0) cat(sprintf('\tWARNING DUPLICATE RAW FILE ENTRIES IN FILE %s:\t%s\n',file, paste(intersect(unique_files_current, unique_files),collapse=',')))
    select_colnames = grep('Raw\ file|Intensity|Proteins|Modifications|Sequence|Modified\ sequence|Charge|Protein\ group\ IDs|id|Retention\ time|Reverse|Contaminant',colnames(tmp), ignore.case = F)
    tmp = tmp[,select_colnames,with=F]
    res = rbind(res, tmp, fill=T)
    unique_files = c(unique_files, unique_files_current)  
  }
  
  cat(sprintf('\tWRITING\t%s\n',output))
  write.table(res, file=output, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
  cat(sprintf('\tWRITING\t%s\n',gsub('.txt','-keys.txt',output)))
  #print(unique_files)
  write.table(unique_files, file=gsub('.txt','-keys.txt',output), eol='\n', sep='\t', quote=F, row.names=F, col.names=F)
}

MQutil.getKeys = function(filename, output){
  library(data.table)
  file = Sys.glob(filename)
  cat(sprintf('\tPROCESSING:\n\t%s\n',paste(file,collapse='\n\t')))
  #tmp = data.table(read.delim(file, stringsAsFactors=F))
  tmp = fread(file, stringsAsFactors=F, integer64 = 'double')
  write.table(unique(tmp$Raw.file),file=output, eol='\n', sep='\t', quote=F, row.names=F, col.names=F)
  cat(sprintf('\tWRITTEN\t%s\n',output))
}


MQutil.ProteinToSiteConversion <- function (maxq_file, ref_proteome_file, output_file, mod_type='ub') {
  
  if(mod_type=='ub'){
    maxq_mod_residue='K\\(gl\\)'
    mod_residue = 'K'
  }else if(mod_type=='ph'){
    maxq_mod_residue='(S|T|Y)\\(ph\\)'  
    mod_residue = 'S|T|Y'
  }
  
  ## read in ref. proteome
  ref_proteome = read.fasta(file = ref_proteome_file, 
                            seqtype = "AA", as.string = T,
                            set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = FALSE)
  
  ######################
  ## make mod-site index
  
  p_seqs = c()
  p_names = c()
  p_annots = c()
  
  for(e in ref_proteome){
    p_seqs = c(p_seqs, e[1])
    p_names = c(p_names, attr(e,'name'))
    p_annots = c(p_annots, attr(e,'Annot'))
  }
  
  ref_table = data.table(names=p_names, annots=p_annots, seqs=p_seqs)
  ref_table[,uniprot_ac:=gsub('([a-z,0-9,A-Z]+\\|{1})([A-Z,0-9,\\_]+)(\\|[A-Z,a-z,0-9,_]+)','\\2',names)]
  
  indices = lapply(ref_table$seqs, function(x) as.vector(str_locate_all(x,pattern=mod_residue)[[1]][,1]))
  ptm_sites = lapply(ref_table$seqs, function(x) as.vector(str_match_all(x,pattern=mod_residue)))
  lengths = unlist(lapply(indices, FUN = length))
  keys = rep(ref_table$uniprot_ac, lengths)
  protein_indices = data.table(uniprot_ac=keys, ptm_site=unlist(ptm_sites), res_index = unlist(indices))
  
  #################################
  ## map mod sites in data to index 
  
  ## read in maxq. data
  maxq_data = fread(maxq_file, integer64 = 'double')
  maxq_data = maxq_data[grep("CON__|REV__",maxq_data$Proteins, invert=T),]
  unique_peptides_in_data = unique(maxq_data[,c('Proteins','Modified sequence'),with=F])
  setnames(unique_peptides_in_data,'Modified sequence','sequence')
  
  mod_sites = c()
  mod_seqs = c()
  
  for(i in 1:nrow(unique_peptides_in_data)){
    entry = unique_peptides_in_data[i,]
    peptide_seq = entry$sequence
    ## cleanup the sequence (removing all modifications) for matching the protein sequence
    peptide_seq_clean = gsub('[a-z,0-9,\\(,\\),_]','', peptide_seq)
    mod_sites_in_peptide = str_locate_all(string = peptide_seq, pattern = maxq_mod_residue)[[1]][,1]
    
    if(length(mod_sites_in_peptide)>0){
      uniprot_acs = entry$Proteins
      uniprot_acs = str_split(string = uniprot_acs, pattern = ';')[[1]]
      
      for(uac in uniprot_acs){
        protein_seq = ref_table[uniprot_ac==uac,]$seqs
        if(length(protein_seq)>0){
          ## get the position of the peptide in the protein sequence
          peptide_index_in_protein = str_locate(protein_seq, peptide_seq_clean)[[1]][1]
          
          for(m in 1:length(mod_sites_in_peptide)){
            mod_site = mod_sites_in_peptide[m]   
            peptide_seq_before_site = str_sub(peptide_seq, 1, mod_site-1)
            ## count all AA (not counting all modifications) before the modification to get the relative position of the modification in the peptide sequence
            residues_before_site = str_count(string = peptide_seq_before_site, pattern = 'A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y')
            mod_site_index_in_protein = peptide_index_in_protein+residues_before_site
            protein_mod_sites = protein_indices[uniprot_ac==uac,]
            if(!is.na(mod_site_index_in_protein)){
              #cat(sprintf('%s\n',mod_site_id))
              mod_res = protein_mod_sites[res_index==mod_site_index_in_protein,ptm_site]
              mod_site_id = sprintf('%s_%s%s', uac, str_sub(protein_seq,mod_site_index_in_protein,mod_site_index_in_protein), mod_site_index_in_protein)
              mod_sites = c(mod_sites, mod_site_id)
              mod_seqs = c(mod_seqs, peptide_seq)
              stopifnot(length(mod_sites)==length(mod_seqs))
            }else{
              cat(sprintf('MISMATCH\t%s\n\tPEPTIDE_SEQ\t%s\n\tMOD_SITE\t%s\n\tPEPTIDE_IDX_IN_PROTEIN\t%s\n\tRESIDUES_BEFORE_SITE\t%s\n\tPROTEIN_SEQ\t%s\n',mod_site_id, peptide_seq, mod_site, peptide_index_in_protein, residues_before_site, protein_seq))
            } 
          }
        }
      }  
    }
  }
  
  mod_site_mapping = data.table(mod_sites, mod_seqs)
  mod_site_mapping_agg = aggregate(mod_sites ~ mod_seqs, mod_site_mapping, FUN=function(x)paste(x,collapse=','))
  
  setnames(maxq_data,'Modified sequence','mod_seqs')
  unmapped_mod_seqs = maxq_data[!(mod_seqs %in% mod_site_mapping_agg$mod_seqs) & grepl('(gl)',mod_seqs) & !grepl('REV__|CON__',Proteins),]
  unmapped_mod_seqs = unique(unmapped_mod_seqs[,c('mod_seqs','Proteins'),with=F])
  cat('UNABLE TO MAP\n')
  print(unmapped_mod_seqs)
  
  final_data = merge(maxq_data, mod_site_mapping_agg, by='mod_seqs')
  setnames(final_data,c('Proteins','mod_sites','mod_seqs'),c('Proteins_ref','Proteins','Modified sequence'))
  write.table(final_data, file = output_file, eol='\n', sep='\t',quote=F, row.names=F, col.names=T)
  
  ## write a mapping table
  protein_seq_mapping = unique(maxq_data[,c('Proteins','mod_seqs'),with=F])
  setnames(protein_seq_mapping,'Proteins','Protein')
  mapping_table = merge(protein_seq_mapping, mod_site_mapping_agg, by='mod_seqs', all=T)
  write.table(mapping_table, file=gsub('.txt','-mapping.txt',output_file), eol='\n', sep='\t',quote=F, row.names=F, col.names=T)
}


MQutil.annotate = function(input_file=opt$input, output_file=opt$output, uniprot_ac_col='Protein', group_sep=';', uniprot_dir = '~/github/kroganlab/source/db/', species='HUMAN'){
  cat(">> ANNOTATING\n")
  results = read.delim(input_file, stringsAsFactors=F, sep='\t')
  
  # remove unnamed proteins that are listed as ""
  results <- results[-which(results$Protein==""),]
  
  # read in all the annotation files from the uniprot_db directory
  species_split = unlist(strsplit(species, "-"))
  Uniprot = NULL
  for(org in species_split){
    cat(sprintf("\tLOADING %s\n",org))
    tmp = read.delim2(sprintf("%s/uniprot_protein_descriptions_%s.txt",uniprot_dir,org), stringsAsFactors=F, quote="")    
    if(is.null(Uniprot)){
      Uniprot = as.data.frame(tmp)
    }else{
      Uniprot = rbind(Uniprot, tmp)  
    }
  }
  
  # get list of all unique prey entries in this file. Keep 'group_sep' in mind.
  preys <- unique(results$Protein)
  preys <- preys.original <- data.frame(prey = preys, idx = 1:length(preys), stringsAsFactors=F)
  # split apart all the preys and index them so we can piece them back together when merging
  preys <- do.call(rbind, apply(preys, 1, function(y){ data.frame(prey = unlist(strsplit(y[1], ";")), idx = as.numeric(y[2]), stringsAsFactors = F) } ))
  
  # annotate all the preys wiht the Uniprot info
  preys <- merge( preys, Uniprot[,c("Entry","Entry.name","Protein.names","Gene.names")], by.x ="prey", by.y="Entry", all.x=T)
  # aggregate all the preys on the indexes so we can merge with the original data
  
  # merge protein name
  tmp <- aggregate(data = preys[,c("prey",'idx','Entry.name')], .~idx, paste, collapse=";")
  names(tmp) <- c('idx','uniprot_ac','Protein_name')
  preys.new <- merge(preys.original, tmp, by='idx', all.x=T)
  
  # merge protein description
  tmp <- aggregate(data = preys[,c('idx','Protein.names')], .~idx, paste, collapse=";")
  names(tmp) <- c('idx', 'Protein_desc')
  preys.new <- merge(preys.new, tmp, by='idx', all.x=T)
  
  # merge protein description
  preys$Gene.names <- gsub(" .*","", preys$Gene.names)
  tmp <- aggregate(data = preys[,c('idx','Gene.names')], .~idx, paste, collapse=";")
  names(tmp) <- c('idx','Gene.names')
  preys.new <- merge(preys.new, tmp, by='idx',all.x=T)
  
  # merge the annotations all back into the original data
  results_out <- merge(results, preys.new, by.x="Protein", by.y="prey", all.x=T)
  results_out$idx = c()
  
  # alert user of any unmapped proteins
  unmapped = unique(results_out[is.na(results_out$uniprot_ac),"Protein"]) 
  cat('UNMAPPED PROTEINS\t', length(unmapped), '\n')
  cat('\t',paste(unmapped,collapse='\n\t'),'\n')
  write.table(results_out, file=output_file, sep='\t', quote=F, row.names=F, col.names=T)
  cat(">> ANNOTATING COMPLETE!\n")
}


MQutil.resultsWide = function(input_file, output_file){
  input = fread(input_file, integer64 = 'double')
  input_l = melt(data = input[,c('Protein', 'Label','log2FC','adj.pvalue'),with=F],id.vars=c('Protein', 'Label'))
  
  ## then cast to get combinations of LFCV/PVAl and Label as columns
  input_w = dcast.data.table( Protein ~ Label+variable, data=input_l, value.var=c('value'))
  write.table(input_w, file=output_file, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
}

MQutil.mapSitesBack = function(input_file, mapping_file, output_file){
  input = fread(input_file, integer64 = 'double')
  setnames(input,'Protein','mod_sites')
  mapping = fread(mapping_file, integer64 = 'double')
  mapping = unique(mapping[!is.na(mod_sites),c('Protein','mod_sites'),with=F])
  mapping = aggregate(Protein ~ mod_sites, data=mapping, FUN=function(x)paste(x,collapse=';'))
  out = merge(input, mapping, by='mod_sites', all.x=T)
  write.table(out[,c(ncol(out),1:(ncol(out)-1)),with=F], file=output_file, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
}


MQutil.plotHeat = function(input_file, output_file, labels=NULL, names='Protein', cluster_cols=F, display='log2FC', row_names='Protein', col_names='Label'){
  heat_data = fread(input_file, integer64 = 'double')
  #heat_data = mss_F[,c('uniprot_id','Label','log2FC')]
  
  ## good
  if(display=='log2FC'){
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
    pheatmap(heat_data_w, scale="none", cellheight=10, cellwidth=10, filename =out_file, color=colors_tot, breaks=seq(from=-extreme_val, to=extreme_val, by=bin_size), cluster_cols=cluster_cols, fontfamily="mono")  
  }else{
    heat_data_w = heat_data_w[,labelOrder]
    pheatmap(heat_data_w, scale="none", cellheight=10, cellwidth=10, filename=out_file, color=colors_tot, breaks=seq(from=-extreme_val, to=extreme_val, by=bin_size), cluster_cols=cluster_cols, fontfamily="mono")
  }
  
  heat_data_w
}

## still need to make fileType, conditions, cluster_cols, display configurable through command line
MQutil.plotHeatmap = function(input_file, output_file, fileType='l', labels='*', cluster_cols=F, display='log2FC', lfc_lower=-2, lfc_upper=2, FDR=0.05){
  ## read input
  input = read.delim(input_file, stringsAsFactors = F)
  
  ## select data points  by LFC & FDR criterium in single condition and adding corresponding data points from the other conditions
  sign_hits = significantHits(input,labels=labels,LFC=c(lfc_lower,lfc_upper),FDR=FDR)
  sign_labels = unique(sign_hits$Label)
  cat(sprintf("\tSELECTED HITS FOR PLOTS WITH LFC BETWEEN %s AND %s AT %s FDR:\t%s\n",lfc_lower, lfc_upper, FDR, nrow(sign_hits)/length(sign_labels))) 
  
  ## REPRESENTING RESULTS AS HEATMAP
  ## plot heat map for all contrasts
  if(any(grepl('uniprot_genename',colnames(sign_hits)))){
    heat_labels = paste(sign_hits$Protein,sign_hits$uniprot_genename,sep=' ')  
  }else{
    heat_labels = sign_hits$Protein
  }
  
  heat_labels = gsub('\\sNA$','',heat_labels)
  heat_data_w = plotHeat(mss_F=sign_hits, out_file=output_file, names=heat_labels, cluster_cols=cluster_cols, display=display)  
}

MQutil.simplify = function(input_file, output_file){
  input = fread(input_file, integer64 = 'double')
  output = simplifyOutput(input)
  write.table(output, file=output_file, sep='\t', quote=F, row.names=F, col.names=T, eol = '\n')
}

MQutil.MaxQToSaint = function(data_file, keys_file, ref_proteome_file, quant_variable='spectral_count'){
  cat(">> CONVERTING TO SAINT FORMAT\n")
  
  data = fread(data_file, integer64 = 'double')
  keys = fread(keys_file, integer64 = 'double')
  
  ## write baits in format
  ## hIP101-10       PV_2C_co_uni    T
  
  saint_baits = keys[,c('BioReplicate','Condition','SAINT'),with=F]
  
  ## write interactions in format
  ## hIP101-10       PV_2C_co_uni    Q9NTG7  1
  
  tryCatch(setnames(data, 'Raw file', 'RawFile'), error=function(e) cat('Raw.file not found\n'))
  tryCatch(setnames(keys, 'Raw.file', 'RawFile'), error=function(e) cat('Raw.file not found\n'))
  
  cat('\tVERIFYING DATA AND KEYS\n')
  if(any(!c('RawFile','IsotopeLabelType','Condition','BioReplicate','Run','SAINT') %in% colnames(keys))){
    stop('COLNAMES IN KEYS NOT CONFORM TO SCHEMA\n\tRawFile\tIsotopeLabelType\tCondition\tBioReplicate\tRun\tSAINT\n')
  } 
  if(!'IsotopeLabelType' %in% colnames(data)) data[,IsotopeLabelType:='L']
  data = mergeMaxQDataWithKeys(data, keys, by = c('RawFile','IsotopeLabelType'))
  data_f = filterMaxqData(data)
  data_f = removeMaxQProteinGroups(data_f) ## do we want this or not?
  
  cat("\tAGGREGATING ON", quant_variable,"VALUES...\n")
  ## aggregate over technical replicates if necessary
  if(quant_variable=='spectral_count'){
    setnames(data_f,'MS/MS Count','spectral_counts')
    data_f_agg = aggregate(spectral_counts ~ BioReplicate+Condition+Proteins+Sequence+Charge,data=data_f,FUN = max)
    data_f_agg = aggregate(spectral_counts ~ BioReplicate+Condition+Proteins,data=data_f_agg,FUN = sum)
  }else if(quant_variable=='ms1'){
    data_f_agg = aggregate(Intensity ~ BioReplicate+Condition+Proteins+Sequence+Charge,data=data_f,FUN = max)
    data_f_agg = aggregate(Intensity ~ BioReplicate+Condition+Proteins,data=data_f_agg,FUN = sum)
  }else{
    stop("ERROR!! Wrong value for variable to quantify. Please use 'spectral_count' or 'ms1'.")
  }
  
  ## IP name, bait name, prey name, and spectral counts or intensity values
  saint_interactions = data_f_agg
  
  ## write preys in format
  ## Q9NTG7  43573.5 Q9NTG7
  
  ref_proteome = read.fasta(file = ref_proteome_file, 
                            seqtype = "AA", as.string = T,
                            set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = FALSE)
  p_lengths = c()
  p_names = c()
  for(e in ref_proteome){
    p_lengths = c(p_lengths, nchar(e[1]))
    p_names = c(p_names, attr(e,'name'))
  }
  ref_table = data.table(names=p_names, lengths=p_lengths)
  ref_table[,uniprot_ac:=gsub('([a-z,0-9,A-Z]+\\|{1})([A-Z,0-9,\\_]+)(\\|[A-Z,a-z,0-9,_]+)','\\2',names)]
  ref_table[,uniprot_id:=gsub('([a-z,0-9,A-Z]+\\|{1})([a-z,0-9,A-Z]+\\|{1})([A-Z,a-z,0-9,_]+)','\\3',names)]
  
  unique_preys = data.table(uniprot_ac=unique(data_f_agg$Proteins))
  saint_preys = ref_table[,c('uniprot_ac','lengths','uniprot_id'),with=F]
  saint_preys = merge(unique_preys, saint_preys, by='uniprot_ac', all.x=T)
  missing_lengths = nrow(saint_preys[is.na(saint_preys$uniprot_id),])
  saint_preys[is.na(saint_preys$uniprot_id),]$uniprot_id=saint_preys[is.na(saint_preys$uniprot_id),]$uniprot_ac  
  if(missing_lengths>0){
    cat(sprintf("\tWARNING!\tCOMPUTING %s MISSING LENGTHS WITH THE MEDIAN LENGTH FROM THE DATASET\n",missing_lengths))
    saint_preys[is.na(saint_preys$lengths),]$lengths=median(saint_preys$lengths,na.rm = T)
  }
  
  ## WRITE
  write.table(saint_baits,file = gsub('.txt','-saint-baits.txt',keys_file), eol='\n',sep='\t', quote=F, row.names=F, col.names=F)
  write.table(saint_preys,file = gsub('.txt','-saint-preys.txt',keys_file), eol='\n',sep='\t', quote=F, row.names=F, col.names=F)
  write.table(saint_interactions,file = gsub('.txt','-saint-interactions.txt',keys_file), eol='\n',sep='\t', quote=F, row.names=F, col.names=F)
}

MQutil.dataPlots = function(input_file, output_file){
  
  data_mss = fread(input_file, integer64 = 'double')
  unique_subjects = unique(data_mss$PROTEIN)
  condition_length = length(unique(data_mss$GROUP_ORIGINAL))
  min_abu = min(data_mss$ABUNDANCE, na.rm = T)
  max_abu = max(data_mss$ABUNDANCE, na.rm=T)
  
  pdf(output_file, width = condition_length*1.5, height = 3)
  
  cat('PRINTING CONDITION PLOTS\n')
  for(subject in unique_subjects){
    subject_data = data_mss[PROTEIN==subject,]
    cat(sprintf('\t%s\n',subject))
    p = ggplot(data = subject_data, aes(x=SUBJECT_ORIGINAL,y=ABUNDANCE, colour=FEATURE))
    p = p + geom_point(size=2) + 
      facet_wrap(facets = ~ GROUP_ORIGINAL, drop = T, scales = 'free_x', ncol = condition_length) + 
      ylim(min_abu,max_abu) +
      theme(axis.text.x=element_text(angle=-90,hjust=1)) +
      guides(colour=FALSE) +
      xlab(NULL) +
      ggtitle(subject)
    print(p)
  }
  dev.off()
}

MQutil.spectralCounts = function(input_file, keys_file, output_file){
  data = fread(input_file, integer64 = 'double')
  keys = fread(keys_file, integer64 = 'double')
  
  tryCatch(setnames(data, 'Raw file', 'RawFile'), error=function(e) cat('Raw.file not found\n'))
  tryCatch(setnames(keys, 'Raw.file', 'RawFile'), error=function(e) cat('Raw.file not found\n'))
  
  cat('\tVERIFYING DATA AND KEYS\n')
  if(!'IsotopeLabelType' %in% colnames(data)) data[,IsotopeLabelType:='L']
  data = mergeMaxQDataWithKeys(data, keys, by = c('RawFile','IsotopeLabelType'))
  data_sel = data[,c('Proteins','Condition','BioReplicate','Run','MS/MS Count'),with=F]
  setnames(data_sel,'MS/MS Count','spectral_counts')
  data_sel = aggregate( spectral_counts ~ Proteins+Condition+BioReplicate+Run, data=data_sel, FUN = sum)
  data_sel = data.frame(data_sel, bait_name=paste(data_sel$Condition, data_sel$BioReplicate, data_sel$Run, sep='_'))
  write.table(data_sel[,c('bait_name','Proteins','spectral_counts')], file=output_file, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
}



# Convert MaxQuant file into a Protein Prospector like format to run through the Mist pipeline
MQutil.MISTformat = function(input_file, keys_file, output_file, species="HUMAN", uniprot_dir='~/github/kroganlab/source/db/'){
  cat('\tREADING IN DATA AND KEYS\n')
  #data = fread(input_file, integer64 = 'double')
  data <- data.table(read.delim(input_file, stringsAsFactors=F))
  keys = data.table(read.delim(keys_file, stringsAsFactors = F))
  
  tryCatch(setnames(data, 'Raw file', 'RawFile'), error=function(e) cat('Raw.file not found\n'))
  tryCatch(setnames(data, 'Raw.file', 'RawFile'), error=function(e) cat('Raw.file not found\n'))
  tryCatch(setnames(keys, 'Raw file', 'RawFile'), error=function(e) cat('Raw.file not found\n'))
  tryCatch(setnames(keys, 'Raw.file', 'RawFile'), error=function(e) cat('Raw.file not found\n'))
  tryCatch(setnames(data,'MS/MS Count','spectral_counts'), error=function(e) cat('MS/MS Count column not found\n'))
  tryCatch(setnames(data,'MS.MS.Count','spectral_counts'), error=function(e) cat('MS/MS Count column not found\n'))
  
  cat('\tVERIFYING DATA AND KEYS\n')
  if(!'IsotopeLabelType' %in% colnames(data)) data[,IsotopeLabelType:='L']
  data = mergeMaxQDataWithKeys(data, keys, by = c('RawFile','IsotopeLabelType'))
  data_sel = data[,c('Proteins','Condition','BioReplicate','Run','RawFile','spectral_counts'),with=F]
  
  data_sel = aggregate( spectral_counts ~ Proteins+Condition+BioReplicate+Run+RawFile, data=data_sel, FUN = sum)
  data_sel = data.frame(data_sel, bait_name=paste(data_sel$Condition, data_sel$BioReplicate, data_sel$Run, sep='_'))
  
  # clean up proteins & annotate
  #~~~~~~~~~~~~~~~~~~~~~~~
  # remove CON's
  if( length(grep("^CON__",data_sel$Proteins))>0 ) data_sel = data_sel[-grep("^CON__",data_sel$Proteins),]
  # remove the party sets
  if( length(grep(";",data_sel$Proteins))>0 ) data_sel = data_sel[-grep(";",data_sel$Proteins),]     # NOTE!!! We lose a lot of entries this way... :\
  # keep only uniprot id
  #data_sel$uniprot_id = gsub("^.*\\|","", data_sel$Proteins)
  data_sel$Proteins = gsub("(^.*\\|)([A-Z0-9]+)(\\|.*$)","\\2",data_sel$Proteins)
  # remove blank protein names
  if(data_sel$Proteins == ""){ data_sel <- data_sel[-which(data_sel$Proteins == ""),]}
  data_sel$ms_unique_pep = ""
  # re-order
  data_sel <- data_sel[,c("RawFile",'Proteins','ms_unique_pep', 'spectral_counts')]
  names(data_sel) = c('id','ms_uniprot_ac','ms_unique_pep','ms_spectral_counts')
  
  # remove interactions with spectral_counts=0
  data_sel <- data_sel[which(data_sel$ms_spectral_counts==0),]
  
  # annotate proteins and add Masses for Mist
  species_split = unlist(strsplit(species, "-"))
  Uniprot = NULL
  for(org in species_split){
    cat(sprintf("\tLOADING %s\n",org))
    tmp = read.delim2(sprintf("%s/uniprot_protein_descriptions_%s.txt",uniprot_dir,org), stringsAsFactors=F, quote="")    
    if(is.null(Uniprot)){
      Uniprot = as.data.frame(tmp)
    }else{
      Uniprot = rbind(Uniprot, tmp)  
    }
  }
  results_annotated = merge(data_sel, Uniprot, all.x=T, by.x='ms_uniprot_ac',by.y='Entry')
  results_annotated$Mass <- as.numeric(gsub(",","",results_annotated$Mass))
  
  write.table(results_annotated, file=output_file, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
  cat('MIST FILES CREATED! Have a nice day :)\n')
}



main <- function(opt){
  if(opt$command %in% ALLOWED_COMMANDS){
    cat(sprintf('>> EXECUTING:\t%s\n',opt$command))
    loadLibs()
    if(opt$command == 'concat'){
      MQutil.concat(filenames=opt$files, output = opt$output)
    }else if(opt$command == 'convert-silac'){
      MQutil.SILACToLong(filename = opt$files, output = opt$output)
    }else if(opt$command == 'keys'){
      MQutil.getKeys(filename = opt$files, output = opt$output)
    }else if(opt$command == 'convert-sites'){
      MQutil.ProteinToSiteConversion (maxq_file = opt$files, output_file = opt$output, ref_proteome_file = opt$proteome, mod_type = opt$mod_type)
    }else if(opt$command == 'annotate'){
      MQutil.annotate(input_file = opt$files, output_file = opt$output, species=opt$species, uniprot_dir=opt$uniprot_dir)
    }else if(opt$command == 'results-wide'){
      MQutil.resultsWide(input_file = opt$files, output_file = opt$output )
    }else if(opt$command == 'mapback-sites'){
      MQutil.mapSitesBack(input_file = opt$files, output_file = opt$output , mapping_file = opt$mapping)
    }else if(opt$command == 'heatmap'){
      MQutil.plotHeatmap(input_file = opt$files, output_file = opt$output, labels = opt$labels, lfc_lower = opt$lfc_lower, lfc_upper=opt$lfc_upper, FDR = opt$q_value)
    }else if(opt$command == 'simplify'){
      MQutil.simplify(input_file = opt$files, output_file = opt$output)
    }else if(opt$command == 'saint-format'){
      MQutil.MaxQToSaint(data_file = opt$files, keys_file =  opt$keys, ref_proteome_file = opt$proteome, quant_variable = opt$identifier_column)
    }else if(opt$command == 'data-plots'){
      MQutil.dataPlots(input_file = opt$files, output_file = opt$output)
    }else if(opt$command == 'spectral-counts'){
      MQutil.spectralCounts(input_file = opt$files, keys_file =  opt$keys, output_file = opt$output)
    }else if(opt$command == 'mist'){
      MQutil.MISTformat(input_file = opt$files, keys_file =  opt$keys, output_file = opt$output, species=opt$species, uniprot_dir=opt$uniprot_dir)
    }
  }else{
    cat(sprintf('COMMAND NOT ALLOWED:\t%s\n',opt$command)) 
    cat(sprintf('ALLOWED COMMANDS:\t%s\n',paste(ALLOWED_COMMANDS,collapse=','))) 
  }
}

# opt$command = 'data-plots'
# opt$files = '~/Projects/HPCKrogan/Data/TBLMSE/data/swissprot/TBLMSE-cox-ub-swissprot-modK-mss-normalized.txt'
# opt$output = '~/Projects/HPCKrogan/Data/TBLMSE/data/swissprot/TBLMSE-cox-ub-swissprot-modK-mss-normalized.pdf'

# opt$command = 'saint-format'
# opt$files = '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-human-exvivo/AP-MS/M2/data/FLU-HUMAN-ALL-APMS-M2-data.txt'
# opt$keys = '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-human-exvivo/AP-MS/M2/data/FLU-HUMAN-ALL-APMS-M2-keys.txt'
# opt$proteome = '~/Projects/HPCKrogan/Data/Uniprot/homo-sapiens-swissprot.fasta'

# opt$command = 'heatmap'
# opt$input = '~/Projects/HPCKrogan/Data/FluOMICS/projects/Flu-human-exvivo/H5N1/ph//results//20141203//FLU-HUMAN-H5N1-PH-results.txt'
# opt$output = '~/Projects/HPCKrogan/Data/FluOMICS/projects/Flu-human-exvivo/H5N1/ph//results//20141203//FLU-HUMAN-H5N1-PH-results.pdf'
# opt$labels=NULL, names='Protein', cluster_cols=F, display='log2FC'

# opt$command = 'results-wide'
# opt$input = '~/Projects/HPCKrogan/Data/HIV-proteomics/results/20141124-ub-proteins/HIV-UB-SILAC-KROGAN-results.txt'
# opt$output = '~/Projects/HPCKrogan/Data/HIV-proteomics/results/20141124-ub-proteins/HIV-UB-SILAC-KROGAN-results.txt'

# opt$command = 'annotate'
# opt$input = '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-human-exvivo/H5N1/ph//results//20141214//FLU-HUMAN-H5N1-PH-results-mapped.txt'
# opt$output = '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-human-exvivo/H5N1/ph//results//20141214//FLU-HUMAN-H5N1-PH-results-mapped-ann.txt'

# opt$command = 'convert-silac'
# opt$files = '~/Code/RMSQ/tests/abundance/ENTERO-LG-ABU-data.txt'
# opt$output = '~/Code/RMSQ/tests/abundance/ENTERO-LG-ABU-data-long.txt'

# opt$command = 'concat'
# opt$files = '~/Projects/HPCKrogan/Data/Mtb/Files/073113*evidence.txt'
# opt$output = '~/Projects/HPCKrogan/Data/HIV-proteomics/Meena/abundance/HIV_vs_MOCK_PROTEIN_evidence_split.txt'

# opt$command = 'convert-sites'
# opt$files =  '~/Projects/HPCKrogan/Data/HIV-proteomics/data//UB-silac//HIV-UB-SILAC-KROGAN-data.txt' 
# opt$output = '~/Projects/HPCKrogan/Data/HIV-proteomics/data//UB-silac//HIV-UB-SILAC-KROGAN-data-modK.txt'
# opt$mod_type = 'ub'
# opt$proteome = '~/Projects/HPCKrogan/Data/Uniprot/homo-sapiens-swissprot.fasta'

# opt$command = 'convert-sites'
# opt$files =  '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-mouse-invivo/Pilot/Mouse-H5N1-Ph//data//FLU-MOUSE-H5N1-PH-PILOT-data.txt' 
# opt$output =  '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-mouse-invivo/Pilot/Mouse-H5N1-Ph//data//FLU-MOUSE-H5N1-PH-PILOT-data-modSTY.txt'
# opt$mod_type = 'ph'
# opt$proteome = '~/Projects/HPCKrogan/Data/Uniprot/mus-musculus-uniprot.fasta'

# opt$command = 'mapback-sites'
# opt$files =  '~/Projects/HPCKrogan/Data/TBLMSE/results/20141124-sites-eqmed-noint/TBLMSE-cox-ub-swissprot-modK-results.txt'
# opt$mapping =  '~/Projects/HPCKrogan/Data/TBLMSE/data/swissprot/TBLMSE-cox-ub-swissprot-modK-mapping.txt'
# opt$output = '~/Projects/HPCKrogan/Data/TBLMSE/results/20141124-sites-eqmed-noint/TBLMSE-cox-ub-swissprot-modK-results.txt'

# opt$command = 'heatmap'
# opt$files =  '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-mouse-invivo/H1N1/ub/results/20150114/FLU-MOUSE-H1N1-UB-results-ann.txt'
# opt$output = '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-mouse-invivo/H1N1/ub/results/20150114/FLU-MOUSE-H1N1-UB-results-ann.pdf'

# opt$command = 'simplify'
# opt$files =  '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-mouse-invivo/H1N1/ub/results/20150114-sites/FLU-MOUSE-H1N1-UB-results.txt'
# opt$output = '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-mouse-invivo/H1N1/ub/results/20150114-sites/FLU-MOUSE-H1N1-UB-results-simplified.txt'

# opt$command = 'heatmap'
# opt$files =  '~/Code/RMSQ/tests/heatmap/031615-lm-1-6-ph-mss-results.txt'
# opt$output = '~/Code/RMSQ/tests/heatmap/031615-lm-1-6-ph-mss-results.pdf'
# opt$labels = '*'
# opt$lfc_lower = -2
# opt$lfc_upper = 2
# opt$q_value = 0.05

# opt$command = 'spectral-counts'
# opt$files =  '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-mouse-invivo/H1N1/ub/data//FLU-MOUSE-H1N1-UB-data.txt'
# opt$keys = '~/Projects/HPCKrogan/Data/FluOMICS/projects/Proteomics/Flu-mouse-invivo/H1N1/ub/data/FLU-MOUSE-H1N1-UB-keys.txt'
# opt$output = '~/Desktop/test.txt'

main(opt)
