#! /usr/bin/Rscript --vanilla

###############################
## FILE AND LIB LOADING #######

cat(">> LOADING EXTERNAL FILES AND LIBRARIES\n")

suppressMessages(library(getopt))
suppressMessages(library(yaml))
suppressMessages(library(biomaRt))

# set source directory
args <- commandArgs(trailingOnly = F) 
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))

## load all external files
if(length(scriptPath)==0){
  source("R/MSstats_functions.R")
}else{
  source(paste(scriptPath,"/MSstats_functions.R",sep=""))  
}

#########################
## MAIN FUNCTIONS #######

filterData = function(data, config){
  cat(">> FILTERING\n")
  if(config$filters$protein_groups == 'remove'){
    cat("\tPROTEIN GROUPS\tREMOVE\n")
    data_f = removeMaxQProteinGroups(data)  
  }else if(config$filters$protein_groups == 'explode'){
    cat("\tPROTEIN GROUPS\tEXPLODE\n")
    data_f = explodeMaxQProteinGroups(data)  
  }else{
    cat("\tPROTEIN GROUPS\tIGNORE\n")
    data_f = data
  }
  if(config$filters$contaminants){
    cat("\tCONTAMINANTS\tREMOVE\n")
    data_f = filterMaxqData(data_f)  
  }
  if(!is.null(config$filters$modification)){
    cat(sprintf("\tMODIFICATIONS\t%s\n",config$filters$modification))
    if(config$filters$modification == 'UB'){
      data_f = data_f[Modifications %like% 'GlyGly']
    }else if(config$filters$modification == 'PH'){
      data_f = data_f[Modifications %like% 'Phospho']
    }
  }
  return(data_f)
}

aggregateData = function(data_w, keys, config, castFun){
  cat(">>AGGREGATING TECHNICAL REPEATS\n")
  ##  we need to convert to long format for more easy aggregation
  data_l = meltMaxQToLong(data_w, na.rm = T)
  keysagg = flattenKeysTechRepeats(keys)
  keysmapping = merge(keys, keysagg[,c('BioReplicate','RawFile')], by='BioReplicate')
  setnames(keysmapping,c('RawFile.x','RawFile.y'),c('RawFile','RawFileCombined'))
  data_l_combined = merge(data_l, keysmapping, by=c('RawFile','IsotopeLabelType'))
  aggregate_fun = tryCatch(eval(parse(text=config$aggregation$aggregate_fun)), 
                           error = function(e) cat("\tWARNING: argument for aggregate_fun not a valid R function - defaulting to:\t 'max'"), 
                           finally=max) 
  
  data_l_combined_agg = data.table(aggregate(Intensity ~ RawFileCombined + Proteins + Sequence + Charge + IsotopeLabelType, FUN=sum, data=data_l_combined))
  setnames(data_l_combined_agg,'RawFileCombined','RawFile')
  data_w_agg = castMaxQToWide(data_l_combined_agg)
  return(list(data_w_agg = data_w_agg, keys_agg = keysagg))
}

## returns data tabel in wide format
normalizeData = function(data_w, config){
  cat(">> NORMALIZING\n")

  if(grepl('scale|quantile|cyclicloess',config$normalization$method)){
    cat(sprintf("\tNORMALIZATION\t%s\n",config$normalization$method))
    data_fn = normalizeSingle(data_w=data_w, NORMALIZATION_METHOD=config$normalization$method)  
  }else if(grepl('reference',config$normalization$method) && !is.null(config$normalization$reference) && config$normalization$reference %in% unique(data_w$Proteins)){
    
    ## CURRENTLY UING MSSTATS METHODS FOR REFERENCE NORMALIZATION
    data_fn = data_w
    
#     cat(sprintf("\tNORMALIZATION\tTO REFERENCE\t%s\n",config$normalization$reference))
#     ref_files = keys[keys$NormalizationGroup == 'REFERENCE', 'RawFile']
#     data_l_ref = data_l[data_l$Raw.file %in% ref_files & data_l$Intensity > 0 & is.finite(data_l$Intensity), ]
#     data_l_nonref = data_l[!(data_l$Raw.file %in% ref_files) & data_l$Intensity > 0 & is.finite(data_l$Intensity), ]
#     data_l_ref_n = normalizeToReference(data_l_ref=data_l_ref, ref_protein = config$normalization$reference, output_file = config$files$output)
#     data_fn= rbind(data_l_ref_n, data_l_nonref)
#     data_fn = castMaxQToWide(data_fn)
    
  }else{
    data_fn = data_w
  }

  return(data_fn)
}

runMSstats = function(dmss, contrasts, config){
  
  if(grepl('before', config$msstats$profilePlots)){
    mssquant = dataProcess(dmss, normalization=F, fillIncompleteRows = T)  
    dataProcessPlots(data=mssquant, type="ProfilePlot", featureName="Peptide", address=gsub('.txt','-before',config$files$output))
    dataProcessPlots(data=mssquant, type="QCPlot", address=gsub('.txt','-before',config$files$output))
  }
  
  if(!is.null(config$msstats$normalization_reference) & config$msstats$normalization_method == 'globalStandards'){
    normalization_refs = unlist(lapply(strsplit(config$msstats$normalization_reference, split = ','), FUN=trim))
    mssquant = dataProcess(dmss, normalization=config$msstats$normalization_method, nameStandards=normalization_refs , fillIncompleteRows=T)
  }else{
    cat(sprintf('>> NORMALIZATION\t%s\n',config$msstats$normalization_method))
    mssquant = dataProcess(dmss, normalization=config$msstats$normalization_method , fillIncompleteRows = F, betweenRunInterferenceScore = F, FeatureSelection = F)
  } 
  
  if(grepl('after', config$msstats$profilePlots)){
    dataProcessPlots(data=mssquant, type="ProfilePlot", featureName="Peptide", address=gsub('.txt','-after',config$files$output))
    dataProcessPlots(data=mssquant, type="QCPlot", address=gsub('.txt','-after',config$files$output))
  }
  
  if(!all(levels(mssquant$GROUP_ORIGINAL) == colnames(contrasts))){
    cat(sprintf('\tERROR IN CONTRAST COMPARISON: GROUP LEVELS DIFFERENT FROM CONTRASTS FILE\n\tGROUP LEVELS\t%s\n\tCONTRASTS FILE\t%s\n', paste(levels(mssquant$GROUP_ORIGINAL),collapse=','),paste(colnames(contrasts),collapse=',')))
    quit()
  } 
  
  cat(sprintf('\tFITTING CONTRASTS:\t%s\n',paste(rownames(contrasts),collapse=',')))
  write.table(mssquant, file=gsub('.txt','-mss-normalized.txt',config$files$output), eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
  results = groupComparison(data = mssquant, contrast.matrix = contrasts, labeled = as.logical(config$msstats$labeled), scopeOfBioReplication = config$msstats$scopeOfBioReplication, scopeOfTechReplication = config$msstats$scopeOfTechReplication, interference = as.logical(config$msstats$interference), equalFeatureVar = as.logical(config$msstats$equalFeatureVar), missing.action = config$msstats$missing_action)$ComparisonResult
  write.table(results, file=config$files$output, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)  
  cat(sprintf(">> WRITTEN\t%s\n",config$files$output))
  return(results)
}


convertDataLongToMss = function(data_w, keys, config){
  cat(">> CONVERTING DATA TO MSSTATS FORMAT\n")
  data_l = meltMaxQToLong(data_w, na.rm = F)
  data_lk = mergeMaxQDataWithKeys(data_l, keys, by=c('RawFile','IsotopeLabelType'))
  dmss = dataToMSSFormat(data_lk)
  ## sanity check for zero's
  if(nrow(dmss[!is.na(dmss$Intensity) & dmss$Intensity == 0,]) > 0){
    dmss[!is.na(dmss$Intensity) & dmss$Intensity == 0,]$Intensity = NA
  } 
  write.table(dmss, file=gsub('.txt','-mss.txt',config$files$data), eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
  return(dmss)  
}

writeExtras = function(results, config){
  cat(">> ANNOTATING\n")
  if(config$output_extras$biomart){
    #mart = useMart(biomart = 'unimart')
    mart = useMart(biomart = 'unimart', dataset = 'uniprot', verbose = T, host="www.ebi.ac.uk", path="/uniprot/biomart/martservice")
    mart_anns = AnnotationDbi::select(mart, keytype='accession', columns=c('accession','name','protein_name','gene_name','ensembl_id'), keys=as.character(unique(results$Protein)))
    if(nrow(mart_anns) > 0){
      mart_anns = aggregate(. ~ accession, data=mart_anns, FUN=function(x)paste(unique(x),collapse=','))
      results_ann = merge(results, mart_anns, by.x='Protein', by.y='accession', all.x=T)
      cat(sprintf('\tCHANGED OUTPUT FILE TO\t%s\n',config$files$output))
      write.table(results_ann, file=gsub('.txt','-ann.txt',config$files$output), eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
      cat(sprintf(">> WRITTEN\t%s\n",config$files$output)) 
    }else{
      results_ann = results
    }
  }else{
    results_ann = results
    config$files$output = config$output_extras$msstats_output
  }
  
  lfc_lower = as.numeric(unlist(strsplit(config$output_extras$LFC,split=" "))[1])
  lfc_upper = as.numeric(unlist(strsplit(config$output_extras$LFC,split=" "))[2])
  ## select subset of labels for heatmap and volcan plots
  selected_labels = config$output_extras$comparisons
  if(is.null(selected_labels) || selected_labels=='all') selected_labels='*'
  
  ## select data points  by LFC & FDR criterium in single condition and adding corresponding data points from the other conditions
  sign_hits = significantHits(results_ann,labels=selected_labels,LFC=c(lfc_lower,lfc_upper),FDR=config$output_extras$FDR)
  sign_labels = unique(sign_hits$Label)
  cat(sprintf("\tSELECTED HITS FOR PLOTS WITH LFC BETWEEN %s AND %s AT %s FDR:\t%s\n",lfc_lower, lfc_upper, config$output_extras$FDR, nrow(sign_hits)/length(sign_labels))) 
	  
	
	cat(paste("output_file: ", config$files$output, "\n"))

  ## REPRESENTING RESULTS AS HEATMAP
  if(config$output_extras$heatmap){
    ## plot heat map for all contrasts
    cat(">>   PLOTTING HEATMAP FOR ALL CONTRASTS\n")
    heat_labels = prettyPrintHeatmapLabels(uniprot_acs=sign_hits$Protein,uniprot_ids=sign_hits$name, gene_names=sign_hits$gene_name)
    heat_data_w = plotHeat(mss_F = sign_hits, out_file =  gsub('.txt','-sign.pdf',config$files$output), names=heat_labels, cluster_cols=config$output_extras$heatmap_cluster_cols, display = config$output_extras$heatmap_display)  
  }
  
  if(config$output_extras$volcano){
  	cat(">>   PLOTTING VOLCANO PLOT\n")
    file_name = gsub('.txt','-volcano.pdf',config$files$output)
    volcanoPlot(results_ann[grep(selected_labels,results_ann$Label),], lfc_upper, lfc_lower, FDR=config$output_extras$FDR, file_name=file_name)  
  }
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)


main <- function(opt){
  cat(">> MSSTATS PIPELINE\n")
  config = tryCatch(yaml.load_file(opt$config), error = function(e) {cat(opt$config);break} )
  
  if(config$data$enabled){
    cat(">> LOADING DATA\n")
    data = fread(config$files$data, stringsAsFactors=F, integer64 = 'double')
    setnames(data, colnames(data),gsub('\\s','.',colnames(data)))
    keys = fread(config$files$keys, stringsAsFactors=F, integer64 = 'double')
    
    ## the following lines were added to integrate the Label with the Filename when using multiple labels (e.g. H/L)
    ## currently we leave this in because the MSstats discinction on labeltype doesn't work 
    ## see ISSUES https://github.com/everschueren/RMSQ/issues/1
    
    tryCatch(setnames(data, 'Raw.file', 'RawFile'), error=function(e) cat('Raw.file not found\n'))
    tryCatch(setnames(keys, 'Raw.file', 'RawFile'), error=function(e) cat('Raw.file not found\n'))
    
    cat('\tVERIFYING DATA AND KEYS\n')
    if(!'IsotopeLabelType' %in% colnames(data)) data[,IsotopeLabelType:='L']

    data = mergeMaxQDataWithKeys(data, keys, by = c('RawFile','IsotopeLabelType'))
    data$RawFile = paste(data$RawFile, data$IsotopeLabelType, sep='')
    keys$RawFile = paste(keys$RawFile, keys$IsotopeLabelType, sep='')
    keys$Run = paste(keys$IsotopeLabelType,keys$Run , sep='')
    data$IsotopeLabelType = 'L'
    keys$IsotopeLabelType = 'L'
    data[Intensity<1,]$Intensity=NA ## fix for weird converted values from fread
    
    ## end hacks for SILAC
    
    ## FILTERING
    if(config$filters$enabled) data_f = filterData(data, config) else data_f=data
    
    ## FORMATTING IN WIDE FORMAT FOR NORMALIZATION PURPOSES
    if(config$files$sequence_type == 'modified') castFun = castMaxQToWidePTM else castFun = castMaxQToWide
    data_w = castFun(data_f)
    
    ## test code for heatmap
    if(!is.null(config$files$sample_plots) && config$files$sample_plots){
      keys_in_data = keys[keys$RawFile %in% unique(data$RawFile),]
      sampleCorrelationHeatmap(data_w = data_w, keys = keys_in_data, config = config) 
      samplePeptideBarplot(data_f, config)
    }
    
    ## AGGREGATION
    if(config$aggregation$enabled){
      res = aggregateData(data_w, keys, config)
      data_w=res$data_w_agg
      keys=res$keys_agg
    } 
    
    ## NORMALIZATION -- we normalize with MSstats now
    #if(config$normalization$enabled) data_fn = normalizeData(data_w, config) else data_fn=data_w
    
  }
  
  ## MSSTATS
  if(config$msstats$enabled){
    
    if(is.null(config$msstats$msstats_input)){
      dmss = data.table(convertDataLongToMss(data_w, keys, config))
      
      ## make sure there are no doubles !!
      ## doubles could arise when protein groups are being kept and the same peptide is assigned to a unique Protein. Not sure how this is possible but it seems to be like this in maxquant output. A possible explanation is that these peptides have different retention times (needs to be looked into)
      dmss = data.frame(dmss[,j=list(ProteinName=paste(ProteinName,collapse=';'),Intensity=median(Intensity, na.rm=T)),by=c('PeptideSequence','ProductCharge','PrecursorCharge','FragmentIon','IsotopeLabelType','Run','BioReplicate','Condition')])
      
    }else{
      cat(sprintf("\tREADING PREPROCESSED\t%s\n", config$msstats$msstats_input)) 
      dmss = fread(config$msstats$msstats_input, stringsAsFactors=F, integer64 = 'double')
    }
    
    cat(sprintf('>>LOADING MSSTATS %s VERSION\n', config$msstats$version))
    if(!is.null(config$msstats$version) & config$msstats$version == 'MSstats.daily'){
      library('MSstats.daily', character.only = T) 
    }else if(!is.null(config$msstats$version) & config$msstats$version == 'MSstats'){
      library('MSstats', character.only = T) 
    }else{
      if(file.exists(config$msstats$version)){
        for (fname in list.files(config$msstats$version, pattern = '\\.[Rr]$')){
          cat(sprintf("\tLOADING %s\n", fname))
          source(file.path(config$msstats$version, fname))
        }
      }else{
        cat(sprintf('COULD NOT LOAD MSSTATS FUNCTIONS AT %s', config$msstats$version)) 
      }
    }
    
    ## solely for debugging purposes
    if(exists("DEBUG")){
      if(DEBUG){
        if(!DEBUGALL){
          set.seed(7)
          protein_sample  = sample(unique(dmss$ProteinName), 100)
          dmss_sample = dmss[dmss$ProteinName %in% protein_sample,]
          if(DEBUGONE){
            dmss_sample_1 = dmss[dmss$Protein==DEBUGSUBJECT,]
            dmss_sample = rbind(dmss_sample, dmss_sample_1)
          }
          dmss = dmss_sample
        }
      }  
    }
    
    contrasts = as.matrix(read.delim(config$files$contrasts, stringsAsFactors=F))
    results = runMSstats(dmss, contrasts, config)
  }
  
  ## ANNOTATING RESULT FILE
  if(config$output_extras$enabled){
    if(! is.null(config$output_extras$msstats_output)) results = read.delim(config$output_extras$msstats_output, stringsAsFactors=F)
    writeExtras(results, config)
  }
}

#########################
## CONFIG LOADING #######

spec = matrix(c(
  'verbose', 'v', 2, "integer", "",
  'help'   , 'h', 0, "logical", "available arguments (this screen)",
  'config'  , 'c', 1, "character", "configuration file in YAML format"),
  byrow=TRUE, ncol=5)

opt = getopt(spec = spec, opt = commandArgs(TRUE), command = get_Rscript_filename(), usage = FALSE, debug = FALSE)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

## TEST WORKS WITH LATEST CODE
# opt = c(opt, config='tests/LabelFree-ub/LabelFree-ub-test.yaml')

if(!exists("DEBUG")){
  cat(">> RUN MODE\n")
  main(opt)
}else{
  cat(">> DEBUG MODE\n")
} 
# rm(list = ls(all = TRUE)) 
