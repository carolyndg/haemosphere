# Return a matrix resembling limma's topTable function.
#
# Example: Suppose x is the data frame of expression values, with columns ('BCell(1)','TCell(1)','TCell(2)','BCell(2)','HSC(1)','HSC(2)').
# haemosphere.topTable(x, c('B','T','T','B','HSC','HSC'), c('B','T','T','B','HSC','HSC'), 'B', 'T') 
# will return topTable of DE genes between 'B' and 'T'. 
# Use microarray=False if x is rna-seq data (filterMinCPM, filterMinExpressedSamples and normalizationMethod are ignored for microarray).
#
# Now we can also do haemosphere.topTable(x, c('B','T','T','B','HSC','HSC'), c('Committed','Committed','Committed','Committed','Prog','Prog'), 'Committed', 'Prog') 
# to find DE genes between progenitor cell types and committed cell types. In this case, the design matrix is
# constructed the same way as before (hence we need replicateGroups), but contrast matrix is constructed with correct
# mixtures of each replicates.
#
# minRows is used to return the minimum number of rows for the returned table, regardless of adjPCutoff. 
# This is useful when there are 0 or hardly any DE genes - at least the user can see that there hasn't been an error.
#
# Parameters
# ----------
# x (data.frame): expression data frame with features as rows, sample ids as columns
# replicateGroups (string vector): sample group values deemed "replicates" designated to each sample id, in the same order and length as columns of x
# groups (string vector): sample group values of interest for analysis, in the same order and length as columns of x
# group1 (string): a member of groups
# group2 (string): another member of groups, different to group1
# microarray (boolean): True if x is microarray expression data, False otherwise
# filterMinCPM (float), filterMinExpressedSamples (integer): 
#		keep only rows where filterMinExpressedSamples or more samples out of group1 and group2 must be greater than filterMinCPM 
# normalizationMethod (string): same as used in the 'method' argument of edgeR's calcNormFactors fuction.
# adjPCutoff (float): consider rows with adj.P.Val (from limma's topTable) less than this value as significant.
# minRows (integer): minimum number of rows to return regardless of significant rows.
#
# Returns
# -------
# data.frame with columns ("features","logFC","adjPVal","AveExpr").
# 		matrix(0,0,0) if something went wrong.
#


haemosphere.topTable = function(x, replicateGroups, groups, group1, group2, microarray=T, 
                                filterMinCPM=0.5, filterMinExpressedSamples=2, normalizationMethod='TMM', adjPCutoff=0.05, minRows=50)

  {	
  if (microarray) library(limma)
  else library(edgeR)
  options(digits=4)
  
  # find indices of group1 and group2 from groups vector
  group1Indices = which(groups==group1)
  group2Indices = which(groups==group2)
  
  # Design matrix is created using replicate groups; it can't have certain characters, so use make.names first.
  # (So do not try to match the elements of these vectors before and after make.names.)
  # make.names doesn't recognise special characters like '+' or '-' properly, so both 'CD8+' and 'CD8-' will become 'CD8.'
  # To avoid this, create new names to both groups and replicateGroups based on their position in unique
  # Note also that if an element of a char vector starts with a 'number', such as '6hrFXII', make.names will prepend 'X' (same for numeric vectors)
  replicateGroups = make.names(match(replicateGroups, unique(replicateGroups)))  # we'll end up with c('X1','X2','X1',...)
  groups = make.names(match(groups, unique(groups)))
  design = model.matrix(~0+as.factor(replicateGroups))
  colnames(design) = levels(as.factor(replicateGroups))
  
  # Since we need to supply the arguments to makeContrasts dynamically, we can't just use the normal form of makeContrasts,
  # because it expects variable names rather than strings. Using parse() after constructing a string works in a normal R session,
  # but it doesn't work through rpy2. Instead, found the trick using do.call from 
  # https://stat.ethz.ch/pipermail/bioconductor/2009-June/027913.html, which works wonders.
  if (identical(replicateGroups, groups)) {  # makeContrasts the normal way
    cont.matrix = do.call(makeContrasts, list(paste(replicateGroups[group1Indices][1], replicateGroups[group2Indices][1], sep="-"), levels=design))
  }
  else {
    # construct strings to be used for makeContrasts; they should look like "(B+T)/4"
    group1String = paste("(", paste(unique(replicateGroups[group1Indices]),collapse="+"), ")/", length(replicateGroups[group1Indices]), sep="")
    group2String = paste("(", paste(unique(replicateGroups[group2Indices]),collapse="+"), ")/", length(replicateGroups[group2Indices]), sep="")
    cont.matrix = do.call(makeContrasts, list(paste(group1String, group2String, sep="-"), levels=design))
  }
  
  # perfrom lmFit
  if (microarray)
    fit = lmFit(x, design)
  else {
    # create a DGEList object and normalize if specified
    x = DGEList(x)
    if (!is.null(normalizationMethod)) x = calcNormFactors(x, method=normalizationMethod)
    
    # keep only rows where filterMinExpressedSamples or more samples out of group1 and group2 must be greater than filterMinCPM
    x.cpm = cpm(x)
    x = x[rowSums(x.cpm[,c(group1Indices, group2Indices)]>filterMinCPM)>=filterMinExpressedSamples,]
    
    # change normalization to normalize.method as voom is updated
    fit = lmFit(voom(x, design=design, normalize.method='none'), design)
  }
  
  eb = eBayes(contrasts.fit(fit, cont.matrix))
  
  # fetch topTable
  m = topTable(eb, adjust='fdr', sort.by='P', number=Inf)   # this should fetch every gene/probe sorted by adj p
  
  sigRows = which(m[,"adj.P.Val"]<adjPCutoff)
  
  if (nrow(m)<minRows)
    minRows = nrow(m)
  
  if (length(sigRows)<minRows) # too few rows, so return at least minRows
    m = m[1:minRows,]
  else	# return all significant rows
    m = m[sigRows,]
  
  # 2015-09-25: After a recent update of rpy2, it seems row index is not preserved after converting from R to pandas.
  # so add the rownames into m
  m['features'] = rownames(m)
  
  # select column subset and remove dots in adj.P.Val colname
  m = m[c("features","logFC","adj.P.Val","AveExpr")]
  colnames(m)[3] = "adjPValue"
  
  if (nrow(m)==0)  # something went wrong
    return(matrix(0,0,0))
  else
    return(data.frame(m))
}

# Wrapper for the main function. This is done to enable object saving. Note that wrapper has one extra argument of 'saveToFile'. 
# If this points to a filename or filepath, all the input objects as well as the main function will be saved to this file. 
# This enables replication of the function in a different environment by using load().
haemosphere.topTableWrapper = function(x, replicateGroups, groups, group1, group2, microarray=T, 
                                       filterMinCPM=0.5, filterMinExpressedSamples=2, normalizationMethod='TMM', adjPCutoff=0.05, minRows=50, saveToFile='')
{
  if (saveToFile!='') {
    save(haemosphere.topTable, x, replicateGroups, groups, group1, group2, microarray, filterMinCPM, 
         filterMinExpressedSamples, normalizationMethod, adjPCutoff, minRows, file=saveToFile)
  }
  return(haemosphere.topTable(x, replicateGroups, groups, group1, group2, microarray=microarray, filterMinCPM=filterMinCPM, 
                              filterMinExpressedSamples=filterMinExpressedSamples, normalizationMethod=normalizationMethod, 
                              adjPCutoff=adjPCutoff, minRows=minRows))
}
