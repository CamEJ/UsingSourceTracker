# Running SourceTracker version 1.0.1
# https://github.com/danknights/sourcetracker
# Download tar.gz and then extract this and set folder as wkdir
# Open example.r and run script the first time to check everything works 
# correctly using their example data. 

# ================== getting data ready ====================================# 

# you will need 1) a metadatafile and 2) an OTU table

# 1. metadata file formate
# this must have:
# a) a column called 'SourceSink' saying if your sample is a source or a sink.
# b) a column called 'Env' saying which environment your sample is from
# c) a column called Description which you can use to describe your samples
# d) any samples you do now want to be analysed as a source or sink should 
# have NA in these 3 columns. 

# head(metadata)
#           SlurryAdd   Phase Treatment SourceSink    Env Description
# Slurry1        <NA>    <NA>      <NA>     source slurry      slurry
# Slurry2        <NA>    <NA>      <NA>     source slurry      slurry
# T0_Tmt1_a     Minus NoFlood   Control     source   soil        soil
# T0_Tmt1_b     Minus NoFlood   Control     source   soil        soil
# T0_Tmt1_c     Minus NoFlood   Control     source   soil        soil
# T0_Tmt2_a      Plus NoFlood    Slurry       sink   tmt2    0_slurry

# 2. OTU table.

# This needs to have same sample names as in your metadata file
# They suggest using trimmed OTU to remove low abundance OTUs
# I ran this with OTUs trimmed to > 0.005% as per Bokulich et al 2013. 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3531572/

# The script below expects that initially your samples will be columns and 
# OTUs will be rows (& they will therefore be transformed). 
# your OTU table should not include any taxo information. 


# ====================== loading data & running scrips ================================

# load sample metadata
metadata <- read.table('metadata2slurry.txt',sep='\t',h=T,row.names=1)

# load OTU table
# This 'read.table' command is designed for a 
# QIIME-formatted OTU table.
# namely, the first line begins with a '#' sign
# and actually _is_ a comment; the second line
# begins with a '#' sign but is actually the header
otus <- read.table('OTUtable_2slurry.txt',sep='\t', header=T,row.names=1,check=F,skip=1,comment='')
head(otus)
dim(otus)
otus2 = otus[,c(1:8,15:90)] # removing extra T0 samples & taxo info as non in eg data
head(otus2) # check
otus <- t(as.matrix(otus2))

# extract only those samples in common between the two tables
common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids,]
metadata <- metadata[common.sample.ids,]
# double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}

# extract the source environments and source/sink indices
train.ix <- which(metadata$SourceSink=='source')
test.ix <- which(metadata$SourceSink=='sink')
envs <- metadata$Env
if(is.element('Description',colnames(metadata))) desc <- metadata$Description


# load SourceTracker package
source('src/SourceTracker.r')

# tune the alpha values using cross-validation (this is slow!)
# tune.results <- tune.st(otus[train.ix,], envs[train.ix])
# alpha1 <- tune.results$best.alpha1
# alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st <- sourcetracker(otus[train.ix,], envs[train.ix])

# Error in sum(x[i, ]) : invalid 'type' (character) of argument
# I initially got this error when I had only one source sample.
# added a second and it worked a charm

# Estimate source proportions in test data
results <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2)

# Estimate leave-one-out source proportions in training data 
results.train <- predict(st, alpha1=alpha1, alpha2=alpha2)

# plot results
labels <- sprintf('%s %s', envs,desc)
plot(results, labels[test.ix], type='pie')

# other plotting functions
# plot(results, labels[test.ix], type='bar')
plot(results, labels[test.ix], type='dist')
# plot(results.train, labels[train.ix], type='pie')
# plot(results.train, labels[train.ix], type='bar')
# plot(results.train, labels[train.ix], type='dist')

# plot results with legend
plot(results, labels[test.ix], type='pie', include.legend=TRUE, env.colors=c('#47697E','#5B7444','#CC6666','#79BEDB','#885588'))

# ===================reading data out to plot more clearly======================

xx <- lapply(results, unlist)
max <- max(sapply(xx, length))
job = do.call(rbind, lapply(xx, function(z)c(z, rep(NA, max-length(z)))))

jobbie = as.data.frame(job)