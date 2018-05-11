
#!/bin/bash

############################################################################
##  This script will create OTU tables via three different methods:
##  dada2, vsearch de novo clustering, and deblur
##  
##  This script will
##  1. Import your data
##  2. Direct it through "clustering" methods
##  3. Generate basic richness measurements
##  4. Generate representative sequences
##
##  This script MUST be customized to your own data
##  Make sure to read through it and adjust accordingly
##  Major points of required customization should be pointed out
##  Before Starting:
##  1. change your original input data folder#
##  2. each data element is prepended by "samA"; change "samA" to something
##  unique to your study (bulk find and replace) and you should be clear to go
##  just remember to activate (remove "#") or deactivate (prepend with "#")
##  sections that you want to run or not run accordingly
##  Also, if you are running this on another computer, change the thread calls
##  Search and replace "24" with your desired thread call
##  note also that this pipeline is designed for the 341f 785r primer set
##  trying with a different primer set will not work
############################################################################


## It is encouraged that the user copies out desired subsections of this
## master script and run them as needed



############################################################################
########################### main data import ###############################
############################################################################

# Only use your main data import once
# deactivated by default

# qiime tools import \
#--type 'SampleData[PairedEndSequencesWithQuality]' \
#--input-path testreads \
#--source-format CasavaOneEightSingleLanePerSampleDirFmt \
#--output-path samAreads.qza

#qiime demux summarize \
#--i-data samAreads.qza \
#--o-visualization samA-demux.qzv \
#--verbose


############################################################################
## the next steps prepare the data for two different processing pipelines ##
############################################################################

#1#
#remove the primers#
#you MUST customize with YOUR primers

qiime cutadapt trim-paired \
--i-demultiplexed-sequences samAreads.qza \
--p-cores 24 \
--p-front-f CCTACGGGNGGCWGCAG \
--p-front-r GACTACHVGGGTATCTAATCC \
--o-trimmed-sequences samA-data-cutadapt \
--verbose

#2#
#join pairs with some strictness and QC filtering#
#this is used only for the VSearch and deblur pipelines
#you MUST customize for YOUR overlap and YOUR amplicon joined size

qiime vsearch join-pairs \
--i-demultiplexed-seqs samA-data-cutadapt.qza \
--p-minovlen 50 \
--p-minmergelen 380 \
--p-maxmergelen 480 \
--o-joined-sequences samA-data-joined \
--verbose

############################################################################
######################### VSearch De Novo Clustering #######################
############################################################################

#1#
#this is a strict quality filtering#
#any sequence with a base with lower than PHRED25 is removed#

qiime quality-filter q-score-joined \
--i-demux samA-data-joined.qza \
--p-min-quality 25 \
--o-filtered-sequences samA-seqs-strict-qtrim \
--o-filter-stats samA-seqs-strict-qtrim-stats \
--verbose

#2#
#this is required, creates a feature table#
#very time consuming step#

qiime vsearch dereplicate-sequences \
  --i-sequences samA-seqs-strict-qtrim.qza \
  --o-dereplicated-table samA-denovo-derep-table \
  --o-dereplicated-sequences samA-denovo-derep-seqs

#4#
# chimera removal #

qiime vsearch uchime-denovo \
  --i-sequences samA-denovo-derep-seqs.qza \
  --i-table samA-denovo-derep-table.qza \
  --o-chimeras samA-denovo-chimera-seqs \
  --o-nonchimeras samA-denovo-nonchimera-seqs \
  --o-stats samA-denovo-chimera-stats

# 6 #
# after separating chimeras, i have to filter the input data, removing chimeras#

qiime feature-table filter-features \
  --i-table samA-denovo-derep-table.qza \
  --m-metadata-file samA-denovo-chimera-seqs.qza \
  --p-exclude-ids \
  --o-filtered-table samA-denovo-filtered-table.qza

qiime feature-table filter-seqs \
  --i-data samA-denovo-derep-seqs.qza \
  --m-metadata-file samA-denovo-chimera-seqs.qza \
  --p-exclude-ids \
  --o-filtered-data samA-denovo-filtered-seqs.qza

qiime feature-table summarize \
  --i-table samA-denovo-filtered-table.qza \
  --o-visualization samA-denovo-filtered-table.qzv

#3#
#this is the clustering but does not actually create a taxonomy table#
#essentially, clusters are created primarily using the ref database as seed

##### I should also run a 98.5% clustering here #####

qiime vsearch cluster-features-de-novo \
--i-sequences samA-denovo-filtered-seqs.qza \
--i-table samA-denovo-filtered-table.qza \
--i-reference-sequences silva_128_97.qza \
--p-perc-identity 0.97 \
--p-threads 24 \
--o-clustered-table samA-denovo-OTU-table \
--o-clustered-sequences samA-denovo-OTU-seqs \
--o-new-reference-sequences samA-denovo-OTU-ref-seqs


# 7 #
# work towards alpha diversity analysis, first build a tree #

qiime feature-table summarize \
  --i-table samA-denovo-filtered-table.qza \
  --o-visualization samA-denovo-filtered-table.qzv \
  --m-sample-metadata-file Sample_metadata2.txt

qiime feature-table tabulate-seqs \
  --i-data samA-denovo-filtered-ref-seqs.qza \
  --o-visualization samA-denovo-filtered-ref-seqs.qzv

qiime alignment mafft \
  --i-sequences samA-denovo-OTU-filtered-ref-seqs.qza \
  --o-alignment samA-denovo-aligned.qza \
  --p-n-threads 24

qiime alignment mask \
  --i-alignment samA-denovo-aligned.qza \
  --o-masked-alignment samA-denovo-aligned-masked.qza

qiime phylogeny fasttree \
  --i-alignment samA-denovo-aligned-masked.qza \
  --o-tree samA-denovo-unrooted-tree.qza

qiime phylogeny midpoint-root \
--i-tree samA-denovo-unrooted-tree.qza \
--o-rooted-tree samA-denovo-rooted-tree.qza

# 8 #
# generate core diversity metrics #
# pay close attention to the sampling depth command #

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny samA-denovo-rooted-tree.qza \
  --i-table samA-denovo-filtered-table.qza \
  --p-sampling-depth 10 \
  --m-metadata-file Sample_metadata2.txt \
  --output-dir samA-denovo-core-metrics-results

## visualization of diversity metrics ##

qiime diversity alpha-rarefaction \
  --i-table samA-denovo-filtered-table.qza \
  --i-phylogeny samA-denovo-rooted-tree.qza \
  --p-min-depth 1 \
  --p-max-depth 20000 \
  --m-metadata-file Sample_metadata2.txt \
  --o-visualization samA-denovo-alpha-rarefaction.qzv

############################################################################
############################# Deblur #######################################
############################################################################

#### don't use this section for now, deblur I don't especially trust #######

# this takes the joined paired ends as input
# first a light quality filter #

#qiime quality-filter q-score-joined \
#--i-demux samA-data-joined.qza \
#--o-filtered-sequences samA-qtrimmed-seqs-deblur \
#--o-filter-stats samA-qtrim-stats-deblur \
#--verbose

# dereplicating, making feature table to reduce computational load #

#qiime vsearch dereplicate-sequences \
#  --i-sequences samA-qtrimmed-seqs-deblur.qza \
#  --o-dereplicated-table samA-deblur-qtrimmed-derep-table \
#  --o-dereplicated-sequences samA-deblur-qtrimmed-derep-seqs

#the deblur denoising is mainly an error correction "true" read determination
#it is supposedly not designed for joined sequences or even for paired end
#it also must have a set of "positive" sequences, so as to remove junk?
#This is problematic for me... I would have to cluster the BOLD database and thus make some 
#reduced data set
#but qiime2 allows it; use with caution I suppose?
#an 88% clustering of the BOLD database would be equivalent to the 16S implementation

###########################################################################################
# CONTINUE EDITING PIPELINE HERE ####################################

#qiime deblur denoise \
#--i-demultiplexed-seqs samA-qtrimmed-seqs-deblur.qza \
#--p-trim-length 380 \
#--p-sample-stats \
#--p-jobs-to-start 24 \
#--o-representative-sequences samA-deblur-rep-seqs \
#--o-table samA-deblur-table \
#--o-stats samA-deblur-stats

#then we can simply classify using the same method

#qiime feature-classifier classify-sklearn \
#  --i-classifier classifier.qza \
#  --p-n-jobs 24 \
#  --i-reads samA-deblur-rep-seqs.qza \
#  --o-classification samA-deblur-taxonomy

#qiime metadata tabulate \
#  --m-input-file samA-deblur-taxonomy.qza \
#  --o-visualization samA-deblur-taxonomy.qzv

#qiime taxa barplot \
#  --i-table samA-deblur-table.qza \
#  --i-taxonomy samA-deblur-taxonomy.qza \
#  --m-metadata-file Sample_metadata2.txt \
#  --o-visualization samA-deblur-taxa-bar-plots.qzv

############################################################################
# the dada2 pipeline, results in an otu table which can then be classified #
############################################################################

## this uses the cutadapt output ##
## it does not use the joined sequences, it joins itself

qiime dada2 denoise-paired \
--i-demultiplexed-seqs samA-data-cutadapt.qza \
--o-table samA-dada2-table \
--o-representative-sequences samA-dada2_seqs \
--p-n-threads 4 \
--p-trunc-len-f 250 \
--p-trunc-len-r 200

qiime feature-table summarize \
--i-table samA-dada2-table.qza \
--o-visualization samA-dada2-seq-stats \
--verbose

#qiime feature-classifier classify-sklearn \
#  --i-classifier classifier.qza \
#  --i-reads samA-dada2_seqs.qza \
#  --p-n-jobs 24 \
#  --o-classification samA-dada2-taxonomy

qiime metadata tabulate \
  --m-input-file samA-dada2-taxonomy.qza \
  --o-visualization samA-dada2-taxonomy

#qiime taxa barplot \
#  --i-table samA-dada2-table.qza \
#  --i-taxonomy samA-dada2-taxonomy.qza \
#  --m-metadata-file Sample_metadata2.txt \
#  --o-visualization samA-dada2-taxa-bar-plots

# this is removing samples at a low frequency count that I calculated
# independently (5% of median)
# this part MUST to be customized to which pipeline results you want to proceed
# with... this uses dada2 currently

# THIS MUST BE CUSTOMIZED AND IS SUBJECTIVE #

#qiime feature-table filter-samples \
#  --i-table samA-dada2-table.qza \
#  --p-min-frequency 1200 \
#  --o-filtered-table samA-dada2-table-lowcountfilter

#qiime feature-table summarize \
#  --i-table samA-dada2-table-lowcountfilter.qza \
#  --o-visualization samA-dada2-table-lowcountfilter.qzv

# rebuilding the barplots with the lowcountfilter #

#qiime taxa barplot \
#  --i-table samA-dada2-table-lowcountfilter.qza \
#  --i-taxonomy samA-dada2-taxonomy.qza \
#  --m-metadata-file Sample_metadata2.txt \
#  --o-visualization samA-dada2-taxa-bar-plots-lowcountfiltered

########### DIVERSITY METRICS FOR DADA2 AFTER LOW COUNT FILTERING ########

# work towards alpha diversity analysis, first build a tree #

#qiime feature-table summarize \
#  --i-table samA-dada2-table-lowcountfilter.qza \
#  --o-visualization samA-dada2-table-lowcountfilter.qzv \
#  --m-sample-metadata-file Sample_metadata2.txt

qiime feature-table tabulate-seqs \
  --i-data samA-dada2_seqs.qza \
  --o-visualization samA-dada2_seqs.qzv

qiime alignment mafft \
  --i-sequences samA-dada2_seqs.qza \
  --o-alignment samA-dada2_seqs-aligned.qza \
  --p-n-threads 24

qiime alignment mask \
  --i-alignment samA-dada2_seqs-aligned.qza \
  --o-masked-alignment samA-dada2_seqs-aligned-masked.qza

qiime phylogeny fasttree \
  --i-alignment samA-dada2_seqs-aligned-masked.qza \
  --o-tree samA-unrooted-tree-dada2.qza

qiime phylogeny midpoint-root \
--i-tree samA-unrooted-tree-dada2.qza \
--o-rooted-tree samA-rooted-tree-dada2.qza


# generate core diversity metrics #
# pay close attention to the sampling depth command #

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny samA-rooted-tree-dada2.qza \
  --i-table samA-dada2-table-lowcountfilter.qza \
  --p-sampling-depth 10 \
  --m-metadata-file Sample_metadata2.txt \
  --output-dir samA-core-metrics-results-dada2

## visualization of diversity metrics ##

qiime diversity alpha-rarefaction \
  --i-table samA-dada2-table-lowcountfilter.qza \
  --i-phylogeny samA-rooted-tree-dada2.qza \
  --p-min-depth 1 \
  --p-max-depth 300000 \
  --m-metadata-file Sample_metadata2.txt \
  --o-visualization samA-open-alpha-rarefaction-dada2.qzv

