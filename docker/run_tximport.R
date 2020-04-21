library(tximport)

# function used to extract the sample name from the full
# path provided as a command line arg
getSampleName <- function(x, suffix){
    b <- basename(x)
    return(substr(b, 1, nchar(b) - nchar(suffix) - 1))
}

args <- commandArgs(TRUE)
argsLen <- length(args)
TRANSCRIPT_TO_GENE_MAPPING_FILE <- args[1] # file mapping the transcript to gene names. Needs column headers of `enst` and `name`
SUFFIX <- args[2] # the ending we strip off to get the sample name.  e.g. for sampleA.abundance.tsv, SUFFIX='abundance.tsv'
OUTPUT_TABLE_BASE <- args[3] # the name of the output file to write the counts table
IS_PDX <- args[4]
HUMAN_TAG <- args[5]
MOUSE_TAG <- args[6]
files <- args[7:argsLen] # the remainder of the args are paths to the files

# cast the commandline arg to a boolean
if(IS_PDX == 'true'){
    IS_PDX = TRUE
} else {
    IS_PDX = FALSE

}

# tximport expects a list where the elements are named by the sample
# TODO: guard against "weird" file names
sample_names <- sapply(files, getSampleName, suffix=SUFFIX)
names(files) <- sample_names

mapping = read.table(TRANSCRIPT_TO_GENE_MAPPING_FILE, sep='\t', header=T, stringsAsFactors=F)
txi.kallisto <- tximport(files, 
    type='kallisto', 
    tx2gene=mapping[c('enst','name')], 
    ignoreTxVersion=T,
    countsFromAbundance="lengthScaledTPM"
)

counts <- as.data.frame(txi.kallisto$counts)

if(IS_PDX){
    human_genes = unique(mapping[mapping['human_or_mouse']=='human', 'name'])
    mouse_genes = unique(mapping[mapping['human_or_mouse']=='mouse', 'name'])

    human_counts <- na.omit(counts[human_genes,])
    mouse_counts <- na.omit(counts[mouse_genes,])

    write.table(
        human_counts,
        file=paste(OUTPUT_TABLE_BASE, HUMAN_TAG, 'tsv', sep='.'),
        quote=F,
        sep='\t'
    )

    write.table(
        mouse_counts,
        file=paste(OUTPUT_TABLE_BASE, MOUSE_TAG, 'tsv', sep='.'),
        quote=F,
        sep='\t'
    )
} else {
    
    write.table(
        counts,
        file=paste(OUTPUT_TABLE_BASE, 'tsv', sep='.'),
        quote=F,
        sep='\t'
    )
}

