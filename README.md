---
title: "T. cruzi Y strain infecting H. sapiens: UTR Length and Alternative Trans-splicing / Poly-adenylation Analysis (72hrs)"
output:
  knitrBootstrap::bootstrap_document:
    theme: flatly
---



<div style='font-size:40px;'>T. cruzi Y strain infecting H. sapiens: UTR Length and Alternative Trans-splicing / Poly-adenylation Analysis (72hrs)</div>

[view source](README.rmd)

<a href='mailto:khughitt@umd.edu'>Keith Hughitt</a> (<time>2016-03-22</time>)

Overview
========

The goal of this analysis is to parse the output from our [UTR analysis
pipeline](https://github.com/elsayed-lab/utr_analysis), analyze spliced leader
acceptor site and poly-adenylation site counts, and determine the most likely
primary UTR boundaries for each gene where information is available.

5'- and 3'-UTR boundaries, along with some other useful information such as UTR
GC- and CT-richness are outputted in a format that is convenient for downstream
analysis.

Finally, we will look for evidence of either alternative trans-splicing or
poly-adenylation and attempt to visualize the prevalence and magnitude of
these events across the different developmental stages.

This analysis is based on a version originally written and applied separately
to *T. cruzi* and *L. major* in 2014. The code has been refactored and cleaned
up to make it easier to apply to new datasets.

Settings
========



```r
# species information
parasite = 'T. cruzi strain Y strain'
host     = 'H. sapiens'
stages  = c('trypomastigote', 'amastigote', 'epimastigote')

# Output verbosity
verbose = FALSE

# Number of individual UTR plots to display for each developmental stage and
# UTR side. The purpose of these plots are to provide a sense of the data
# distribution and to visualize the effect of primary # site selection with and
# without smoothing.
max_plots = 10

# output file prefix
outfile_prefix = 'tcruzi_ystrain'

# fasta/gff
input_gff   = file.path(Sys.getenv('REF'), 'tcruzi_clbrener_esmeraldo-like',
                        '/annotation/TriTrypDB-8.1_TcruziCLBrenerEsmeraldo-like.gff')
input_fasta = file.path(Sys.getenv("REF"), 
                        'tcruzi_clbrener_esmeraldo-like/genome/TriTrypDB-8.1_TcruziCLBrenerEsmeraldo-like_Genome.fasta')

# spliced leader / polyadeynlation sites
# Taken from T. cruzi UTR analysis output from 2015/07
sl_input = list(
    amastigote='input/tcruzi_infecting_hsapiens_amastigote_sl_sorted.gff.gz',
    trypomastigote='input/tcruzi_infecting_hsapiens_trypomastigote_sl_sorted.gff.gz',
    epimastigote='input/tcruzi_infecting_hsapiens_epimastigote_sl_sorted.gff.gz'
)

polya_input = list(
    amastigote='input/tcruzi_infecting_hsapiens_amastigote_polya_sorted.gff.gz',
    trypomastigote='input/tcruzi_infecting_hsapiens_trypomastigote_polya_sorted.gff.gz',
    epimastigote='input/tcruzi_infecting_hsapiens_epimastigote_polya_sorted.gff.gz'
)

# Non-coding RNA filters
id_filter_string   = 'rRNA|5SRRNA|snRNA|snoRNA|SLRNA|TRNA|SRP'
type_filter_string = 'rRNAencoding|snRNAencoding|snoRNAencoding|tRNAencoding'
```



Methods
=======

## Libraries



## Settings


## Load sequence and annotations


```r
# Load genome sequence and annotations
gff = import.gff3(input_gff)
fasta = readDNAStringSet(input_fasta)

chromosomes = gff[gff$type %in% c('contig', 'chromosome')]
genes       = gff[gff$type == 'gene']

# Copy input GFF to output directory; we will create a modified version of
# the GFF which includes UTR boundary coordinates
file.copy(input_gff, outdir, overwrite=TRUE)
```

```
## [1] TRUE
```

```r
Sys.setenv(combined_gff=file.path(outdir, basename(input_gff)))

# Fix FASTA names
# e.g. "LmjF.24 | organism=Leishmania_major_strain_Friedlin |..." -> "LmjF.24"
names(fasta) = sapply(strsplit(names(fasta), ' | '), function(x) {x[1]})

# L. major: load unannotated ORFs detected from ribosome profiling data
if (parasite == 'L. major') {
    orfs      = gff[gff$type == 'ORF']
    orfs$Name = orfs$ID
    orfs$description = 'Unannotated ORF'

    # 2015/06/30 A few of the unannotated ORFs appear to have multiple conflicting
    # entries -- removing these for now...
    orfs = orfs[!duplicated(orfs$Name)]

    # Drop GFF columns not shared between TriTrypDB GFF and uORF GFF
    keep_cols = intersect(colnames(mcols(genes)), colnames(mcols(orfs)))
    genes = genes[,keep_cols]
    orfs = orfs[,keep_cols]

    genes = append(genes, orfs)

    # Fix names (L. major chromosome identifiers)
    names(fasta) = substring(names(fasta), 0, 7)
}
```

## Remove noncoding RNAs


```r
# filter noncoding RNAs
noncoding_ids = genes$ID[grepl(id_filter_string, genes$ID)]

ncrna_mask = rep(FALSE, length(genes))

if (!is.null(id_filter_string)) {
    ncrna_mask = ncrna_mask | grepl(id_filter_string, genes$ID)
}

if (!is.null(type_filter_string)) {
    ncrna_mask = ncrna_mask | grepl(type_filter_string, as.character(genes$type))
}

# remove noncoding RNAs
genes = genes[!ncrna_mask,]

# remaining gene ids
gene_ids = genes$Name
```

## Load spliced leader and poly(A) site coordinates


```r
sl_gffs = list()
polya_gffs = list()

# names of developmental stages being analyzed
stages = names(sl_input)

for (stage in stages) {
    # SL acceptor sites
    sl_gffs[stage] = import.gff3(gzfile(sl_input[[stage]]))
    sl_gffs[stage] = sl_gffs[[stage]][sl_gffs[[stage]]$Name %in% gene_ids]

    # Poly(A) sites
    polya_gffs[stage] = import.gff3(gzfile(polya_input[[stage]]))
    polya_gffs[stage] = polya_gffs[[stage]][polya_gffs[[stage]]$Name %in% gene_ids]
}

# Normalize seqinfo (e.g. when reads cover differing sets of contigs)
sl_seqinfo = seqinfo(sl_gffs[[stages[1]]])
polya_seqinfo = seqinfo(polya_gffs[[stages[1]]])

if (length(stages) > 1) {
    for (stage in stages[2:length(stages)]) {
        sl_seqinfo = suppressWarnings(merge(sl_seqinfo, seqinfo(sl_gffs[[stage]])))
        polya_seqinfo = suppressWarnings(merge(polya_seqinfo, seqinfo(polya_gffs[[stage]])))
    }
}

for (stage in stages) {
    # Add levels that aren't already represented
    sl_levels = seqlevels(sl_seqinfo)[!seqlevels(sl_seqinfo) 
                                        %in% seqlevels(sl_gffs[[stage]])]
    seqlevels(sl_gffs[[stage]]) = c(seqlevels(sl_gffs[[stage]]), sl_levels)

    # Reorder levels
    seqlevels(sl_gffs[[stage]]) = sort(seqlevels(sl_gffs[[stage]]))

    # Add levels that aren't already represented
    polya_levels = seqlevels(polya_seqinfo)[!seqlevels(polya_seqinfo) 
                                        %in% seqlevels(polya_gffs[[stage]])]
    seqlevels(polya_gffs[[stage]]) = c(seqlevels(polya_gffs[[stage]]), 
                                       polya_levels)

    # Reorder levels
    seqlevels(polya_gffs[[stage]]) = sort(seqlevels(polya_gffs[[stage]]))
}
```

## Methods and Results

### 5'UTR
























