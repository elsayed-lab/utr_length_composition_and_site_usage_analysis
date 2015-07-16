#
# find_peak
#
find_peak = function(gene, acceptor_sites, gene_strand, feature_type,
                     smoothed=FALSE, include_plot=FALSE, secondary=FALSE) {
    # Determine orientation of feature relative to CDS
    if ((feature_type == 'sl'    && gene_strand == '+') ||
        (feature_type == 'polya' && gene_strand == '-')) {
        feature_side = 'left'
    } else {
        feature_side = 'right'
    }

    # Create a vector from the furthest SL site to the CDS boundary
    if (feature_side == 'left') {
        cds_boundary = start(gene)
        raw_scores = rep(0, cds_boundary - min(start(acceptor_sites)))
        rel_start =  cds_boundary - start(acceptor_sites) + 1
    } else {
        cds_boundary = end(gene)
        raw_scores = rep(0, max(start(acceptor_sites)) - cds_boundary)
        rel_start = start(acceptor_sites) - cds_boundary + 1
    }

    # add scores at indices where acceptor sites were detected
    raw_scores[rel_start] = score(acceptor_sites)
    x = 1:length(raw_scores)

    # create smooted version of scores (optional)
    if (smoothed) {
        input_scores = ksmooth(x, raw_scores, kernel="normal", bandwidth=6,
                               n.points=length(x))$y
    } else {
        input_scores = raw_scores
    }

    # find acceptor site peak
    if (secondary && (length(acceptor_sites) >= 2)) {
        # secondary site
        j = 2
    } else {
        # primary site
        j = 1
    }
    sorted_scores = sort(input_scores, decreasing=TRUE)
    raw_idx       = which(input_scores == sorted_scores[j])

    if (length(raw_idx) > 1) {
        if (feature_side == 'left') {
            raw_idx = head(raw_idx, 1)
        } else {
            raw_idx = tail(raw_idx, 1)
        }
    }

    # find raw score (number of reads) for the desired primary or secondary
    # site peak
    raw_score = raw_scores[raw_idx]

    # If smoothing was used, the actual location of the smoothed peak may not
    # have any support. In this case we will find the nearest site to the peak
    # with a non-zero score in the raw scores vector
    while(raw_score == 0) {
        # find next highest peak
        j = j + 1 
        raw_idx = which(input_scores == sorted_scores[j])

        # if more than one peak, choose the furthest one from CDS
        if (length(raw_idx) > 1) {
            if (feature_side == 'left') {
                raw_idx = head(raw_idx, 1)
            } else {
                raw_idx = tail(raw_idx, 1)
            }
        }
        raw_score = raw_scores[raw_idx]
    }

    # plot raw and smoothed peaks
    if (include_plot) {
        if (smoothed) {
            df = melt(data.frame(x, raw=raw_scores, smoothed=input_scores),
                      id=c("x"), variable.name='type')
            print(qplot(x, value, data=df, color=type, geom='line') 
                  + ggtitle(gene$ID))
        } else {
            df = melt(data.frame(x, raw=raw_scores), id=c("x"),
                    variable.name='type')
            print(qplot(x, value, data=df, geom='line')
                  + ggtitle(gene$ID))
        }
    }

    # find index of the site with the highest score (or second highest score,
    # in the case of secondary site selection)
    idx = which(score(acceptor_sites) == raw_score)

    # if there is a tie, use the furthest site from CDS
    if (length(idx) > 1) {
        if (feature_side == 'left') {
            idx = head(idx, 1)
        } else {
            idx = tail(idx, 1)
        }
    }
    return(idx)
}

#
# find_primary_site
#
# Determines the primary site among a list of acceptor sites and their
# associated score (number of reads mapped).
#
# Returns a list containing the site and score for the primary site.
#
find_primary_site = function(acceptor_sites, feature, gene, gene_strand, smoothed=FALSE) {
    if (length(acceptor_sites) == 0) {
        return(list(location=NA, num_reads=NA))
    } else if (length(acceptor_sites) == 1) {
        # if only one site found, use it
        return(list(location=start(acceptor_sites)[1],
                    num_reads=score(acceptor_sites)[1]))
    } else {
        # two or more sites

        # find highest smoothed peak which has non-zero coverage in the raw
        # data as well
        max_index = find_peak(gene, acceptor_sites, gene_strand, feature,
                              smoothed=smoothed)
        max_index_smoothed = find_peak(gene, acceptor_sites, gene_strand,
                                       feature, smoothed=TRUE)

        # if the smoothing would result in a different primary UTR choice,
        # plot the raw and smoothed versions of the data
        if (abs(max_index - max_index_smoothed) > 15) {
            num_diff = num_diff + 1
            if (plot_num <= max_plots) {
                sprintf("PLOTTING %d/%d", plot_num, max_plots)
                #tmp = find_peak(gene, acceptor_sites, gene_strand, feature,
                #                smoothed=TRUE, include_plot=TRUE)
                plot_num = plot_num + 1
            }
        }

        return(list(location=start(acceptor_sites)[max_index], 
                    num_reads=score(acceptor_sites)[max_index]))
    }
}

#
# find_secondary_site
#
# Determines the secondary site among a list of acceptor sites and their
# associated score (number of reads mapped).
#
# Returns a list containing the site and score for the secondary site.
#
find_secondary_site = function(acceptor_sites, feature, gene, gene_strand,
                               smoothed=FALSE) {
    if (length(acceptor_sites) == 0) {
        return(list(location=NA, num_reads=NA))
    } else if (length(acceptor_sites) == 1) {
        # if only one site found, use it
        return(list(location=start(acceptor_sites)[1],
                    num_reads=score(acceptor_sites)[1]))
    } else {
        # two or more sites

        # find highest smoothed peak which has non-zero coverage in the raw
        # data as well
        secondary_site_index = find_peak(gene, acceptor_sites, gene_strand,
                                         feature, smoothed=smoothed,
                                         secondary=TRUE)
        secondary_site_index_smoothed = find_peak(gene, acceptor_sites,
                                                  gene_strand, feature,
                                                  smoothed=TRUE,
                                                  secondary=TRUE)

        # if the smoothing would result in a different secondary UTR choice,
        # plot the raw and smoothed versions of the data
        #if (abs(secondary_site_index - secondary_site_index_smoothed) > 15) {
        #    num_diff = num_diff + 1
        #    if (plot_num <= max_plots) {
        #        tmp = find_peak(gene, acceptor_sites, gene_strand, feature,
        #                        smoothed=TRUE, include_plot=TRUE,
        #                        secondary=TRUE)
        #        plot_num = plot_num + 1
        #    }
        #}

        return(list(location=start(acceptor_sites)[secondary_site_index], 
                    num_reads=score(acceptor_sites)[secondary_site_index]))
    }
}

#
# get_utr_sequences
#
# Returns a vector of Biostrings instances containing the UTR sequence for each
# input gene.
#
get_utr_sequences = function(genes, fasta, utr_lengths, default_width,
                             utr5=TRUE) {
    # retrieve utr length if known
    widths = data.frame(
        id=genes$ID,
        width=NA)

    widths$width = utr_lengths[match(genes$ID, utr_lengths$name),]$length

    # For UTRs of unknown length, set to 0 for now; will set back to NA later
    unknown_utrs = is.na(widths$width)
    widths$width[unknown_utrs] = 0

    # get positive and negative strand genes
    if (utr5) { 
        utr = flank(genes, widths$width)
    } else {
        utr = flank(genes, widths$width, start=FALSE)
    }
    pos_strand = utr[as.character(strand(utr)) == "+"]
    neg_strand = utr[as.character(strand(utr)) == "-"]

    # for genes that were assigned the default UTR size, make sure that the
    # assigned boundaries fall within the chromosome
    start(pos_strand) = pmax(1, start(pos_strand))
    start(neg_strand) = pmax(1, start(neg_strand))

    end(pos_strand) = pmin(end(pos_strand), width(fasta[seqnames(pos_strand)]))
    end(neg_strand) = pmin(end(neg_strand), width(fasta[seqnames(neg_strand)]))

    seqs = fasta[pos_strand]
    seqs = append(seqs, reverseComplement(fasta[neg_strand]))
    names(seqs) = c(pos_strand$ID, neg_strand$ID)

    return(seqs)
}

#
# get_num_reads
#
# Returns the number of reads mapped to a given position for the specified
# stage, or 0 if none are found
#
get_num_reads = function(utr_reads, site) {
    if (is.na(site)) {
        return(NA)
    }
    num_reads = score(utr_reads[start(utr_reads) == site])
    return(ifelse(length(num_reads) == 0, 0, num_reads))
}

#
# get_utr_composition
#
# Creates a dataframe containing the GC richness and CT richness for each
# input sequence.
#
get_utr_composition = function (seqs) {
    freqs = alphabetFrequency(seqs)

    df = data.frame(
        name=names(seqs),
        gc=round((freqs[,'G'] + freqs[,'C']) / rowSums(freqs[,1:4]), 3),
        ct=round((freqs[,'C'] + freqs[,'T']) / rowSums(freqs[,1:4]), 3) 
    )
    
    # If length is unknown or sequence contains mostly N's, make NA
    na_ind = (freqs[,'N'] / width(seqs)) > 0.5
    na_ind[is.na(na_ind)] = TRUE

    df[na_ind, 'gc'] = NA
    df[na_ind, 'ct'] = NA
    df[width(seqs) == 0, 'gc'] = NA
    df[width(seqs) == 0, 'ct'] = NA

    return(df)
}

# rescale to range 0-1
# if trimmed=TRUE, outliers will be trimmed first to reduce their influence on
# the resulting colorscale
rescale = function (x, trimmed=TRUE, trim_quantile=0.99) {
    if (trimmed) {
        lower_bound = quantile(x, 1 - trim_quantile)
        upper_bound = quantile(x, trim_quantile)
   
        x = pmax(pmin(x, upper_bound), lower_bound)

        return(pmax(0, ((x-min(x)) / (max(x) - min(x)))))
    }
    # no trimming
    return(pmax(0, ((x-min(x)) / (max(x) - min(x)))))
}
