# The file contains functions that will be used throughout all analyses

#----01_PR1_identification----
# Get AA sequences from Uniprot and save them into an object of class AAStringSet
uniprot2stringset <- function(ID=NULL) {
    seqs <- Rcpi::getSeqFromUniProt(ID)
    seqs_char <- NULL
    for (i in 1:length(seqs)) {
        seqs_char[i] <- seqs[[i]][[1]]
    }
    stringset <- Biostrings::AAStringSet(seqs_char)
    return(stringset)
}

# Perform BLAST search
blast_search <- function(query=NULL, db=NULL, dbtype="prot", type="blastp",
                         BLAST_args = "-evalue 10e-10", 
                         dbname="Gmax_proteome_db") {
    # Write FASTA
    if(!dir.exists(dbname)) {
        dir.create(dbname)
        Biostrings::writeXStringSet(db, filepath = file.path(dbname, "proteome.fa"))
        # Make db
        rBLAST::makeblastdb(file.path(dbname, "proteome.fa"), dbtype = dbtype)
    }
    # Open db
    database <- rBLAST::blast(file.path(dbname, "proteome.fa"), type = type)
    # Perform search
    search <- predict(database, query, BLAST_args = BLAST_args)
    return(search)
}

# Filter by PFAM/SMART IDs
find_domains <- function(ids=NULL, annotation=NULL, pfam=NULL, smart=NULL) {
    filtannot <- annotation[annotation[,1] %in% ids, c("Gene.Name", "PFAM.ID", "SMART.ID")]
    filtannot <- filtannot[filtannot$PFAM.ID %in% pfam | filtannot$SMART.ID %in% smart, ]
    filt_annot <- filtannot[!duplicated(filtannot$Gene.Name), ]
    return(filt_annot)
}


#----03_evolutionary_analyses----

# Process the output from universalmotif::run_meme and create a data frame of motif ranges
meme2ranges <- function(meme_output) {
    mlist <- lapply(1:length(meme_output$sites), function(x) {
        return(as.data.frame(meme_output$sites[[x]]@ranges))
    })
    names(mlist) <- names(meme_output$sites)
    
    # Add a column with motif name to each element of the list
    mlist2 <- lapply(1:length(mlist), function(x) {
        y <- cbind(mlist[[x]], Motif_name=names(mlist)[x])
        return(y)
    })
    
    # Reduce list to a data.frame
    final.ranges <- Reduce(rbind, mlist2)
    final.ranges <- final.ranges[, c(2,3,5,6)]
}


#----05. Transcriptomic analyses----

fTau <- function(x) {
    if(all(!is.na(x))) {
        if(min(x, na.rm=TRUE) >= 0) {
            if(max(x)!=0) {
                x <- (1-(x/max(x)))
                res <- sum(x, na.rm=TRUE)
                res <- res/(length(x)-1)
            } else {
                res <- 0
            }
        } else {
            res <- NA
        }
    } else {
        res <- NA
    }
    return(res)
}



#' Perform differential expression analysis using DESeq2
#'
#' @param counts Data frame containing read counts per gene/transcript.
#' @param metadata Data frame of two columns containing sample names in the first column and sample information in the second column.
#' @param alpha Significance level to be used for differential expression. Default is 0.05.
#' @param fc_cutoff Minimum fold-change cutoff to be used for differential expression. Default is 1, which means no fold change cut-off.
#' @param comparison Character vector of length 2 describing the level for the case condition and the level for the control condition. The default is \code{c("case", "control")}, which indicates the levels "case" and "control" of the data frame metadata.
#' @param plot_MA Logical indicating whether to save an MA plot or not.
#'
#' @return List containing: \itemize{
#'   \item Up-regulated genes
#'   \item Down-regulated genes
#'   \item All differentially expressed genes
#' }
#'
DGE_deseq2 <- function(counts, metadata, alpha = 0.05, fc_cutoff = 1,
                       comparison=c("treatment", "control"), 
                       plotMA = FALSE) {
    metadata <- as.data.frame(metadata)
    counts <- counts[, colnames(counts) %in% metadata[,1]]
    colData <- data.frame(sampleName = metadata[,1], condition = metadata[,2],
                          stringsAsFactors = TRUE)
    
    #Create a DESeq2 object
    dds_object <- DESeq2::DESeqDataSetFromMatrix(countData = counts, 
                                                 colData = colData, 
                                                 design = ~ condition)
    dds_object$condition <- relevel(dds_object$condition, ref = comparison[2])
    
    #Differential expression analysis at 95% confidence level
    dds <- DESeq2::DESeq(dds_object)
    contrast <- c("condition", comparison)
    res <- as.data.frame(DESeq2::results(dds, contrast = contrast))
    res$log2FoldChange[res$padj > 0.05] <- NA
    res_df <- data.frame(Gene=rownames(res),
                         log2FC=res[, 2])
    return(res_df)
}

#'
#' @param counts Data frame of raw read counts, with gene IDs in row names and sample names in column names
#' @param contrast_list List of contrast data frames. For each data frame, first column must contain control level and second column must contain case level

parallel_DGE <- function(counts, contrast_list, samplelist) {
    
    dge_list <- lapply(seq_along(samplelist), function(class) {
        message("Starting DGE on class: ", names(samplelist)[class], "...")
        fcount <- counts[, colnames(counts) %in% samplelist[[class]]$BioSample]
        meta <- samplelist[[class]][, c(3,2)]
        
        dge <- lapply(1:nrow(contrast_list[[class]]), function(x) {
            message("Performing DGE with contrast number ", x)
            comparison <- DGE_deseq2(counts = fcount,
                                     metadata = as.data.frame(meta),
                                     comparison = as.character(rev(contrast_list[[class]][x,])))
            colnames(comparison)[2] <- paste0(names(samplelist)[class], 
                                              "-C", 
                                              x)
            return(comparison)
        })
        
        dge <- Reduce(merge, dge)
        return(dge)
    })
    
    final_dge_df <- Reduce(merge, dge_list)
    return(final_dge_df)
}








