##########################################################################################################################
## Functions to be called from algorithm-specific scripts to produce and deal with result DGE lists
##########################################################################################################################

##-------------------------------------------------------------------------------------------------------------------------
extractCapturedSubstrings <- function(pattern, string) {
    tmp <- attributes(regexpr(pattern,string,perl=TRUE))[c("capture.start","capture.length")]
    substring(string,first=as.vector(tmp$capture.start),last=as.vector(tmp$capture.start+tmp$capture.length-1))
}
##-------------------------------------------------------------------------------------------------------------------------
import_kallisto <- function(kallisto_output_files,
                            design_info,
                            design,
                            tx2gene,
                            aggregation_column="ens_gene",
                            downstream_DGE_algo="sleuth") {
    
    n <- attributes(design)$dimnames[[1]]
    sample_to_covariates <- data.frame(sample=n,
                                    design_info[n,names(attributes(design)$contrast)],
                                    path=kallisto_output_files[n])
    
    if        (downstream_DGE_algo=="sleuth") {
        out <- sleuth_prep(sample_to_covariates,
                           target_mapping = tx2gene_versioned,
                           aggregation_column = aggregation_column,
                           extra_bootstrap_summary = TRUE)
        
    } else if (downstream_DGE_algo=="limma") {
         txi <- tximport(kallisto_output_files[n],
                         type="kallisto",
                         tx2gene=tx2gene,
                         countsFromAbundance="lengthScaledTPM")
         out <- DGEList(txi$counts)
         
    } else if (downstream_DGE_algo=="edgeR") {
        txi <- tximport(kallisto_output_files[n],
                        type = "kallisto",
                        tx2gene=tx2ensembl_versioned,
                        countsFromAbundance="no")
        
        ## use the gene-level estimated counts from the quantification tool (here: kallisto)
        ## (the summarized count is simply the sum of the per-tx counts)
        ## then additionally use the transcript-level abundance estimates to calculate a gene-level offset
        ## that corrects for changes to the average transcript length across samples:
        
        cts_gene <- txi$counts
        normMat_gene <- txi$length
        normMat_gene <- normMat_gene /exp(rowMeans(log(normMat_gene ))) ## sample-specific deviation from the average summarized length(over samples)
        normCts_gene  <- cts_gene/normMat_gene 
            
        ## Computing effective library sizes from scaled counts, to account for composition biases between samples.
        eff.lib_gene <- calcNormFactors(normCts_gene) * colSums(normCts_gene)
            
        ## Combining effective library sizes with the length factors, and calculating offsets for a log-link GLM.
        normMat_gene <- sweep(normMat_gene, 2, eff.lib_gene, "*")
        ## this command multiplies column i of normMat_gene by  eff.lib_gene[i]
            
        normMat_gene <- log(normMat_gene)
            
        out <- DGEList(cts_gene)
        out <- scaleOffset(out, normMat_gene)
        ## scaleOffset ensures that the scale of offsets are consistent with library sizes.
        ## This is done by ensuring that the mean offset for each gene is the same as the mean log-library size.
        ## The length or dimensions of offset should be consistent with the number of libraries in y.
        

    } else {
        cat("Unknown downstream_DGE_algo in import_kallisto: returning NULL!\n")
        out <- NULL
    }
    out
    
}
##-------------------------------------------------------------------------------------------------------------------------
design_matrix <- function(tbl,            ## a data.frame, expected to minimally contain:
                          R1_col, R2_col, ## names of input columns holding the fastq(.gz) filenames of Read1 and Read2
                                          ## (Read2=NA -> single end)
                          fqname_pattern, ## a regular expression, which
                                          ## -- describes the structure of the fastq filenames in R1_col and R2_col,
                                          ## -- identifies, in a single sub-expression captured by (), a substring in R1_col (and R2_col, if not NA), which
                                          ## ----- must be the same in R1_col and R2_col of a given input row (unless R2_col=NA)
                                          ## ----- will be used as rowname in the output (to be interpreted as a sample identifier)
                                          ## if is.na(fqname_pattern) --> no rownames on output
                        
                          label_col,      ## name of a column holding an experiment-specific sample id (not derived from the fastq filenames)
                                          ## if is.na(label_col) --> output[,"label"] := 1:(nrow(input))]
                          factor_baselevels ## a named character vector, where
                                            ## -- element names :=  input column names holding the factors of the the experimental design
                                            ## -- elements      :=  the desired baseline level for each factor
                          ) {

    ## Extract and name the output columns
    if(is.na(label_col)) {smpl <- 1:nrow(tbl)} else {smpl <- tbl[,label_col]} 
    out <- data.frame(Sample=smpl,tbl[,colnames(tbl) %in% names(factor_baselevels)])
    colnames(out)[2:ncol(out)] <- colnames(tbl)[colnames(tbl) %in% names(factor_baselevels)]

    ## Set the output rownames
    if(!is.na(fqname_pattern)) {
        v <- extractCapturedSubstrings(fqname_pattern,tbl[,R1_col])
        if(!is.na(R2_col)) { ## paired-end data -- check whether the captured substring identifies both R1 and R2
            if(identical(v, extractCapturedSubstrings(fqname_pattern,tbl[,R2_col]))) {
                rownames(out) <- v
            } else {
                rownames(out) <- NULL ## if not, don't use it as a sample id (rowname)!
            }
        } else {             ## single-end data
            rownames(out) <- v
        }
    } 

    ## Set the factor reference levels as specified by the "factor_baselevels" parameter
    for(cl in names(factor_baselevels)) {
        out[,cl] <- factor(out[,cl],
                           levels=c(factor_baselevels[cl], setdiff(out[,cl],factor_baselevels[cl])) )
    }

    ## Return result
    out
}
## test
##tbl <- read.table("/data/public/ugoebel/Analysis/Kononenko/Project_Overhoff/overhoff_design.csv",
##                  sep="\t",quote="",header=TRUE)
##design_matrix(tbl=tbl,
##              R1_col="Filename1", R2_col="Filename2",
##              fqname_pattern="/data/public/ugoebel/Data/Kononenko/Project_Overhoff/(A\\d+_\\d+_S\\d+)_R[12]\\.fq\\.gz$",
##              label_col="SampleLabel",
##              factor_baselevels=c(Genotype="WT",Batch="A006200061"))
            
            
##-------------------------------------------------------------------------------------------------------------------------
DGEoutput2Enrichr <- function(d, rule="QVAL", qval_col="Adjusted.P.value", t_col="t", logFC_col="logFC", gene_col="__ROWNAMES",
                              gene_map=NULL, N_t=500, maxQ=0.05, cmp=">",geneset_dbs=biojupies_dbs) {
    ## d: the result data.frame of a differential gene expression program
    ## qval_col:  the column holding the adjusted p-value (or any other measure to be handled by the QVAL   rule) 
    ## t_col:     the column holding the t-value          (or any other measure to be handled by the TRULE  rule)
    ## logFC_col: the column holding the effect strength AND EFFECT DIRECTION (logFC, or the estimated model coefficient)
    ## rule==QVAL:   return genes (with d[,qval_col]<=maxQ) & get(cmp)(d[,logFC_col],0)]
    ##                                                      ## The get() function takes a string as input and returns a reference
    ##                                                      ## to the function corresponding to that string, if it exists.
    ##                                                      ## Used here to pass the comparison operator (to 0) as a parameter.
    ## rule==TRULE:  return the N_t genes with the highest (cmp==">") or lowest (cmp=="<") values in d[,t_col]
    ##               (used by BioJupies on the limma result tables, to generate input for Enrichr)
    
    ## gene_col:  the column holding the gene names to be used in the output
    ##            If gene_col=="__ROWNAMES" --> use rownames(d)
    ##
    ## gene_map:  a vector v with v[gene_col]=external which maps gene names in d to mouse MGI names. If NULL, use as is.

    if(gene_col=="__ROWNAMES") {
        genes <- rownames(d)
    } else {
        genes <- d[,gene_col]
    }
    if(!is.null(gene_map)) {
        genes <- gene_map[genes]
    }
    if("THIS_GENE_NAMES" %in% colnames(d)) {
        cat("*** Column name conflict in bioJupiesList2enrichR -- returning NULL!")
    } else {
        d$THIS_GENE_NAMES <- genes  ## append a column with the output gene names
    }
    
    if(rule=="QVAL") {
        l <- unique(d$THIS_GENE_NAMES[(d[,qval_col]<=maxQ) & get(cmp)(d[,logFC_col],0)])
    } else if (rule=="TRULE") { 
        l <- head(unique(d[order(d[,t_col],decreasing=(cmp==">")),]$THIS_GENE_NAMES),
                  n=500)
    }

    list(e=enrichr(l,geneset_dbs),l=l)
}

##-------------------------------------------------------------------------------------------------------------------------
make_featureCount_script <- function(bam_base="/data/public/ugoebel/Analysis/Nephrolab/Project_Fabretti/Subread-align",
                                     bam_pattern="^A\\d+_\\d+_S\\d+_L004_trimmed_plusERCC92.bam$",
                                     sampleName_pattern="A\\d+_\\d+_S\\d+",
                                     annotation_base="/data/public/ugoebel/Data/Human",
                                     annotation_file="Homo_sapiens.GRCh38.101.ERCC.gtf.gz",
                                     samtools_params=c(q=30,f=3,F=256), ## if NULL, no pre-filtering of the bamfiles by samtools is done
                                     outfile_prefix="featureCounts",
                                     scriptfile_prefix="run_featureCounts",
                                     NCORES=3
                                     ) {
    bamfiles <- list.files(path=bam_base,pattern=bam_pattern)
    names(bamfiles) <- sapply(bamfiles,
                              extractCapturedSubstrings,
                              pattern=paste("(",sampleName_pattern,")",sep=""))
    
    infiles <- paste(bam_base,bamfiles,sep="/")
    if(is.null(samtools_params)) {
        samtools_calls <- ""
    } else {
        outfiles <- paste("FILTERED_", paste(names(samtools_params),samtools_params,sep="",collapse="_"),
                          "_",bamfiles,sep="")
        samtools_calls <- paste("samtools view -hb ",
                                paste(paste("-",names(samtools_params), " ", samtools_params, sep=""),collapse=" "), " ",
                                "-o ", outfiles, " ", infiles,"\n",
                                "samtools index ", outfiles, sep="")
        infiles <- outfiles
    }
    outfiles <- paste(outfile_prefix,"_",sub("\\.bam$","",infiles),".csv",sep="")
    
    cmds <- paste(samtools_calls,"\n",
                  "featureCounts -a ", annotation_base, "/", annotation_file, " ",
                  "-o ", outfiles, " ",
                  "-F \"GTF\" -t \"exon\" -g \"gene_id\" ",
                  "--minOverlap 20 ",
                  "-M --primary ", ## Multi-mapping reads will also be counted; (but) Count primary alignments only
                  "-O ", ## Assign reads to all their overlapping meta-features
                  "-J ", ## Count number of reads supporting each exon-exon junction
                  ##"-G " ## FASTA-format file that contains the reference sequences used in read mapping (improves -J; optional)
                  "-T ",NCORES, " --verbose ",
                  infiles, "\n",
                  sep="")
    
    fn <- paste(scriptfile_prefix, "_",
                "NCORES",NCORES,"_",
                ifelse(is.null(samtools_params),"",paste(names(samtools_params),samtools_params,sep="",collapse="_")),
                ".csh",sep="")
    write('#!/usr/bin/csh',fn)
    write(cmds,fn, append=TRUE)
}

##-------------------------------------------------------------------------------------------------------------------------
make_featureCount_table <- function(table_dir,table_pattern,file_pattern="A\\d+_\\d+_S\\d+",
                                    meta_feature=TRUE) {

    
    #####first read all tables and take the union of all referenced (meta-)features
    files <- list.files(path=table_dir,pattern=table_pattern,full.names=TRUE)

    L <- NULL
    for(fn in files) {
        tbl <- read.table(fn, sep="\t",quote="",header=TRUE)
        ## NOTE: Length is the total number of non-overlapping bases in the focal feature or meta-feature;
        ##       when counting at the meta-feature (gene) level, then the semicolon-separated [Start..End] ranges
        ##       of the constituent features may overlap, so sum(end[i]-start[i]+1) may be != Length.
        ##       *** For now, always count at the meta-feature level ***, so the the relevant value is Length.

        if(is.null(L)) {
            L <- tbl$Length
            names(L) <- tbl$Geneid
        } else {
            i <- which(  tbl$Geneid %in% names(L))
            if(!all(tbl$Length[i] == L[tbl$Geneid[i]])) {
                cat("*** Incongruent feature lengths in make_featureCount_table -- returning Null!\n")
                return(NULL)
            }

            i <- which(!(tbl$Geneid %in% names(L)))
            if(length(i)>0) {
                cat("make_featureCount_table: merging (meta-)feature sets of different tables!\n")
                l <- tbl$Length[i]
                names(l) <- tbl$Geneid[i]
                L <- c(L,l)
                rm(l); gc()
            }
        }
    }
    

    #####then record the counts of each (meta-)feature in each input table
    output <- data.frame(Length=L)
    rownames(output) <- names(L)  
    
    for(fn in files) {
        tbl <- read.table(fn, sep="\t",quote="",header=TRUE)
        cols <- colnames(tbl)[grep(file_pattern, colnames(tbl),perl=TRUE)]
        ptn <- paste("(",file_pattern,")",sep="")
        names(cols) <- sapply(cols,function(x)extractCapturedSubstrings(ptn,x))
        
        for(n in names(cols)) {
            v <- tbl[,cols[n]]; names(v) <- tbl$Geneid
            output[,n] <- 0; output[names(v),n] <- v
        }
    }
    ##### re-order to have the Length column at the end
    output <- cbind(output[,-1], Length=output[,1])

    ##### for testing: write the table, only to read it in later from the main script:
    ##write.table(output,file=paste(table_dir,"featureCounts_subread-align_q30__trimmed_plusERCC92.csv",sep="/"),
    ##            sep="\t", quote=FALSE,row.names=TRUE,col.names=TRUE)

    output
}
