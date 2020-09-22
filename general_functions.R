##########################################################################################################################
## Functions to be called from algorithm-specific scripts to produce and deal with result DGE lists
##########################################################################################################################

##-------------------------------------------------------------------------------------------------------------------------
extractCapturedSubstrings <- function(pattern, string, global_search=FALSE) {
    if(global_search) {
        tmp <- attributes(gregexpr(pattern,string,perl=TRUE)[[1]])[c("capture.start","capture.length")] ## only 1 input string -->> [[1]]
    } else {
        tmp <- attributes( regexpr(pattern,string,perl=TRUE)     )[c("capture.start","capture.length")]
    }
    substring(string,first=as.vector(tmp$capture.start),last=as.vector(tmp$capture.start+tmp$capture.length-1))
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
