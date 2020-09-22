options(width=200)

source("~/CECAD/Programming/DGEanalysis/dge_functions.R")
library(tximport)
library(limma)
library(edgeR) ## also needed by limma for calcNormFactors()
library(csaw)
library(sleuth)
library(enrichR)
library(pcaExplorer)
library(DESeq2) ## for preparing input for pcaExplorer


do_counts <- TRUE

####--------------- Read the design information ----------------------------------------------------------------------------------------
DGE_designfile <- "/data/public/ugoebel/Analysis/Nephrolab/Project_Fabretti/Info_on_Data_and_Experiment/fabretti_design2.csv"

design_info <- design_matrix(read.table(DGE_designfile,sep="\t",quote="",header=TRUE),
                             R1_col="FileName1", R2_col="FileName2",
                             fqname_pattern="^(A\\d+_\\d+_S\\d+)_L004_R[12]_001.fastq.gz",
                             label_col="Sample_Name",
                             factor_baselevels=c(Condition="healthy",Individual=1,Urine_Sample=1,Passage=1))

####--------------- Specify one or more specific designs for DGE tests using a linear model  -------------------------------------------
design <- model.matrix(~ Condition, design_info)

####--------------- Mappings between ensembl and MGI, for genes and transcripts  -------------------------------------------------------
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
descriptors <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version", "ensembl_gene_id",
                                             "external_gene_name", "description"), mart = mart)
descriptors <- cbind(versioned_tx=paste(descriptors$ensembl_transcript_id,
                                        descriptors$transcript_version,sep="."), descriptors)

tx2ensembl_versioned <- data.frame(TXNAME=descriptors$versioned_tx,
                                   GENEID=descriptors$ensembl_gene_id) ## expected by tximport
tx2gene_versioned    <- data.frame(target_id=descriptors$versioned_tx,
                                   ens_gene=descriptors$ensembl_gene_id,
                                   ext_gene=descriptors$external_gene_name) ## expected by sleuth

tmp <- unique(data.frame(GENEID=descriptors$ensembl_gene_id,
                         D     =descriptors$description,
                         E     =descriptors$external_gene_name)) ## info pertaining to genes, not individual tx
ensembl2description <- tmp$D; names(ensembl2description) <- tmp$GENEID
ensembl2external    <- tmp$E; names(ensembl2external)    <- tmp$GENEID
rm(tmp); gc()

## new 18.9.20
tmp <- unique(data.frame(TXID=descriptors$versioned_tx,
                         D     =descriptors$description,
                         E     =descriptors$external_gene_name)) ## info pertaining to genes, not individual tx
tx2description <- tmp$D; names(tx2description) <- tmp$TXID
tx2external    <- tmp$E; names(tx2external)    <- tmp$TXID
rm(tmp); gc()

####------------------------------------------------------------------------------------------------------------------------------------
KALLISTO_dir <- "/data/public/ugoebel/Analysis/Nephrolab/Project_Fabretti/Kallisto/"
KALLISTO_filepattern <- "(A\\d+_\\d+_S\\d+)_L004_kallisto_cutadapt_trimmed.out$"

##TABLE_dir <- "/home/goebel/CECAD/Work/PIs/Nephrolab/Project_Fabretti/Analysis/Subread-align/FeatureCounts"
TABLE_dir <- "/data/public/ugoebel/Analysis/Nephrolab/Project_Fabretti/Subread-align/FeatureCounts/"
TABLE_file <- "featureCounts_q30_f3_trimmed_plusERCC92.csv"

if(do_counts) {
    count_tbl <- read.table(paste(TABLE_dir,TABLE_file,sep="/"),sep="\t",quote="",comment.char="#",header=TRUE)

    lengths_of_counted_tags <- count_tbl[,ncol(count_tbl)]
    names(lengths_of_counted_tags) <- rownames(count_tbl)
    
    count_tbl <- count_tbl[,rownames(design_info)] ## make sure they are in the same order! -- and remove the length column
                                                   ## ASSUMES that all(colnames(count_tbl) %in% rownames(design_info)) !

    files <- NA
} else {
    count_tbl <- NA

    files <- list.files(path=KALLISTO_dir,pattern=KALLISTO_filepattern,full.names=TRUE)
    n <-  extractCapturedSubstrings(KALLISTO_filepattern,files)
    files <- paste(files,"abundance.h5",sep="/")
    names(files) <- n

    files <- files[rownames(design_info)] ## make sure they are in the same order!
}

####------------------------------------------------------------------------------------------------------------------------------------
DGE_algorithm_data <- list()

DGE_algorithm_data$columns <- list(
    limma=c(qval_col  = "adj.P.Val",
            pval_col  = "P.Value",
            t_col     = "t",
            logFC_col = "logFC",
            gene_col  = "__ROWNAMES"),
    
    edgeR=c(qval_col  = "FDR",
            pval_col  = "PValue", 
            t_col     =  NULL,
            logFC_col = "logFC",
            gene_col  = "__ROWNAMES"),

    sleuth=c(qval_col  = "qval",
             pval_col  = "pval",
             t_col     =  NULL,
             logFC_col = "b",
             gene_col  = "target_id")) ## changed 18.9.20: use versioned tx ids instead of gene ids!
             ##gene_col  = "ens_gene"))

DGE_algorithm_data$possible_selection_rules <- list(
    limma= c("QVAL","TRULE"),
    edgeR= "QVAL",
    sleuth="QVAL")

DGE_algorithm_data$gene_map <- list(
    limma= ensembl2external,
    edgeR= ensembl2external,
    sleuth=tx2external) ## changed 18.9.20: use versioned tx ids instead of gene ids!
    ##sleuth=ensembl2external)

DGE_algorithm_data$description_map <- list( ## added 18.9.20
    limma= ensembl2description,
    edgeR= ensembl2description,
    sleuth=tx2description) 

####------------------------------------------------------------------------------------------------------------------------------------

DGE_algorithms_used <- c("limma", "sleuth", "edgeR")

if(do_counts) {
    use_algos <- setdiff(DGE_algorithms_used,"sleuth")
} else {
    use_algos <- DGE_algorithms_used
}
#####-----------------------------------------------------------------------------------------------------------------------------------

DGE_output <- list()
for(algo in use_algos) {
    DGE_output[[algo]] <- list()

    if(do_counts) {
        y <- DGEList(counts=count_tbl[,rownames(design)],
                     genes=lengths_of_counted_tags,
                     group=design[,2]) ## design[,2] is meaningful only for a one-covariate desifn!
            
    } else {
        ## import_kallisto() returns
        ## for downstream_DGE_algo=="sleuth": a ‘sleuth’ object with sleuth_prep() already called upon
        ## for downstream_DGE_algo %in% c("limma", "edgeR"): a DGEList with only the "counts" field set
        y <- import_kallisto(kallisto_output_files=files,
                             design_info=design_info,
                             design=design,
                             tx2gene=tx2gene_versioned,
                             aggregation_column="ens_gene",
                             downstream_DGE_algo=algo)
    }
    
    if(algo == "limma") {
        keep <- filterByExpr(y,design); print(table(keep))
        y <- y[keep,]
        y <- calcNormFactors(y)
        
        v <- voom(y, design)
        
        fit <- lmFit(v, design)
        fit <- eBayes(fit)
        
        DGE_output[[algo]][[i_design]] <- sapply(colnames(design)[-1],
                                                 function(cf) {
                                                     topTable(fit,
                                                              number=nrow(y$counts),coef=cf)
                                                 },simplify=FALSE)
    } else if (algo == "edgeR") {
        keep <- filterByExpr(y,design); print(table(keep))
        y <- y[keep, ]
        ##################### y <- calcNormFactors(y)
        ## See https://support.bioconductor.org/p/121087/:
        ##If you have an offsets matrix in your DGEList then you won't use the norm.factors anyway,
        ##so it wouldn't matter if you did something with them or not.
        ##Put a different way, the offsets are supposed to be better than simple normalization factors,
        ##and are preferentially used by glmFit. 
        y <- estimateDisp(y, design)
        
        fit <- glmQLFit(y, design)
        
        DGE_output[[algo]][[i_design]] <- sapply(colnames(design)[-1],
                                                 function(cf) {
                                                     as.data.frame(topTags(glmQLFTest(fit,coef=cf),
                                                                           n=nrow(y$counts)))
                                                 },simplify=FALSE)
    } else if (algo == "sleuth") {
        ## NOTE that there is no low-expression filtering with sleuth!!
        
        y  <- sleuth_fit(y,  design, "dummy_design_name")
        
        DGE_output[[algo]][[i_design]] <- sapply(colnames(design)[-1],
                                                 function(cf) {
                                                     y  <- sleuth_wt(y,  which_beta=cf, which_model="dummy_design_name")
                                                     sleuth_results(y,  cf, 'wt', which_model="dummy_design_name",
                                                                    show_all = FALSE,pval_aggregate = FALSE)
                                                 },simplify=FALSE)
        
        
    }
}


#####-----------------------------------------------------------------------------------------------------------------------------------
#####============== (2) Gene Set Analysis ==============================================================================================
#####-----------------------------------------------------------------------------------------------------------------------------------

GSEA_base    <- "/data/public/ugoebel/Analysis/Niessen/Project_Persa/BioJupies" ## also expected to hold the GSEA results!

##----------------- parsers for the individual GSEA databases -------------------------------------------------------------
##----------------- NOTE: the database names are implicitly specified here, as the names of the parsers!! -----------------
##-----------------       For each parser name, a file with this name is expected in directory GSEA_base --------------------
GSEA_parsers <- list()
GSEA_parsers[["ChEA_2016"]]                      <- function(x) {
                                                       y <- toupper(strsplit(x[1],"\\s+")[[1]])
                                                       ## split the leading part on whitespece
                                                       has_organism <- y[length(y)] %in% c("HUMAN","MOUSE","RAT")
                                                       c(term=y[1],
                                                         id=y[2],
                                                         method=y[3],
                                                         cell_type=paste(y[4:(length(y)-as.numeric(has_organism))],collapse="_"),
                                                         genome=ifelse(has_organism,y[length(y)],NA))
                                                    }
##-------------------------------------------------------------------------------------------------------------------------
GSEA_parsers[["ENCODE_TF_ChIP-seq_2015"]]        <- function(x) {
                                                       y <- toupper(strsplit(x[1],"\\s+")[[1]]) 
                                                       c(term=y[1],
                                                         cell_type=paste(y[-c(1,length(y))],collapse="_"), ## pure guess ..
                                                         method="ChIP-seq",
                                                         genome=y[length(y)])
                                                    }
##-------------------------------------------------------------------------------------------------------------------------
GSEA_parsers[["GO_Biological_Process_2018"]] <-
    GSEA_parsers[["GO_Molecular_Function_2018"]] <-
    GSEA_parsers[["GO_Cellular_Component_2018"]] <- function(x) {
                                                       y <- strsplit(x[1],"\\s+")[[1]]
                                                       c(term=paste(y[1:(length(y)-1)],collapse="_"),
                                                         id=extractCapturedSubstrings("\\(([^()]+)\\)$",x[1]))
                                                    }
##-------------------------------------------------------------------------------------------------------------------------
##GSEA_parsers[["KEGG_2019_Mouse"]]                <- function(x) {
GSEA_parsers[["KEGG_2019_Human"]]                <- function(x) {
                                                       y <- strsplit(x[1],"\\s+")[[1]]
                                                       c(term=paste(y,collapse="_"))
                                                    }
##-------------------------------------------------------------------------------------------------------------------------
GSEA_parsers[["Reactome_2016"]]                  <- function(x) {
                                                       y <- strsplit(x[1],"\\s+")[[1]]
                                                       c(term=paste(y[1:(length(y)-3)],collapse="_"),
                                                         genome=paste(y[(length(y)-2):(length(y)-1)],collapse="_"), ## always?
                                                         id=y[length(y)])
                                                    }
##-------------------------------------------------------------------------------------------------------------------------
##GSEA_parsers[["WikiPathways_2019_Mouse"]]       <- function(x) {
GSEA_parsers[["WikiPathways_2019_Human"]]       <- function(x) {
                                                      y <- strsplit(x[1],"\\s+")[[1]]
                                                      c(term=paste(y[1:(length(y)-1)],collapse="_"),
                                                        id=extractCapturedSubstrings("(WP\\d+)$",x[1]))
                                                    }
##-------------------------------------------------------------------------------------------------------------------------
## results from a given group of databases will be merged
GSEA_groups <- list(TFs=c("ChEA_2016", "ENCODE_TF_ChIP-seq_2015"),
                    GO=c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018"),
                    ##pathways=c("KEGG_2019_Mouse", "Reactome_2016", "WikiPathways_2019_Mouse"))
                    pathways=c("KEGG_2019_Human", "Reactome_2016", "WikiPathways_2019_Human"))
##-------------------------------------------------------------------------------------------------------------------------

#####-----------------------------------------------------------------------------------------------------------------------------------
#####============== (2a) Parse the GSEA database files, extract gene sets  =============================================================

GSEA_dbs <- list()

for(n in names(GSEA_parsers)) {
    cat(n,"\n")
    
    v <- readLines(paste(GSEA_base, n,sep="/"))
    if(!all(sapply(v,function(x)grepl("\t\t",x)))) {
        cat("*** Unexpected database format for", db, "-- skipping! ***\n")
        next
    }
    GSEA_dbs[[n]] <- list()
    
    tmp <- strsplit(v,'\t\t')
    tag <- sapply(tmp,function(x)gsub("\\s+","_",x[1]))              ## x[1], "as is", "_"-separated

    GSEA_dbs[[n]]$info <- sapply(tmp, GSEA_parsers[[n]])             ## x[1], parsed
    if(!is.null(dim(GSEA_dbs[[n]]$info))) {
        GSEA_dbs[[n]]$info <- t(GSEA_dbs[[n]]$info)
        
    } else {
        cn <- unique(names(GSEA_dbs[[n]]$info))
        cn <- ifelse(length(cn)==1,cn,"term")
        GSEA_dbs[[n]]$info <- matrix(GSEA_dbs[[n]]$info,ncol=1)
        colnames(GSEA_dbs[[n]]$info) <- cn
    }
    rownames(GSEA_dbs[[n]]$info) <- tag
    
    GSEA_dbs[[n]]$genes <- sapply(tmp,function(x)Reduce(union,sapply(x[2:length(x)],function(y)strsplit(y,"\\t"))))
    names(GSEA_dbs[[n]]$genes) <- tag                                ## x[2:(length(x)],each split into genes on "\t",
                                                                     ## then merged.
                                                                     ## The GO dbs have >1 "\t\t" in some rows
                                                                     ## (?? indicating hierarchical levels??).
                                                                     ## Simply ignoring the additional "\t\t"s
                                                                     ## yields correct set sizes as reported by Enrichr.
                                             ##-------------------->>   But: if possible, clarify the meaning of the format!!
 
}

                                                      
#####-----------------------------------------------------------------------------------------------------------------------------------
#####============== (2b) Run GSEA on the DGE outputs  ==================================================================================
#####============== (2c) Intersect GSEA and DGE gene lists; output results  ============================================================
#####-----------------------------------------------------------------------------------------------------------------------------------

relevant_info <- list(TFs     =c("term","method","genome","cell_type"),
                      GO      =c("term","id"),
                      pathways=c("term","id")) ## NA if missing!

maxq_GSEA <- 0.05
maxq_DGE <- 0.05 ##1e-3 ## 1e-5 ##0.05
rule <- "QVAL" ## use DGE genes with adjusted p-value <= 0.05

for(g in names(GSEA_groups)) {
    for(algo in use_algos) {
        cf <- names(DGE_output[[algo]])[1] ## ASSUMing that there is only one coefficient in the model -- otherwise need a loop over cf!

        GSEA2genes <- c()
        
        for(effect in c("up","down")) {
            cat("***", g, algo, effect,"\n")
            
            tmp <- DGEoutput2Enrichr(DGE_output[[algo]][[cf]],
                                     rule=rule,
                                     maxQ=maxq_DGE,
                                     qval_col =DGE_algorithm_data$columns[[algo]]["qval_col"],
                                     logFC_col=DGE_algorithm_data$columns[[algo]]["logFC_col"],
                                     gene_col =DGE_algorithm_data$columns[[algo]]["gene_col"],
                                     gene_map =DGE_algorithm_data$gene_map[[algo]],
                                     cmp=ifelse(effect=="up",">","<"),
                                     geneset_dbs=GSEA_groups[[g]])
            genes <- tmp$l ## the differentially expressed query genes
            for(db in GSEA_groups[[g]]) {
                d <- tmp$e[[db]]; if(is.null(d)) next
                d <- d[d[,"Adjusted.P.value"]<=maxq_GSEA,]; if(nrow(d)==0) next
                
                dbtag <- gsub("\\s+","_",d[,1])
                missing <- c()
                for(colmn in relevant_info[[g]]) {
                    if(colmn %in% colnames(GSEA_dbs[[db]]$info)) {
                        missing <- c(missing,setdiff(dbtag,rownames(GSEA_dbs[[db]]$info)))
                    }
                }
                if(length(missing) > 0) { ## CHECK why this can happen -- may have to do with duplicated entries! (over-truncated?)
                    ##cat("*** ", fct, subst, algo, db, ": some tags are not in GSEA_dbs[[db]]$info! ****\n")
                    good <- which(!(dbtag %in% missing))
                } else {
                    good <- 1:length(dbtag)
                }
                
                this_G2g <- data.frame(GSEAdb=db, GSEAset=dbtag[good],
                                       DGEalgorithm=algo,
                                       DGEeffect=effect,
                                       GSEAqval=d[good,"Adjusted.P.value"])
                for(colmn in relevant_info[[g]]) {
                    if(colmn %in% colnames(GSEA_dbs[[db]]$info)) {
                        this_G2g[,paste("GSEA",colmn,sep="")] <- GSEA_dbs[[db]]$info[dbtag[good],colmn]
                    } else {
                        this_G2g[,paste("GSEA",colmn,sep="")] <- NA
                    }
                }
                
                this_G2g[,"DGEgenes_in_GSEAset"] <- sapply(dbtag[good], 
                                                           function(tag) paste(genes[which(toupper(genes) %in% GSEA_dbs[[db]]$genes[[tag]])],
                                                                               collapse=","))
                
                GSEA2genes <- rbind(GSEA2genes,this_G2g)
                cat(nrow(GSEA2genes),"\n")
            }
        }
        if(!is.null(GSEA2genes)) {
            write.table(GSEA2genes,file=paste(g,"_", algo,"_",
                                       ifelse(do_counts,"counts","kallisto"),
                                       ".csv",sep=""),
                        sep="\t",quote=FALSE,
                        row.names=FALSE,col.names=TRUE)
        }
    }
}



