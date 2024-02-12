#' @title get_all_GO
#' @description extracts all the GO IDs that belongs to a specific GO branch
#' @param ont one of the 'BP', 'MF', or 'CC'
#' @return a vector that contains all GO terms of the specific 'ont'
#' @importFrom AnnotationDbi Ontology
#' @importFrom GO.db GOTERM
#' @export
get_all_GO <- function(ont = "BP") {
    goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
    goterms <- goterms[goterms == ont]
    return(names(goterms))
}

#' @title get_all_genes
#' @description extracts all genes that have GO annotation
#' @param ont one of the 'BP', 'MF', or 'CC'
#' @return a vector that contains all genes that have GO annotation
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @export
get_all_genes <- function(ont = "BP") {
    goterms <- get_all_GO(ont)
    go2gene <- AnnotationDbi::mapIds(org.Hs.eg.db, 
                    keys=goterms, column='ENTREZID',
                    keytype="GOALL", multiVals='list')
    all_genes <- unique(unlist(go2gene))
    all_genes <- all_genes[!is.na(all_genes)]
    return(all_genes)
}

#' @title simpleGO
#' @description Provides a simple method for GO ORA analysis
#' @param gene a gene vector (e.g., DE genes)
#' @param GOID selected GO ID to test with
#' @param background background genes (e.g., all genes that have GO annotation)
#' @return a data frame that contains a summary of the analysis
#' @examples 
#' data(geneList, package = "DOSE")
#' de <- names(geneList)[1:100]
#' background <- get_all_genes('BP')
#' simpleGO(de, "GO:0007059", background)
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom clusterProfiler bitr
#' @importFrom stats phyper
#' @export
simpleGO <- function(gene, GOID, background) {
    go2gene <- bitr(GOID, OrgDb = org.Hs.eg.db, 
                fromType = "GOALL", toType = "ENTREZID")
    # all genes that are annotated by the GOID
    all_genes_in_category <- unique(go2gene$ENTREZID)
    # query genes that are annotated by the GOID
    genes_in_category <- gene[gene %in% all_genes_in_category]

    k <- length(gene[gene %in% background])
    m <- length(all_genes_in_category)
    n <- length(background) - m
    q <- length(genes_in_category)
    pvalue <- phyper(q-1, m, n, k, lower.tail = FALSE)
    GeneRatio <- paste0(q, '/', k)
    BgRatio <- paste0(m, '/', m+n)
    res <- data.frame(GOID = GOID, 
                    GeneRatio = GeneRatio, 
                    BgRatio = BgRatio,
                    pvalue = pvalue)
    return(res)
}


