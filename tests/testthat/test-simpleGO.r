context('GO-ORA')

test_that('simpleGO', {
    data(geneList, package = "DOSE")
    de <- names(geneList)[1:100]
    background <- get_all_genes('BP')
    GOID <- "GO:0007059"
    res <- simpleGO(de, GOID, background)
    res2 <- clusterProfiler::enrichGO(de, OrgDb = 'org.Hs.eg.db', ont = "BP")

    expect_equal(
        res$pvalue, res2[GOID, 'pvalue']
    )
})
