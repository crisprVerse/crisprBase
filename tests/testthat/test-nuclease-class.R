context('Testing Nuclease class')

data(SpCas9, package="crisprBase")
data(AsCas12a, package="crisprBase")
data(enAsCas12a, package="crisprBase")

test_that('Motif representations', {
    expect_error(crisprBase:::.checkRebaseMotif("(9/10)AHT^HCT(9/10)"))
    expect_error(crisprBase:::.checkRebaseMotif("AHT^^HCT"))
    expect_error(crisprBase:::.checkRebaseMotif("(9/10)AHT^HCT"))
    expect_error(crisprBase:::.checkRebaseMotif("AHTU"))
    expect_null(crisprBase:::.checkRebaseMotif("(9/10)AHTHCT"))
    expect_null(crisprBase:::.checkRebaseMotif("AHTHCT(9/10)"))
    expect_null(crisprBase:::.checkRebaseMotif("AHTHCT"))
    expect_null(crisprBase:::.checkRebaseMotif("AHT^HCT"))
})


test_that('Motif sequences', {
    cas9_pams <- c("AGG", "CGG", "GGG", "TGG")
    cas12a_pams <- c("TTTA", "TTTC", "TTTG")
    names(cas9_pams)   <- rep("NGG",4)
    names(cas12a_pams) <- rep("TTTV",3)
    expect_equal(motifs(SpCas9, primary=TRUE, expand=TRUE, as.character=TRUE), cas9_pams)
    expect_equal(motifs(AsCas12a, primary=TRUE, expand=TRUE, as.character=TRUE), cas12a_pams)
})



 
# # Outside
# (3/3) A^AAA[NGG]CGG -3
# (2/2) AA^AA[NGG]CGG -2
# (1/1) AAA^A[NGG]CGG -1
# # Within
# AAAA^[NGG]CGG 0
# AAAA[N^GG]CGG 1
# AAAA[NG^G]CGG 2
# AAAA[NGG]^CGG 3
# # Outside
# (1/1) AAAA[NGG]C^GG 4 (pam_length+1)
# (2/2) AAAA[NGG]CG^G 5 (pam_length+2)
# (3/3) AAAA[NGG]CGG^ 6 (pam_length+3)





test_that('Extracting motif cut sites (upstream)', {

    motifs <- c("(3/3)AGG",
                "(2/2)AGG",
                "(1/1)AGG",
                "(3/1)AGG",
                "(3/2)AGG",
                "(3/5)AGG")
    nucleases <- lapply(motifs,
                        function(motif){
        Nuclease(nucleaseName="", metadata="", motifs=motif)
    })
    expectedCutSitesFwdStrand <- c(-3,-2,-1,-3,-3,-3)
    expectedCutSitesRevStrand <- c(-3,-2,-1,-1,-2,-5)

    expect_equal(expectedCutSitesFwdStrand,
                 vapply(nucleases,
                        cutSites,
                        combine=FALSE,
                        middle=FALSE,
                        FUN.VALUE=0))
    expect_equal(expectedCutSitesRevStrand,
                 vapply(nucleases,
                        cutSites,
                        combine=FALSE,
                        strand="-",
                        middle=FALSE,
                        FUN.VALUE=0))
})



test_that('Extracting motif cut sites (downstream)', {

    motifs <- c("AGG(1/1)",
                "AGG(2/2)",
                "AGG(3/3)",
                "AGG(3/2)",
                "AGG(3/1)",
                "AGG(3/5)")
    nucleases <- lapply(motifs,
                        function(motif){
        Nuclease(nucleaseName="", metadata="", motifs=motif)
    })
    expectedCutSitesFwdStrand <- c(4,5,6,6,6,6)
    expectedCutSitesRevStrand <- c(4,5,6,5,4,8)

    expect_equal(expectedCutSitesFwdStrand,
                 vapply(nucleases,
                        cutSites,
                        combine=FALSE,
                        middle=FALSE,
                        FUN.VALUE=0))
    expect_equal(expectedCutSitesRevStrand,
                 vapply(nucleases,
                        cutSites,
                        combine=FALSE,
                        middle=FALSE,
                        strand="-",
                        FUN.VALUE=0))
})



test_that('Extracting motif cut sites (within)', {

    motifs <- c("^AGG",
                "A^GG",
                "AA^G",
                "AAG^")
    nucleases <- lapply(motifs,
                        function(motif){
        Nuclease(nucleaseName="", metadata="", motifs=motif)
    })
    expectedCutSitesFwdStrand <- c(0,1,2,3)
    expectedCutSitesRevStrand <- c(3,2,1,0)

    expect_equal(expectedCutSitesFwdStrand,
                 vapply(nucleases,
                        cutSites,
                        combine=FALSE,
                        middle=FALSE,
                        FUN.VALUE=0))
    expect_equal(expectedCutSitesRevStrand,
                 vapply(nucleases,
                        cutSites,
                        combine=FALSE,
                        middle=FALSE,
                        strand="-",
                        FUN.VALUE=0))
})





