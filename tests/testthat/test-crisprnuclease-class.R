context('Testing CrisprNuclease class')

test_that('PAM', {
    pams_cas9 <- c("AGG", "CGG", "GGG", "TGG")
    pams_cas9_extended <- c("AGG", "CGG", "GGG", "TGG",
                            "AAG", "CAG", "GAG", "TAG",
                            "AGA", "CGA", "GGA", "TGA")
    names(pams_cas9) <- rep("NGG", 4)
    names(pams_cas9_extended) <- c(rep("NGG", 4),
                                   rep("NAG", 4),
                                   rep("NGA", 4))
    expect_true(pamLength(SpCas9)==3)
    expect_true(pamLength(AsCas12a)==4)
    expect_equal(sort(pams(SpCas9, as.character=TRUE)), sort(pams_cas9))
    expect_equal(sort(pams(SpCas9, primary=FALSE, as.character=TRUE)), sort(pams_cas9_extended))
    expect_equal(pamIndices(SpCas9), 21:23)
    expect_equal(pamIndices(AsCas12a), 1:4)
    expect_true(pamSide(SpCas9)=="3prime")
    expect_true(pamSide(AsCas12a)=="5prime")
})


test_that('Spacer', {
    expect_true(spacerLength(SpCas9)==20)
    expect_true(spacerLength(AsCas12a)==23)
    expect_equal(spacerIndices(SpCas9), 1:20)
    expect_equal(spacerIndices(AsCas12a), 5:27)
    expect_true(spacerSide(SpCas9)=="5prime")
    expect_true(spacerSide(AsCas12a)=="3prime")
})

test_that('Protospacer', {
    expect_true(protospacerLength(SpCas9)==23)
    expect_true(protospacerLength(AsCas12a)==27)
})


