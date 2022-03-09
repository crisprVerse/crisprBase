context('Testing arithmetic functions')

data(SpCas9, package="crisprBase")
data(AsCas12a, package="crisprBase")
data(enAsCas12a, package="crisprBase")

targets_cas9   <- c("ATGGTGGTAGTCTCTCGATACGG",
                         "ATGGGTTTAGTCTAAACGATAGG")
targets_cas12a <- c("TTTAATGGTGGTAGTCTCTCGATACGG",
                         "TTTCATGGGTTTAGTCTAAACGATAGG")
spacers_cas9   <- c("ATGGTGGTAGTCTCTCGATA",
                    "ATGGGTTTAGTCTAAACGAT")
spacers_cas12a <- c("ATGGTGGTAGTCTCTCGATACGG",
                    "ATGGGTTTAGTCTAAACGATAGG")
pams_cas9   <- c("CGG", "AGG")
pams_cas12a <- c("TTTA", "TTTC")

# Grangesstuff:
chr <- c("chr10", "chr11")                    
pam_site <- c(1998467, 220000)
strand <- c("+", "-")
gr <- GenomicRanges::GRanges(chr,
                             IRanges::IRanges(start=pam_site,
                                              width=1),
                             strand=strand)

gr_pam_cas9 <- GenomicRanges::GRanges(chr,
                       strand=strand,
                       IRanges::IRanges(start=c(1998467, 219998),
                               end=c(1998469,220000)))
gr_pam_cas12a <- GenomicRanges::GRanges(chr,
                       strand=strand,
                       IRanges::IRanges(start=c(1998467, 219997),
                               end=c(1998470,220000)))

gr_spacer_cas9 <- GenomicRanges::GRanges(chr,
                       strand=strand,
                       IRanges::IRanges(start=c(1998447, 220001),
                               end=c(1998466,220020)))
gr_spacer_cas12a <- GenomicRanges::GRanges(chr,
                       strand=strand,
                       IRanges::IRanges(start=c(1998471, 219974),
                               end=c(1998493,219996)))
gr_target_cas9 <- GenomicRanges::GRanges(chr,
                       strand=strand,
                       IRanges::IRanges(start=c(1998447, 219998),
                               end=c(1998469,220020)))
gr_target_cas12a <- GenomicRanges::GRanges(chr,
                       strand=strand,
                       IRanges::IRanges(start=c(1998467, 219974),
                               end=c(1998493,220000)))

gr_spacer_cas9_long <- GenomicRanges::GRanges(chr,
                       strand=strand,
                       IRanges::IRanges(start=c(1998367, 220001),
                               end=c(1998466,220100)))
gr_spacer_cas12a_long <- GenomicRanges::GRanges(chr,
                       strand=strand,
                       IRanges::IRanges(start=c(1998471, 219897),
                               end=c(1998570,219996)))

gr_target_cas9_long <- GenomicRanges::GRanges(chr,
                       strand=strand,
                       IRanges::IRanges(start=c(1998367, 219998),
                               end=c(1998469,220100)))
gr_target_cas12a_long <- GenomicRanges::GRanges(chr,
                       strand=strand,
                       IRanges::IRanges(start=c(1998467, 219897),
                               end=c(1998570,220000)))



test_that('Extraction functions', {
    expect_equal(extractSpacerFromTarget(targets_cas9, SpCas9), spacers_cas9)
    expect_equal(extractSpacerFromTarget(targets_cas12a, AsCas12a), spacers_cas12a)
    expect_equal(extractPamFromTarget(targets_cas9, SpCas9), pams_cas9)
    expect_equal(extractPamFromTarget(targets_cas12a, AsCas12a), pams_cas12a)
})

test_that('Genomic ranges functions', {
    # Using GRanges as inputs:
    expect_equal(getPamRanges(gr, nuclease=SpCas9), gr_pam_cas9)
    expect_equal(getPamRanges(gr, nuclease=AsCas12a), gr_pam_cas12a)
    expect_equal(getProtospacerRanges(gr, nuclease=SpCas9), gr_spacer_cas9)
    expect_equal(getProtospacerRanges(gr, nuclease=AsCas12a), gr_spacer_cas12a)
    expect_equal(getTargetRanges(gr, nuclease=SpCas9), gr_target_cas9)
    expect_equal(getTargetRanges(gr, nuclease=AsCas12a), gr_target_cas12a)

    # Using coordinates as inputs:
    expect_equal(getPamRanges(seqnames=chr, pam_site=pam_site, strand=strand, nuclease=SpCas9), gr_pam_cas9)
    expect_equal(getPamRanges(seqnames=chr, pam_site=pam_site, strand=strand, nuclease=AsCas12a), gr_pam_cas12a)
    expect_equal(getProtospacerRanges(seqnames=chr, pam_site=pam_site, strand=strand, nuclease=SpCas9), gr_spacer_cas9)
    expect_equal(getProtospacerRanges(seqnames=chr, pam_site=pam_site, strand=strand, nuclease=AsCas12a), gr_spacer_cas12a)
    expect_equal(getTargetRanges(seqnames=chr, pam_site=pam_site, strand=strand, nuclease=SpCas9), gr_target_cas9)
    expect_equal(getTargetRanges(seqnames=chr, pam_site=pam_site, strand=strand, nuclease=AsCas12a), gr_target_cas12a)

    # Different spacer length
    expect_equal(getProtospacerRanges(gr, nuclease=SpCas9, spacer_len=100), gr_spacer_cas9_long)
    expect_equal(getProtospacerRanges(gr, nuclease=AsCas12a, spacer_len=100), gr_spacer_cas12a_long)
    expect_equal(getTargetRanges(gr, nuclease=SpCas9, spacer_len=100), gr_target_cas9_long)
    expect_equal(getTargetRanges(gr, nuclease=AsCas12a, spacer_len=100), gr_target_cas12a_long)
})




