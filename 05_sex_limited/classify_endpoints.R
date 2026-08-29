#!/usr/bin/env Rscript
## Rebuild the four-class endpoint labelling described in sci7-28.tex Methods
## ("Sex-limited disease and the two fertility axes"). The original labelling was
## produced inline and never written to disk, so this reconstructs it from the
## stated rule and checks it against the published counts
## (69 gynaecological / 13 male-limited / 98 fertility-bound / 1,606 sex-shared
## among the 1,786 endpoints having both rho_GE and r_g).
suppressMessages(library(data.table))

MAN <- "${BASEDIR}/finngen_r13/R13_manifest.tsv"
RG  <- "${BASEDIR}/finngen_r13/fg_rg_results.tsv"
RHO <- Sys.getenv("RHO_FILE", "${BASEDIR}/rhoGE_50tissue.tsv")
OUT <- Sys.getenv("CLASS_OUT", "${BASEDIR}/endpoint_class.tsv")

man <- fread(MAN, select = c("phenocode","phenotype","category"))
setnames(man, c("phenocode","pheno","cat"))

## ---- fertility-bound: consequences of having reproduced, circular with a
## number-of-children trait. ICD-10 chapters XV + XVI, plus the three named
## exceptions from the Methods.
is_ch15 <- grepl("^XV Pregnancy", man$cat)
## Chapter XVI is the NEWBORN's illness, not the mother's. Its circularity runs
## through the offspring rather than through maternal reproductive history, and
## empirically it behaves the opposite way (delta rho +0.005, indistinguishable
## from sex-shared, P=0.96) to the maternal endpoints (-0.011, P=1e-8). Pooling
## the two only dilutes the control, so perinatal endpoints go to sex-shared.
is_ch16 <- rep(FALSE, nrow(man))
## Results text (sci7-28.tex:148) states the class as "pregnancy, childbirth,
## perinatal, infertility and ovulation-induction endpoints" -- infertility and
## puerperium belong here, not with the gynaecological diseases.
named_fb <- grepl("gestational diabetes", man$pheno, ignore.case = TRUE) |
            grepl("GEST_DIABETES|GESTDM", man$phenocode) |
            (grepl("endometriosis", man$pheno, ignore.case = TRUE) &
             grepl("infertil", man$pheno, ignore.case = TRUE)) |
            grepl("ovulation induction|induction of ovulation", man$pheno, ignore.case = TRUE) |
            grepl("infertil|sterilit", man$pheno, ignore.case = TRUE) |
            grepl("puerper|maternal care|childbirth|obstetric|labour and delivery",
                  man$pheno, ignore.case = TRUE)
man[, fertbound := is_ch15 | is_ch16 | named_fb]
## The keyword rule can pull a chapter XVI endpoint back in when its description
## mentions labour or delivery -- but it is still the newborn who is ill.
man[grepl("^XVI Certain conditions originating in the perinatal", cat), fertbound := FALSE]

## ---- female reproductive tract. Word-boundary matching on the phenotype
## description, not the phenocode: "cervical" in M13_CERVICALGIA is a neck, and
## "salpingitis" in H8_EUSTSALP is an ear.
## Stems are deliberately open-ended (and the leading \\b dropped where a prefix
## is normal): "Amenorrhoea", "Menorrhagia", "Postmenopausal", "Salpingitis",
## "Vulvovaginal" and "oophoritis" all failed a stricter word-boundary form.
fem_rx <- paste0("(\\buterus\\b|\\buterine\\b|\\buteri\\b|\\bovar[a-z]*|cervix|",
                 "cervical (cancer|intraepithelial|dysplasia|smear)|vulv[a-z]*|",
                 "vagin[a-z]*|endometri[a-z]*|myometri[a-z]*|fallopian|oviduct|",
                 "female pelvic|female genital|menstrua[a-z]*|menopaus[a-z]*|",
                 "menorrh[a-z]*|premenstrual|leiomyoma|myoma|salping[a-z]*|",
                 "oophor[a-z]*|\\badnexa[a-z]*|parametri[a-z]*|",
                 "habitual abort[a-z]*|hirsutism in women)")
mal_rx <- paste0("\\b(prostate|prostatic|testis|testes|testicular|scrotum|scrotal|",
                 "penis|penile|epididym[a-z]*|spermato[a-z]*|seminal vesicle|",
                 "erectile|impoten[a-z]*|phimosis|paraphimosis|balanitis|priapism|",
                 "hydrocele|varicocele|prepuce|foreskin|cryptorchid[a-z]*|",
                 "undescended test[a-z]*|hypospadias|epispadias|gynaecomastia|gynecomastia|",
                 "male genital|male breast)\\b")

man[, fem_hit := grepl(fem_rx, pheno, ignore.case = TRUE)]
man[, mal_hit := grepl(mal_rx, pheno, ignore.case = TRUE)]

## Methods: endpoints defined explicitly in both sexes go to sex-shared.
man[grepl("women and men|men and women|in both sexes", pheno, ignore.case = TRUE),
    `:=`(fem_hit = FALSE, mal_hit = FALSE)]
## Breast is hormone-driven and >99% of cases are female, so breast endpoints are
## grouped with the female tract even though FinnGen's definition admits male cases.
## "Male breast" endpoints stay male; congenital and accessory breast stay shared.
man[grepl("breast|mammary", pheno, ignore.case = TRUE) &
    !grepl("male breast|congenital|accessory", pheno, ignore.case = TRUE),
    `:=`(fem_hit = TRUE, mal_hit = FALSE)]
## Eustachian salpingitis is an ear condition despite the "salping" stem.
man[grepl("eustachian", pheno, ignore.case = TRUE), fem_hit := FALSE]
## Cervical spine / neck endpoints live in the musculoskeletal chapter.
man[grepl("^XIII Diseases of the musculoskeletal", cat), fem_hit := FALSE]
## "adnexa" is the EYE adnexa as often as the uterine one. Guard on the organ word
## rather than the chapter: C3_EYE_ADNEXA_WIDE is a cancer in chapter II and slips
## past a chapter-VII-only rule.
man[grepl("\\beye\\b|ocular|conjunctiv|orbit", pheno, ignore.case = TRUE), fem_hit := FALSE]

man[, class := fifelse(fertbound, "fertility_bound",
              fifelse(fem_hit & !mal_hit, "gynaecological",
              fifelse(mal_hit & !fem_hit, "male_limited", "sex_shared")))]

## ---- restrict to endpoints with both estimates, as the paper does
rg  <- fread(RG)[!is.na(rg_father) & !is.na(rg_mother), .(phenocode)]
rho <- fread(RHO)[!is.na(rho_father) & !is.na(rho_mother), .(phenocode)]
keep <- intersect(rg$phenocode, rho$phenocode)
cl <- man[phenocode %in% keep, .(phenocode, pheno, cat, class)]

cat(sprintf("endpoints with both rho_GE and r_g: %d  (paper: 1,786)\n", nrow(cl)))
print(cl[, .N, by = class][order(-N)])
cat("paper: gynaecological 69 | male_limited 13 | fertility_bound 98 | sex_shared 1,606\n\n")

fwrite(cl, OUT, sep = "\t")
cat("written", OUT, "\n")
cat("\n--- gynaecological ---\n"); print(cl[class == "gynaecological", .(phenocode, pheno)])
cat("\n--- male_limited ---\n");   print(cl[class == "male_limited",   .(phenocode, pheno)])
