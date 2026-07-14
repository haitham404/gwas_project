# scripts/phase10_heritability/heritability.R

library(sommer)


load_grm <- function(grm_file, grm_id_file) {

  grm <- as.matrix(
    read.table(
      grm_file,
      header = FALSE
    )
  )

  id_data <- read.table(
    grm_id_file,
    header = FALSE
  )

  grm_ids <- as.character(id_data$V2)

  if (nrow(grm) != length(grm_ids)) {
    stop("GRM size does not match number of sample IDs")
  }

  rownames(grm) <- grm_ids
  colnames(grm) <- grm_ids

  return(
    list(
      grm = grm,
      ids = grm_ids
    )
  )
}


estimate_h2 <- function(
  phenotype,
  sample_ids,
  grm,
  grm_ids
) {

  pheno <- data.frame(
    IID = as.character(sample_ids),
    phenotype = as.numeric(phenotype)
  )

  # Remove missing phenotype values
  pheno <- pheno[
    !is.na(pheno$phenotype),
  ]

  # Match phenotype samples with GRM samples
  common_ids <- intersect(
    grm_ids,
    pheno$IID
  )

  if (length(common_ids) < 2) {
    stop("Not enough matched samples")
  }

  pheno <- pheno[
    match(common_ids, pheno$IID),
  ]

  grm_index <- match(
    common_ids,
    grm_ids
  )

  G <- grm[
    grm_index,
    grm_index,
    drop = FALSE
  ]

  rownames(G) <- common_ids
  colnames(G) <- common_ids

  pheno$IID <- factor(
    pheno$IID,
    levels = common_ids
  )

  fit <- mmes(
    phenotype ~ 1,

    random = ~ vsm(
      ism(IID),
      Gu = G
    ),

    rcov = ~ units,

    data = pheno,

    verbose = FALSE
  )

  h2 <- vpredict(
    fit,
    h2 ~ V1 / (V1 + V2)
  )

  return(
    list(
      h2 = h2,
      sample_size = nrow(pheno),
      model = fit
    )
  )
}