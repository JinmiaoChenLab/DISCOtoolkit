#' Cluster-level cell type prediction
#'
#' @param rna A data.frame or matrix with genes listed as row names
#' @param ref.data Reference gene expression. If not specified, reference gene
#' expression data will be downloaded from the DISCO database.
#' @param ref.deg DEG list.
#' @param atlas Specify the atlas. If not specified, CELLiD will use all data.
#' @param n.predict Number of predicted cell type. Maximum value is 3
#' @param ref.path Specify the file path to save the downloaded reference data.
#' @param ncores Number of core used for CELLiD
#' @export
#' @importFrom pbmcapply pbmclapply
#' @importFrom stringr str_match
CELLiDCluster <- function(rna, ref.data = NULL, ref.deg = NULL, atlas = NULL, n.predict = 1, verbose = T, ref.path = NULL, ncores = 10) {

  # check input data
  if (!(is.data.frame(rna) | is.matrix(rna))) {
    stop("The rna must be a data frame or a matrix")
  }

  if (n.predict > 3) {
    message("Any value of n.predict that exceeds 3 will be automatically adjusted to 3.")
    n.predict = 3
  }

  # Download reference data if missing
  if (is.null(ref.data) | is.null(ref.deg)) {
    op = options(timeout=600)
    on.exit(options(op), add = TRUE)
    if (is.null(ref.path)) {
      dir.create("DISCOtmp")
      ref.path = "DISCOtmp"
    }
    if (is.null(ref.data)) {
      download.file("http://www.immunesinglecell.org/toolkitapi/getRef",
                    paste0(ref.path, "/ref_data.rds"), method = "curl")
      ref.data = readRDS(paste0(ref.path, "/ref_data.rds"))
    }
    if (is.null(ref.deg)) {
      download.file("http://www.immunesinglecell.org/toolkitapi/getRefDeg",
                    paste0(ref.path, "/ref_deg.rds"), method = "curl")
      ref.deg = readRDS(paste0(ref.path, "/ref_deg.rds"))
    }
  }

  if (!is.null(atlas)) {
    ref.data = ref.data[,which(str_match(colnames(ref.data), "--(.*)")[,2] == atlas)]
  }
  genes = intersect(rownames(rna), rownames(ref.data))

  if (length(genes) <= 2000) {
    stop("Less than 2000 genes are shared between the input data and the reference dataset!")
  }

  if (length(genes) <= 5000) {
    warning("The input data and reference dataset have a limited number of overlapping genes, which may potentially impact the accuracy of the CELLiD.")
  }

  rna = rna[genes, ,drop=F]
  ref.data = ref.data[genes, ]

  # Initial round of prediction
  predicted.cell = pbmclapply(
    1:ncol(rna), function(j) {
      predicted = as.numeric(
        apply(ref.data, 2, function(i) {
          cor(as.numeric(i)[which(rna[,j] > 0)], as.numeric(rna[which(rna[,j] > 0), j]), method = "spearman", use="complete.obs")
        })
      )
      return(predicted)
    }, mc.cores = ncores
  )

  predicted.cell = do.call(cbind, predicted.cell)
  rownames(predicted.cell) = colnames(ref.data)

  ct = lapply(1:ncol(predicted.cell), function(i) {
    x = predicted.cell[,i]
    return(as.numeric(which(rank(-x) <=5)))
  })

  # Second round of prediction
  predicted.cell = pbmclapply(
    1:ncol(rna), function(i) {
      ref = ref.data[,ct[[i]]]
      g = unique(ref.deg$gene[which(ref.deg$group %in% colnames(ref))])
      g = intersect(rownames(ref), g)
      ref = ref[g,]
      input = rna[g,i]
      predict = apply(ref, 2, function(i) {
        cor(as.numeric(i), input, method = "spearman", use="complete.obs")
      })
      return(c(str_match(names(sort(predict, decreasing = T)[1:n.predict]), "^(.*)--")[,2],
               str_match(names(sort(predict, decreasing = T)[1:n.predict]), "--(.*)$")[,2],
               as.numeric(sort(predict, decreasing = T)[1:n.predict])))
    }, mc.cores = ncores
  )

  predicted.cell = do.call(rbind, predicted.cell)
  predicted.cell = data.frame(predicted.cell)
  colnames(predicted.cell)[1:(ncol(predicted.cell)/3)] = paste0("predict_cell_type_", 1:n.predict)
  colnames(predicted.cell)[(ncol(predicted.cell)/3 + 1):(ncol(predicted.cell)/3*2)] = paste0("source_atlas_", 1:n.predict)
  colnames(predicted.cell)[(ncol(predicted.cell)/3*2 + 1):ncol(predicted.cell)] = paste0("score_", 1:n.predict)

  for (i in (ncol(predicted.cell)/3*2 + 1):ncol(predicted.cell)) {
    predicted.cell[,i] = round(as.numeric(predicted.cell[,i]), 3)
  }

  rownames(predicted.cell) = colnames(rna)



  return(predicted.cell)
}


#' Geneset enrichment analysis using DISCO data
#'
#' @param input A one or two columns data.frame. The first column is gene and the second column is log2FC (optional)
#' @param reference Reference geneset If not specified, reference geneset will be downloaded from the DISCO database.
#' @param ref.path Specify the file path to save the downloaded reference geneset data.
#' @param ncores Number of core used for CELLiD
#' @export
CELLiDEnrichment <- function(input, reference = NULL, ref.path = NULL, ncores = 10){
  # check input data
  if (!(is.data.frame(input))) {
    stop("The input must be a data frame")
  }
  if (ncol(input) > 2) {
    stop("The input must be one or two columns")
  }
  if (ncol(input) == 2) {
    colnames(input) = c("gene", "fc")
  }else {
    colnames(input) = c("gene")
  }


  if (is.null(reference)) {
    #TODO reference url
    op = options(timeout=600)
    on.exit(options(op), add = TRUE)
    if (is.null(ref.path)) {
      dir.create("DISCOtmp")
      ref.path = "DISCOtmp"
    }
    download.file("http://www.immunesinglecell.org/toolkitapi/getGeneSet",
                  paste0(ref.path, "/ref_geneset.rds"), method = "curl")
    reference = readRDS(paste0(ref.path, "/ref_geneset.rds"))
  }
  reference$name = paste0(reference$name, " in ", reference$atlas)

  if(ncol(input) == 2) {
    input$fc = 2^input$fc
    rownames(input) = input$gene
    input = input[intersect(reference$gene, input$gene),]

    res = pbmclapply(
      unique(reference$name), function(i) {
        atlas = str_match(i, " in (.*?$)")[,2]
        reference.filter = reference[which(reference$name == i),]
        reference.full = reference[which(reference$atlas == atlas),]
        rownames(reference.filter) = reference.filter$gene
        input.filter = input[intersect(reference.full$gene, input$gene),]

        if (length(input.filter) == 0) {
          return(NULL)
        }

        a = sum(reference.filter[intersect(reference.filter$gene, input.filter$gene),1] * input.filter[intersect(reference.filter$gene, input.filter$gene),2]) + 1
        b = sum(input.filter[setdiff(input.filter$gene, reference.filter$gene),2]) + 1
        c = sum(reference.filter[setdiff(reference.filter$gene,input.filter$gene),1]) + 1
        d = length(unique(reference.full$gene)) - length(unique(c(reference.filter$gene,input.filter$gene)))
        r = fisher.test(matrix(c(a,b,c,d), nrow = 2))
        if (r$p.value < 0.01) {
          return(c(r$p.value, as.numeric(r$estimate), i,
                   paste0(intersect(reference.filter$gene, input.filter$gene), collapse  = ","),
                   length(unique(reference.full$gene)), length(intersect(reference.filter$gene,input.filter$gene)), length(reference.filter$gene)))
        } else {
          return(NULL)
        }
      }, mc.cores = ncores
    )
  }else {
    input = toupper(input$gene)
    input = intersect(reference$gene, input)

    res = pbmclapply(
      unique(reference$name), function(i) {
        atlas = str_match(i, " in (.*?$)")[,2]
        reference.filter = reference[which(reference$name == i),]
        reference.full = reference[which(reference$atlas == atlas),]

        input.filter = intersect(reference.full$gene, input)

        if (length(input.filter) == 0) {
          return(NULL)
        }

        a = length(intersect(reference.filter$gene, input.filter)) + 1
        b = length(setdiff(input.filter, reference.filter$gene)) + 1
        c = length(setdiff(reference.filter$gene, input.filter)) + 1
        d = length(unique(reference.full$gene)) - a - b -c
        r = fisher.test(matrix(c(a,b,c,d), nrow = 2))
        if (r$p.value < 0.01) {
          return(c(r$p.value, as.numeric(r$estimate), i,
                   paste0(intersect(reference.filter$gene, input.filter), collapse  = ","),
                   length(unique(reference.full$gene)), length(intersect(reference.filter$gene,input.filter)), length(reference.filter$gene)))
        } else {
          return(NULL)
        }
      }, mc.cores = ncores
    )
  }

  res = do.call(rbind, res)
  res = data.frame(res)
  res$X1 = as.numeric(res$X1)
  res$X2 = as.numeric(res$X2)
  res = res[order(res$X1,-res$X2),]
  colnames(res) = c("pval", "or", "name", "gene", "background", "overlap", "geneset")
  res$or = round(as.numeric(res$or), 2)
  res$pval = signif(as.numeric(res$pval), digits=2)
  res = res[1:min(50, nrow(res)),,drop=F]
  return(res)
}


