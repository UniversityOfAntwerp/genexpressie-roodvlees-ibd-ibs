# Functions for differential expression analysis:

# Run differential expression tests for all genes
diffex.test.all <- function(form, data, meta, var=NULL) {
  
  #' Test differential expression for each gene
  #'
  #' Fits a linear mixed model for each gene in a normalized expression matrix.
  #'
  #' @param form Model formula.
  #' @param data Data frame or matrix with genes as rows and samples as columns.
  #' @param meta Metadata data frame containing sample information.
  #' @param var Optional coefficient name or index to extract from the model.
  #'
  #' @return Data frame with model results and FDR-adjusted p-values.
  
  updated.form <- update.formula(form, gene ~ .)
  meta.gene <- meta
  
  pb <- txtProgressBar(
    min = 0,
    max = nrow(data),
    initial = 0,
    style = 3
  )
  
  R <- Reduce(
    rbind,
    apply(data, 1, function(expr) {
      
      tryCatch({
        
        meta.gene$gene <- expr
        fit <<- lmer(updated.form, data = meta.gene)
        
        res <- if (is.null(var)) {
          as.data.frame(summary(fit)$coefficients)[2, ]
        } else {
          as.data.frame(summary(fit)$coefficients)[var, ]
        }
        
        res$singular <- isSingular(fit)
        return(res)
        
      }, error = function(cond) {
        
        missing <- as.data.frame(list(NA, NA, NA, NA, NA, NA))
        colnames(missing) <- c("Estimate", "Std. Error", "df", "t value", "Pr(>|t|)", "singular")
        return(missing)
        
      }, finally = {
        setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
      })
    })
  )
  
  rownames(R) <- rownames(data)
  R$qvalue <- p.adjust(R$`Pr(>|t|)`, method = "fdr")
  
  return(R)
}
# Create a volcano plot from differential expression results
volcano <- function(diffex.res, q.thresh=0.05, fc.thresh=1, only_sig=T) {
  
  p.thresh <- max(diffex.res[diffex.res$qvalue < q.thresh,]$`Pr(>|t|)`)
  
  significant <- diffex.res[
    (diffex.res$qvalue < q.thresh) &
      (abs(diffex.res$Estimate) >= fc.thresh),
  ]
  
  pv <- -log2(significant$`Pr(>|t|)`)
  
  insignificant <- diffex.res[
    !((diffex.res$qvalue < q.thresh) &
        (abs(diffex.res$Estimate) >= fc.thresh)),
  ]
  
  pv <- pv[is.finite(pv)]
  ylim <- c(0, max(pv)) * 1.1
  
  xlim <- if (only_sig & (nrow(significant) > 0)) {
    c(min(significant$Estimate), max(significant$Estimate)) * 1.1
  } else {
    c(min(diffex.res$Estimate), max(diffex.res$Estimate)) * 1.1
  }
  
  plot(
    insignificant$Estimate,
    -log2(insignificant$`Pr(>|t|)`),
    xlim = xlim,
    ylim = ylim,
    xlab = "Log Fold Change",
    ylab = "-log p-value"
  )
  
  points(
    significant$Estimate,
    -log2(significant$`Pr(>|t|)`),
    col = sign(significant$Estimate) + 3
  )
  
  lines(xlim, -log2(c(p.thresh, p.thresh)), col = "black")
  lines(-c(fc.thresh, fc.thresh), ylim, col = "black")
  lines(c(fc.thresh, fc.thresh), ylim, col = "black")
}
# Test annotation-term enrichment among significant genes
diffex.enrich <- function(diffex.res, annotations, direction="all") {
  
  # Use only genes present in both the results and the annotation table
  universe <- intersect(rownames(diffex.res), annotations$gene)
  diffex.rel <- diffex.res[universe,]
  
  # Select significant genes in the requested direction
  diffex.sig <- rownames(diffex.rel[
    (diffex.rel$qvalue < 0.05) &
      if (direction == "up") {
        diffex.rel$Estimate > 0
      } else if (direction == "down") {
        diffex.rel$Estimate < 0
      } else {
        TRUE
      },
  ])
  
  # Keep only annotation terms linked to significant genes
  annot.rel.terms <- unique(annotations[annotations$gene %in% diffex.sig,]$term)
  
  annot.rel <- annotations[
    (annotations$term %in% annot.rel.terms) &
      (annotations$gene %in% universe),
  ]
  
  # Keep term descriptions if they are available
  annot.description <- if ("description" %in% colnames(annot.rel)) {
    ad <- annot.rel[!duplicated(annot.rel$term), ]
    ad[ad$term %in% annot.rel.terms, c("term", "description")]
  } else {
    NULL
  }
  
  # Test each annotation term
  terms <- unique(annot.rel.terms)
  pb <- txtProgressBar(min=1, max=length(terms), initial=0, style=3)
  
  R <- Reduce(
    rbind,
    lapply(terms, function(t) {
      
      setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
      
      annot.term <- annot.rel[annot.rel$term == t,]$gene
      
      # Build the 2x2 contingency table
      a <- length(intersect(diffex.sig, annot.term))
      b <- length(setdiff(annot.term, diffex.sig))
      c <- length(setdiff(diffex.sig, annot.term))
      d <- length(setdiff(universe, union(diffex.sig, annot.term)))
      
      contingency <- matrix(c(a, b, c, d), nrow=2)
      
      # Use Fisher's exact test for small counts, otherwise use chi-square test
      p.value <- if (min(contingency) < 5) {
        fisher.test(contingency)$p.value
      } else {
        chisq.test(contingency)$p.value
      }
      
      odds.ratio <- (a/c) / (b/d)
      
      as.data.frame(
        list(t, a, b, c, d, odds.ratio, p.value),
        col.names = c("term", "a", "b", "c", "d", "odds.ratio", "p.value")
      )
    })
  )
  
  R$qvalue <- p.adjust(R$p.value, "fdr")
  
  if (!is.null(annot.description)) {
    R <- merge(R, annot.description)
  }
  
  R
}
# Run a GSEA-like enrichment test using ranked differential expression results
diffex.gsea <- function(diffex.res, annotations, direction="all") {
  
  # Use only genes present in both the results and the annotation table
  universe <- intersect(rownames(diffex.res), annotations$gene)
  diffex.rel <- diffex.res[universe,]
  
  # Rank genes so the strongest results are placed at the top
  diffex.rel$qvalue <- 1 - diffex.rel$qvalue
  
  # Adjust ranking based on the requested expression direction
  diffex.rel$qvalue <- if (direction == "up") {
    diffex.rel$qvalue * sign(diffex.rel$Estimate)
  } else if (direction == "down") {
    diffex.rel$qvalue * -1 * sign(diffex.rel$Estimate)
  } else {
    diffex.rel$qvalue
  }
  
  # Keep only relevant annotation terms
  annot.rel.terms <- unique(annotations[annotations$gene %in% universe,]$term)
  
  annot.rel <- annotations[
    (annotations$term %in% annot.rel.terms) &
      (annotations$gene %in% universe),
  ]
  
  # Keep term descriptions if they are available
  annot.description <- if ("description" %in% colnames(annot.rel)) {
    ad <- annot.rel[!duplicated(annot.rel$term), ]
    ad[ad$term %in% annot.rel.terms, c("term", "description")]
  } else {
    NULL
  }
  
  # Test each annotation term
  terms <- unique(annot.rel.terms)
  pb <- txtProgressBar(min=1, max=length(terms), initial=0, style=3)
  
  R <- Reduce(
    rbind,
    lapply(terms, function(t) {
      
      setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
      
      annot.term <- annot.rel[annot.rel$term == t,]$gene
      
      # Compare ranked values for annotated and non-annotated genes
      q.annot <- diffex.rel$qvalue[rownames(diffex.rel) %in% annot.term]
      q.noannot <- diffex.rel$qvalue[!(rownames(diffex.rel) %in% annot.term)]
      
      r <- wilcox.test(q.annot, q.noannot, alternative="greater")
      
      as.data.frame(
        list(t, r$statistic, r$p.value),
        col.names = c("term", "statistic", "p.value")
      )
    })
  )
  
  R$qvalue <- p.adjust(R$p.value, "fdr")
  
  if (!is.null(annot.description)) {
    R <- merge(R, annot.description)
  }
  
  R
}
