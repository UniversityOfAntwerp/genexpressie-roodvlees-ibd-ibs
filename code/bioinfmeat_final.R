if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") 

BiocManager::install("GEOquery") 


library(GEOquery) 

library(lmerTest)

gse <- getGEO("GSE25220", GSEMatrix = TRUE) 
expr <- exprs(gse[[1]]) 
meta <- pData(gse[[1]]) 

dim(expr) 
head(expr[,1:3]) 
head(meta) 

gpl <- getGEO(annotation(gse[[1]]))
annot <- Table(gpl) 


#kolomnamen van expressietabel vervangen met samplenamen

map2 <- setNames(meta$title, rownames(meta)) 

new_names2 <- map2[colnames(expr)] 
new_names2 
colnames(expr) <- ifelse(is.na(new_names2), colnames(expr), new_names2) 

#mensen die niet meededen excluderen 

remove_samples <- meta$title[meta$description == "No intervention"] 

expr <- expr[, !colnames(expr) %in% remove_samples] 
meta <- meta[meta$description == "Participated in red meat intervention", ]



#histogram van alle expressiewaarden (ze zijn al log-getransformeerd), ter illustratie

values2 <- as.vector(expr[1001, ]) 

hist(values2, breaks = 100, main = "Expression value distribution", xlab = "Log2 expression intensity") 



#diff. experession code
fit <- NULL
diffex.test.all <- function(form, data, meta, var=NULL) {
  
  #' Differential expression testing function
  #'
  #' A differential expression test for each gene in a dataframe of normalised gene expression counts.
  #' The test is performed using the negative binomial distribution.
  #'
  #' @param form A formula
  #' @param data A data frame of normalised expression counts.
  #'        Rows are genes, columns are samples.
  #' @param meta A metadata dataframe
  #'
  #' @return A data frame of test results with coefficients and FDR-adjusted p-values
  
  require(MASS)
  library(lmerTest)
  
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
        
        res <- if (is.null(var)){as.data.frame(summary(fit)$coefficients)[2, ]} else {as.data.frame(summary(fit)$coefficients)[var, ]}
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

#gemakkelijkere kolomnamen installeren
colnames(meta) <- make.names(colnames(meta))

#asfactoren maken waar toepasselijk
meta$intervention.status.ch1 <- as.factor(meta$intervention.status.ch1)
meta$subject.number.ch1 <- as.factor(meta$subject.number.ch1)
meta$disease.status.ch1 <- as.factor(meta$disease.status.ch1)
meta$gender.ch1 <- as.factor(meta$gender.ch1)

#leeftijd als numerisch encoderen
typeof(meta$age.ch1)
meta$age.ch1 <- as.numeric(meta$age.ch1)


#IBS als baseline zetten ivg m. IBD
meta$disease.status.ch1 <- relevel(
  meta$disease.status.ch1,
  ref = "IBS"
)



#"Before" als baseline zetten ivg m "After"
meta$intervention.status.ch1 <- relevel(
  meta$intervention.status.ch1,
  ref = "before"
)


#volledig model (voor DiseaseEffect, intervention effect, en interaction effect respectievelijk). Voor de rest van de model wordt enkel met DiseaseEffect gewerkt, maar de andere modellen werken analoog
R8_3_DiseaseEffect <- diffex.test.all(~intervention.status.ch1:disease.status.ch1 + disease.status.ch1 + intervention.status.ch1 + gender.ch1 + age.ch1 + (1| subject.number.ch1),
                       expr, 
                       meta,
                       var="disease.status.ch1IBD") 

#R8_3_InterventionEffect <- diffex.test.all(~intervention.status.ch1:disease.status.ch1 + disease.status.ch1 + intervention.status.ch1 + gender.ch1 + age.ch1 + (1| subject.number.ch1),
                                      #expr, 
                                      #meta,
                                      #var="intervention.status.ch1after") 

#R8_3_InteractionEffect <- diffex.test.all(~intervention.status.ch1:disease.status.ch1 + disease.status.ch1 + intervention.status.ch1 + gender.ch1 + age.ch1 + (1| subject.number.ch1),
                                      #expr, 
                                      # meta,
                                      # var="disease.status.ch1IBD:intervention.status.ch1after") 








#singularities en NA's eruitfilteren

R8_3_DiseaseEffect<- R8_3_DiseaseEffect[R8_3_DiseaseEffect$singular == "FALSE",]

R8_3_DiseaseEffect <- R8_3_DiseaseEffect[!is.na(R8_3_DiseaseEffect$singular), ]
R8_3_DiseaseEffect$singular <- NULL


#Voor het InteractionModel deze lijn ook runnen
#R8_3_InteractionEffect <- R8_3_InteractionEffect[!is.na(R8_3_InteractionEffect$Estimate), ]



#volcano plots maken
volcano <- function(diffex.res, q.thresh=0.05, fc.thresh=1, only_sig=T){
  p.thresh <- max(diffex.res[diffex.res$qvalue < q.thresh,]$`Pr(>|t|)`)
  significant <- diffex.res[(diffex.res$qvalue < q.thresh) & (abs(diffex.res$Estimate) >= fc.thresh),]
           pv <- -log2(significant$`Pr(>|t|)`)
              
  insignificant <- diffex.res[!((diffex.res$qvalue < q.thresh) & (abs(diffex.res$Estimate)
                                                                  >= fc.thresh)),]
  pv <- pv[is.finite(pv)]
  ylim <- c(0, max(pv))*1.1
  xlim <- if (only_sig & (nrow(significant) > 0)){
    c(min(significant$Estimate), max(significant$Estimate))*1.1
  } else {
    c(min(diffex.res$Estimate), max(diffex.res$Estimate))*1.1
  }
  plot(insignificant$Estimate, -log2(insignificant$`Pr(>|t|)`),
       xlim=xlim, ylim=ylim,
       xlab="Log Fold Change",
       ylab="-log p-value")
  points(significant$Estimate, -log2(significant$`Pr(>|t|)`),
         col=sign(significant$Estimate)+3)
  lines(xlim, -log2(c(p.thresh, p.thresh)), col='black')
  lines(-c(fc.thresh, fc.thresh), ylim, col='black')
  lines(c(fc.thresh, fc.thresh), ylim, col='black')
}
volcano(R8_3_DiseaseEffect)





#Extra kolom met gennamen maken voor expressietabel
map <- setNames(annot$GENE_SYMBOL, annot$ID)
map
R8_3_DiseaseEffect$gene <- map[rownames(R8_3_DiseaseEffect)]
R8_3_DiseaseEffect$gene[is.na(R8_3_DiseaseEffect$gene)] <- rownames(R8_3_DiseaseEffect)[is.na(R8_3_DiseaseEffect$gene)]


#Beste probe per gen behouden (slechte duplicaten wegwerken)
R8_3_DiseaseEffect <- R8_3_DiseaseEffect[order(R8_3_DiseaseEffect$qvalue), ]
R8_3_DiseaseEffect <- R8_3_DiseaseEffect[!duplicated(R8_3_DiseaseEffect$gene), ]
R8_3_DiseaseEffect = R8_3_DiseaseEffect[ (R8_3_DiseaseEffect$gene !=""),  ]


#Row names veranderen door Gennamen
rownames(R8_3_DiseaseEffect) <- R8_3_DiseaseEffect$gene






#ENSG's voor overeenkomstige genen mappen [DOWNLOAD EERSTE de gPROFILER dataset (en importeer ze als "geneiddata")!]
map <- setNames(geneiddata$converted_alias, geneiddata$initial_alias)

map[R8_3_DiseaseEffect$gene]
R8_3_DiseaseEffect$ENSG_ID <- map[rownames(R8_3_DiseaseEffect)]

R8_3_DiseaseEffect <- R8_3_DiseaseEffect[(R8_3_DiseaseEffect$ENSG_ID != ""), ] #

R8_3_DiseaseEffect <- R8_3_DiseaseEffect[!is.na(R8_3_DiseaseEffect$ENSG_ID), ]



#Download de EnsemblReactome-dataset en noem ze 'reactome'
colnames(reactome)[colnames(reactome) == "V1"] <- "gene"
colnames(reactome)[colnames(reactome) == "V2"] <- "term"
colnames(reactome)[colnames(reactome) == "V3"] <- "description"


#Map ENSG's met overeenkomstige beschrijvingen (indien deze aanwezig zijn).
map <- setNames(reactome$description, reactome$gene)

map[R8_3_DiseaseEffect$ENSG_ID]
R8_3_DiseaseEffect$Description <- map[R8_3_DiseaseEffect$ENSG_ID]

R8_3_DiseaseEffect <- R8_3_DiseaseEffect[!is.na(R8_3_DiseaseEffect$Description), ]

R8_3_DiseaseEffect <- R8_3_DiseaseEffect[!duplicated(R8_3_DiseaseEffect$ENSG_ID), ]

#Vervang rownames door ENSG id's
rownames(R8_3_DiseaseEffect) <- R8_3_DiseaseEffect$ENSG_ID







# Define a function which performs the Fisher exact test
# For a differential expression result.
diffex.enrich <- function(diffex.res, annotations, direction="all"){
  # Define the universe
  universe <- intersect(rownames(diffex.res), annotations$gene)
  diffex.rel <- diffex.res[universe,]
  # Select the set of significant genes, depending on the direction
  # we want ti investigate up/down/all
  diffex.sig <- rownames(diffex.rel[(diffex.rel$qvalue < 0.05) &
                                      if(direction == "up") {
                                        (diffex.rel$Estimate > 0)
                                      } else if(direction == "down"){
                                        (diffex.rel$Estimate < 0)
                                      } else {TRUE},])
  # to speed up computation, select only relevant genes and annotations
  annot.rel.terms <- unique(annotations[annotations$gene %in% diffex.sig,]$term)
  annot.rel <- annotations[(annotations$term %in% annot.rel.terms) &
                             (annotations$gene %in% universe),]
  annot.description <- if ("description" %in% colnames(annot.rel)){
    ad <- annot.rel[!duplicated(annot.rel$term), ]
    ad[ad$term %in% annot.rel.terms, c("term","description")]
  } else {NULL}
  # Which terms to test?
  terms <- unique(annot.rel.terms)
  pb <- txtProgressBar(min=1, max=length(terms), initial=0, style=3)
  R <- Reduce(rbind,
              lapply(terms, function(t){
                setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
                annot.term <- annot.rel[annot.rel$term == t,]$gene
                # Determine the contingency matrix
                a <- length(intersect(diffex.sig, annot.term))
                b <- length(setdiff(annot.term, diffex.sig))
                c <- length(setdiff(diffex.sig, annot.term))
                d <- length(setdiff(universe, union(diffex.sig, annot.term)))
                contingency <- matrix(c(a,b,c,d), nrow=2)
                # Calculate the p-value. If we have enough data, use the chi2
                # approximation to calculate the p-value
                p.value <- if(min(contingency) < 5){
                  fisher.test(contingency)$p.value
                } else {
                  chisq.test(contingency)$p.value
                }
                odds.ratio <- (a/c)/(b/d)
                as.data.frame(list(t, a, b, c, d, odds.ratio, p.value),
                              col.names=c("term", "a", "b", "c", "d", "odds.ratio", "p.value"))
              })
  )
  R$qvalue <- p.adjust(R$p.value, "fdr")
  if (! is.null(annot.description)){
    R <- merge(R, annot.description)
  }
  R
}



#Functional enrichment m.b.v. Fisher's exact test uitvoeren 

enrich.resR8_3_DiseaseEffect<- diffex.enrich(R8_3_DiseaseEffect, reactome)

enrich.res.upR8_3_DiseaseEffect <- diffex.enrich(R8_3_DiseaseEffect, reactome, direction="up")

enrich.res.downR8_3_DiseaseEffect <- diffex.enrich(R8_3_DiseaseEffect, reactome, direction="down")







# Define a function which performs a Gene Set Enrichment Analysis (GSEA)
# For a differential expression result.
diffex.gsea <- function(diffex.res, annotations, direction="all"){
  # Define the universe
  universe <- intersect(rownames(diffex.res), annotations$gene)
  diffex.rel <- diffex.res[universe,]
  # Order the genes by p-value.
  # Depending on direction up/down/all, change the sign so that
  # the ordering puts the correct direction on top
  diffex.rel$qvalue <- 1 - diffex.rel$qvalue # put best results at top
  diffex.rel$qvalue <- if(direction == "up") {
    diffex.rel$qvalue * sign(diffex.rel$Estimate)
  } else if(direction == "down"){
    diffex.rel$qvalue * -1 * sign(diffex.rel$Estimate)
  } else {diffex.rel$qvalue}
  # to speed up computation, select only relevant genes and annotations
  annot.rel.terms <- unique(annotations[annotations$gene %in% universe,]$term)
  annot.rel <- annotations[(annotations$term %in% annot.rel.terms) &
                             (annotations$gene %in% universe),]
  annot.description <- if ("description" %in% colnames(annot.rel)){
    ad <- annot.rel[!duplicated(annot.rel$term), ]
    ad[ad$term %in% annot.rel.terms, c("term","description")]
  } else {NULL}
  # which terms to test?
  terms <- unique(annot.rel.terms)
  pb <- txtProgressBar(min=1, max=length(terms), initial=0, style=3)
  R <- Reduce(rbind,
              lapply(terms, function(t){
                setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
                annot.term <- annot.rel[annot.rel$term == t,]$gene
                # Make two lists with q-values for annotated and unannotated genes
                q.annot <- diffex.rel$qvalue[rownames(diffex.rel) %in% annot.term]
                q.noannot <- diffex.rel$qvalue[!(rownames(diffex.rel) %in% annot.term)]
                r <- wilcox.test(q.annot, q.noannot, alternative="greater")
                as.data.frame(list(t, r$statistic, r$p.value),
                              col.names=c("term", "statistic", "p.value"))
              })
  )
  R$qvalue <- p.adjust(R$p.value, "fdr")
  if (! is.null(annot.description)){
    R <- merge(R, annot.description)
  }
  R
}


#GSEA uitvoeren

#Ongericht
gsea.resR8_3_DiseaseEffect <- diffex.gsea(R8_3_DiseaseEffect, reactome)
#Opreguleerd
gsea.res.upR8_3_DiseaseEffect <- diffex.gsea(R8_3_DiseaseEffect, reactome, direction="up")
#Downreguleerd
gsea.res.downR8_3_DiseaseEffect <- diffex.gsea(R8_3_DiseaseEffect, reactome, direction="down")




