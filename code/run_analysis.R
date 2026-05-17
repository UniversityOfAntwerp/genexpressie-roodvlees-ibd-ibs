source("code/libraries.R")
source("code/function.R")

#matrices uit file halen
gse <- getGEO("GSE25220", GSEMatrix = TRUE) 
expr <- exprs(gse[[1]]) 
meta <- pData(gse[[1]]) 
gpl <- getGEO(annotation(gse[[1]]))
annot <- Table(gpl) 
#rijnamen van expressietabel wat duidelijker maken 
map1 <- setNames(annot$GENE_SYMBOL, annot$ID) 
new_names1 <- map1[rownames(expr)] 
rownames(expr) <- ifelse(is.na(new_names1), rownames(expr), new_names1) 
#rijen met lege gennamen deleten 
expr <- expr[rownames(expr) != "" & !is.na(rownames(expr)), ] 
#kolomnamen van expressietabel vervangen met samplenamen
map2 <- setNames(meta$title, rownames(meta)) 
new_names2 <- map2[colnames(expr)] 
colnames(expr) <- ifelse(is.na(new_names2), colnames(expr), new_names2) 
#mensen die niet meededen excluderen 
remove_samples <- meta$title[meta$description == "No intervention"] 
expr <- expr[, !colnames(expr) %in% remove_samples] 
meta <- meta[meta$description == "Participated in red meat intervention", ]
#histogram van alle expressiewaarden (ze zijn al log-getransformeerd), ter illustratie
values2 <- as.vector(expr[1001, ]) 
hist(values2, breaks = 100, main = "Expression value distribution", xlab = "Log2 expression intensity") 

#Lmertest diffex.test.all
fit <- NULL

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

#Voor het InteractionModel deze lijn ook runnen!
#R8_3_InteractionEffect <- R8_3_InteractionEffect[!is.na(R8_3_InteractionEffect$Estimate), ]

#volcano plots maken
volcano(R8_3_DiseaseEffect)

#Extra kolom met gennamen maken voor expressietabel
map <- setNames(annot$GENE_SYMBOL, annot$ID)
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

#Functional enrichment m.b.v. Fisher's exact test uitvoeren 
enrich.resR8_3_DiseaseEffect<- diffex.enrich(R8_3_DiseaseEffect, reactome)
enrich.res.upR8_3_DiseaseEffect <- diffex.enrich(R8_3_DiseaseEffect, reactome, direction="up")
enrich.res.downR8_3_DiseaseEffect <- diffex.enrich(R8_3_DiseaseEffect, reactome, direction="down")
                                    
#GSEA uitvoeren
#Ongericht
gsea.resR8_3_DiseaseEffect <- diffex.gsea(R8_3_DiseaseEffect, reactome)
#Opreguleerd
gsea.res.upR8_3_DiseaseEffect <- diffex.gsea(R8_3_DiseaseEffect, reactome, direction="up")
#Downreguleerd
gsea.res.downR8_3_DiseaseEffect <- diffex.gsea(R8_3_DiseaseEffect, reactome, direction="down")




