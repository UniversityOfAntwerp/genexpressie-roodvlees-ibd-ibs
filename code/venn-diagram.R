#installeer als nog niet geinstalleerd!
#install.packages("RColorBrewer")
#install.packages("eulerr")

library(eulerr)
library(RColorBrewer)

metabolism_termen <- c(
  "Metabolism",
  "The citric acid (TCA) cycle and respiratory electron transport",
  "Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins.",
  "Respiratory electron transport",
  "Citric acid cycle (TCA cycle)",
  "Pyruvate metabolism and Citric Acid (TCA) cycle",
  "Metabolism of amino acids and derivatives",
  "Selenoamino acid metabolism",
  "Selenocysteine synthesis",
  "Metabolism of proteins",
  "Metabolism of RNA",
  "Metabolism of non-coding RNA",
  "Mitochondrial biogenesis",
  "Mitochondrial protein import",
  "Mitochondrial translation",
  "Mitochondrial translation initiation",
  "Mitochondrial translation elongation",
  "Mitochondrial translation termination",
  "Complex I biogenesis",
  "Cristae formation",
  "Peroxisomal protein import",
  "Class I peroxisomal membrane protein import",
  "Metal ion SLC transporters",
  "Response to metal ions",
  "Metallothioneins bind metals",
  "Cytoprotection by HMOX1",
  "Response of EIF2AK4 (GCN2) to amino acid deficiency",
  "Cellular response to starvation",
  "Macroautophagy",
  "Autophagy",
  "Translation",
  "Eukaryotic Translation Initiation",
  "Eukaryotic Translation Elongation",
  "Eukaryotic Translation Termination",
  "Peptide chain elongation",
  "Translation initiation complex formation",
  "Cap-dependent Translation Initiation",
  "Viral mRNA Translation",
  "SRP-dependent cotranslational protein targeting to membrane",
  "Nonsense-Mediated Decay (NMD)",
  "Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC)",
  "Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC)",
  "rRNA processing",
  "rRNA processing in the nucleus and cytosol",
  "Major pathway of rRNA processing in the nucleolus and cytosol",
  "rRNA modification in the nucleus and cytosol",
  "RNA Polymerase III Transcription"
)
Diffexdown_Disease <- paste(gsea.res.downR8_3_DiseaseEffect$description[gsea.res.downR8_3_DiseaseEffect$qvalue < 0.05])
Diffexdown_Intervention <- paste(gsea.res.downR8_3_InterventionEffect$description[gsea.res.downR8_3_InterventionEffect$qvalue < 0.05])
Metabolism_Linked <- paste(metabolism_termen)


# kleuren
myCol <- brewer.pal(3, "Pastel2")

# de venn-diagram 
fit <- list(
  "DF-Disease" = Diffexdown_Disease,
  "DF-Intervention" = Diffexdown_Intervention,
  "Metabolism-associated" = Metabolism_Linked
)

# plotten
plot(euler(fit),
     fills = myCol,
     quantities = TRUE
)

#fisher's exact test om te zien voor significante overlap tussen
# downreguleerde functionaliteiten in disease vs intervention

A_sig <-  paste(gsea.res.downR8_3_DiseaseEffect$description[gsea.res.downR8_3_DiseaseEffect$qvalue < 0.05])  # significante termen voor disease-effect (waarbij we de overlap willen testen)
B_sig <-  paste(gsea.res.downR8_3_InterventionEffect$description[gsea.res.downR8_3_InterventionEffect$qvalue < 0.05])
universe <- gsea.res.downR8_3_DiseaseEffect$description #alle geteste downgereguleerde functionaliteiten bij de GSEA

a <- length(intersect(A_sig, B_sig))
b <- length(setdiff(A_sig, B_sig))
c <- length(setdiff(B_sig, A_sig))
d <- length(setdiff(universe, union(A_sig, B_sig)))
contingency <- matrix(c(a,b,c,d), nrow=2)

colnames(contingency) <- c("In B_sig", "Niet in B_sig")
rownames(contingency) <- c("In A_sig", "Niet in A_sig")

fisher.test(contingency)$p.value





#fisher's exact test om te testen of metabolism-associated functions overgepresenteerd zijn
  #in de significant downregulated functionalities

B <- metabolism_termen
A <- unique(c(
  gsea.res.downR8_3_DiseaseEffect$description[
    gsea.res.downR8_3_DiseaseEffect$qvalue < 0.05
  ],
  gsea.res.downR8_3_InterventionEffect$description[
    gsea.res.downR8_3_InterventionEffect$qvalue < 0.05
  ]
))

universe <- gsea.res.downR8_3_DiseaseEffect$description)  #totaal aantal termen na GSEA zijn evenveel en hetzelfde voor DiseaseEffect en InterventionEffect, dus arbitrair DiseaseEffect gekozen.

a <- length(intersect(A, B))
b <- length(setdiff(A, B))
c <- length(setdiff(B, A))
d <- length(setdiff(universe, union(A, B)))

contingency <- matrix(c(a,b,c,d), nrow = 2)

colnames(contingency) <- c("In B", "Niet in B")
rownames(contingency) <- c("In A", "Niet in A")

fisher.test(contingency)$p.value
