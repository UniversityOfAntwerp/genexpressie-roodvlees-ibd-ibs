# Roodvlees genexpressie analysise in IBD and IBS patiënten.

Deze repository bevat een R-analyse van genexpressie in darmbiopten van personen met IBD en IBS. De samples zijn genomen voor en na een korte interventie met rood vlees. We gebruiken hiervoor Study 2 van de publieke GEO-dataset [GSE25220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25220).

De data komen uit de paper van Hebels et al. over rood vlees, fecal water genotoxicity en veranderingen in genexpressie in de menselijke colon:

**Red meat intake-induced increases in fecal water genotoxicity correlate with pro-carcinogenic gene expression changes in the human colon**  
DOI: [10.1016/j.fct.2011.10.038](https://doi.org/10.1016/j.fct.2011.10.038)

We hebben deze dataset vooral gebruikt om te kijken of er met een mixed model duidelijke verschillen terug te vinden zijn tussen IBD en IBS, en of de roodvleesinterventie daar nog een extra effect op heeft. De analyse is dus niet bedoeld als volledige reproductie van de originele paper, maar eerder als een eigen uitwerking op basis van dezelfde publieke data.

De vragen waar we vooral naar kijken zijn:                                      

- verschillen IBD en IBS samples in genexpressie?
- verandert genexpressie na de roodvleesinterventie?
- is het effect van de interventie anders bij IBD dan bij IBS?

## Repository-structuur

```text
.
├── README.md                                          #uitleg over het project
├── code                                               #Rscripts voor de analyse
│   ├── functions.R                                    #eigen functies
│   ├── import_data.R                                  #data inladen
│   ├── venn-diagram.R                                 #venn-diagram weergeven    
│   ├── libraries.R                                    #packages inladen
│   └── run_analysis.R                                 #analyse uitvoeren
├── data                                               #gebruikte data
│   ├── Ensembl2Reactome_All_Levels.human.txt          #genen koppelen aan pathways
│   ├── GSE25220_series_matrix.txt.gz                  #genexpressie data
│   └── gProfiler_hsapiens_2026-05-06_12-58-07.csv     #extra geninformatie
└── results                                            #resultaten van de analyse
    ├── R8_3_disease_fullmodel                         #verschil tussen ibd en ibs
    ├── R8_3_interaction_fullmodel                     #interactie met de interventie
    └── R8_3_intervention_fullmodel                    #effect van rood vlees
```

## Opbouw van de code

Het script waarmee de analyse start is:

```text
code/run_analysis.R
```

Dat script laadt eerst de andere bestanden in:

```r
source("code/libraries.R")
source("code/functions.R")
source("code/import_data.R")
```

De code opgesplitst in aparte bestanden omdat het anders snel onoverzichtelijk wordt. De algemene structuur is:

- `libraries.R`: laadt de packages.
- `functions.R`: bevat de functies voor de modellen, volcano plots, enrichment en GSEA.
- `import_data.R`: leest de extra databestanden uit de map `data` in.
- `run_analysis.R`: voert de analyse zelf uit.
- `venn-diagram.R`: illustratie functionaliteiten m.b.v. venn-diagram, evenals fischer's exact testen voor significantie overlappen

De belangrijkste packages zijn `GEOquery`, `lmerTest`.

## Data inladen

De expressiedata en metadata worden opgehaald via GEO:

```r
gse <- getGEO("GSE25220", GSEMatrix = TRUE)
expr <- exprs(gse[[1]])
meta <- pData(gse[[1]])
```

Hierbij is `expr` de matrix met expressiewaarden en `meta` de tabel met sample-informatie.

De platformannotatie wordt ook ingeladen:

```r
gpl <- getGEO(annotation(gse[[1]]))
annot <- Table(gpl)
```

Die annotatie is later nodig om de probe-ID’s te koppelen aan gennamen. Zonder die stap zijn de resultaten moeilijker te interpreteren, omdat je dan vooral met probe-ID’s blijft zitten.

## Data voorbereiden

Eerst worden rijen zonder bruikbare probe-ID verwijderd:

```r
expr <- expr[rownames(expr) != "" & !is.na(rownames(expr)), ]
```

Daarna worden de kolomnamen aangepast op basis van de sampletitels. Dat maakt het makkelijker om de expressiematrix en de metadata naast elkaar te gebruiken.

Samples die niet deelnamen aan de roodvleesinterventie worden verwijderd:

```r
remove_samples <- meta$title[meta$description == "No intervention"]
expr <- expr[, !colnames(expr) %in% remove_samples]
meta <- meta[meta$description == "Participated in red meat intervention", ]
```

Die samples eruit gehaald omdat de analyse vooral draait rond het verschil voor en na de interventie. Samples zonder interventie zouden in dit geval vooral extra ruis geven.

## Metadata voorbereiden

Voor het model moeten enkele variabelen als factor worden behandeld:

```r
meta$intervention.status.ch1 <- as.factor(meta$intervention.status.ch1)
meta$subject.number.ch1 <- as.factor(meta$subject.number.ch1)
meta$disease.status.ch1 <- as.factor(meta$disease.status.ch1)
meta$gender.ch1 <- as.factor(meta$gender.ch1)
```

Leeftijd wordt numeriek gemaakt:

```r
meta$age.ch1 <- as.numeric(meta$age.ch1)
```

Daarna worden de referentiegroepen ingesteld:

```r
meta$disease.status.ch1 <- relevel(meta$disease.status.ch1, ref = "IBS")
meta$intervention.status.ch1 <- relevel(meta$intervention.status.ch1, ref = "before")
```

We gebruiken IBS als referentie, zodat het IBD-effect geïnterpreteerd kan worden ten opzichte van IBS. Voor de interventie gebruik ik `before` als referentie, waardoor het effect van `after` overeenkomt met de verandering na de interventie.

## Mixed model

Voor elk gen wordt een mixed model gefit. Het model is:

```r
~ intervention.status.ch1:disease.status.ch1 +
  disease.status.ch1 +
  intervention.status.ch1 +
  gender.ch1 +
  age.ch1 +
  (1 | subject.number.ch1)
```

De random intercept voor `subject.number.ch1` is belangrijk omdat sommige personen meerdere metingen hebben. Zonder die term zouden die metingen te veel als volledig onafhankelijke observaties behandeld worden.

Met dit model kunnen drie effecten bekeken worden:

| Effect | Coëfficiënt |
|---|---|
| Verschil tussen IBD en IBS | `disease.status.ch1IBD` |
| Effect na de interventie | `intervention.status.ch1after` |
| Verschillend interventie-effect bij IBD | `disease.status.ch1IBD:intervention.status.ch1after` |

In de huidige versie van de code wordt vooral het disease-effect actief uitgevoerd. De intervention- en interaction-analyses staan ook in het script, maar zijn uitgecommentarieerd. De enige veranderingen vereist hierna is de "R8_3_DiseaseEffect" namen in opvbolgende lijnen veranderen naar respectievelijke bestandnamen.

## Belangrijkste functies

De functie `diffex.test.all()` voert het mixed model uit voor alle genen. Per gen worden onder andere de estimate, p-waarde, q-waarde en informatie over singular fits opgeslagen.

De q-waarde wordt berekend met FDR-correctie:

```r
R$qvalue <- p.adjust(R$`Pr(>|t|)`, method = "fdr")
```

Verder staan er functies in voor visualisatie en pathway-analyse:

- `volcano()` maakt een volcano plot.
- `diffex.enrich()` voert Reactome enrichment uit m.b.v. de Fischer's exact test of Chi-kwadraat test.
- `diffex.gsea()`  voert een GSEA analyse uit.

De volcano plot gebruikt de estimate op de x-as en `-log2(p-value)` op de y-as. Dat geeft snel een eerste beeld van genen met een groter effect en een lage p-waarde.

## Resultaten verwerken

Na het fitten van de modellen worden resultaten met problemen verwijderd, bijvoorbeeld ontbrekende waarden of singular fits. Die stap is nodig omdat sommige genen geen betrouwbaar modelresultaat opleveren.

Daarna worden probe-ID’s gekoppeld aan gennamen:

```r
map <- setNames(annot$GENE_SYMBOL, annot$ID)
R8_3_DiseaseEffect$gene <- map[rownames(R8_3_DiseaseEffect)]
```

Wanneer meerdere probes naar hetzelfde gen verwijzen, wordt de probe met de laagste q-waarde behouden. Dat is een praktische keuze om per gen maar één resultaat over te houden.

Daarna worden gennamen gekoppeld aan Ensembl-ID’s:

```r
map <- setNames(geneiddata$converted_alias, geneiddata$initial_alias)
R8_3_DiseaseEffect$ENSG_ID <- map[rownames(R8_3_DiseaseEffect)]
```

Tot slot worden Reactome-beschrijvingen toegevoegd:

```r
map <- setNames(reactome$description, reactome$gene)
R8_3_DiseaseEffect$Description <- map[R8_3_DiseaseEffect$ENSG_ID]
```

Die extra annotatie maakt het makkelijker om de resultaten biologisch te interpreteren, vooral bij de pathway-analyse.

## Resultaatbestanden gebruiken

In de map `results/` staan drie resultaatbestanden:

```text
R8_3_disease_fullmodel
R8_3_interaction_fullmodel
R8_3_intervention_fullmodel
```

Deze bestanden kunnen direct worden ingelezen. Dat is handig als je niet telkens de volledige analyse opnieuw wilt uitvoeren.

```r
disease_results <- read.csv(
  "results/R8_3_disease_fullmodel",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

interaction_results <- read.csv(
  "results/R8_3_interaction_fullmodel",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

intervention_results <- read.csv(
  "results/R8_3_intervention_fullmodel",
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
```
