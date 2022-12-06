load("F6_outputData.RData")

load("../My_code_KEGG_GOBP_functional_analysis/geneNameData/Result_pollenData_geneNames.RData")

################################################################################
# ------------------------------------------------------------------------------
get_modules_entrezID = function(biclusters, geneNames.symbols.entrezID){
  geneNames_entrezID = geneNames.symbols.entrezID$NCBI.entrezID
  geneNames_symbols = geneNames.symbols.entrezID$symbols
  modules_entrezID = list()
  for(i in 1:length(biclusters)){
    ids = which(geneNames_symbols%in%biclusters[[i]]$genes==T)
    modules_entrezID[[i]] = geneNames_entrezID[ids] 
  }
  return(modules_entrezID)
}
# ------------------------------------------------------------------------------
biclusters = modres$modules
modules = get_modules_entrezID(biclusters, geneNames.symbols.entrezID)


write.table(modules[[1]], file = "Result/F6_m1_geneList_entrezID.txt", append = FALSE, quote = F, sep = "",
            eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F, qmethod = "double")

################################################################################
# KEGG_GOBP_functional_analysis
load("../My_code_KEGG_GOBP_functional_analysis/InputData_EnrichmentAnalysis.RData")
source('../My_code_KEGG_GOBP_functional_analysis/F2_enrichment_analysis_20201203.R')

geneNames_entrezID = geneNames.symbols.entrezID$NCBI.entrezID
modules_information = list(Gene.Modules=modules, 
                           AllConsideredGenes.entrezId=geneNames_entrezID)

Res = modules_functional_analysis(modules_information,known_functional_sets)

fileName1 = "Number of Enriched terms"
fileName2 = "Enriched KEGG Terms"
fileName3 = "Enriched GOBP terms"
fileName4 = "Enriched reactome terms"

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, fileName1)
addWorksheet(wb, fileName2)
addWorksheet(wb, fileName3)
addWorksheet(wb, fileName4)

writeData(wb, sheet = fileName1, Res$stat_result,rowNames=T)
writeData(wb, sheet = fileName2, Res$KEGG)
writeData(wb, sheet = fileName3, Res$GOBP)
writeData(wb, sheet = fileName4, Res$reactome)

saveWorkbook(wb,paste("Result/Enriched_anlysis_pollenData_NMF_l20_v2.xlsx"),overwrite=T)
################################################################################
################################################################################




