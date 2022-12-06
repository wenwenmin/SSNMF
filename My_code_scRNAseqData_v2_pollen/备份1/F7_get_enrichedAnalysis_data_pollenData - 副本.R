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
source('../My_code_KEGG_GOBP_functional_analysis/EnrichmentAnalysisFunctions.R')
load("../My_code_KEGG_GOBP_functional_analysis/InputData_EnrichmentAnalysis.RData")

geneNames_entrezID = geneNames.symbols.entrezID$NCBI.entrezID
modules_information = list(Gene.Modules=modules, 
                           AllConsideredGenes.entrezId=geneNames_entrezID)

Res = GeneModules.EnrichmentAnalysis(modules_information,known_functional_sets)

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

writeData(wb, sheet = fileName1, Res$modules.StatisticalResults)
writeData(wb, sheet = fileName2, Res$KEGGTerm.Table)
writeData(wb, sheet = fileName3, Res$GOBPTerm.Table)
writeData(wb, sheet = fileName4, Res$reactomeTerm.Table)

saveWorkbook(wb,paste("Result/Enriched_anlysis_pollenData_NMF_l20.xlsx"),overwrite=T)
################################################################################
################################################################################




