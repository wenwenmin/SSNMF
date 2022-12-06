# -------------------------------------------------------
# -------------------------------------------------------
# 1_pollenData
# Samples: 301
# Genes: 8747
# Classes: 4
#
#     SONMF      ONMF      SNMF       NMF 
# 0.8495811 0.8313297 0.8251795 0.8495811 
load("../My_data/1_pollenData/pollenData.RData")
dataSet = pollenData
Exp = dataSet$Exp
X = log2(Exp+1)
true = dataSet$sample.group
ncl = length(unique(true))
A = as.matrix(dataSet$ppi_Net)
# -------------------------------------------------------
# -------------------------------------------------------
# 2_yanData
# 1_pollenData
# Samples: 90
# Genes: 9307
# Classes: 7
#
# SONMF      ONMF      SNMF       NMF 
# 0.8034293 0.7966970 0.6730864 0.6616088 
load("../My_data/2_yanData/yanData.RData")
dataSet = yanData

Exp = dataSet$Exp
X = log2(Exp+1)
true = dataSet$sample.group
ncl = length(unique(true))

A = as.matrix(dataSet$ppi_Net)
# -------------------------------------------------------
# -------------------------------------------------------
# 3_camp1Data
# Samples: 425
# Genes: 8058
# Classes: 5
#
# SONMF      ONMF      SNMF       NMF 
# 0.8198361 0.9215709 0.4555782 0.8106942 
load("../My_data/3_camp1Data/camp1Data.RData")
dataSet = camp1Data

Exp = dataSet$Exp
X = log2(Exp+1)
true = dataSet$sample.group
ncl = length(unique(true))

A = as.matrix(dataSet$ppi_Net)
# -------------------------------------------------------
# -------------------------------------------------------
# 4_Filbin_gliomas_science2018
# Samples: 2458
# Genes: 6719
# Classes: 6
#
# SONMF      ONMF      SNMF       NMF 
# 0.5729846 0.5392929 0.4598188 0.4664106 
load("../My_data/4_Filbin_gliomas_science2018/filbinData.RData")
dataSet = filbinData

Exp = dataSet$Exp
X = log2(Exp+1)
true = dataSet$sample.group
ncl = length(unique(true))

A = as.matrix(dataSet$ppi_Net)
# -------------------------------------------------------
# -------------------------------------------------------
# 5_Lake_Brain_science2016
# Samples: 3042
# Genes: 11207
# Classes: 16
#
# SONMF      ONMF      SNMF       NMF 
# 0.6312378 0.6685399 0.5941330 0.5941330 

load("../My_data/5_Lake_Brain_science2016/lakeData.RData")
dataSet = lakeData
Exp = dataSet$Exp
X = log2(Exp+1)
true = dataSet$sample.group
ncl = length(unique(true))

A = as.matrix(dataSet$ppi_Net)
# -------------------------------------------------------
# -------------------------------------------------------