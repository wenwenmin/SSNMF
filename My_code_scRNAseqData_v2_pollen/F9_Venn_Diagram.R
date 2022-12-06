# Load library
library(VennDiagram)
library(RColorBrewer)

get_venn_4sets = function(set_list,
                          myfilename,
                          category_names){
  myCol <- brewer.pal(4, "Pastel2")
  # Chart
  fig = venn.diagram(
    x = set_list,
    category.names = category_names,
    filename = myfilename,
    output=TRUE,
    
    # Output features
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 1,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .5,
    # fontface = "bold",
    # fontfamily = "sans",
    
    #Set names
    cat.cex = 0.6,
    #cat.fontface = "bold",
    cat.default.pos = "outer",
    # cat.pos = c(-27, 27, 135),
    # cat.dist = c(0.055, 0.055, 0.085),
    # cat.fontfamily = "sans",
    # rotation = 1
  )
  return(fig)
}

# set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
# set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
# set3 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
# set4 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
# 
# set_list = list(set1, set2, set3, set4)
# category_names = c("Set 1" , "Set 2 " , "Set 3", "Set 4")
# myfilename = 'F9_venn_diagramm2.png'
# 
# fig = get_venn_4sets(set_list,myfilename, category_names)


# set_list = list(biclusters[[1]]$genes, 
#                 biclusters[[2]]$genes, 
#                 biclusters[[3]]$genes, 
#                 biclusters[[4]]$genes)

set_list = list(modules[[1]], 
                modules[[2]], 
                modules[[3]], 
                modules[[4]])

category_names = c("M1" , "M2 " , "M3", "M4")
myfilename = 'F9_venn_diagramm.png'
fig = get_venn_4sets(set_list,myfilename, category_names)
#grid.draw(fig)

