library(ComplexHeatmap)
library(circlize)

L <- list.files(path = ".", pattern = ".out.mapping.ARG", all.files = FALSE,  
    full.names = FALSE, recursive = FALSE,
    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

cnames <-  c('ARG',	'query_start', 'query_end', 'read_id', 'predicted_ARG_class' ,'best_hit',
          'probability', 'identity', 'alignment_length', 'alignment_bitscore', 'alignment_evalue', 'counts')

for (i in 1:length(L) ){
    if (nrow(read.delim(L[i])) > 0) {
    x <- read.table(L[i], head=TRUE)
    colnames(x) <- cnames 
    x$sample <- strsplit( L[i], split= ".out.mapping.ARG"  ) [[1]]
    x$isolate <- strsplit( strsplit( L[i], split= ".out.mapping.ARG"  ) [[1]] , split="_")[[1]][2] 
    if(i==1 | i==2 ){ ggdata <- x }
    if(i!=1 & i!=2) { ggdata <- rbind(ggdata, x ) }
    }
}

Iso  <- unique(ggdata$isolate)
Anti <- as.character(unique(ggdata$predicted_ARG_class) )
#ARGs matrix
M <- matrix(data=NA, ncol=length(Iso), nrow=length(Anti) )

for( j in 1:length(Iso) ){
    for( k in 1:length(Anti) ){    
        temp <- ggdata [ which(as.character( ggdata$predicted_ARG_class ) == Anti[k]  )  , ]
        M[k,j] <- length(    which(  temp$isolate == Iso[j]  ) )
        rm(temp)
    }
} 
rownames(M) <- Anti
colnames(M) <- Iso 

# read metadata
m1 <- rep(NA, length (Iso))
m2 <- rep(NA, length (Iso))
meta1 <- read.csv('blast.csv', header = TRUE )
meta2 <- read.csv('MT-1_metadata.csv', header = TRUE)
for( j in 1:length(Iso) ){
    temp  <- which( as.character(meta1$sample) == paste0("ORF_", Iso[j]) ) 
    temp2 <- which( as.character(meta2$sample) == Iso[j]) 
    if (length(temp)>0){
        m1[j] <- as.character( meta1$Species[temp] )
    }
    if (length(temp2)>0){
        m2[j] <- as.character( meta2$DESCRIPTION[temp2] )
    }
    rm(temp,temp2)
}
colnames(M) <- paste0( m1," (", m2,")")


column_ha = HeatmapAnnotation( bar1 = anno_barplot( colSums(M) ))  
column_haBottom = HeatmapAnnotation( species = m1   , annotation_legend_param = list( species = list( labels=unique(sort(m1)), ncol = 1 ))  )  
row_ha = rowAnnotation(foo2 = anno_barplot(rowSums(M), which = "row")   )

pdf('fig3.pdf', width=30, height = 9)
Heatmap(M, top_annotation = column_ha, 
    raster_quality=2,
    col = colorRamp2(c(0, 10, 60), c("white", "cyan", "darkblue")),
    bottom_annotation = column_haBottom, name='ARG count', column_names_max_height =unit(140, "mm"), top_annotation_height =unit(20, "mm")) + row_ha
dev.off()

pdf('fig3_2.pdf', width=12, height = 6)
Heatmap(M, show_column_names = FALSE,
    col = colorRamp2(c(0, 10, 60), c("white", "cyan", "darkblue")),
    raster_quality=2,
    top_annotation = column_ha, name='ARG count', column_names_max_height =unit(140, "mm"), top_annotation_height =unit(10, "mm")) + row_ha
dev.off()
