library(ggplot2)
library(randomcoloR)
library(wesanderson)

L <- list.files(path = ".", pattern = "_FLT_.out.mapping.ARG", all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

cnames <-  c('ARG',	'query_start',	'query_end',	'read_id',	'predicted_ARG_class'	,'best_hit'	,
          'probability',	'identity',	'alignment_length',	'alignment_bitscore',	'alignment_evalue',	'counts')

for (i in 1:length(L) ){
  x <- read.table(L[i], head=TRUE)
  colnames(x) <- cnames 
  head(x)
  dim(x)
  x$sample <- strsplit( L[i], split= "_FLT_.out.mapping.ARG"  ) [[1]]
  x$Flight <- strsplit( strsplit( L[i], split= "_FLT_.out.mapping.ARG"  ) [[1]] , split="_")[[1]][1] 
  x$location <- strsplit( strsplit( L[i], split= "_FLT_.out.mapping.ARG"  ) [[1]] , split="_")[[1]][2] 
  
  if(i==1){ ggdata <- x}
  if(i!=1){ ggdata <- rbind(ggdata, x ) }
  
}

# remove miltidrug class
ggdata <- ggdata [which(ggdata$predicted_ARG_class != 'multidrug' ) , ]


#barplot 1
ggdata2 <- as.data.frame(table(ggdata$predicted_ARG_class))
ggdata2$Var1 <- factor(ggdata2$Var1)

ggdata2 <- ggdata2[which(ggdata2$Freq>0),]

ggplot(ggdata, aes(x=predicted_ARG_class, fill=Flight) )+
  geom_bar(stat="count", width=0.7 ) +
  theme_classic() +
  coord_flip() + 
  scale_fill_manual(values=wes_palette(n=3, name="Moonrise2")) +
  theme(legend.position="right")  + ylab('Read counts') + xlab("antibiotic class") + 
  scale_x_discrete (limits = ggdata2$Var1[ sort(ggdata2$Freq, decreasing=T, index.return=T)$ix])

ggsave('barplot1_without_multidrug.pdf', width=6, height = 4.5)

#barplot 2
COL <-  distinctColorPalette(k = length(unique( ggdata$predicted_ARG_class )  ), altCol = FALSE, runTsne = FALSE)
ggplot(ggdata, aes(x=factor(sample), fill= predicted_ARG_class))+
  geom_bar(stat="count", width=0.7, color="black", size=0.15 )+
  theme_classic()  + #coord_flip() + 
  theme(legend.position="right")   + scale_fill_manual(values=COL ) + ylab('Read counts') + xlab("")+
  theme(  axis.text.y = element_text(face="bold") )+
  theme(axis.text.x=element_text(angle = -90, hjust = 0,  colour = c( rep("#798E87",7 ),rep("#C27D38",8 ),rep("#CCC591",7 ) )  ) ) +
  guides(fill=guide_legend(title="Predicted ARG class"))
        
ggsave('barplot2_without_multidrug.pdf', width=9, height = 5)

#save color scale
D <- data.frame(COL, sort(unique(ggdata$predicted_ARG_class )) )
barplot(1:nrow(D), col=COL, names.arg = D$sort.unique.ggdata.predicted_ARG_class..,las=2)
write.csv(D, file='colorscale.csv', row.names = F)

