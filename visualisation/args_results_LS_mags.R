library(ggplot2)
library(randomcoloR)
library(wesanderson)

L <- list.files(path = ".", pattern = ".out.mapping.ARG", all.files = FALSE,  
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

cnames <-  c('ARG', 'query_start', 'query_end', 'read_id', 'predicted_ARG_class' ,'best_hit',	
             'probability', 'identity', 'alignment_length', 'alignment_bitscore', 'alignment_evalue', 'counts')

for (i in 1:length(L) ){
  print(i)
  if ( nrow(read.delim(L[i])) > 0 ) {
  x <- read.table(L[i], head = TRUE )
  colnames(x) <- cnames 
  x$sample <- strsplit( L[i], split= ".out.mapping.ARG" ) [[1]]
  x$Flight <- paste0("F",substr( strsplit( strsplit( L[i], split= ".out.mapping.ARG" ) [[1]] , split="_")[[1]][2] , 1,1) )
  x$location <- strsplit( strsplit( L[i], split= ".out.mapping.ARG"  ) [[1]] , split="_")[[1]][3] 
  if(i==1){ ggdata <- x}
  if(i!=1){ ggdata <- rbind(ggdata, x ) }
  }
}

# remove multidrug class
ggdata <- ggdata [which(ggdata$predicted_ARG_class != 'multidrug' ) , ]

# barplot 1
ggdata2 <- as.data.frame(table(ggdata$predicted_ARG_class))
ggdata2$Var1 <- factor(ggdata2$Var1)
ggdata2 <- ggdata2[which(ggdata2$Freq>0),]

ggplot(ggdata, aes(x=predicted_ARG_class, fill=Flight) ) +
  geom_bar(stat="count", width=0.7 ) +
  theme_classic() + coord_flip() + 
  theme(legend.position="right") + ylab('ARG genes predicted in ORFs') + xlab("antibiotic class" ) + 
  scale_x_discrete (limits = ggdata2$Var1[ sort(ggdata2$Freq, decreasing=T, index.return=T)$ix] ) +
  scale_fill_manual(values=wes_palette(n=3, name="Moonrise2")) 

ggsave('ls_mags_barplot1_without_multidrug.pdf', width=6, height = 3.5)

ggplot(ggdata )+
  geom_bar(aes(y = Flight, fill=location), width=0.7 ) +
  theme_classic() + coord_flip() + 
  theme(legend.position="left")  + ylab('Flight') + xlab("ARG genes predicted in ORFs") + 
  scale_fill_brewer(palette="Accent") 

ggsave('ls_mags_barplot2_without_multidrug.pdf', width=3, height = 3.5)

# read color scale
COL <- read.csv('colorscale.csv')
COL2 <- as.vector(COL$COL)
names(COL2) <- COL$sort.unique.ggdata.predicted_ARG_class..

# barplot 2
ggplot(ggdata, aes(x=factor(sample), fill= predicted_ARG_class))+
  geom_bar(stat="count", width=0.7, color="black", size=0.15 ) +
  theme_classic() + 
  theme(legend.position="right")   + scale_fill_manual(values=COL2 ) + ylab('ARGs predicted in ORFs') + xlab("") +
  theme(axis.text.y = element_text(face="bold") ) +
  theme(axis.text.y=element_text(face="bold", angle = 0, hjust = 0,  colour = c( rep("#798E87",3 ),rep("#C27D38",8 ),rep("#CCC591",8 ) )  ) ) +
  guides(fill=guide_legend(title="Predicted ARG class")) + ggtitle('MAGs with at least 1 predicted ARG') + coord_flip()
        
ggsave('ls_mags_barplot3_without_multidrug.pdf', width=9, height = 5 )

# scatterplot
ggplot(ggdata, aes(identity, probability, color=location, shape=Flight))+ 
  geom_point(alpha=0.7, size=2.5)  +  theme_classic() +  scale_color_manual(values=COL2 ) +
  scale_color_brewer(palette="Accent") +
  guides(color = FALSE) 

ggsave('ls_mags_scatterplot_without_multidrug.pdf', width=4.7, height = 4)

#plot(ggdata$probability, ggdata$identity)
