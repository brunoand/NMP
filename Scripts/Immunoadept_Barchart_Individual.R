cription     :new service
# version         :1.0
# copyright       :NSilico Life Science Ltd.
# author          :Cintia Palu


#========================================================================================
#                         Required libraries, functions, info
#========================================================================================

individualMicrobiomeBarplot = function(zz){
  source("custom.color.R")
  
  require(ggplot2)
  require(ggfortify)
  require(stringr)
   
  
  inputfilename <- zz
  df.input <- read.delim(file=inputfilename, as.is=TRUE,row.names=1)#skip=1
  
  df=df.input$Metaphlan2_Analysis[order(rownames(df.input))]
  names(df)=rownames(df.input)[order(rownames(df.input))]
  
  
  for(taxon in c('|t__','|s__','|g__','|f__')){
    df.input=df[grep(paste('\\',taxon,sep=''),names(df))]

    switch(taxon,
         '|t__'={taxa='Strain'},
         '|s__'={taxa='Species'},
         '|g__'={taxa='Genus'},
         '|f__'={taxa='Family'},
         '|o__'={taxa='Order'},
         '|c__'={taxa='Class'},
         '|p__'={taxa='Phylum'},
         'k__'={taxa='Kingdom'}
    
    )
    temp.name=NULL
    switch(taxa,
           'Phylum'={temp.name=str_split(names(df.input),'p\\|_*')},
           'Class'={temp.name=str_split(names(df.input),'p\\|_*')},
           'Order'={temp.name=str_split(names(df.input),'\\|c_*')},
           'Family'={temp.name=str_split(names(df.input),'\\|o_*')},
           'Genus'={temp.name=str_split(names(df.input),'\\|f_*')},
           'Species'={temp.name=str_split(names(df.input),'\\|s_*')},#the species name already contains the genus
           'Strain'={temp.name=str_split(grep(taxon,names(df.input),value = TRUE),'\\|s_*')}
    )
  
    if(!is.null(temp.name)){
      for(i in 1:length(names(df.input))){
        if(is.na(temp.name[[i]][2])){
          names(df.input)[i]=temp.name[[i]][1]
        }else{
          names(df.input)[i]=temp.name[[i]][2]
        }
      }
    }else{
      names(df.input)=sub('k_*','',names(df.input))
    }
    names(df.input)=gsub('_',' ',names(df.input))
    
    
    ##########################################################################################
    ##################################### Bar Chart ##########################################
    ##########################################################################################
    
    ##########################
    # Defining the file name
    
    file=gsub("profile.txt",paste(taxa,"_BarchartByPhylo.png",sep=''), zz)
    
    
    sampleid=gsub("_profile.txt","", zz)
    nplots=ceiling(length(sampleid)/22)
    nsample = ceiling(length(sampleid)/nplots)
    
    ######################
    ##     Plotting     ##
    ######################
    # Stablishing color palette
    cl= custom.color(length(df.input))
    
    H=2+length(temp.name)*0.125
    png(file, res=720,units="in", width=H*2.15,#ceiling(max(nchar(names(df.input)))/5),#7 com 17
        height=H)
    # Defining plot area
    
    par(mfrow=c(nplots,1), oma = c(3,ceiling(max(nchar(names(df.input)))/2.125),0,0), 
        mar = c(1.5, 0, 1.2, 0.5),lwd = 0.4,cex=0.8)
    
    barplot(df.input, col=cl, mgp = c(0.5, 0.5, 0),horiz =TRUE,
            main = "", las=2, border="gray40")
      
    title( main = paste('Microbiome diversity by', taxa), 
           outer = TRUE, line = -0.8 )
    title( xlab = 'Organism Relative Presence (%)', 
           outer = TRUE, line = 1 )
    
    dev.off()
    
    df=df[-grep(paste('\\',taxon,sep=''),names(df))]
  }#for(taxon
}

