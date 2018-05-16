#!/usr/local/bin/Rscript

# description     :new service

# version         :2.0

# copyright       :NSilico Life Science Ltd.

# author          :Cintia Palu





#========================================================================================

#                         Required libraries, functions, info

#========================================================================================

source("Immunoadept_Barchart_Individual.R")



# library(ggplot2)

# library(ggfortify)

# library(stringr)



toDate = function(ymd){

  ymd=as.character(ymd)

  return(paste(substr(ymd,1,4),substr(ymd,5,6),substr(ymd,7,8),sep='.'))

}

# The arguments will always be inputed in pairs of Metaphlan output file + the data of collection

# All files must belong to the same patient

arguments <- commandArgs(trailingOnly = TRUE)

# Organising dates and files

arg=data.frame(Files=arguments[(2*1:(length(arguments)/2))-1],Date=arguments[2*1:(length(arguments)/2)])

# Sorting from most recent 

arg=arg[order(as.character(arg$Date),decreasing = TRUE),]

arg$Files=as.character(arg$Files)

arg$Date=toDate(arg$Date)



#Plotting individual profile for the most recent data

individualMicrobiomeBarplot(arg$Files[1])



# If there is data from previous analysis

if (length(arg$Files)>1){

  

  

  for(zz in 1:length(arg$Files)){

    #zz=1

    df.input <- read.delim(file=arg$Files[zz], as.is=TRUE)#skip=1

    rownames(df.input)=df.input[,1]

    colnames(df.input)[2]=arg$Date[zz]

    

    #df=rbind(df,df.input)

    if(zz==1){

      df=df.input

    }else{

      df=merge(df,df.input,all=TRUE)

    }

  } #for(zz in 1:length(arg$Files)

  

  rownames(df)=df[,1]

  df=df[,-1]

  

  df=df[order(rownames(df)),]

  

    

    # T=c('|t__','|s__','|g__','|f__','|o__','|c__','|p__','k__')

    # taxon=T[1]

  zz=arg$Files[1]

    for(taxon in c('|t__','|s__','|g__','|f__')){

        df.input=df[grep(paste('\\',taxon,sep=''),rownames(df)),]

    

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

               'Phylum'={temp.name=str_split(rownames(df.input),'p\\|_*')},

               'Class'={temp.name=str_split(rownames(df.input),'p\\|_*')},

               'Order'={temp.name=str_split(rownames(df.input),'\\|c_*')},

               'Family'={temp.name=str_split(rownames(df.input),'\\|o_*')},

               'Genus'={temp.name=str_split(rownames(df.input),'\\|f_*')},

               'Species'={temp.name=str_split(rownames(df.input),'\\|s_*')},#the species name already contains the genus

               'Strain'={temp.name=str_split(grep(taxon,rownames(df.input),value = TRUE),'\\|s_*')}

        )

      

        if(!is.null(temp.name)){

          for(i in 1:length(rownames(df.input))){

            if(is.na(temp.name[[i]][2])){

              rownames(df.input)[i]=temp.name[[i]][1]

            }else{

              rownames(df.input)[i]=temp.name[[i]][2]

            }

          }

        }else{

          rownames(df.input)=sub('k_*','',rownames(df.input))

        }

        rownames(df.input)=gsub('_',' ',rownames(df.input))

        

        

        ##########################################################################################

        ##################################### Bar Chart ##########################################

        ##########################################################################################

        

        ##########################

        # Defining the file name

        

        file=gsub("profile.txt",paste(taxa,"_Barchart",dim(df)[2],"Series0.png",sep=''), zz)

        

        

        sampleid=gsub("_profile.txt","", zz)

        nplots=ceiling(length(sampleid)/22)

        nsample = ceiling(length(sampleid)/nplots)

        

        ######################

        ##     Plotting     ##

        ######################

        # Stablishing color palette

        cl= custom.color(dim(df.input)[2])

        

        H=2+length(temp.name)*0.125

        png(file, res=720,units="in", width=H*2.15+(0.125*dim(df.input)[2]),#ceiling(max(nchar(rownames(df.input)))/5),#7 com 17

            height=H+(0.125*dim(df.input)[2]))

        

        # Defining plot area

        par(mfrow=c(nplots,1), oma = c(1,ceiling(max(nchar(rownames(df.input)))/2.125),0,0), 

            mar = c(1.5, 0, 1.2, 10.5),lwd = 0.4,cex=0.8)

        

        barplot(t(as.matrix(df.input)), beside = TRUE,col=cl,#legend.text=colnames(df.input),

                mgp = c(0.5, 0.5, 0),horiz =TRUE,main = "", las=2, border="gray40")

                #args.legend = list(x ='topright', bty='n', inset=c(1.25,0)))

          

        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)

        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")







        legend("topright", colnames(df.input)[length(cl):1], xpd = TRUE, horiz = FALSE,bty = 'n',

               inset = c(0, 0), col = cl[length(cl):1], cex = 1.2, pch=15)

        title( main = paste('Microbiome diversity by', taxa), 

               outer = TRUE, line = -0.8 )

        title( xlab = 'Organism Relative Presence (%)', 

               outer = TRUE, line = 1 )

        

        dev.off()

        

      #}#if(taxon!='|t__')

      df=df[-grep(paste('\\',taxon,sep=''),rownames(df)),]

    }#for(taxon



}#length(arg$Files)>1


