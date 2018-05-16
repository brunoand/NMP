custom.color<-function(n){
# This function provides a vector size n with RGB colors
# with hexadecimal values. It has a spectrum smilar to 
# rainbow colors. It was successfully tested for n=1:1000
# Author: Cintia C Palu
# NSilico & UCC
  
  palette=array(NA,n)
  
  if(n>54){
    f=as.hexmode(floor(240/(n/7)))
  }else if(n>46){
    f=as.hexmode(33)
  }else if(n>38){
    f=as.hexmode(39)
  }else if(n>30){
    f=as.hexmode(45)
  }else if(n>19){
    f=as.hexmode(55)
  }else {f=as.hexmode(80)}
   
  i=1
 
  g=as.hexmode(16)
  r=as.hexmode(15+f)
  b=as.hexmode(16)
  while(r<256){
    palette[i]=paste("#",r,g,b,sep="")
    i=i+1
    r=r+f
  }
  max=r-f
  g=b=as.hexmode(15+f)
  r=r-f
  while(g<r){
    palette[i]=paste("#",r,g,b,sep="")
    i=i+1
    g=g+f
  }
  
  g=r
  b=as.hexmode(16)
  
  while(b<g){
    palette[i]=paste("#",r,g,b,sep="")
    i=i+1
    b=b+f
  }
  
  r=g=g-f
  b=b-2*f
  
  while(b>15){
    palette[i]=paste("#",r,g,b,sep="")
    i=i+1
    b=b-f
  }
  
  r=r-f
  b=b+f
  
  while(r>15){
    palette[i]=paste("#",r,g,b,sep="")
    i=i+1
    r=r-f
  }
 
   r=b
  
  while(g<256){
    palette[i]=paste("#",r,g,b,sep="")
    i=i+1
    g=g+f
  }

  g=g-f
  b=b+f
    while(b<256){
      palette[i]=paste("#",r,g,b,sep="")
      i=i+1
      b=b+f
    }
  g=r
  r=as.hexmode(15)+f
  
  while(b<256){
    palette[i]=paste("#",r,g,b,sep="")
    i=i+1
    b=b+f
  }  
  
  r=as.hexmode(16)
  b=b-f
  g=b
  while(g>16){
    palette[i]=paste("#",r,g,b,sep="")
    i=i+1
    g=g-f
  }  
  
  r=r+f
  while(r<256){
    palette[i]=paste("#",r,g,b,sep="")
    i=i+1
    r=r+f
  } 


# Checking for error 
  if(is.na(palette[n])){print(paste('Error on n=',n,', i=',i))}

# Palette is usually longer than n. Here we pick random n values
# within palette.
  palette=palette[sort(sample(1:length(palette),n,replace=FALSE))]
  
 return(as.vector(palette)) 
}


