
sim <- function(gens,start,s) {
  ss <- sum(start)
  generation <- function(n) {
    p <- c(2*n[1]+n[2],n[2]+2*n[3])
    sel <- c(s[1]*p[1]*p[1],s[2]*2*p[1]*p[2],s[3]*p[2]*p[2])
    sel <- sel/sum(sel)
    x <- rmultinom(1,ss,prob=sel)
    x[,1] 
  }
  st <- start
  pp <- numeric(gens)
  for (i in 1:gens) {
    xx <- generation(st)
    pp[i] <- c(xx[2]+2*xx[3])/(2*sum(xx))
    st <- xx
    if (st[1]==ss) {
      break
    } else if (st[3]==ss) {
      pp[i:gens] <- 1
      break
    }
  }
  pp
}

getFirstAbsorption <- function(xx) {
  ## Find the row for which p is first zero.  
  ## 0 for no absorbtion, positive for 0 negative for 1
  lastValues <- xx[nrow(xx),]
  nc <- ncol(xx)
  ret <- numeric(nc)     ### Set to zero
  for (col in 1:nc) {
    if (lastValues[col]==0) {
      ret[col] <- min(which(xx[,col]==0))
    } else if (lastValues[col]==1) {
      ret[col] <- -min(which(xx[,col]==1))
    }
  }
return(ret)
}


library(shiny)

# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {

  output$ui <- renderUI({    
    if (input$selection == "het_input") 
      return(numericInput(inputId="s" , label = "Selection coefficient, s", value=0.02, min=0, max=0.2,0.01))
    if (input$selection == "positive")
      return(numericInput(inputId="s" , label = "Selection coefficient, s", value=0.02, min=0, max=0.2,0.01))
    if (input$selection == "homo_recess") 
      return(numericInput(inputId="s" , label = "Selection coefficient, s", value=0.02, min=0, max=0.2,0.01))
    return()
  })
  
  output$distPlot <- renderPlot({
    ## Generates the allele frequency paths and plot    
    N <- as.integer(isolate(input$N))
    alpha <- 80
    myselection <- isolate(input$selection)
    pops <- isolate(input$pops)
    generations <- isolate(input$gens)
    s <- isolate(input$s)
    startingP <- isolate(input$startingP)
    input$run
    if (myselection=="none") {
      sels <- c(1,1,1)
    } else if (myselection=="het_input") {
      sels <- c(1,1+s,1)
    } else if (myselection=="homo_recess") {
      sels <- c(1,1,1-s)
    } else {
      sels <- c(1,1+s,1+2*s)
    }   

    if (startingP==1) {
      start <- c(N-1,1,0) 
    } else if (startingP==2) {
        start <- rmultinom(1,N,prob=c(0.99*0.99,2*0.01*0.99,0.01*0.01))
    } else if (startingP==3) {
        start <- rmultinom(1,N,prob=c(0.95*0.95,2*0.05*0.95,0.05*0.05))
    } else {
        start <- rmultinom(1,N,prob=c(0.25,0.5,0.25))
    }

   p <-  replicate(pops,
                   sim(generations, start=start, s=sels))
   breaks <- seq(0,1,by=0.02)
   hs <- hist(p[generations,], breaks=breaks, plot=FALSE)
   
   old.par <- par(no.readonly=TRUE)
   mar.default <- par('mar')
   mar.left <- mar.default
   mar.right <- mar.default
   mar.left[4] <- 0
   mar.right[2] <- 0
   
   # Main plot 
   par (fig=c(0,0.8,0,1.0), mar=mar.left)
   plot(NULL, type="l", xlim=c(0,input$gens)
        ,ylim=c(0,1)
        ,ylab="Frequency of Variant Allele"
        ,xlab="Generation", axes=FALSE
        ,cex.lab=1.6, cex.axis=1.6
   )

    apply(p, 2, function(x)  {
      succ <- 1+(x[generations]>0) + (x[generations]==1.0);
      lines(x,col=c(
        rgb(255,0,0, max=255,alpha=alpha),
        rgb(0, 255, 0, max=255, alpha=alpha),
        rgb(0, 0, 255, max=255, alpha=alpha))[succ],lwd=2)
       })
   text(generations, 0.97, paste(sum(p[generations,]==1),"fixed"),col="black")
   text(generations, 0.03, paste(sum(p[generations,]==0),"lost"),col="black")
 
   firstAbsorption <- getFirstAbsorption(p)
   ones <- sum(firstAbsorption<0)
   segregating <- sum(firstAbsorption==0)
   
   if (sum(firstAbsorption>0) > 0) {
    rug(jitter(firstAbsorption[firstAbsorption>0]),col="red",ticksize=-0.03,lwd=2)
    text(median(firstAbsorption[firstAbsorption>0]), -0.02, "*", pch=1, col="red", cex=2 )
   }
   if (sum(firstAbsorption<0) > 0) {
     rug(jitter(-firstAbsorption[firstAbsorption<0]),side=3,col="blue",ticksize=-0.03,lwd=2)
     text(median(-firstAbsorption[firstAbsorption<0]), 1.02, "*", pch=1, col="blue", cex=2 )
   }
   
    axis(1);axis(2)
 
 par (fig=c(0.8,1.0,0.0,1.0), mar=mar.right, new=TRUE)
 plot (NA, type='n', axes=FALSE, yaxt='n',
       xlab='Frequency', ylab=NA, main=NA,
       xlim=c(0,max(hs$counts)),
       ylim=c(1,length(hs$counts)))
 arrows(rep(0,length(hs$counts)), 1:length(hs$counts),
        hs$counts, 1:length(hs$counts),
        length=0,angle=0)
 axis(1)
 par(old.par)

  },height=600,width=900)
})
