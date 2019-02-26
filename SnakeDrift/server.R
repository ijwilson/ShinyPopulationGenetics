
library(shiny)
library(rje)    ## for colours

plotFreqs <- function(flist) {
  gens <- length(flist[[1]]$freq)

  plotpoints <- function(col_freq_pair) {
    colours <- rep(1, gens)
    colours[col_freq_pair[[2]]==0] <- 2
    colours[col_freq_pair[[2]]==1] <- 3
    points((1:gens)-0.5, col_freq_pair[[2]], col=colours, cex=2, pch=18)
    points((1:gens)-0.5, col_freq_pair[[2]], type="c",col=col_freq_pair[[1]])
  }
  
  par(mar=c(2, 2, 0, 0), mgp=c(1, 0.5, 0))      ## fix margins
  plot(NULL, xlim=c(0.2, gens), ylim=c(0,1), axes=FALSE,xlab="",ylab="")
  lapply(c(0,0.5,1), function(y) segments(0, y, gens, y, col="lightgrey"))

  lapply(flist, plotpoints)

  axis(1, at=1:gens - 0.5, labels=1:gens, tick=FALSE, cex=1.2)
  axis(2, at=c(0,0.5,1))
}

plotSnakes <- function(n, gens, colours) {
  freq <- numeric(gens)
  plot(1, type="n", ylim=c(1-0.5, n+0.5), xlim=c(0.2, gens), axes=FALSE, ylab="", xlab="")
  ybottom <- seq((1:n))-0.3
  ytop <- ybottom + 0.6
  ## initial population
  rect(0.2, ybottom, 0.8, ytop, col=colours)
  prev <- 1:n
  freq[1] <- sum(prev <= n/2)/n
  for (i in 2:gens) {
    samp <- sort(sample(1:n, replace=TRUE))
    curr <- sort(prev[samp])
    if (gens<=20)
      arrows(i-1.18, samp , i-0.82, 1:n, lwd=2, length=0.1, col=grey(0.4))  
    rect(i-0.8, ybottom, i-0.2, ytop, col=colours[curr])
    prev <- curr
    freq[i] <- sum(curr<=n/2)/n
  }
  return(freq)
}

shinyServer(function(input, output) {
  
  v <- reactiveValues()
  v$frequency_list <- NULL
  v$current_generations <- 1

  output$showSnakes <- renderPlot({
    ## get parameters
    n <- as.integer(isolate(input$n))                 
    gens <- as.integer(isolate(input$gens))
    linecolour <- as.integer(isolate(input$line_colour))
    
    input$run

    if (isolate(input$two_colours)) {
      cols <- sort(rep(cubeHelix(10)[c(3,9)], length.out=n))
    } else {
      cols <- cubeHelix(n)
    }
    ## blank plots 
    layout(matrix(c(1,2), ncol=1), heights=c(3,1))
    par(mar=c(0,2,0,0), mgp=c(0,0,0))      ## fix margins
    
    if (gens != isolate(v$current_generations) ) {
      v$frequency_list <- list()    ## empty the list
      v$current_generations <- gens
    }
    
    f <- plotSnakes(n, gens, cols)
    isolate(v$frequency_list[[length(isolate(v$frequency_list))+1]] <- list(col=linecolour, freq=f))
    plotFreqs(isolate(v$frequency_list))
  })
})
