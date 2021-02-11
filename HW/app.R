#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

hwe_geno<-function(mm, mn, nn)
{
    obs <- c(mm, mn, nn)
    n <- sum(obs)
    f <- (2*mm+mn)/(2*n)
    expected <- c(n*f^2, 2*n*f*(1-f), n*(1-f)^2)
    xx <- sum((obs-expected)*(obs-expected)/expected)
    p <- pchisq(xx, df=1, lower.tail=FALSE)
    
    expected_d <- setNames(as.data.frame(rbind(expected, 100*expected/n), row.names=c("genotypes", "percentage")), c("MM", "MN", "NN"))

    
    al_dist <- setNames(as.data.frame(rbind(c(2*mm+mn, mn+2*nn), 100*c(f, 1-f)), row.names=c("genotypes", "percentage")), c("M", "N")) 
  
    list(
        allele.dist = al_dist,
        observed.dist = data.frame(MM=c(mm, 100*mm/n), MN = c(mn, 100*mn/n), NN=c(nn, 100*nn/n)),
        expected.dist = expected_d,
        chi.value=xx,
        p.value=p,
        df=1
    )
}

library(shiny)

# Define UI for application that draws a histogram

ui <- fluidPage(
    
    # Application title
    titlePanel("Chi-square test of Hardy-Weinberg Equilibrium "),
    
    sidebarLayout(sidebarPanel(
        h4("Observed genotype distribution"),
        numericInput(
            inputId = "mm",
            label = "MM:",
            value = NA
        ),
        numericInput(
            inputId = "mn",
            label = "MN:",
            value = NA
        ),
        numericInput(
            inputId = "nn",
            label = "NN:",
            value = NA
        )
    ),
    
        
        mainPanel(
            tabsetPanel(
                tabPanel("Summary",
                         h3(textOutput("obs.dist", container = span)),
                         tableOutput("obs.tbl"),
                         
                         h3(textOutput("exp.dist", container = span)),
                         tableOutput("exp.tbl"),
                         
                         h3(textOutput("allele.dist", container = span)),
                         htmlOutput("allele.tbl", container = span),
                         value=1),
                
                tabPanel("Calculations",
                         h3(textOutput("chi", container = span)),
                         htmlOutput("chi.val", container = span),
                         h3(textOutput("p", container = span)),
                         htmlOutput("p.val", container = span),
                         value=2),
                
                selected= 2, type = "tabs")
        )
    )
)



server <- function(input, output) {
    cale <- reactive({
        as.numeric(input$ale)
    })
    
    dat <- reactive({
        df <- data.frame(
            lbls = c("MM", "MN", "NN"),
            value = rbind(
                input$mm,
                input$mn,
                input$nn
            ),
            stringsAsFactors = FALSE
        )
    #    print(df)
       # df
    })
    
    cmm <- reactive({
        as.numeric(input$mm)
    })
    
    cmn <- reactive({
        as.numeric(input$mn)
    })
    
    cnn <- reactive({
        as.numeric(input$nn)
    })

    
    hwe_p <-
        function()
            ({
                hwe_geno(cmm(), cmn(), cnn())
            })
    
    output$allele.tbl <- renderTable({
        hwe_p()$allele.dist
    }, rownames = TRUE)
    
    output$obs.tbl <- renderTable({
        hwe_p()$observed.dist
    }, rownames = TRUE)
    
    output$exp.tbl <- renderTable({
        hwe_p()$expected.dist
    }, rownames = TRUE)
    
    output$chi.val <- renderTable({
        hwe_p()$chi.value
    })
    
    output$p.val <- renderTable({
        hwe_p()$p.value
    })
    
    output$allele.dist <- renderText({"Allele distribution"})
    output$obs.dist <- renderText({"Observed distribution" })
    
    output$exp.dist <- renderText({"Expected distribution"    })
    
    output$chi <- renderText({"Chi square value"    })
    
    output$p <- renderText({        "P value"    })
    
}

# Run the application
shinyApp(ui, server)


