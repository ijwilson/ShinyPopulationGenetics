
library(shiny)
# Define UI for application that plots simulations
shinyUI(fluidPage(
  # Application title
  # headerPanel("Drift Simulation"),
  titlePanel("Bottlenecks and Selection"),
  ## First row
  fluidRow(
    column(4,
           selectInput("selection", "Type of Selection", 
                       choices=list("None (pure drift)"="none", 
                                    "Heterozygote Advantage  (1, 1+s, 1)"="het_input",
                                    "Positive  (1, s, 1+2s)"="positive",
                                    "Homozygote Recessive  (1, 1, 1-s)"="homo_recess")
                       , selected="none"),
           selectInput(inputId = "N", label="Standard Population Size",
                       choices = list(100, 200, 500, 1000, 5000)
                       , selected=1000),

            uiOutput("ui")  ## the dynamic ui component
    ),
    column(4,

           numericInput("initial", "Before Bottleneck  (generations)", value=20 , min=1, max=100, step=10),
           
           numericInput(inputId = "duration", "Bottleneck Duration",
                        value=10, min=2, max=100),
           numericInput(inputId = "fraction", 
                        label="bottleneck Severity %",
                        value=20, min=1, max=100)
    ),
   column(4,
          selectInput("startingP",
                      "Starting variant frequency",
                      choices=list("One copy"=1, "1%"=2, "5%"=3, "50%"=4)
                      , selected=4),
          numericInput("pops", "Replicates", value=10, min=1, max=1000, step=2),
          actionButton("run","Go!")
  )),
  plotOutput("distPlot")
))
