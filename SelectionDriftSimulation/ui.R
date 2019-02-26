library(shiny)
# Define UI for application that plots simulations
shinyUI(fluidPage(
  # Application title
  # headerPanel("Drift Simulation"),
  titlePanel("Selection and Drift"),
  ## First row
  fluidRow(
    column(4,
           selectInput("selection", "Type of Selection", 
                       choices=list("None (pure drift)"="none", 
                                    "Heterozygote Advantage  (1, 1+s, 1)"="het_input",
                                    "Positive  (1, s, 1+2s)"="positive",
                                    "Homozygote Recessive  (1, 1, 1-s)"="homo_recess")
                       , selected="none"),
            uiOutput("ui")  ## the dynamic ui component
    ),
    column(4,
           selectInput(inputId = "N", label="Population Size",
                       choices = list(10,20,50,100,200,500, 1000, 5000)
                       , selected=1000),
           selectInput("startingP",
                       "Starting variant frequency",
                       choices=list("One copy"=1, "1%"=2, "5%"=3, "50%"=4)
                       , selected=4)
    ),
   column(4,
          numericInput("gens", "Number of Generations", value=100 , min=100, max=1000, step=100),
          numericInput("pops", "Replicates", value=100, min=1, max=1000, step=100),
          actionButton("run","Go!")
  )),
  plotOutput("distPlot")
))
