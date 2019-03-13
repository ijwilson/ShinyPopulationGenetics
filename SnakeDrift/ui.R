
library(shiny)

shinyUI(
  fluidPage(
  ## Application title
    title = "Snake Drift",
    fluidRow(
      h1("Snake Drift"),
      column(4,
      sliderInput("n",
                    "Population Size:",
                    min = 5,
                    max = 40,
                    value = 10) 
      ),
      column(4,
        sliderInput("gens",
                    "Generations:",
                    min = 5,
                    max = 40,
                    value = 10) 
      ),
      column(2,
             selectInput("line_colour", "Line Colour", 
                         choices=list("Black"=1,  "Red"=2, "Blue"=4, "Green"=3, "Cyan"=5, "Magenta"=6)
                         , selected=1)
      ),
      column(2,
             checkboxInput("two_colours", "Two Colours", value = FALSE),
             actionButton("run", "Go!")#,
      )
    ),
    plotOutput("showSnakes", height="500px")
  )
)

