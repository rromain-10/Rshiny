#load shiny 
library(shiny)

#inputs and outputs
ui<- fluidPage("ARPE19 cancer cell progression through the Ras pathway",
               fileInput(inputId = "Counts", label= "Upload data file",
                         multiple = , accept = NULL, width = NULL,
                         buttonLabel = "Browse..."),
               checkboxGroupInput(inputId = "Gene", label = " Choose Gene of Interest",
                           c("PPP1CB", "PPP2R5C", "BUB1B","AURKB" ,"SKA1")),
               checkboxGroupInput(inputId = "Sample", label = " Choose sample",
                           c("ARPE19", "T53D4", "RasV12","Aktmyr" ,"MekDD")),
               submitButton(text = "Apply Changes", icon = NULL, width = NULL),
               plotOutput(outputId = "plot"),
              )


#server
server<- function(input, output) {
  output$plot<- renderPlot({input$counts}),
  read.csv(inFile$datapath, header = input$header))
 
#
shinyApp(ui=ui, server =server)

