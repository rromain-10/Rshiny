#load shiny 
library(shiny)
library(reshape2)

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
               mainPanel(
                 plotOutput(outputId = "boxplot")
               )
              )


#server
server<- function(input, output) {
  output$boxplot<- renderPlot({ 
    req(input$Counts)
    tryCatch(
      {
        df<- read.csv(input$Counts$datapath
                      )
      },
      error = function (e){
        stop(safeError(e))
      }
    )
   #return(df)
    melted = melt(df)
    boxplot(melted$value ~ melted$Gene)
    #boxplot(log(df[[2]]))
    })
}

 
#
shinyApp(ui=ui, server =server)

