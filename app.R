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
               mainPanel(
                 tableOutput(outputId = "contents")
               )
              )


#server
server<- function(input, output) {
  output$contents<- renderTable({ req(input$counts)
    tryCatch(
      {
        df<- read.csv(input$counts$datapath,
                      header = input$header,
                      sep= input$sep,
                      quote= input$quote)
      },
      error = function (e){
        stop(safeError(e))
      }
    )
   return(df)
    })
}

 
#
shinyApp(ui=ui, server =server)

