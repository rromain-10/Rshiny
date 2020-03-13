#load shiny 
library(shiny)
library(reshape2)

#inputs and outputs
ui<- fluidPage("ARPE19 cancer cell progression through the Ras pathway",
               fileInput(inputId = "counts", label= "Upload data file",
                         multiple = , accept = NULL, width = NULL,
                         buttonLabel = "Browse..."),
               checkboxGroupInput(inputId = "gene", 
                                  label = " Choose Gene of Interest",
                                  unlist(strsplit('AURKB BUB1 BUB1B NUF2 PPP1CA PPP1CB PPP1CC PPP2R2A PPP2R2A.1 PPP2R5A PPP2R5B PPP2R5C PPP2R5D PPP2R5E SGO1 SKA1 SKA2',
                                                   split=' ')
                                         )
                                  ),
               checkboxGroupInput(inputId = "sample", label = " Choose sample",
                           c("ARPE19", "T53D4", "RasV12","Aktmyr" ,"MekDD")),
               submitButton(text = "Apply Changes", icon = NULL, width = NULL),
               mainPanel(
                 plotOutput(outputId = "boxplot")
               )
              )


#server
server<- function(input, output) {
  output$boxplot<- renderPlot({ 
    req(input$counts)
    tryCatch(
      {
        df<- read.csv("DEseq2 means_goi_ncounts.csv")
      },
      error = function (e){
        stop(safeError(e))
      }
    )
    # the selected genes are: in input$gene
    # figure out which are selected first
    melted_df = melt(df[df$Gene == input$gene,])
    ggplot(melted_df, 
           aes(x=variable,y=log(value),
               group=Gene,
               col=Gene, 
               linetype=Gene)) + 
      geom_line(size = 1.5)  +
      scale_linetype_manual(values = c(rep("solid", 10), 
                                       rep("dotted", 10)
                                       )
                            ) +
      scale_color_manual(values = c(brewer.pal(10,"Spectral"),
                                    brewer.pal(10,"Spectral")
                                    )
                         ) +
      theme_bw()
    })
}

 
#
shinyApp(ui=ui, server =server)

