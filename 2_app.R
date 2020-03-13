library(shiny)
library(ggplot2)
library(reshape2)
library(datasets)

# goi_counts<- read.csv("DEseq2 means_goi_ncounts.csv")
# rownames(goi_counts)<-goi_counts$Gene
#                          )
#goi_counts


ui<- fluidPage(
  titlePanel("ARPE19 progression"),
  sidebarLayout(
    sidebarPanel(
      fileInput(inputId = "counts",
                label= "Upload data file",
                multiple = F,
                accept = NULL,
                width = NULL,
                buttonLabel = "Browse..."),
      # checkboxGroupInput(inputId = "Gene",
      #                    label = " Choose Gene of Interest(s)",
      #                    choices =  c("PPP1CB", "PPP2R5C", "BUB1B","AURKB" ,"SKA1"),
      #                    selected = "PPP1CB" ),
      selectizeInput(inputId = "Sample", 
                         label = " Choose sample(s)",
                         choices= c("ARPE19", "T53D4", "RasV12","Aktmyr" ,"MekDD"),
                         selected = "ARPE19" ),
      selectizeInput(inputId = "Gene",
                         label = " Choose Gene of Interest(s)",
                         choices =  c("PPP1CB", "PPP2R5C", "BUB1B","AURKB" ,"SKA1"),
                         selected = "PPP1CB" )
      # submitButton(text = "Apply Changes", 
      #              icon = NULL,
      #              width = NULL),
      # conditionalPanel(condition = "output.nrows")
    ),
    mainPanel(
      plotOutput("boxplot1"),
      tableOutput("contents"),
    )
  )
)


server<- function(input,output, session) {
  session$onSessionEnded(stopApp)
    output$contents<- renderTable({
      in_counts<-input$counts
      
      df<-read.csv("DEseq2 means_goi_ncounts.csv")
    })
    output$boxplot1<- renderPlot({
      req(input$counts)
    melted = melt(df)
    boxplot(melted$value ~ melted$Gene, xaxt="n")
    text(x = 1:length(levels(melted$Gene)),
         ## Move labels to just below bottom of chart.
         y = par("usr")[3] - 0.90,
         ## Use names from the data list.
         labels = levels(melted$Gene),
         ## Change the clipping region.
         xpd = NA,
         ## Rotate the labels by 35 degrees.
         srt = 35,
         ## Adjust the labels to almost 100% right-justified.
         adj = 0.965,
         ## Increase label size.
         cex = 1)
    })
    # output$ARPE19plot<- renderPlot({
    #   boxplot(as.formula(formulaText)),
    #   data= df
    # })
  }

shinyApp (ui, server)

