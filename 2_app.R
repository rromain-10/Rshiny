library(shiny)
library(ggplot2)
library(reshape2)
library(datasets)
library(dplyr)

# goi_counts<- read.csv("DEseq2 means_goi_ncounts.csv")
# rownames(goi_counts)<-goi_counts$Gene
#                          )
#goi_counts

df<- read.csv("DEseq2 means_goi_ncounts.csv")
print(str(df))

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
      checkboxGroupInput(inputId = "gene", 
                         label = " Choose Gene of Interest",
                         unlist(strsplit('AURKB BUB1 BUB1B NUF2 PPP1CA PPP1CB PPP1CC PPP2R2A PPP2R2A.1 PPP2R5A PPP2R5B PPP2R5C PPP2R5D PPP2R5E SGO1 SKA1 SKA2',
                                         split=' ')
                         )
      ),
      checkboxGroupInput(inputId = "sample", label = " Choose sample",
                         c("ARPE19", "T53D4", "RasV12","Aktmyr" ,"MekDD")
      ),
      submitButton(text = "Apply Changes",
                   icon = NULL,
                   width = NULL),
      conditionalPanel(condition = "output.nrows")
    ),
    mainPanel(
      plotOutput("boxplot1"),
      tableOutput("contents")
    )
  )
)


server<- function(input,output, session) {
  session$onSessionEnded(stopApp)
  ({
    output$contents<- renderTable({
      in_counts<-input$counts 
      df<-read.csv(input$counts$datapath
    )
    })
    output$boxplot1<- renderPlot({
      filtered <-
        df %>%
        filter(Gene == df$Gene,
                 ARPE19 == input$ARPE19,
                 T53D4 == input$T53D4,
                 RasV12 == input$RasV12,
                 MekDD == input$MekDD,
                 Aktmyr == input$Aktmyr
        )
      ggplot(filtered, aes(x= Gene, y= sample, counts))+
        geom_boxplot()
    })
    # output$boxplot1<- renderPlot({
    #   req(input$counts)
    # melted = melt(df)
    # boxplot(melted$value ~ melted$Gene, xaxt="n")
    # text(x = 1:length(levels(melted$Gene)),
    #      ## Move labels to just below bottom of chart.
    #      y = par("usr")[3] - 0.90,
    #      ## Use names from the data list.
    #      labels = levels(melted$Gene),
    #      ## Change the clipping region.
    #      xpd = NA,
    #      ## Rotate the labels by 35 degrees.
    #      srt = 35,
    #      ## Adjust the labels to almost 100% right-justified.
    #      adj = 0.965,
    #      ## Increase label size.
    #      cex = 1)
    # })
    # output$ARPE19plot<- renderPlot({
    #   boxplot(as.formula(formulaText)),
    #   data= df
    # })
     })
    }

shinyApp (ui=ui, server=server)

