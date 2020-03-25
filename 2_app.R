library(shiny)
library(ggplot2)
library(reshape2)
library(datasets)
library(dplyr)

# goi_counts<- read.csv("DEseq2 means_goi_ncounts.csv")
# rownames(goi_counts)<-goi_counts$Gene
#                          )
#goi_counts


#print(str(df))

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
      checkboxGroupInput(inputId = "Gene", 
                         label = " Choose Gene of Interest",
                         selected=unlist(strsplit('AURKB BUB1 BUB1B NUF2 PPP1CA PPP1CB PPP1CC PPP2R2A PPP2R2A.1 PPP2R5A PPP2R5B PPP2R5C PPP2R5D PPP2R5E SGO1 SKA1 SKA2',
                                         split=' ')),
                         unlist(strsplit('AURKB BUB1 BUB1B NUF2 PPP1CA PPP1CB PPP1CC PPP2R2A PPP2R2A.1 PPP2R5A PPP2R5B PPP2R5C PPP2R5D PPP2R5E SGO1 SKA1 SKA2',
                                         split=' ')
                         )
      ),
      checkboxGroupInput(inputId = "sample", label = " Choose sample", selected=c("ARPE19", "T53D4", "RasV12","Aktmyr" ,"MekDD"),
                         c("ARPE19", "T53D4", "RasV12","Aktmyr" ,"MekDD")
      ),

      conditionalPanel(condition = "output.nrows")
    ),
    mainPanel(
      plotOutput("boxplot1"),
      tableOutput("contents")
    )
  )
)


server<- function(input,output, session) {
  session$onSessionEnded(stopApp);
  goi_ncounts<- read.csv("DEseq2_means_goi_ncounts.csv");
  ({
    # output$contents<- renderTable({
    #    in_counts<-input$counts 
    #   
    # )
    # });
    output$boxplot1<- renderPlot({
      print("rendering plot")
      filtered = goi_ncounts
      melted = melt(goi_ncounts)
      colnames(melted) <- c("Gene","sample","counts")
      print(input$sample)
      filtered = melted[melted$Gene %in% input$Gene & melted$sample %in% input$sample,]
      print(filtered)
      ggplot(filtered, 
             aes(x=sample,y=log(counts),
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

