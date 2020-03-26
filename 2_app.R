library(shiny)
library(ggplot2)
library(reshape2)
library(datasets)
library(dplyr)
library(RColorBrewer)
library(d3heatmap)

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
      tableOutput("contents"),
      plotOutput("d3heatmap")
    )
  )
)


server<- function(input,output, session) {
  session$onSessionEnded(stopApp);
  goi_ncounts<- read.csv("DEseq2_means_goi_ncounts.csv");
  changing_lrt_rdl<- read.csv("changing_lrt_rdl.csv");
  genes_of_interest_means<- read.csv("genes_of_interest_means.csv")
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
        ))
      })
    output$d3heatmap<- renderPlot({
      print("rendering plot")
      filtered2 = genes_of_interest_means
      melted2 = melt(genes_of_interest_means)
      colnames(melted2) <- c("Gene","sample","counts")
      print(input$sample)
      filtered2 = melted2[melted2$Gene %in% input$Gene & melted$sample %in% input$sample,]
      print(filtered2)
      ggplot(filtered2, 
             aes(x=sample,y=Gene,
                 group=counts,
                 col=Gene, 
                 # linetype=Gene
                 )) + 
        ggtitle ("Means counts of genes of interest") +
        geom_tile(aes(fill=counts))
      })
    
     })
    }

shinyApp (ui=ui, server=server)

