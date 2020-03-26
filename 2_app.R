library(shiny)
library(ggplot2)
library(reshape2)
library(datasets)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
# install.packages('rsconnect')
# rsconnect::setAccountInfo(name='rromain-10', 
#                           token='708AE86DD216401F0304DE58972E2A49', 
#                           secret='B+AAnRaJVGVFzu98RU3u6LXF7cVqeSSVmk4YJqOE')
# library(rsconnect)
# getwd()
#rsconnect::deployApp('/Users/romarioromain/OneDrive - Colostate/RR_ARPE_DELUCA_COLLAB/Rshiny/2_app.R')

## Inputs and Outputs
ui<- fluidPage(
  titlePanel("ARPE19 progression data view"),
  sidebarLayout(
    sidebarPanel(
      checkboxGroupInput(inputId = "Gene", 
                         label = " Choose Gene of Interest",
                         selected=unlist(strsplit('AURKB BUB1 BUB1B NUF2 PPP1CA PPP1CB PPP1CC PPP2R2A PPP2R2A.1 PPP2R5A PPP2R5B PPP2R5C PPP2R5D PPP2R5E SGO1 SKA1 SKA2',
                                         split=' ')),
                         unlist(strsplit('AURKB BUB1 BUB1B NUF2 PPP1CA PPP1CB PPP1CC PPP2R2A PPP2R2A.1 PPP2R5A PPP2R5B PPP2R5C PPP2R5D PPP2R5E SGO1 SKA1 SKA2',
                                         split=' ')
                         )
      ),
      checkboxGroupInput(inputId = "sample", label = " Choose sample", selected=c("ARPE19"),
                         c("ARPE19", "T53D4", "RasV12","Aktmyr" ,"MekDD")
      ),
      conditionalPanel(condition = "output.nrows"),
      selectInput(inputId = "All",
                  label = "Select genes from all signficantly changing genes",
                  choices = means_of_changing_genes$X
      )
    ),
    mainPanel(
      plotOutput("linegraph"),
      tableOutput("contents"),
      plotOutput("pheatmap1"),
      plotOutput("pheatmap2"),
      plotOutput("pheatmap3")
    )
  )
)



server<- function(input,output, session) {
  session$onSessionEnded(stopApp);
  goi_ncounts<- read.csv("DEseq2_means_goi_ncounts.csv");
  changing_lrt_rdl<- read.csv("changing_lrt_rdl.csv");
  genes_of_interest_means<- read.csv("genes_of_interest_means.csv");
  means_of_changing_genes<- read.csv("means_of_changing_genes.csv");
  ncounts_goi_no_genecol<- read.csv("ncounts_goi_no_genecol.csv")
  ({
    #linegraph showing normalized counts of the TOP 20 genes of interest
    output$linegraph<- renderPlot({
      print("rendering plot")
      filtered = goi_ncounts
      melted = melt(goi_ncounts)
      colnames(melted) <- c("Gene","sample","counts")
      print(input$sample)
      filtered = melted[melted$Gene %in% input$Gene & melted$sample %in% input$sample,]
      print(filtered)
      ggplot(filtered, 
             aes(x=sample,y=counts,
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
        )) +
        ggtitle("Normalized counts of the Top 20 genes of interest")
      })
    
    #heatmap showing mean counts of significantly changing genes of interest
    output$pheatmap1<- renderPlot({
      print("rendering plot")
      filtered2 = genes_of_interest_means
      melted2 = melt(genes_of_interest_means)
      colnames(melted2) <- c("Gene","sample","counts")
      print(input$sample)
      filtered2 = melted2[melted2$Gene %in% input$Gene & melted2$sample %in% input$sample,]
      print(filtered2)
      ggplot(filtered2, 
             aes(x=sample,y=Gene,
                 group=counts,
                 col=Gene)) + 
        ggtitle ("Mean counts of significantly changing genes of interest") +
        geom_raster(aes(fill=log10(counts))) +
        scale_fill_gradient(low = "green", high = "red")
      })
    
    #heatmap showing mean counts for ALL significantly changing genes of interest
    output$pheatmap2<- renderPlot({
      print("rendering plot")
      filtered3 = means_of_changing_genes
      melted3 = melt(means_of_changing_genes)
      colnames(melted3) <- c("Gene","sample","counts")
      print(input$sample)
      filtered3 = melted3[melted3$Gene %in% input$All & melted3$sample %in% input$sample,]
      print(filtered3)
      ggplot(filtered3, 
             aes(x=sample,y=Gene,
                 group=counts,
                 col=Gene)) + 
        geom_raster(aes(fill=log(counts))) +
        ggtitle ("Means counts for all significantly changing genes")
      })
      
      #heatmap showing mean counts for ALL significantly changing genes
      output$pheatmap2<- renderPlot({
        print("rendering plot")
        filtered3 = means_of_changing_genes
        melted3 = melt(means_of_changing_genes)
        colnames(melted3) <- c("Gene","sample","counts")
        print(input$sample)
        filtered3 = melted3[melted3$Gene %in% input$All & melted3$sample %in% input$sample,]
        print(filtered3)
        ggplot(filtered3, 
               aes(x=sample,y=Gene,
                   group=counts,
                   col=Gene)) + 
          geom_raster(aes(fill=log(counts))) +
          ggtitle ("Mean counts for all significantly changing genes")
        })
      })
    }


shinyApp(ui=ui, server=server)

