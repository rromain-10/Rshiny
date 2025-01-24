## ----include=F--------------------------------------------------------------------------------------
knitr::opts_chunk$set (echo=F)


## ----packages---------------------------------------------------------------------------------------
#install.packages('devtools')
##install.packages('viridisLite')

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(version = "3.11")
# BiocManager::install("DESeq2")
# BiocManager::install("biomaRt")
# BiocManager::install("ensembldb")

library(BiocManager)
options(repos = BiocManager::repositories()) # to locate XVector and who knows what else
library(rmarkdown)
library(rsconnect)
library(biomaRt)
library(devtools)
library(shiny)
library(ggplot2)
library(reshape2)
library(datasets)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(DESeq2)
library(shinydashboard)
library(viridis)
library(tidyr)
library(DT)
library(plotly)
library(heatmaply)

  ##files needed
  gene_divisions<- read.csv('20200923_gene_divisions_means.csv');
  goi_ncounts<- read.csv("DEseq2_means_goi_ncounts.csv");
  changing_lrt_rdl<- read.csv("20200923_changing_lrt_rdl.csv");
  genes_of_interest_means<- read.csv("20200923_genes_of_interest_means.csv");
  means_of_changing_genes<- read.csv("20200923_means_of_changing_genes.csv");
  ncounts_goi_no_genecol<- read.csv("ncounts_goi_no_genecol.csv");
  normalized_genecounts<- read.csv("20200923_normalized_genecounts.csv");
  dds<- readRDS("dds.RDS")
  dds$sample <- factor(dds$sample, levels=c("ARPE19", "T53D4", "RasV12", "Aktmyr", "MekDD"))
  dds$sample
  cts<- readRDS("cts.RDS")
  #res_T53D4_vs_MekDD<- readRDS('res/res_T53D4_vs_MekDD.RDS')
  res_ARPE19_vs_T53D4<-readRDS ('res/res_ARPE19_vs_T53D4')
  res_ARPE19_vs_RasV12<- readRDS ('res/res_ARPE19_vs_RasV12')
  res_ARPE19_vs_MekDD<- readRDS('res/res_ARPE19_vs_MekDD')
  res_ARPE19_vs_Aktmyr<- readRDS('res/res_ARPE19_vs_Aktmyr')
  selectedGenes = c('BUB1B','AURKB','BUB1','SKA1','SKA2','SKA3','SGO1','PPP2R5C',
                    'PPP2R5B','PPP2R2A','PPP2R2A','PPP2R5D','PPP2R5A','PPP2R5E',
                    'PPP1CA','PPP1CC','PPP1CB','HEC1','NUF2','SPC24','SPC25')
  #Gene associated with MEK pathway
  Mek_genes=c('MAP2K1', 'MAP2K2','MAP2K3','MAP2K4','MAP2K5','MAP2K6','MAP2K7',
              'PTPN11', 'GBR2','SOS1')
  
  #Gens associated with the RAS pathway
  Ras_genes=c('KRAS','NRAS','BRAF','HRAS','RAF1')
  
  #Gene associated with AKT pathway
  Akt_genes =c('AKT1','AKT2','AKT3','AKT8','PDK1','PDK2','MTOR')
  
  #Genes assocaited with T53D4
  T53D4_genes=c('CDK4', 'TP53')
  
  means_of_changing_genes2<- read.csv('means_of_changing_genes2.csv', header = T)
      row.names(means_of_changing_genes2)<- means_of_changing_genes2$X
      counts<-means_of_changing_genes2[,2:length(means_of_changing_genes2[1,])]
      counts_matrix<-data.matrix(counts)

bigheatmap_img = base64enc::dataURI(file="BigHeatmap.png", mime="image/png")

## ----Inputs and outputs-----------------------------------------------------------------------------
ui<- fluidPage(
  titlePanel("ARPE (human retinal cells) progression data view"),
  dashboardBody( "
                 
                 This app shows RNA sequencing data analysis for the Deluca lab's
                 ARPE19 cell line progression model.
                 
                 "
  ),
  
    dashboardBody(
      
      "______________________________________________________________________________________________________________________________________________________________________________________________________",
      
    "MAplots displayed below shows significant up or down regulation of expressed genes, where the comparsions are between plot titles listed (lfc=0.01)  "
  ),
  
  fluidRow(
    
      column(3, 
             plotOutput('MAplot1')
            ),
      
      column (3,
              plotOutput('MAplot2')
              ),
     
     column(3,
            plotOutput('MAplot3')
            ),
     
     column(3,
            plotOutput('MAplot4')
            )
  ),
  
  dashboardBody(
    
  "__________________________________________________________________________________________________________________________________________________________________________________________________________",
    
    "The figure below displays nomralized counts in boxplot format showing interquartiles ranges, outliers, median and mean for the gene of choice across the 5 sample progressive cel line.  "
  ),
  
  #boxplot for significantly changing genes 
  dashboardBody(
      sidebarLayout(
      sidebarPanel(
        
        conditionalPanel(condition = "output.nrows"),
        selectInput(inputId = "All",
                    label = "Select genes from all signficantly changing genes",
                    choices = means_of_changing_genes$X
                   )
        ),
  mainPanel(
    plotOutput('boxplot1')
      )
    )
  ),
  
   dashboardBody(
     
       "_____________________________________________________________________________________________________________________________________________________________________________________________________",
    
    "The figure below displays a heatmap shwoing gene expression for ONLY genes of interest that are significantly changing (lfc=0.01 & padj<0.1).  "
    ),  
  #heatmap for Genes of interest
  dashboardBody(
    sidebarLayout(
      sidebarPanel(
  
        #heatmap for Genes of interest
        checkboxGroupInput(inputId = "Gene",
                           label = " Choose Gene of Interest",
                           selected=unlist(strsplit('AURKB NUF2 PPP1CB PPP1CC PPP2R5B PPP2R5D SGO1 SKA1',
                                                    split=' ')),
                           unlist(strsplit('AURKB NUF2 PPP1CB PPP1CC PPP2R5B PPP2R5D SGO1 SKA1',
                                           split=' ')
                           )
        ),
        
        checkboxGroupInput(inputId = "sample", label = " Choose sample", selected=c("ARPE19", "T53D4", "RasV12","Aktmyr" ,"MekDD"),
                           c("ARPE19", "T53D4", "RasV12","Aktmyr" ,"MekDD")
                          )
        ),
      mainPanel(
        plotOutput('pheatmap1')
            )
    )
   ),
  
  dashboardBody(
    
    "________________________________________________________________________________________________________________________________________________________________________________________________________",
    
    "The figure below displays a heatmap shwoing gene expression for ALL genes that are significantly changing in the entire dataset (lfc=0.01 & padj<0.1).  "
    ),   

  #heatmap for Genes of interest
  dashboardBody(
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput( inputId = 'none',
                            label = 'Heatmap'
                          )
        ),
      mainPanel(
        img(src=bigheatmap_img, width="732px", heigt="497px")
      )
    )
   ),
  
  dashboardBody(
    
    "________________________________________________________________________________________________________________________________________________________________________________________________________",
    
    "The data tables shows the cluster each gene is associated with from the heatmap above. Cluster 1 is the first cluster of the heatmap and counts downwards to 7 "
    ),   

  dashboardBody(
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput( inputId = 'none',
                            label = 'Type gene name into search bar to see associated cluster'
                          )
        ),
      mainPanel(DT::dataTableOutput('gene_div'),
                            label = 'Interactive Heatmap'
                          )
        ),
      mainPanel(
        plotOutput('whole_heatmaply')
            )
   ),
  
  #boxplot for significantly changing genes 
  dashboardBody(
    
      "______________________________________________________________________________________________________________________________________________________________________________________________________",
    
    sidebarLayout(
      sidebarPanel(
                    
        conditionalPanel(condition = "output.nrows"),
        selectInput(inputId = "sig_chan",
                    label = "Select genes from all signficantly changing genes",
                    choices = changing_lrt_rdl$X
          ),
        
     checkboxGroupInput(inputId = "sample1", label = " Choose sample", selected=c("ARPE19", "T53D4", "RasV12","Aktmyr" ,"MekDD"),
                           c("ARPE19", "T53D4", "RasV12","Aktmyr" ,"MekDD")
                          )
        ),
  mainPanel(
    plotOutput('heatmap2')
          )
       )
     )
   )
   



## ----figures and displays---------------------------------------------------------------------------
server<- function(input,output, session) {
  session$onSessionEnded(stopApp);
   ({
  
    #MAplots ARPE19 vs T53D4 
    output$MAplot1 <- renderPlot ({
      print("rendering plot: MAplot ARPE19 vs T53D4")
      #ARPE19 vs T53D4 plots##
      par(mfrow=c(1,1))
      #print(input$sample1)
      plotMA(res_ARPE19_vs_T53D4, main="ARPE19 vs T53D4", ylim = c(-20,20),
             ylab = "log fold change (ratio of normalized ARPE19 / T53D4)",
             xlab = "means of normalized counts",
             cex.main= 2,
             cex= 0.5,
             alpha=0.5)
      # with(res_ARPE19_vs_T53D4 [selectedGenes, ],
      #      { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
      #               lwd=2)
      #        text(baseMean,log2FoldChange,
      #             selectedGenes, pos=2, col="dodgerblue")})
      # with(res_ARPE19_vs_T53D4 [T53D4_genes, ],
      #      { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
      #               lwd=2)
      #        text(baseMean,log2FoldChange, T53D4_genes, pos=2, col="dodgerblue")})
    })
       
    #MAplots ARPE19 vs RasV12    
    output$MAplot2 <- renderPlot ({
        print("rendering plot: MAplot ARPE19 vs RasV12")
         ##ARPE19 vs RasV12 plots##
         par(mfrow=c(1,1))
         #print(input$sample1)
         plotMA(res_ARPE19_vs_RasV12, main="ARPE19 vs RasV12", ylim = c(-20,20),
                ylab = "log fold change (ratio of normalized ARPE19 / RasV12)",
                xlab = "means of normalized counts",
                cex.main= 2,
                cex= 0.5,
                alpha= 0.5)
         # with(res_ARPE19_vs_RasV12 [selectedGenes, ],
         #      { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
         #               lwd=2)
         #        text(baseMean,log2FoldChange,
         #             selectedGenes, pos=2, col="dodgerblue")})
         # with(res_ARPE19_vs_RasV12 [Ras_genes, ],
         #      { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
         #               lwd=2)
         #        text(baseMean,log2FoldChange,
         #             Ras_genes, pos=2, col="dodgerblue")})
       })
       
    #MAplots ARPE19 vs MekDD 
    output$MAplot3 <- renderPlot ({
       print("rendering plot: MAplot ARPE19 vs MekDD")
       ##ARPE19 vs MekDD plots##
       par(mfrow=c(1,1))
       #print(input$sample1)
       plotMA(res_ARPE19_vs_MekDD, main="ARPE19 vs MekDD", ylim = c(-20,20),
              ylab = "log fold change (ratio of normalized ARPE19 / MekDD)",
              xlab = "means of normalized counts",
              cex.main= 2,
              cex= 0.5,
              alpha= 0.5)
       # with(res_ARPE19_vs_MekDD [selectedGenes, ],
       #      { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
       #               lwd=2)
       #        text(baseMean,log2FoldChange,
       #             selectedGenes, pos=2, col="dodgerblue")})
       # with(res_ARPE19_vs_MekDD [Mek_genes, ],
       #      { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
       #               lwd=2)
       #        text(baseMean,log2FoldChange, Mek_genes, pos=2, col="dodgerblue")})
       
       })
      
    #MAplots ARPE19 vs Aktmyr  
    output$MAplot4 <- renderPlot ({
      print("rendering plot: MAplot ARPE19 vs Akymyr")
       ##ARPE19 vs Aktmyr plots##
       par(mfrow=c(1,1))
       #print(input$sample1)
       plotMA(res_ARPE19_vs_Aktmyr, main="ARPE19 vs Aktmyr ", ylim = c(-20,20),
              ylab = "log fold change (ratio of normalized ARPE19 / Aktmyr)",
              xlab = "means of normalized counts",
              cex.main= 2,
              cex= 0.5,
              alpha= 0.5)
       # with(res_ARPE19_vs_Aktmyr [selectedGenes, ],
       #      { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
       #               lwd=2)
       #        text(baseMean,log2FoldChange, selectedGenes, pos=2, col="dodgerblue")})
       # with(res_ARPE19_vs_Aktmyr [Akt_genes, ],
       #      { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
       #               lwd=2)
       #        text(baseMean,log2FoldChange, Akt_genes, pos=2, col="dodgerblue")})


    })
    
        #heatmap showing mean counts for ALL significantly changing genes of interest
    output$boxplot1<- renderPlot({
      print("rendering plot")
      print(input$All)
      #Reorder samples to match experimental design
      dds$sample <- factor(dds$sample, levels=c("ARPE19", "T53D4", "RasV12", "Aktmyr", "MekDD"))
      dds$sample
      #AKT1
      plotCounts(dds, gene=print(input$All) ,intgroup = c("sample"),cex.main=2, cex=1.5, xlab= "Sample", returnData = T)
      pc<-plotCounts(dds, gene=print(input$All), intgroup =c("sample"), cex.main=2, cex=1.5, xlab= "Sample", returnData=T)
      ggplot (pc,aes(x=sample,y=count)) +
        geom_boxplot(aes(fill=sample, alpha=0.9))+ stat_summary( geom="point", shape=20, size=5,color="red", fill="red") +
        geom_jitter()+scale_y_log10() + 
        ggtitle("Gene expression")
      
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
                 col=Gene,
                 fill =counts)) + 
        ggtitle ("Mean counts of significantly changing genes of interest") +
        geom_tile() +
        scale_fill_gradient(low = 'Blue',high = 'Red')+
        theme_minimal()
      })
    
    # output$whole_heatmap<- renderPlot ({
    #   print("rendering pheatmap")
    #   pheatmap(counts_matrix,
    #                      scale="row",
    #                      color = colorRampPalette(c("blue", "white", "red"), space = "Lab")(300),
    #                      border_color = TRUE,
    #                      treeheight_row = 100,
    #                      cluster_rows=TRUE,
    #                      cluster_cols=F,
    #                      cutree_rows = 7,
    #                      cutree_cols = 5,
    #                      clustering_distance_rows = "euclidean",
    #                      clustering_method = "complete",
    #                      show_rownames = T)
    # 
    #   
    # 
    #   })
    
      
    output$gene_div = DT::renderDataTable ({
      print("rendering gene_divs")
      colnames(gene_divisions)<- c('Gene', 'Cluster')
      gene_divisions
      
    })
      
    # output$whole_heatmaply<- renderPlot ({
    #   print("rendering heatmaply")
    #   # heatmaply(counts_matrix,
    #   #           k_col = 5,
    #   #           k_row = 7,
    #   #           scale = 'row')
    # })
    
    #Heatmap
     output$heatmap2<- renderPlot ({
      print("rendering plot")
      filtered3 = counts_matrix
      melted3 = melt(counts_matrix)
      colnames(melted3) <- c("Gene","sample","counts")
      print(input$sample1)
      print(input$sig_chan)
      filtered3 = melted3[melted3$Gene %in% input$sig_chan & melted3$sample %in% input$sample1,]
      print(filtered3)
      ggplot(filtered3, 
             aes(x=sample,y=Gene,
                 group=Gene,
                 col=sample,
                 fill =counts)) + 
        ggtitle ("Mean counts of significantly changing genes of interest") +
        geom_tile() +
        scale_fill_viridis(na.value = "transparent") +
        theme(legend.position = "bottom") +
        theme_minimal()+
        scale_fill_gradient(low = 'Blue',high = 'Red')
     })
    })
  }
   
  
     


## ----App--------------------------------------------------------------------------------------------
shinyApp(ui=ui, server=server)


## ---------------------------------------------------------------------------------------------------
# library(rsconnect)
# rsconnect::setAccountInfo(name='onishlab',
#                           token='C1E0A73570EF891C8ACF5D1C00144C1C',
#                           secret='00jn0wfb//WoESIvtONyd/w3qqaHr1la8lb+b2Wa')
# 
# 
# getwd()
# rsconnect::deployApp("/Users/romarioromain/OneDrive - Colostate/RR_ARPE_DELUCA_COLLAB/Rshiny", appTitle = 'ARPE19 cancer cell progression app', account = 'onishlab', appPrimaryDoc = "Deulca_collab_app.R")
