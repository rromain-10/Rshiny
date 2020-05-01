install.packages('shinydashboard')
# install.packages('rsconnect')
# rsconnect::setAccountInfo(name='rromain-10', 
#                           token='708AE86DD216401F0304DE58972E2A49', 
#                           secret='B+AAnRaJVGVFzu98RU3u6LXF7cVqeSSVmk4YJqOE')
# library(rsconnect)
# getwd()
#rsconnect::deployApp('/Users/romarioromain/OneDrive - Colostate/RR_ARPE_DELUCA_COLLAB/Rshiny/2_app.R')


library(shiny)
library(ggplot2)
library(reshape2)
library(datasets)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(DESeq2)
library(shinydashboard)


## Inputs and Outputs
ui<- fluidPage(
  titlePanel("ARPE19 progression data view"),
  dashboardBody(
    sidebarLayout(
    sidebarPanel(
      
      # selectInput(inputId = "sample1",
      #             label = " Choose sample",
      #             choices=c("ARPE19vsT53D4"="res_ARPE19_vs_T53D4","ARPE19vsRasV12"="res_ARPE19_vs_RasV12", "ARPE19vsAktmyr"="res_ARPE19_vs_Aktmyr","ARPE19vsMekDD"="res_ARPE19_vs_MekDD"
      #                      )
      #             ),

      checkboxGroupInput(inputId = "Gene", 
                         label = " Choose Gene of Interest",
                         selected=unlist(strsplit('AURKB BUB1 BUB1B NUF2 PPP1CA PPP1CB PPP1CC PPP2R2A PPP2R2A.1 PPP2R5A PPP2R5B PPP2R5C PPP2R5D PPP2R5E SGO1 SKA1 SKA2',
                                                  split=' ')),
                         unlist(strsplit('AURKB BUB1 BUB1B NUF2 PPP1CA PPP1CB PPP1CC PPP2R2A PPP2R2A.1 PPP2R5A PPP2R5B PPP2R5C PPP2R5D PPP2R5E SGO1 SKA1 SKA2',
                                         split=' ')
                         )
      ),
      checkboxGroupInput(inputId = "sample", label = " Choose sample", selected=c("ARPE19","RasV12"),
                         c("ARPE19", "T53D4", "RasV12","Aktmyr" ,"MekDD")
      ),
      conditionalPanel(condition = "output.nrows"),
      selectInput(inputId = "All",
                  label = "Select genes from all signficantly changing genes",
                  choices = means_of_changing_genes$X
      )
      # selectInput(inputId = "reps",
      #             label = "Select genes from all signficantly changing genes",
      #             choices = changing_lrt_rdl
      #             )
    ),
    mainPanel(
      plotOutput("MAplot"),
      plotOutput("plotcounts"),
      tableOutput("contents"),
      plotOutput("pheatmap1"),
      plotOutput("pheatmap2"),
      plotOutput("pheatmap3"),
      plotOutput("boxplot1")
      )
    )
  )
)




server<- function(input,output, session) {
  session$onSessionEnded(stopApp);
  ##files needed
  goi_ncounts<- read.csv("DEseq2_means_goi_ncounts.csv");
  changing_lrt_rdl<- read.csv("changing_lrt_rdl.csv");
  genes_of_interest_means<- read.csv("genes_of_interest_means.csv");
  means_of_changing_genes<- read.csv("means_of_changing_genes.csv");
  ncounts_goi_no_genecol<- read.csv("ncounts_goi_no_genecol.csv");
  normalized_genecounts<- read.csv("normalized_genecounts.csv")
  dds<- readRDS("dds.RDS")
  cts<- readRDS("cts.RDS")
  res_T53D4_vs_MekDD<- readRDS('res/res_T53D4_vs_MekDD.RDS')
  res_ARPE19_vs_T53D4<-readRDS ('res/res_ARPE19_vs_T53D4')
  res_ARPE19_vs_Rasv12<- readRDS ('res/res_ARPE19_vs_Rasv12')
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
  
  ({
    
    #MAplot 
    output$MAplot <- renderPlot ({
      print("rendering plot")
      par(mfrow=c(1,1))
      #print(input$sample1)
      plotMA(res_T53D4_vs_MekDD, main="T53D4 vs MekDD", ylim = c(-20,20),
             ylab = "log fold change (ratio of normalized T53D4 / MekDD)",
             xlab = "means of normalized counts",
             cex.main= 3,
             cex= 0.5,
             alpha= 0.5)
       with(res_T53D4_vs_MekDD [selectedGenes, ],
            { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
                     lwd=2)
              text(baseMean,log2FoldChange,
                   selectedGenes, pos=2, col="dodgerblue")})
       with(res_T53D4_vs_MekDD [Mek_genes, ],
            { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
                     lwd=2)
              text(baseMean,log2FoldChange, Ras_genes, pos=2, col="dodgerblue")})

       plotMA(res_T53D4_vs_MekDD, main="T53D4 vs MekDD\nshrunken", ylim = c(-7,7),
             ylab = "log fold change (ratio of normalized T53D4 / MekDD)",
            xlab = "means of normalized counts")
       
       # #Identify genes on the plot ARPE19 vs Aktmyr
       # #  Step1 -> execute idx code line below. 
       # #  Step2 -> Click on a dot in the plot. 
       # # Step3 -> To finish, click on "finish" in the upper right hand corner of the plot
       #  idx <- identify(res_T53D4_vs_MekDD$baseMean,res_T53D4_vs_MekDD$log2FoldChange, labels = gene_names$external_gene_name)
       #  #  Step4 -> click here to see what you got!
       #  rownames(res_T53D4_vs_MekDD)[idx]
       
    })
       
       output$MAplot <- renderPlot ({
      #ARPE19 vs T53D4 plots##
       par(mfrow=c(1,1))
       #print(input$sample1)
       plotMA(res_ARPE19_vs_T53D4, main="ARPE19 vs T53D4", ylim = c(-20,20),
              ylab = "log fold change (ratio of normalized ARPE19 / T53D4)",
              xlab = "means of normalized counts",
              cex.main= 3,
              cex= 0.5,
              alpha=0.5)
       with(res_ARPE19_vs_T53D4 [selectedGenes, ],
            { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
                     lwd=2)
              text(baseMean,log2FoldChange,
                   selectedGenes, pos=2, col="dodgerblue")})
       with(res_ARPE19_vs_T53D4 [T53D4_genes, ],
            { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
                     lwd=2)
              text(baseMean,log2FoldChange, T53D4_genes, pos=2, col="dodgerblue")})
       })
       
    #    ##ARPE19 vs RasV12 plots##
    #    par(mfrow=c(1,1))
    #    #print(input$sample1)
    #    plotMA(res_ARPE19_vs_RasV12, main="ARPE19 vs RasV12", ylim = c(-20,20),
    #           ylab = "log fold change (ratio of normalized ARPE19 / RasV12)",
    #           xlab = "means of normalized counts",
    #           cex.main= 3,
    #           cex= 0.5,
    #           alpha= 0.5)
    #    with(res_ARPE19_vs_RasV12 [selectedGenes, ],
    #         { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
    #                  lwd=2)
    #           text(baseMean,log2FoldChange,
    #                selectedGenes, pos=2, col="dodgerblue")})
    #    with(res_ARPE19_vs_RasV12 [Ras_genes, ],
    #         { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
    #                  lwd=2)
    #           text(baseMean,log2FoldChange,
    #                Ras_genes, pos=2, col="dodgerblue")})
    #    
    #    ##ARPE19 vs MekDD plots##
    #    par(mfrow=c(1,1))
    #    #print(input$sample1)
    #    plotMA(res_ARPE19_vs_MekDD, main="ARPE19 vs MekDD", ylim = c(-20,20),
    #           ylab = "log fold change (ratio of normalized ARPE19 / MekDD)",
    #           xlab = "means of normalized counts",
    #           cex.main= 3,
    #           cex= 0.5,
    #           alpha= 0.5)
    #    with(res_ARPE19_vs_MekDD [selectedGenes, ],
    #         { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
    #                  lwd=2)
    #           text(baseMean,log2FoldChange,
    #                selectedGenes, pos=2, col="dodgerblue")})
    #    with(res_ARPE19_vs_MekDD [Mek_genes, ],
    #         { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
    #                  lwd=2)
    #           text(baseMean,log2FoldChange, Mek_genes, pos=2, col="dodgerblue")})
    #    
    #    ##ARPE19 vs Aktmyr plots##
    #    par(mfrow=c(1,1))
    #    #print(input$sample1)
    #    plotMA(res_ARPE19_vs_Aktmyr, main="ARPE19 vs Aktmyr ", ylim = c(-20,20),
    #           ylab = "log fold change (ratio of normalized ARPE19 / Aktmyr)",
    #           xlab = "means of normalized counts",
    #           cex.main= 3,
    #           cex= 0.5,
    #           alpha= 0.5)
    #    with(res_ARPE19_vs_Aktmyr [selectedGenes, ],
    #         { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
    #                  lwd=2)
    #           text(baseMean,log2FoldChange, selectedGenes, pos=2, col="dodgerblue")})
    #    with(res_ARPE19_vs_Aktmyr [Akt_genes, ],
    #         { points(baseMean,log2FoldChange, col= "dodgerblue", cex=1,
    #                  lwd=2)
    #           text(baseMean,log2FoldChange, Akt_genes, pos=2, col="dodgerblue")})
    #    
    #    
    # })
    
    #Plotcounts
    output$plotcounts<- renderPlot ({
      print("rendering plot")
      #Reorder samples to match experimental design
      dds$sample <- factor(dds$sample, levels=c("ARPE19", "T53D4", "RasV12", "Aktmyr", "MekDD"))
      dds$sample
            #HRAS
      plotCounts(dds, gene=which(rownames(normalized_genecounts)=="HRAS"),intgroup = c("sample"), cex.main=2, cex=1.5, xlab= "Sample", returnData = T)
      pc_HRAS1<-plotCounts(dds, gene=which(rownames(normalized_genecounts)=="HRAS"),intgroup = c("sample"), cex.main=2, cex=1.5, xlab= "Sample", returnData = T)
      ggplot (pc_BUB1,aes(x=sample,y=count)) +
        geom_boxplot(aes(fill=sample, alpha=0.9))+ geom_jitter() +
        scale_y_log10() + stat_summary( geom="point", shape=20, size=5, color="red", fill="red") +
        ggtitle("HRAS expression")
    })

    #heatmap showing mean counts for ALL significantly changing genes
    output$boxplot1<- renderPlot({
      print("rendering plot")
      filtered4 = normalized_genecounts
      melted4 = melt(normalized_genecounts)
      colnames(melted4) <- c("Gene","sample","counts")
      print(input$sample)
      filtered4 = melted4[melted4$Gene %in% input$All & melted4$sample %in% input$sample,]
      print(filtered4)
      plotCounts(dds, gene=which(rownames(normalized_genecounts)=="AKT1"),intgroup = c("sample"), cex.main=2, cex=1.5, xlab= "Sample", returnData = T)
      pc_AKT1<-plotCounts(dds, gene=which(rownames(normalized_genecounts)=="AKT1"),intgroup = c("sample"), cex.main=2, cex=1.5, xlab= "Sample", returnData = T)
      ggplot (pc_BUB1,aes(x=sample,y=count)) +
        geom_boxplot(aes(fill=sample, alpha=0.9))+ stat_summary( geom="point", shape=20, size=5, color="red", fill="red") +
        geom_jitter()+scale_y_log10()
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
        scale_fill_gradient(low = "green", high = "red")+
        theme_classic()
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
        geom_raster(aes(fill=log10(counts))) +
        ggtitle ("Means counts for all significantly changing genes")+
        theme_classic()
      })

      })
    }


shinyApp(ui=ui, server=server)

