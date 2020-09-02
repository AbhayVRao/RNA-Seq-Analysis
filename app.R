
library(shiny)
library(tidyverse)
library(dplyr)
library(limma)
library(edgeR)
library(org.Hs.eg.db)
library(ggrepel)
library(readr)
library(plot3D)
library(goseq)
library(fgsea)

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("iPSC Derived Cardiomyocyte RNA Seq. Data"),
    
    # Sidebar with a slider input for number of bins
    sidebarLayout(sidebarPanel(
        selectInput("pathway",
                    "Select a Pathway",
                    choices = gsea_res.top$pathway),
        sliderInput("slider", "Select a Number of Labels:",
                    min = 0, max = 750,
                    value = 0, step = 50,
                    animate =
                        animationOptions(interval = 5000, loop = TRUE)),
    ),
    
    mainPanel(tabsetPanel(
        tabPanel("Genes", plotOutput("genes")),
        tabPanel("Enrichment Score", plotOutput("ES"))
    )))
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$genes <-
        renderPlot({
            x<- which(startsWith(gsea_res.top$pathway, input$pathway))
            genes_a <- c(gsea_res.top$leadingEdge[x]) 
            results.table.1 <- results.table.1 %>% 
                mutate(rank = rank(adj.P.Val, ties.method = "random"))
            results.table.1 %>%
                ggplot(aes(x=logFC, y=-log10(adj.P.Val)))+
                geom_point(color="grey",size=1) +
                geom_point(data=subset(results.table.1, Symbol %in% genes_a[[1]] & adj.P.Val < 0.05), 
                           aes(color=2^(AveExpr)), pch=8) +
                geom_vline(xintercept = 1, linetype=3) +
                geom_vline(xintercept = -1, linetype=3) +
                geom_hline(yintercept = 1.3, linetype=3) +
                geom_text_repel(data=subset(results.table.1, Symbol %in% genes_a[[1]] & rank <= input$slider),
                                aes(label = Symbol)) +
                scale_color_gradient(low="green",high="red") +
                ggtitle("iPSC Derived CM Genes") +
                labs(x="log2 Fold Change", y="Stat. Significance", 
                     color="Avg. Expr. of Gene") +
                theme_minimal() + 
                theme(axis.text = element_text(size = 10))
            
        })
    
    output$ES <- renderPlot({
        plotEnrichment(Pathways[[input$pathway]], gsea.input)
        
    })
}

# Run the application
shinyApp(ui = ui, server = server)

