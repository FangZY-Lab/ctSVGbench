library(spacexr) #C-SIDE
library(CTSV)
library(spVC)
library(sp)
library(BPST)
library(Triangulation)
library(MGLM)
library(CELINA)
library(STANCE)
library(ctSVG)
library(SpatialExperiment)
library(spacexr)
library(Matrix)
library(devtools)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(ggrepel)
library(here)
##get position
dir.create("real/coords", recursive = TRUE, showWarnings = FALSE)
dir.create("real/pos", recursive = TRUE, showWarnings = FALSE)
dir.create("real/bounary", recursive = TRUE, showWarnings = FALSE)

files=list.files(here('real','puck'))

lapply(files,function(file){
dataset=gsub('myRCTD_|.rds',"",file)
puck<- readRDS(here('real','puck',file))
pos <- puck@coords
  ggplot(pos,aes(x=x,y=y))+
    geom_point()
  ggsave(here('real','coord',gsub('.rds','.pdf',file)))
saveRDS(pos,file=here('real','pos',file))
})


for(file in files){
  dataset=gsub('myRCTD_|.rds',"",file)
puck<- readRDS(here('real','puck',file))
pos <- puck@coords
  ggplot(pos,aes(x=x,y=y))+
    geom_point()
  ggsave(here('real','coord',gsub('.rds','.pdf',file)))
  saveRDS(pos,file=here('real','pos',file))
}  


##get boundary
library(spVC)
library(BPST)
library(MGLM)
library(ggplot2)
library(Triangulation)
library(FNN) 

setwd("F:/ctSVGbench/real/pos/")

library(shiny)
library(ggplot2)
library(plotly)

ui <- fluidPage(
  plotlyOutput("plot"),
  tableOutput("click_table")  # Show matched points
)

server <- function(input, output, session) {
  clicked_idx <- reactiveValues(index = integer(0))
  
  output$plot <- renderPlotly({
    ggplotly(
      ggplot(S, aes(x = x, y = y)) +
        geom_point()
    )
  })
  
  # Handle click event
  observeEvent(event_data("plotly_click"), {
    d <- event_data("plotly_click")
    if (!is.null(d)) {
      clicked_point <- c(d$x, d$y)
      
      # Compute Euclidean distance
      distances <- sqrt((S$x - clicked_point[1])^2 + (S$y - clicked_point[2])^2)
      nearest_idx <- which.min(distances)
      
      # Accumulate indices
      clicked_idx$index <- c(clicked_idx$index, nearest_idx)
      
      # Save globally
      clicked_row_indices <<- clicked_idx$index
    }
  })
  
  # Show matched rows
  output$click_table <- renderTable({
    S[clicked_idx$index, , drop = FALSE]
  })
}

shinyApp(ui, server)


i=1
file=list.files('./')[i]
S=readRDS(file)
plot(S[, 1], S[, 2], pch = ".",cex=5)
shinyApp(ui, server)
boundary <- S[clicked_row_indices, ]
points(boundary, type = "l", col = "blue")
Tr.cell <- TriMesh(boundary, n = 2) # n : triangulation fineness
saveRDS(boundary,sprintf("../boundary/%s",file))
i=i+1

