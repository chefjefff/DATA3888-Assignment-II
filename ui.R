#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(GEOquery)
library(R.utils)
library(reshape2)
library(ggplot2)
library(limma)
library(dplyr)


# Define UI for application that draws a histogram
shinyUI(
  pageWithSidebar(
    
    headerPanel('How do hyper-parameters affect classifier accuracy for Kidney Transplant Outcomes'),
    
    sidebarPanel(
      

      sliderInput("R", "Please select the number of repetitions: ", 
                  min = 1, max = 50, value = 25, step = 1),
      
      sliderInput("cvK", "Please select the number of cross-validation folds: ",
                  min = 3, max = 10, value = 5, step = 1),
      
      sliderInput("K", "For K-Nearest Neighbour, please choose the number of neighbours: ",
                  min = 3, max = 20, value = 5, step = 1),

    ),
    
    mainPanel(
      plotOutput("myPlot")
    )
  )
)  