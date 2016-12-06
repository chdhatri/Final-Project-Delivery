#############################CMPE 239 PROJECT#####################
#TEAM VOLTAGE
#Sucharu Gupta
#Dhatri CHennavajula
#Andrew Wong
#Sprush Ujjwals

library("shiny")    
shinyUI(pageWithSidebar(
  headerPanel("CARCINOGENICITY AND MUTAGENICIT PREDICTOR"),
  sidebarPanel(
    numericInput("dose", "Enter the dose:", 0),
    numericInput("count", "Enter the count:", 95),
    selectInput("strain", "Choose a strain:", 
                choices = c("TA100","eColi pKM101", "TA98")),
    selectInput("TrailResult", "Choose a Trail Result:", 
                choices = c("Negative","Equivocal","Positive","Weakly Positive")),
   
   
    
    submitButton("Calculate")
  ),
  
  
  
  mainPanel(
    h4("Test Data"),
    tableOutput("playerinfo"),
    
    h4("Output From Logistic Regression"),
    verbatimTextOutput("result"),
    
    verbatimTextOutput("studyConclusion")
    
   
    
  )
))