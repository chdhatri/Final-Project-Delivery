#############################CMPE 239 PROJECT#####################
#TEAM VOLTAGE
#Sucharu Gupta
#Dhatri CHennavajula
#Andrew Wong
#Sprush Ujjwal

#UNCOMMENT THESE WHEN YOU ARE RUNNING FIRST TIME TO INSTALL PACKAGES
#install.packages("shiny")
#install.packages("e1071")
#install.packages("randomForest")

library(shiny)
library(e1071)
library(randomForest)
shinyServer(function(input, output) {
 setwd("/Users/ssarma/Google Drive/CMPE 239/Project")
  
  output$playerinfo <- renderTable({
    
    genes.data <- read.csv(file='test.data14_15.csv', sep=',', h=T)
    gene.data<-subset(genes.data, 
                      DOSE == input$dose & 
                        newCount == input$count &
                        STRAIN ==input$strain &
                        TRIAL_RESULT == input$TrailResult)
    
    output.data <- gene.data [c("DOSE", "newCount", "STRAIN", "MICROSOMAL_ACTIVATION_USED",
                                "TRIAL_RESULT","RESULT","STUDY_CONCLUSION")]
    output.data
    
  })
  
  output$result <- renderPrint({
    
    
   
    
    ###########################################################
    ## build test set data frame
    ###########################################################
   #test.data <- data.frame(DOSE = as.numeric(gene.data$DOSE),
  
    test.raw.data <- data.frame(DOSE = as.numeric(input$dose),
                                newCount =  as.integer(input$count),
                                STRAIN = as.factor(input$strain),
                                TRIAL_RESULT = as.factor(input$TrailResult))
    ##########################################################
    ## load training data
    ###########################################################
    
    train.raw.data  <- read.table(file='gene.train.csv', sep=',', h=T)
    model.data <- train.raw.data [c("DOSE","newCount","STRAIN","TRIAL_RESULT","RESULT")]
    
    model.data$RESULT <- as.factor(model.data$RESULT)
    
    ########################################################
    # BUILD RF Prediction
    ########################################################
    
   # str(model.data)
    levels(test.raw.data$STRAIN) = levels(model.data$STRAIN)
    levels(test.raw.data$TRIAL_RESULT) = levels(model.data$TRIAL_RESULT)
    # make predictions
  
    rf1.out <- randomForest(RESULT ~ DOSE + newCount +STRAIN+ TRIAL_RESULT, 
                        data = gene.train, importance=T, ntree=500)
    
   
     # make predictions
    glm.pred<-predict(rf1.out, test.raw.data)
    paste("The Chemical with the given inputs is  ", toString(glm.pred))
    
   
  })
  
  
  output$studyConclusion <- renderPrint({
    
    ###########################################################
    ## build test set data frame
    ###########################################################
    #test.data <- data.frame(DOSE = as.numeric(gene.data$DOSE),
    
    test.raw.data <- data.frame(DOSE = as.numeric(input$dose),
                                newCount =  as.integer(input$count),
                                STRAIN = as.factor(input$strain),
                                TRIAL_RESULT = as.factor(input$TrailResult))
    ##########################################################
    ## load training data
    ###########################################################
    
    train.raw.data  <- read.table(file='gene.train.csv', sep=',', h=T)
    model.data <- train.raw.data [c("DOSE","newCount","STRAIN","TRIAL_RESULT","RESULT","STUDY_CONCLUSION")]
    
    model.data$RESULT <- as.factor(model.data$RESULT)
    
    ########################################################
    # BUILD RF Prediction
    ########################################################
    
    test.raw.data <- data.frame(DOSE = as.numeric(5000),
                                newCount =  as.integer(0),
                                STRAIN = as.factor("TA98"),
                                TRIAL_RESULT = as.factor("Negative"))
    
    levels(test.raw.data$STRAIN) = levels(model.data$STRAIN)
    levels(test.raw.data$TRIAL_RESULT) = levels(model.data$TRIAL_RESULT)
   
    
    rf1.out <- randomForest(RESULT ~ DOSE + newCount +STRAIN+ TRIAL_RESULT, 
                            data = gene.train, importance=T, ntree=500)
    
    
    # make predictions
    rf1.pred<-predict(rf1.out, test.raw.data)
    
    test.raw.data$RESULT <- rf1.pred
    levels(test.raw.data$RESULT) = levels(model.data$RESULT)
   # str(test.raw.data)
    paste("The Chemical with the given inputs is  ",  test.raw.data$RESULT)
   
   # levels(as.factor(test.raw.data$RESULT)) = levels(model.data$RESULT)
    rf2.two <- randomForest(STUDY_CONCLUSION ~ DOSE +newCount+STRAIN+TRIAL_RESULT +RESULT, 
                            data = model.data, importance=T, ntree=500)
    #make predictions
    rf.pred<-predict(rf2.two, test.raw.data)
    paste(rf.pred)
    
  })
  
 
  
  
  })