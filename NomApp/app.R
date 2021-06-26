library(shiny)
library(shinyWidgets)
library(openxlsx)
library(data.table)
library(readr)
library(dplyr)

#Reading in Data
df.predict <- read_csv("point_risks.csv")
df.outcome <- read_csv("mat_fin.csv")

#Dropping erraneous first column
df.predict <- df.predict %>% select(-"X1")
df.outcome <- df.outcome %>% select(-"X1")



#Collecting the age values and Phom Scores (Slider Values) by sex
age.values <- subset(df.predict, predictor == "Age")$value %>% as.numeric()

phom.values <- subset(df.predict, predictor == "PHOM Score")$value %>% as.numeric()


parameter_tabs <- tabsetPanel(
  id = "params",
  type = "hidden",
  tabPanel("Ref", #Here, list the inputs for
           
           shinyWidgets::sliderTextInput("PHOM Score", 
                                         "PHOM Score", 
                                         choices = phom.values,
                                         grid = T),
           
           shinyWidgets::sliderTextInput("Age", #Input ID, this gets referenced in the server
                                         "Age at Diagnosis", #Input Label, this is what the user views
                                         choices = age.values,
                                         grid = T),
           
           selectInput('CCI', 
                       'Charlson Comorbidity Index',
                       c("0-4", "5-8", "> 8")),
           
           selectInput('KPS', 
                       'Karnofsky Performance Status',
                       c("KPS ≤ 70", "KPS > 70")),
           
           selectInput('Overall Stage', 
                       'Overall Stage',
                       c("IA","IB", "IIA", "IIB")),
           
           selectInput('Post Radiation Chemo', 
                       'Post Radiation Chemotherapy',
                       c("No", "Yes")),
           
           selectInput('Sex', 
                       'Sex',
                       c("Female","Male"))

  )
)




ui <- fluidPage(

  sidebarLayout(
    sidebarPanel(
      parameter_tabs,
      
    ),
    mainPanel(
      br(),
      h4("Expected Survival Probabilities", align ="center"),
      fluidRow(
        column(12,align="center",
               h5(tableOutput("Predicted_Surv_Prob")))),
      fixedRow(
        #img(src=img.src, height = 500, width = 700, style="display: block; margin-left: auto; margin-right: auto;")
        imageOutput("myImage", hover = FALSE, inline = TRUE, width = "80%")
      ),
      uiOutput(outputId = "text")
    )
    ))
  

vars.set <- c("PHOM Score", "Age", "CCI", "KPS", "Overall Stage",
          "Post Radiation Chemo", "Sex")

server <- function(input, output, session) {

  
  #Creating the data frame with appropriate output variables
  final.df <- reactive({
    
      
    #Swtich function creates the apprporiate value vector for the corresponding variables
    vals <- c(as.character(input$`PHOM Score`), 
              as.character(input$Age),
              as.character(input$`CCI`), 
              as.character(input$`KPS`), 
              as.character(input$`Overall Stage`), 
              as.character(input$`Post Radiation Chemo`), 
              as.character(input$`Sex`))
    
    #creating a list to store the chosen variabl value and associated points
    #For each variable, subsets the variable and the chosen value
    #Goal is to add points at the end
    test.list <- list()
    for (i in 1:length(vals)) {
      test.list[[i]] <- subset(df.predict, predictor == vars.set[i] & 
                                                     value == vals[i])
    }
    
    #Converting to dataframe for compatability with sum
    test.df <- do.call(rbind, test.list)
    
    #TOtal points
    total.points <- sum(test.df$points)
      
      
      data.output <-setDT(as.data.frame(t(subset(df.outcome, 
                                          Points == total.points))), 
                          keep.rownames = TRUE)
      
      data.output$rn <- gsub("[.]", " ", data.output$rn)
      
      return(data.output)
    
  })
  
  #Image output
  output$myImage <- renderImage({
    #Switch the nomogram image source string depending on the input
    img.src <- "./www/nomogram.png"
    
    # Return a list containing the filename
    list(src = eval(img.src),
         contentType = 'image/png',
         width = 700,
         height = 500,
         alt = "nomogram")
  }, deleteFile = FALSE)

  
  
  output$Predicted_Surv_Prob <- renderTable(final.df(), colnames = FALSE)
  
  output$text <- renderText({
    HTML(paste0('placeholder')
    )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
