#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("EM Algorithm (Mixture Model (Gaussian variety))"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          fileInput("file", "Choose a CSV file:",
                    placeholder = "Choose a file now or a calc 1 student does (a + b)^2 = a^2 + b^2"),
          
          # column number in dataset
          numericInput("column",
                       "Column number in dataset",
                       value = 1,
                       min = 1),
          
          # Input for histogram
          sliderInput("bins",
                      "Number of bins:",
                      min = 1,
                      max = 50,
                      value = 30),
          
          # Input for number of mixture components
          numericInput("components",
                       "# of mixture components",
                       value = 2,
                       min = 1,
                       max = 100),
          # Input for how many EM steps to run
          numericInput("n",
                       "Number of iterations of EM to run",
                       value = 50,
                       min = 1,
                       max = 100),
          # Choose size of y axis?
          checkboxInput("manual_y_axis", "Choose y maximum?", F),
          
          # Input for size of y axis
          conditionalPanel(
            condition = "input.manual_y_axis == true",
            numericInput("y_max",
                         "y maximum",
                         value = 1)
          ),
      
          # Input for EM step shown
          numericInput("EM_step",
                       "EM step",
                       value = 1,
                       min = 1,
                       max = 50),
          
          # checkboxes for initial, final, and EM_step
          checkboxInput("initial_plot", "Initial guesstimation plot", F),
          checkboxInput("final_plot", "Final EM algorithm output plot", T),
          checkboxInput("EM_step_plot", "Your EM step plot", T),
          checkboxInput("legend", "Show the legend??????", T),
          
          width = 5
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel(
            tabPanel("Plot", plotOutput("distPlot", width = "600px"))
          ), width = 5
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    #### Observe ####
    observe({
      if(is.null(input$n) || is.na(input$n)){n = 50}else{n=input$n} #Prevents Crash
        # make the max EM_step the number of iterations
        updateNumericInput(session, "EM_step", max = n)
        
        # make sure maxes are not exceeded
        if(is.null(input$EM_step) || is.na(input$EM_step)){EM = 1}else{EM=input$EM_step} #Prevents Crash
        if (EM > n){
          updateNumericInput(session, "EM_step", value = n)
        } 
      })
  
    input_file <- reactive({
      if (is.null(input$file)) {
        return("")
      }
      
      # actually read the file
      read.csv(file = input$file$datapath)
    })
    #### EM Function ####
    EM_function <- function(X, m, n = input$n){
       if(is.null(n) || missing(n) || is.na(n)) n = 50  #Prevents Crash
      # performs the Expectation-maximization algorithm on data X
      # with m mixing components, and n = 50 iterations at most
      
      X = as.numeric(X); N = length(X)
      
      if (m != 1){
        q=(1:(m-1))*1/m
        quantiles=c(unname(quantile(X,q)), max(X))
        tempmu=numeric(m); tempsigma=numeric(m)
        
        # first part
        tempmu[1]=mean(X[X<=quantiles[1]]); tempsigma[1]=sd(X[X<=quantiles[1]])
        for (i in 2:m){
          tempmu[i]=mean(X[X>quantiles[i-1] & X<=quantiles[i]]); tempsigma[i]=sd(X[X>quantiles[i-1] & X<=quantiles[i]])
        }
      } else {
        tempmu=c(mean(X)); tempsigma=c(sd(X))
      }
      
      # initialize the parameters
      pi = matrix(0,n,m); mu = matrix(0,n,m); sigma = matrix(0,n,m)
      pi[1,] = rep(1/m,m) # weight them equally
      # mu[1,] = seq(from = min(X), to = max(X), length.out = m) # space them out evenly
      # sigma[1,] = rep(diff(range(X))/(6*m),m) # spread them out evenly
      mu[1,] = tempmu
      sigma[1,] = tempsigma
      
      theta = list(pi = pi, mu = mu, sigma = sigma)
      
      iter = 1
      while (iter < n){
        # calculate the conditional probabilities
        p = matrix(0,N,m)
        for (j in 1:m){
          p[,j] = theta$pi[iter,j]*dnorm(X, mean = theta$mu[iter,j], sd = theta$sigma[iter,j])
        }
        if (any(is.infinite(p))){
          theta$pi=theta$pi[1:iter,]
          theta$mu=theta$mu[1:iter,]
          theta$sigma=theta$sigma[1:iter,]
          return(theta)
        }
        
        if (iter == 1){
          LogLikelihood = c(log(sum(p[iter,])))
        } else {
          LogLikelihood = c(LogLikelihood, log(sum(p[iter,])))
        }
        
        p.hat = t(apply(p, 1, function(x){x/sum(x)}))
        if (m == 1){p.hat = t(p.hat)} # R doesn't like mathematicians
        
        # calculate the new parameters
        
        new_pi = apply(p.hat,2,sum)/N
        new_mu = apply(p.hat,2,function(x){sum(x*X)/sum(x)})
        new_sigma = sqrt(apply(p.hat,2,function(x){sum(x*(X-sum(x*X)/sum(x))^2)/sum(x)}))
        
        theta$pi[iter+1,]=new_pi
        theta$mu[iter+1,]=new_mu
        theta$sigma[iter+1,]=new_sigma

        if (iter > 1){
          if(abs((LogLikelihood[iter] - LogLikelihood[iter-1])/LogLikelihood[iter-1]) < 1e-5){
            theta$pi=theta$pi[1:(iter+1),]
            theta$mu=theta$mu[1:(iter+1),]
            theta$sigma=theta$sigma[1:(iter+1),]
            break
          }
        }
        
        iter = iter + 1
      }
      
      return(theta)
    }

    output$distPlot <- renderPlot({
      
        if (is.null(input$file)) {
          return("")
        }
      
        # generate bins based on input$bins from ui.R
        user_column = input$column
        if (user_column < 1){
          column = 1
        } else {
          max_column = ncol(input_file())
          column = min(floor(user_column), max_column)
        }
        
        observe({
          # update the user's error in column number
          updateNumericInput(session, "column", value = column)
        })
        
        x    <- input_file()[,column]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)
        
        m = round(input$components)
        
        EM_fit = EM_function(x, m)
        theta = EM_fit
        final=nrow(theta$pi)
        
        observe({
          # make the max EM_step the number of iterations
          updateNumericInput(session, "EM_step", max = final)
          
          # make sure maxes are not exceeded
          if (input$EM_step > final){
            updateNumericInput(session, "EM_step", value = final)
          } 
        })
        
        x_values = seq(min(x),max(x),by=0.05)
        y_values = matrix(0, length(x_values), m)
        
        # initial guesstimate
        for (j in 1:m){
          y_values[,j] = theta$pi[1,j]*dnorm(x_values, mean = theta$mu[1,j], sd = theta$sigma[1,j])
        }
        y_values_initial = rowSums(y_values)
        
        # output at user defined step
        EM_step=round(input$EM_step)
        for (j in 1:m){
          y_values[,j] = theta$pi[EM_step,j]*dnorm(x_values, mean = theta$mu[EM_step,j], sd = theta$sigma[EM_step,j])
        }
        y_values_EM_step = rowSums(y_values)
        
        # final output
        final=nrow(theta$pi)
        for (j in 1:m){
          y_values[,j] = theta$pi[final,j]*dnorm(x_values, mean = theta$mu[final,j], sd = theta$sigma[final,j])
        }
        y_values = rowSums(y_values)
        
        fake_h <- hist(x, breaks = bins, freq = F, plot = F)
        max_hist = fake_h$density[which.max(fake_h$density)]
        max_EM_step = max(y_values_EM_step)
        max_final = max(y_values)
        max_y = 0.1*ceiling(10*max(max_hist, max_EM_step, max_final))
        
        if (input$manual_y_axis){
          y_max = input$y_max
        } else {
          y_max = max_y
        }
        
        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, freq = F, col = 'darkgray', border = 'white',
             xlab = 'Data (data)', ylab = 'Density',
             main = 'Histogram and Gaussian Mixture Model EM fit',
             xlim = c(floor(min(x)), ceiling(max(x))), 
             ylim = c(0, y_max))
        if (input$initial_plot){
          lines(x_values, y_values_initial, col = "green", lwd=5, lty=3)
        }
        if (input$final_plot){
          lines(x_values, y_values, col = "red", lwd=5)
        }
        if (input$EM_step_plot){
          lines(x_values, y_values_EM_step, col = "blue", lwd=5)
        }
        legend_elements <- list(
          initial = list(label = "initial guesstimate", color = "green", type = 3),
          final = list(label = "final guesstimate", color = "red", type = 1),
          EM_step = list(label = sprintf("EM step %d guesstimate", EM_step), color = "blue", type = 1)
        )
        active_elements <- legend_elements[c(input$initial_plot, input$final_plot, input$EM_step_plot)]

        # Construct and display the legend
        if (input$legend && length(active_elements) > 0) {
          legend_labels <- sapply(active_elements, `[[`, "label")
          legend_colors <- sapply(active_elements, `[[`, "color")
          legend_types <- sapply(active_elements, `[[`, "type")
          legend(x = "topright", 
                lty = legend_types, 
                col = legend_colors,
                legend = legend_labels)
        }
        
    })

}

# Run the application 
shinyApp(ui = ui, server = server)
