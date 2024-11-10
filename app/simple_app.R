library(shiny)
library(shinydashboard)
library(tidyverse)
library(minpack.lm)  # Non-linear least squares fitting
library(shinyWidgets)  # For enhanced UI elements

# Load data
df_pr <- read.csv("df_spot_rates.csv")

# Define the Nelson-Siegel model functions
f2 <- function(beta1, theta, tau) {
  beta1 * (1 - exp(-theta / tau)) / (theta / tau)
}

f3 <- function(beta2, theta, tau) {
  beta2 * ((1 - exp(-theta / tau)) / (theta / tau) - exp(-theta / tau))
}

f4 <- function(beta3, theta, tau2) {
  beta3 * ((1 - exp(-theta / tau2)) / (theta / tau2) - exp(-theta / tau2))
}

extendns <- function(theta, beta0, beta1, beta2, beta3, tau, tau2) {
  beta0 + f2(beta1, theta, tau) + f3(beta2, theta, tau) + f4(beta3, theta, tau2)
}

ns_solve_fn <- function(x, params) {
  extendns(x, params[1], params[2], params[3], params[4], params[5], params[6])
}

# Define the residual function for optimization
residuals_ns <- function(params, MAT, VAL) {
  res <- VAL - sapply(MAT, ns_solve_fn, params = params)
  return(res)
}

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Nelson-Siegel Spot Rate Model"),
  
  dashboardSidebar(
    dateInput(
      inputId = "selected_date",
      label = "Select a Date",
      min = min(df_pr$TIME_PERIOD),
      max = max(df_pr$TIME_PERIOD),
      value = min(df_pr$TIME_PERIOD)
    ),
    actionButton("run_analysis", "Run Analysis")
  ),
  
  dashboardBody(
    fluidRow(
      box(
        title = "Fitted Parameters:", status = "primary", solidHeader = TRUE,
        width = 12,
        uiOutput("params_display"),
        tags$hr(),
        verbatimTextOutput("rss_value")
      )
    ),
    fluidRow(
      box(
        title = "Fitted vs Observed Plot", status = "primary", solidHeader = TRUE,
        width = 12,
        plotOutput("yield_plot")
      )
    ),
    fluidRow(
      box(
        title = "Outliers Excluded", status = "warning", solidHeader = TRUE,
        width = 12,
        tableOutput("excluded_points")  # Use `tableOutput()` for the table
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive to filter data based on selected date
  data_filtered <- reactive({
    req(input$selected_date)
    df_pr %>% filter(TIME_PERIOD == input$selected_date)
  })
  
  # Reactive to run the non-linear least squares fitting
  fit_model <- eventReactive(input$run_analysis, {
    data_rel <- data_filtered()
    init_params <- c(1, -1, 1, 1, 1, 1)  # Initial guesses
    
    fit <- nls.lm(
      par = init_params,
      fn = function(x) residuals_ns(x, MAT = data_rel$MAT, VAL = data_rel$OBS_VALUE)
    )
    
    params_est <- fit$par  # Extract parameters
    data_rel <- data_rel %>%
      mutate(
        Predicted_Yield = sapply(MAT, ns_solve_fn, params = params_est),
        Residual = OBS_VALUE - Predicted_Yield
      )
    
    # Identify outliers (absolute residual > 2 * std dev)
    residual_threshold <- 2 * sd(data_rel$Residual)
    data_clean <- data_rel %>% filter(abs(Residual) <= residual_threshold)
    
    list(params = params_est, data_clean = data_clean, data_rel = data_rel)
  })
  
  # Display fitted parameters with LaTeX equations
  output$params_display <- renderUI({
    req(fit_model())
    params <- fit_model()$params
    withMathJax(
      helpText(
        paste0("\\( \\beta_0 = ", round(params[1], 4), "\\)"),
        paste0("\\( \\beta_1 = ", round(params[2], 4), "\\)"),
        paste0("\\( \\beta_2 = ", round(params[3], 4), "\\)"),
        paste0("\\( \\beta_3 = ", round(params[4], 4), "\\)"),
        paste0("\\( \\tau = ", round(params[5], 4), "\\)"),
        paste0("\\( \\tau_2 = ", round(params[6], 4), "\\)")
      )
    )
  })
  
  # Display RSS value
  output$rss_value <- renderText({
    req(fit_model())
    data_clean <- fit_model()$data_clean
    rss <- sum((data_clean$OBS_VALUE - data_clean$Predicted_Yield)^2)
    paste("Residual Sum of Squares (RSS):", round(rss, 4))
  })
  
  # Plot fitted vs observed values
  output$yield_plot <- renderPlot({
    req(fit_model())
    data_clean <- fit_model()$data_clean
    
    data_long <- data_clean %>%
      pivot_longer(cols = c(OBS_VALUE, Predicted_Yield), 
                   names_to = "Type", values_to = "Value")
    
    ggplot(data_long, aes(x = MAT, y = Value, color = Type, linetype = Type)) +
      geom_point(data = filter(data_long, Type == "OBS_VALUE"), size = 1) +
      geom_line(size = 1) +
      labs(
        title = paste("Fitted vs Observed Values on", input$selected_date),
        x = "Maturity (MAT)", y = "Spot Rates"
      ) +
      theme_minimal() +
      scale_color_manual(values = c("OBS_VALUE" = "blue", "Predicted_Yield" = "red"))
  })
  
  # Display excluded points as a table
  output$excluded_points <- renderTable({
    req(fit_model())
    data_rel <- fit_model()$data_rel
    data_clean <- fit_model()$data_clean
    # Identify excluded points by finding the rows in `data_rel` not present in `data_clean`
    excluded <- setdiff(data_rel, data_clean)
    excluded  # Return the full table of excluded points
  })
}

# Run the Shiny app
shinyApp(ui, server)