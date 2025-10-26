# Shiny Dashboard for Drug Discovery Results
# Interactive visualization and exploration of drug discovery models

library(shiny)
library(shinydashboard)
library(shinycssloaders)
library(shinyWidgets)
library(DT)
library(plotly)
library(ggplot2)
library(dplyr)
library(tidyr)
library(fresh)
library(bslib)

# Custom theme
custom_theme <- create_theme(
  adminlte_color(
    light_blue = "#0073b7",
    blue = "#3c8dbc"
  ),
  adminlte_sidebar(
    width = "250px",
    dark_bg = "#222d32",
    dark_hover_bg = "#1e282c"
  ),
  adminlte_global(
    content_bg = "#f4f6f9"
  )
)

# UI Definition
ui <- dashboardPage(
  skin = "blue",
  
  # Header
  dashboardHeader(
    title = span(
      icon("pills"),
      "Drug Discovery ML"
    ),
    titleWidth = 250
  ),
  
  # Sidebar
  dashboardSidebar(
    width = 250,
    sidebarMenu(
      id = "tabs",
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Data Explorer", tabName = "data", icon = icon("database")),
      menuItem("Chemical Space", tabName = "chemspace", icon = icon("atom")),
      menuItem("Model Performance", tabName = "models", icon = icon("chart-line")),
      menuItem("Predictions", tabName = "predictions", icon = icon("magic")),
      menuItem("Feature Importance", tabName = "features", icon = icon("star")),
      menuItem("System Monitor", tabName = "system", icon = icon("microchip")),
      menuItem("Help", tabName = "help", icon = icon("question-circle"))
    )
  ),
  
  # Body
  dashboardBody(
    use_theme(custom_theme),
    
    tabItems(
      # Dashboard Tab
      tabItem(
        tabName = "dashboard",
        fluidRow(
          valueBoxOutput("total_compounds", width = 3),
          valueBoxOutput("active_compounds", width = 3),
          valueBoxOutput("best_model", width = 3),
          valueBoxOutput("best_r2", width = 3)
        ),
        
        fluidRow(
          box(
            title = "Project Overview",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            
            p("Welcome to the Bioinformatics Drug Discovery Dashboard!"),
            hr(),
            h4("Project Summary:"),
            tags$ul(
              tags$li("Target: SARS-CoV-2 3C-like Proteinase"),
              tags$li("Data Source: ChEMBL Database"),
              tags$li("Machine Learning Models: Random Forest, XGBoost, Neural Network"),
              tags$li("Features: Lipinski Descriptors + Molecular Fingerprints")
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Bioactivity Distribution",
            status = "info",
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            withSpinner(plotlyOutput("bioactivity_dist"))
          ),
          
          box(
            title = "Model Comparison",
            status = "success",
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            withSpinner(plotlyOutput("model_comparison_plot"))
          )
        )
      ),
      
      # Data Explorer Tab
      tabItem(
        tabName = "data",
        fluidRow(
          box(
            title = "Dataset Filters",
            status = "primary",
            solidHeader = TRUE,
            width = 3,
            
            selectInput(
              "bioactivity_filter",
              "Bioactivity Class:",
              choices = c("All", "Highly Active", "Active", "Moderately Active", "Inactive"),
              selected = "All"
            ),
            
            sliderInput(
              "pic50_range",
              "pIC50 Range:",
              min = 4,
              max = 10,
              value = c(4, 10),
              step = 0.1
            ),
            
            checkboxInput("druglike_only", "Drug-like compounds only", FALSE),
            
            actionButton("apply_filters", "Apply Filters", class = "btn-primary")
          ),
          
          box(
            title = "Compound Data Table",
            status = "info",
            solidHeader = TRUE,
            width = 9,
            collapsible = TRUE,
            
            withSpinner(DTOutput("data_table"))
          )
        )
      ),
      
      # Chemical Space Tab
      tabItem(
        tabName = "chemspace",
        fluidRow(
          box(
            title = "Chemical Space Analysis",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            fluidRow(
              column(
                width = 4,
                selectInput(
                  "x_axis",
                  "X-Axis:",
                  choices = c("MW", "LogP", "HBD", "HBA", "TPSA", "nRotB"),
                  selected = "MW"
                )
              ),
              column(
                width = 4,
                selectInput(
                  "y_axis",
                  "Y-Axis:",
                  choices = c("LogP", "MW", "HBD", "HBA", "TPSA", "nRotB"),
                  selected = "LogP"
                )
              ),
              column(
                width = 4,
                selectInput(
                  "color_by",
                  "Color By:",
                  choices = c("Bioactivity Class" = "bioactivity_class", 
                             "pIC50" = "pIC50",
                             "Drug-like" = "druglike"),
                  selected = "bioactivity_class"
                )
              )
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Chemical Space Scatter Plot",
            status = "info",
            solidHeader = TRUE,
            width = 8,
            collapsible = TRUE,
            withSpinner(plotlyOutput("chemical_space_plot", height = "600px"))
          ),
          
          box(
            title = "Lipinski Rule of Five",
            status = "warning",
            solidHeader = TRUE,
            width = 4,
            
            h4("Lipinski's Rules:"),
            tags$ul(
              tags$li("Molecular Weight ≤ 500 Da"),
              tags$li("LogP ≤ 5"),
              tags$li("H-bond Donors ≤ 5"),
              tags$li("H-bond Acceptors ≤ 10")
            ),
            hr(),
            withSpinner(plotOutput("lipinski_violations"))
          )
        )
      ),
      
      # Model Performance Tab
      tabItem(
        tabName = "models",
        fluidRow(
          box(
            title = "Model Selection",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            selectInput(
              "selected_model",
              "Choose Model:",
              choices = c("Random Forest", "XGBoost", "Neural Network"),
              selected = "Random Forest"
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Predicted vs Actual",
            status = "info",
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            withSpinner(plotlyOutput("pred_vs_actual"))
          ),
          
          box(
            title = "Residual Plot",
            status = "warning",
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            withSpinner(plotlyOutput("residual_plot"))
          )
        ),
        
        fluidRow(
          box(
            title = "Performance Metrics",
            status = "success",
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            withSpinner(DTOutput("metrics_table"))
          ),
          
          box(
            title = "Training History (Neural Network)",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            withSpinner(plotlyOutput("training_history"))
          )
        )
      ),
      
      # Predictions Tab
      tabItem(
        tabName = "predictions",
        fluidRow(
          box(
            title = "Make New Predictions",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            p("Upload a CSV file with SMILES strings or use example compounds:"),
            
            fluidRow(
              column(
                width = 6,
                fileInput(
                  "smiles_file",
                  "Upload SMILES CSV:",
                  accept = c(".csv")
                )
              ),
              column(
                width = 6,
                actionButton("use_examples", "Use Example Compounds", class = "btn-info"),
                br(), br(),
                actionButton("predict_btn", "Generate Predictions", class = "btn-success")
              )
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Prediction Results",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            
            withSpinner(DTOutput("prediction_results")),
            br(),
            downloadButton("download_predictions", "Download Results", class = "btn-primary")
          )
        )
      ),
      
      # Feature Importance Tab
      tabItem(
        tabName = "features",
        fluidRow(
          box(
            title = "Feature Importance Analysis",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            sliderInput(
              "top_n_features",
              "Number of Top Features to Display:",
              min = 10,
              max = 50,
              value = 20,
              step = 5
            )
          )
        ),
        
        fluidRow(
          box(
            title = "Random Forest Feature Importance",
            status = "info",
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            withSpinner(plotlyOutput("rf_importance", height = "600px"))
          ),
          
          box(
            title = "XGBoost Feature Importance",
            status = "success",
            solidHeader = TRUE,
            width = 6,
            collapsible = TRUE,
            withSpinner(plotlyOutput("xgb_importance", height = "600px"))
          )
        )
      ),
      
      # System Monitor Tab
      tabItem(
        tabName = "system",
        fluidRow(
          valueBoxOutput("cpu_usage", width = 3),
          valueBoxOutput("memory_usage", width = 3),
          valueBoxOutput("gpu_status", width = 3),
          valueBoxOutput("disk_usage", width = 3)
        ),
        
        fluidRow(
          box(
            title = "System Resources",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            verbatimTextOutput("system_info")
          )
        )
      ),
      
      # Help Tab
      tabItem(
        tabName = "help",
        fluidRow(
          box(
            title = "User Guide",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            
            h3("How to Use This Dashboard"),
            hr(),
            
            h4("1. Dashboard"),
            p("Overview of the project with key metrics and visualizations."),
            
            h4("2. Data Explorer"),
            p("Browse and filter the compound dataset. Use filters to focus on specific bioactivity classes or pIC50 ranges."),
            
            h4("3. Chemical Space"),
            p("Visualize molecular properties and chemical space. Check Lipinski rule compliance."),
            
            h4("4. Model Performance"),
            p("Compare different machine learning models. View prediction accuracy and residual plots."),
            
            h4("5. Predictions"),
            p("Upload new SMILES strings to predict bioactivity using trained models."),
            
            h4("6. Feature Importance"),
            p("Analyze which molecular descriptors are most important for predictions."),
            
            h4("7. System Monitor"),
            p("Monitor system resource usage during computation."),
            
            hr(),
            h4("References:"),
            tags$ul(
              tags$li(a("ChEMBL Database", href = "https://www.ebi.ac.uk/chembl/", target = "_blank")),
              tags$li(a("Lipinski's Rule of Five", href = "https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five", target = "_blank")),
              tags$li(a("Random Forest Algorithm", href = "https://github.com/imbs-hl/ranger", target = "_blank"))
            )
          )
        )
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  
  # Reactive data loading
  full_data <- reactive({
    req(file.exists("data/full_dataset_with_descriptors.csv"))
    read_csv("data/full_dataset_with_descriptors.csv", show_col_types = FALSE)
  })
  
  rf_predictions <- reactive({
    req(file.exists("results/rf_test_predictions.csv"))
    read_csv("results/rf_test_predictions.csv", show_col_types = FALSE)
  })
  
  xgb_predictions <- reactive({
    req(file.exists("results/xgb_test_predictions.csv"))
    read_csv("results/xgb_test_predictions.csv", show_col_types = FALSE)
  })
  
  nn_predictions <- reactive({
    req(file.exists("results/nn_test_predictions.csv"))
    read_csv("results/nn_test_predictions.csv", show_col_types = FALSE)
  })
  
  model_comparison <- reactive({
    req(file.exists("results/model_comparison.csv"))
    read_csv("results/model_comparison.csv", show_col_types = FALSE)
  })
  
  # Value boxes
  output$total_compounds <- renderValueBox({
    data <- full_data()
    valueBox(
      nrow(data),
      "Total Compounds",
      icon = icon("flask"),
      color = "blue"
    )
  })
  
  output$active_compounds <- renderValueBox({
    data <- full_data()
    active <- sum(data$bioactivity_class %in% c("Highly Active", "Active"))
    valueBox(
      active,
      "Active Compounds",
      icon = icon("check-circle"),
      color = "green"
    )
  })
  
  output$best_model <- renderValueBox({
    comparison <- model_comparison()
    best <- comparison %>% arrange(RMSE) %>% slice(1)
    valueBox(
      best$Model,
      "Best Model (RMSE)",
      icon = icon("trophy"),
      color = "yellow"
    )
  })
  
  output$best_r2 <- renderValueBox({
    comparison <- model_comparison()
    best_r2 <- max(comparison$R_squared)
    valueBox(
      round(best_r2, 4),
      "Best R² Score",
      icon = icon("star"),
      color = "purple"
    )
  })
  
  # Bioactivity distribution
  output$bioactivity_dist <- renderPlotly({
    data <- full_data()
    
    p <- ggplot(data, aes(x = pIC50, fill = bioactivity_class)) +
      geom_histogram(bins = 40, alpha = 0.7, color = "black") +
      labs(
        title = "pIC50 Distribution",
        x = "pIC50",
        y = "Count",
        fill = "Bioactivity Class"
      ) +
      theme_minimal()
    
    ggplotly(p)
  })
  
  # Model comparison plot
  output$model_comparison_plot <- renderPlotly({
    comparison <- model_comparison()
    
    p <- ggplot(comparison, aes(x = reorder(Model, -R_squared), y = R_squared, fill = Model)) +
      geom_col(alpha = 0.8) +
      geom_text(aes(label = round(R_squared, 3)), vjust = -0.5) +
      labs(
        title = "Model R² Comparison",
        x = "Model",
        y = "R² Score"
      ) +
      theme_minimal() +
      theme(legend.position = "none")
    
    ggplotly(p)
  })
  
  # Data table with filters
  filtered_data <- eventReactive(input$apply_filters, {
    data <- full_data()
    
    if(input$bioactivity_filter != "All") {
      data <- data %>% filter(bioactivity_class == input$bioactivity_filter)
    }
    
    data <- data %>%
      filter(pIC50 >= input$pic50_range[1], pIC50 <= input$pic50_range[2])
    
    if(input$druglike_only) {
      data <- data %>% filter(druglike == TRUE)
    }
    
    data
  }, ignoreNULL = FALSE)
  
  output$data_table <- renderDT({
    data <- filtered_data() %>%
      select(molecule_chembl_id, canonical_smiles, pIC50, bioactivity_class, 
             MW, LogP, HBD, HBA, druglike)
    
    datatable(
      data,
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      filter = "top",
      rownames = FALSE
    )
  })
  
  # Chemical space plot
  output$chemical_space_plot <- renderPlotly({
    data <- full_data()
    
    p <- ggplot(data, aes_string(x = input$x_axis, y = input$y_axis, color = input$color_by)) +
      geom_point(alpha = 0.6, size = 2) +
      labs(
        title = "Chemical Space Analysis",
        x = input$x_axis,
        y = input$y_axis
      ) +
      theme_minimal()
    
    if(input$x_axis == "MW") {
      p <- p + geom_vline(xintercept = 500, linetype = "dashed", color = "red")
    }
    if(input$y_axis == "LogP") {
      p <- p + geom_hline(yintercept = 5, linetype = "dashed", color = "red")
    }
    
    ggplotly(p)
  })
  
  # Lipinski violations
  output$lipinski_violations <- renderPlot({
    data <- full_data()
    
    violations <- data %>%
      count(lipinski_violations) %>%
      mutate(lipinski_violations = as.factor(lipinski_violations))
    
    ggplot(violations, aes(x = lipinski_violations, y = n, fill = lipinski_violations)) +
      geom_col() +
      geom_text(aes(label = n), vjust = -0.5) +
      labs(
        title = "Lipinski Rule Violations",
        x = "Number of Violations",
        y = "Count"
      ) +
      theme_minimal() +
      theme(legend.position = "none")
  })
  
  # Model-specific predictions
  current_predictions <- reactive({
    switch(input$selected_model,
           "Random Forest" = rf_predictions(),
           "XGBoost" = xgb_predictions(),
           "Neural Network" = nn_predictions())
  })
  
  output$pred_vs_actual <- renderPlotly({
    preds <- current_predictions()
    
    p <- ggplot(preds, aes(x = Actual, y = Predicted)) +
      geom_point(alpha = 0.6, color = "steelblue") +
      geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
      geom_smooth(method = "lm", se = TRUE, color = "blue") +
      labs(
        title = paste(input$selected_model, "- Predicted vs Actual"),
        x = "Actual pIC50",
        y = "Predicted pIC50"
      ) +
      theme_minimal()
    
    ggplotly(p)
  })
  
  output$residual_plot <- renderPlotly({
    preds <- current_predictions()
    
    p <- ggplot(preds, aes(x = Predicted, y = Residuals)) +
      geom_point(alpha = 0.6, color = "steelblue") +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
      geom_smooth(se = TRUE, color = "blue") +
      labs(
        title = "Residual Plot",
        x = "Predicted pIC50",
        y = "Residuals"
      ) +
      theme_minimal()
    
    ggplotly(p)
  })
  
  output$metrics_table <- renderDT({
    comparison <- model_comparison()
    
    datatable(
      comparison,
      options = list(pageLength = 10),
      rownames = FALSE
    ) %>%
      formatRound(columns = c("MSE", "RMSE", "MAE", "R_squared"), digits = 4)
  })
  
  # System monitoring
  output$cpu_usage <- renderValueBox({
    cores <- parallel::detectCores()
    valueBox(
      cores,
      "CPU Cores",
      icon = icon("microchip"),
      color = "blue"
    )
  })
  
  output$memory_usage <- renderValueBox({
    mem <- Sys.meminfo()
    used_pct <- round((1 - as.numeric(mem$freeram) / as.numeric(mem$totalram)) * 100, 1)
    
    valueBox(
      paste0(used_pct, "%"),
      "Memory Usage",
      icon = icon("memory"),
      color = if(used_pct > 80) "red" else if(used_pct > 60) "yellow" else "green"
    )
  })
  
  output$gpu_status <- renderValueBox({
    gpu_avail <- torch::cuda_is_available()
    
    valueBox(
      if(gpu_avail) "Available" else "Not Available",
      "GPU Status",
      icon = icon("desktop"),
      color = if(gpu_avail) "green" else "red"
    )
  })
  
  output$system_info <- renderPrint({
    cat("System Information\n")
    cat("==================\n\n")
    cat("R Version:", R.version.string, "\n")
    cat("Platform:", R.version$platform, "\n")
    cat("CPU Cores:", parallel::detectCores(), "\n\n")
    
    mem <- Sys.meminfo()
    cat("Memory Information\n")
    cat("------------------\n")
    cat("Total RAM:", mem$totalram, "\n")
    cat("Free RAM:", mem$freeram, "\n")
    
    if(torch::cuda_is_available()) {
      cat("\nGPU Information\n")
      cat("---------------\n")
      cat("CUDA Available: YES\n")
      cat("GPU:", torch::cuda_get_device_name(0), "\n")
    }
  })
}

# Run the app
shinyApp(ui = ui, server = server)
