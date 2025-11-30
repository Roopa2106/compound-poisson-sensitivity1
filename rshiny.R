# Parameters
lambda <- 3   # arrival rate
mu     <- 2   # rate of exponential jumps
nsim   <- 10000  # number of simulations per t
times  <- c(10, 100, 1000, 10000)

# Function to simulate S(t)
simulate_S <- function(t, lambda, mu, nsim) {
  N <- rpois(nsim, lambda * t)      # number of jumps in each path
  S <- numeric(nsim)
  for (i in seq_len(nsim)) {
    n_i <- N[i]
    if (n_i > 0) {
      S[i] <- sum(rexp(n_i, rate = mu))
    } else {
      S[i] <- 0
    }
  }
  S
}

# Plot histograms in a 2x2 layout
par(mfrow = c(2, 2))
for (t in times) {
  S_t <- simulate_S(t, lambda, mu, nsim)
  hist(S_t,
       breaks = 50,
       main  = paste("Histogram of S(t) at t =", t),
       xlab  = "S(t)",
       col   = "lightblue",
       border = "white")
}
par(mfrow = c(1, 1))

# app.R
library(shiny)
library(ggplot2)
library(dplyr)

# Function to simulate S(t) for a grid of times
simulate_S_grid <- function(times, lambda, mu, nsim) {
  all <- lapply(times, function(t) {
    N <- rpois(nsim, lambda * t)
    S <- numeric(nsim)
    for (i in seq_len(nsim)) {
      n_i <- N[i]
      if (n_i > 0) S[i] <- sum(rexp(n_i, rate = mu))
    }
    data.frame(t = t, S = S)
  })
  bind_rows(all)
}

ui <- fluidPage(
  titlePanel("Compound Poisson Process: S(t) Sensitivity"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Model: S(t) = sum_{i=1}^{N(t)} X_i, N(t) ~ Poisson(lambda * t), X_i ~ Exp(mu)."),
      
      sliderInput("lambda", "Arrival rate lambda:", 
                  min = 0.1, max = 10, value = 3, step = 0.1),
      sliderInput("mu", "Jump rate mu:", 
                  min = 0.1, max = 10, value = 2, step = 0.1),
      
      sliderInput("tmax", "Maximum time t:", 
                  min = 10, max = 10000, value = 1000, step = 10),
      sliderInput("ntimes", "Number of time points:", 
                  min = 5, max = 50, value = 20),
      
      sliderInput("nsim", "Simulations per time:", 
                  min = 500, max = 10000, value = 2000, step = 500),
      
      sliderInput("t_hist", "Time for histogram:", 
                  min = 10, max = 10000, value = 100, step = 10),
      
      actionButton("run", "Run simulation")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Summary vs Time", plotOutput("summaryPlot", height = "450px")),
        tabPanel("Histogram",      plotOutput("histPlot",    height = "450px")),
        tabPanel("Table",          dataTableOutput("summaryTable"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  sim_res <- eventReactive(input$run, {
    times <- seq(10, input$tmax, length.out = input$ntimes)
    times <- sort(unique(round(times)))
    dat   <- simulate_S_grid(times, input$lambda, input$mu, input$nsim)
    
    summary <- dat %>%
      group_by(t) %>%
      summarise(
        mean = mean(S),
        sd   = sd(S),
        p10  = quantile(S, 0.10),
        p50  = median(S),
        p90  = quantile(S, 0.90),
        .groups = "drop"
      )
    
    list(raw = dat, summary = summary)
  })
  
  output$summaryPlot <- renderPlot({
    s <- sim_res()$summary
    ggplot(s, aes(x = t)) +
      geom_ribbon(aes(ymin = p10, ymax = p90), alpha = 0.25) +
      geom_line(aes(y = mean), colour = "blue", size = 1) +
      geom_line(aes(y = p50),  colour = "darkgreen", linetype = "dashed") +
      labs(title = "E[S(t)] and 10–90% band vs t",
           x = "t", y = "S(t)") +
      theme_minimal()
  })
  
  output$histPlot <- renderPlot({
    dat <- sim_res()$raw
    # pick the simulated t closest to requested t_hist
    t_vals <- sort(unique(dat$t))
    t_sel  <- t_vals[which.min(abs(t_vals - input$t_hist))]
    d_sub  <- dplyr::filter(dat, t == t_sel)
    
    ggplot(d_sub, aes(x = S)) +
      geom_histogram(bins = 50, colour = "white") +
      labs(title = paste("Histogram of S(t) at t =", t_sel),
           x = "S(t)", y = "Frequency") +
      theme_minimal()
  })
  
  output$summaryTable <- renderDataTable({
    sim_res()$summary
  })
}

shinyApp(ui, server)

