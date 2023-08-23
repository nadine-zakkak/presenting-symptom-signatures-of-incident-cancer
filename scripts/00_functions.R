# Nadine Zakkak
# Script that contains helper functions required for analysis and output of figures and tables

format_numb <- function(numb) {
  # Format number into 2 decimal places if <10, 1 decimal place if <100, 0 decimal places if >=100
  numb = ifelse(numb<10,
                sprintf(numb, fmt='%.2f'),
                ifelse(numb<100,
                       sprintf(numb, fmt = '%.1f'),
                       sprintf(numb, fmt = '%1.0f')))
}

calculate_prop <- function(n, N){
  # Calculates proportion
  # n = sample size, N = total population 
  n/N
}

calculate_ci <- function(n=NULL, prop=NULL, N){
  # Calculate 95% confidence intervals, using Wilson's method
  #  n = sample size, prop = proportion, N = total population 
  # proportion not in percentage format
  # Returns a list with lower and upper boundaries of confidence interval
  if(is.null(prop) & !is.null(n)) {
    prop = n/N
  }
  lower = (1/(1+(qnorm(0.975)^2)/N))*((prop)+(qnorm(0.975)^2)/(2*N)) - (qnorm(0.975)/(1+((qnorm(0.975)^2)/N))*sqrt((prop)*(1-prop)/N + (qnorm(0.975)^2)/(4*(N^2))))
  upper = (1/(1+(qnorm(0.975)^2)/N))*((prop)+(qnorm(0.975)^2)/(2*N)) + (qnorm(0.975)/(1+((qnorm(0.975)^2)/N))*sqrt((prop)*(1-prop)/N + (qnorm(0.975)^2)/(4*(N^2))))
  ls <- list(lower, upper)
  names(ls) <- c("lower", "upper")
  return(ls)
}
