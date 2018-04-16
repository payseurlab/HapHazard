
# A library of functions for analysing haphazard data

#FUNCTION: calculate the mean block length and 95% CI for a given chromosome, generation, and deme
conf_interval <- function(x,c)
{
  y <- c()
  
  x <- sort(x)
  n <- length(x)
  upper_index <- floor( n - ( n * (1 - c) / 2 ) )
  lower_index <- floor( n * ( 1 - c ) / 2 )
  
  y <- c(x[lower_index], x[upper_index])
  
  return(y)
}
