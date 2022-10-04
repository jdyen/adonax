# function to rebase a factor variable
rebase_factor <- function(x, ...) {
  as.integer(as.factor(x))
}
