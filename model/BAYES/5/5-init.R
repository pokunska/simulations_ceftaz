# Create Stan initial values
#
# This function must return something that can be passed to the `init` argument
#   of `cmdstanr::sample()`. There are several options; see `?cmdstanr::sample`
#   for details.
#
# `.data` represents the list returned from `make_standata()` for this model.
#   This is provided in case any of your initial values are dependent on some
#   aspect of the data (e.g. the number of rows).
#
# `.args` represents the list of attached arguments that will be passed through to
#   cmdstanr::sample(). This is provided in case any of your initial values are
#   dependent on any of these arguments (e.g. the number of chains).
#
# Note: you _don't_ need to pass anything to either of these arguments, you only
#   use it within the function. `bbr` will pass in the correct objects when it calls
#   `make_init()` under the hood.
#
make_init <- function(.data, .args) {
  function(){
    list(CLCHat = exp(rnorm(1, log(5.882),0.25)),
         QCHat  = exp(rnorm(1,  log(2.545),0.25)),
         V1CHat = exp(rnorm(1, log(10.64),0.25)),
         V2CHat = exp(rnorm(1, log(4.227),0.25)), 
         CLTHat = exp(rnorm(1, log(20.8),0.25)),
         QTHat  = exp(rnorm(1,  log(4.06),0.25)),
         V1THat = exp(rnorm(1, log(12.9),0.25)),
         V2THat = exp(rnorm(1, log(5.06),0.25)),
         nu = c(max(3,rgamma(1,2,0.1)),max(3,rgamma(1,2,0.1))),
         omega = rep(0.5,8) * exp(rnorm(8, 0, 0.25)),
         L = diag(8),
         sigma = exp(rnorm(2,log(0.2),0.25)),
         etaStd = matrix(0L, 8, .data$nSubjects))
  }
}
