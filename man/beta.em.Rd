\name{beta.em}
\alias{beta.em}
\title{beta.em}
\usage{
  beta.em(df.in, theta, lrr, tol = 0.01, verbose = FALSE,
    maxit = 10000, eps = 0.01, use.deriv = FALSE)
}
\arguments{
  \item{df.in}{clusterdef object}

  \item{theta}{theta}

  \item{lrr}{LRR}

  \item{tol}{how small a change in group membership
  probabilities prompts continued optimization}

  \item{verbose}{print messages if TRUE}

  \item{maxit}{maximum number of iterations}

  \item{eps}{theta cannot be exactly 0 or 1, and values
  less than eps from 0 or 1 are set to eos and 1-eps
  respectively}

  \item{use.deriv}{if TRUE, use gradients to optimize
  quicker.  Currently off by default because it DOESN'T
  WORK!}
}
\value{
  a clusterfit object
}
\description{
  Fit a beta mixture distribution to theta, assuming LRR
  follows an established Gaussian mixture distribution
}
\author{
  Chris Wallace
}

