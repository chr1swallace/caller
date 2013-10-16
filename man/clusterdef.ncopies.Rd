\name{clusterdef.ncopies}
\alias{clusterdef.ncopies}
\title{clusterdef.ncopies}
\usage{
  clusterdef.ncopies(n,
    R.mean = c(`0` = -3.527211, `1` = -0.664184, `2` = 0, `3` = 0.395621, `4` = 0.678345),
    R.sd = c(`0` = 1.329152, `1` = 0.284338, `2` = 0.159645, `3` = 0.209089, `4` = 0.191579),
    eps = 0.01)
}
\arguments{
  \item{n}{maximum number of copies}

  \item{R.mean}{default values of LRR means for 0..4
  copies}

  \item{R.sd}{default values of LRR sds for 0..4 copies}

  \item{eps}{distance from 0 or 1 to assign theta values
  with theoretical values of 0 or 1}
}
\value{
  a clusterdef object
}
\description{
  Create a clusterdef object, specifying only the maximum
  number of total copies observed
}
\details{
  Default values for LRR means and sd are taken from ???.
  The LRR_SD is quite dependent on data quality and
  normalisation but the LRR-Mean should be fairly good.

  LRR_mean: Copy 0 1 2 3 4+ -3.527211 -0.664184 0.000000
  0.395621 0.678345 LRR_sd: Copy 0 1 2 3 4+ 1.329152
  0.284338 0.159645 0.209089 0.191579
}
\author{
  Chris Wallace
}

