\name{clusterdef}
\alias{clusterdef}
\alias{clusterdef-class}
\title{clusterdef class}
\description{
  Defines parameters used to describe location and spread
  of clusters.
}
\details{
  The definition is split into two parts.  This separation
  allows us to define n clusters which depend on 4xm
  parameters (m is less than n), so that clusters with
  equal total copy number can share the same mean R, for
  example.
}

  \section{Cluster specific slots}{ \describe{

  \item{\code{a1}:}{Count of allele 1}

  \item{\code{a2}:}{Count of allele 2}

  \item{\code{R.index}:}{ index of the R parameters that
  define the location of this cluster}

  \item{\code{theta.index}:}{ index of the theta parameters
  that define the location of this cluster}

  } }

  \section{Parameter slots}{

  \describe{ \item{\code{R.mean}:}{ numeric vector, element
  i defines the mean of the gaussian distribution from
  which R values for clusters with R.index==i are selected}

  \item{\code{R.sd}:}{ numeric vector, element i defines
  the sd of the gaussian distribution from which R values
  for clusters with R.index==i are selected}

  \item{\code{theta.a}:}{ numeric vector, element i defines
  the shape 1 parameter of the beta distribution from which
  theta values for clusters with theta.index==i are
  selected}

  \item{\code{theta.b}:}{ numeric vector, element i defines
  the shape2 parameter of the beta distribution from which
  theta values for clusters with theta.index==i are
  selected} } }


