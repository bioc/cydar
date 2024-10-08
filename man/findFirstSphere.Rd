\name{findFirstSphere}
\alias{findFirstSphere}

\title{Identifies the first non-redundant hyperspheres}
\description{Tests whether each hypersphere is not redundant to (i.e., lies more than a threshold distance away from) another hypersphere with a lower p-value.}

\usage{
findFirstSphere(x, pvalues, threshold=1, block=NULL, num.threads=1)
}

\arguments{
\item{x}{A numeric matrix of hypersphere coordinates (median locations for all markers), where rows correspond to hyperspheres and columns correspond to markers.

Alternatively, a CyData object containing median intensities for groups of cells, such as that produced by \code{\link{countCells}}.}
\item{pvalues}{A numeric vector of p-values, one for each row (i.e., hypersphere) of \code{x}.}
\item{threshold}{A numeric scalar specifying the maximum distance between the locations of two redundant hyperspheres.}
\item{block}{A factor specifying which hyperspheres belong to which block, where non-redundant hyperspheres are identified within each block.}
\item{num.threads}{Integer scalar specifying the number of threads to use.}
}

\details{
This function iterates across the set of hyperspheres, typically ordered by decreasing significance.
It will tag a hypersphere as being redundant if its location lies within \code{threshold} of the location of a higher-ranking hypersphere in all dimensions.
In this manner, the set of all DA hyperspheres can be filtered down to a non-redundant subset that is easier to interpret.

Note that the criterion for redundancy mentioned above is equivalent to a Chebyshev distance, rather than Euclidean.
This is easier to interpret, especially given that the median intensity is defined separately for each marker.
Unlike in \code{\link{countCells}}, the threshold is not scaled by the number of markers because each hypersphere location is computed as an average across cells.
This means that there is generally no need to account for extra distance due to noise between cells.

The default \code{threshold} of unity assumes that the intensities have been transformed to or near a log10 scale.
It means that one hypersphere must vary from another by at least one log10-unit (i.e., a 10-fold change in intensity) in at least one marker to be considered non-redundant.
This avoids reporting many hyperspheres that differ from each other by relatively small, uninteresting shifts in intensity.
Greater resolution can be obtained by decreasing this value, e.g., to 0.5.

If \code{block} is set, non-redundant hyperspheres are only identified within each block (i.e., a hypersphere cannot be redundant to hyperspheres in different blocks).
For example, one can set \code{block} to the sign of the log-fold change.
This ensures that hyperspheres changing in one direction are not considered redundant to those changing in another direction.
By default, all hyperspheres are considered to be part of the same block.

% Does the spatial FDR control the FDR across the reported non-redundant locations?
% In general, using spatialFDR() is fine; the non-redundant locations can be considered as centres of non-overlapping hypercubic partitions of side length "threshold".
% Thus, controlling the spatial FDR should be similar to controlling the FDR across these partitions.
%       Of course, there are some edge cases where this reasoning could fail, e.g., with diffuse DA subpopulations and tight non-DA subpopulations.
% If both subpopulations are detected, the spatial FDR may be low, because the former supplies a lot of true positive volume.
% However, the FDR across the non-redundant hyperspheres might be high if each subpopulation is represented by a single non-redundant hypersphere.
%       Ideally, we would like to control the FDR across the non-redundant hyperspheres directly.
% However, this is not possible because the definition of non-redundancy depends on the p-value.
% Applying the BH method to the non-redundant p-values would fail to account for the implicit multiple testing in the surrounding space.
}

\value{
A logical vector indicating whether each of the hyperspheres in \code{x} is non-redundant.
}

\author{
    Aaron Lun
}

\examples{
# Mocking up some data.
coords <- matrix(rnorm(10000, 2, sd=0.3), nrow=1000)
pval <- runif(1000)
logfc <- rnorm(1000)

# Keep most significant non-redundant ("first") hyperspheres.
findFirstSphere(coords, pval)

# Block on the sign of the log-fold change.
findFirstSphere(coords, pval, block=sign(logfc))
}
