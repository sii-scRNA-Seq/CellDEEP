\name{Subset.seurat.object}
\alias{Subset.seurat.object}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subset seurat object to test
}
\description{
This function allow you to easily subset every object to test them. First specify the ident to subset, then specify which one.
(For cluster 5, Ident.yo.subset = "cluster_id",Specific.Idents.to.subset = 5.)
}
\usage{
Subset.seurat.object(Seurat.object, Ident.to.subset, Specific.Idents.to.subset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Seurat.object}{
    A seurat object
}

  \item{Ident.to.subset}{
    Ident user want to split, eg: cluster_id.
}

  \item{Specific.Idents.to.subset}{
   Specify the ident to split, eg: 5.
}

}
\details{
  It's an optional function to use when prepare data(prepare_data).
}
\examples{
target.object <- Subset.seurat.object(Seurat.object, "cluster_id", 5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
