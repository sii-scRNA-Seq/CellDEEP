\name{prepare_data}
\alias{prepare_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjust input data format
}
\description{
This function adjust input data metadata to add sample_id, cluster_id and group_id, can do normalize,scale and UMAP if needed. If data is simulated by muscat package, will create a metadata with ground truth.
}
\usage{
prepare_data <- function(obj, assay = "RNA",
                         sample_id, group_id, cluster_id,
                         Ident.to.subset = NULL,
                         Need.to.Norm_Transf = FALSE,
                         Need.to.Scale = FALSE,
                         Need.UMAP.generation = FALSE,
                         Ground.truth = FALSE，
                         Specific.Idents.to.subset = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a seurat object}
  \item{assay}{which assay to be used }
  \item{sample_id}{the metadata that include sample(replicate) id information, necessary}
  \item{group_id}{the metadata that include group information, necessary}
  \item{cluster_id}{the metadata that include cluster information, necessary,eg: seurat_cluster }
  \item{Ident.to.subset}{The ident to subset when user want to subset the object}
  \item{Specific.Idents.to.subset}{The ident to subset when user want to subset the object}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
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
