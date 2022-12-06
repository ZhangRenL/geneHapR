\name{windows/GUI_geneHapR}
\alias{startGUI.geneHapR}
\alias{build_app_geneHapR}
\usage{
build_app_geneHapR(...)
startGUI.geneHapR(...)
}
\title{the GUI interface of geneHapR}
\arguments{
\item{...}{not used}
}
\description{
the GUI interface of geneHapR
}
\examples{
if (interactive() & tools:::.OStype() == "windows") {
   options(device.ask.default = FALSE)
   app <- build_app_geneHapR()
   startGUI.geneHapR()
}
}
