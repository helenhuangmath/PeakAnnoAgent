\name{save_figure}
\alias{save_figure}
\title{Write a Vector Figure and PNG Preview}
\description{Render a plot to PDF and to a PNG preview in the .preview subdirectory.}
\usage{save_figure(draw, out_dir, name, width = 7, height = 4.5)}
\arguments{
  \item{draw}{A ggplot object or a zero-argument function drawing on the current graphics device.}
  \item{out_dir}{Output directory, created if necessary.}
  \item{name}{File basename without an extension.}
  \item{width, height}{Figure dimensions in inches.}
}
\value{A named character vector containing pdf and png paths.}
