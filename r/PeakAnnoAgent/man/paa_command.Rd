\name{paa_command}
\alias{paa_command}
\alias{paa_check}
\alias{paa_version}
\title{Find and Check the peakanno-agent Program}
\description{\code{paa_command()} returns the command used to run the Python program,
\code{paa_check()} prints its installation check, and \code{paa_version()} returns its version.}
\usage{
paa_command()
paa_check(offline = TRUE)
paa_version()
}
\arguments{
  \item{offline}{Check local dependencies only, without contacting public services.}
}
\value{\code{paa_command()}: a character vector, the program followed by any leading arguments.
\code{paa_check()}: \code{TRUE} invisibly if all required dependencies are present.
\code{paa_version()}: the version string.}
\details{The command is looked up in this order: the R option \code{peakannoagent.command}
(a character vector such as \code{c("/path/to/python", "-m", "peakannoagent")}), the environment
variable \code{PEAKANNOAGENT_COMMAND}, \code{peakanno-agent} on the PATH, and
\code{python3 -m peakannoagent}.}
\examples{
\dontrun{
options(peakannoagent.command = c("/opt/venv/bin/python", "-m", "peakannoagent"))
paa_version()
paa_check()
}
}
