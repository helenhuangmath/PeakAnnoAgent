#!/usr/bin/env Rscript
# JSON bridge between the Python agent and the R analysis engine.
#
#   Rscript scripts/peak_tools.R <tool> <params.json> <result.json>
#
# Always writes <result.json> as {"status": "ok"|"error", "result": ..., "error": ...}
# and exits 0 unless the arguments themselves are invalid.

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "scripts/peak_tools.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

for (f in c("utils.R", "figures.R", "annotation.R", "enrichment.R", "qc.R", "chippeakanno.R",
  "motifs.R", "networks.R", "agent_tools.R")) {
  source(file.path(repo_root, "R", f))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  cat("Usage: Rscript scripts/peak_tools.R <tool> <params.json> <result.json>\n")
  cat("Tools:", paste(agent_tool_names(), collapse = ", "), "\n")
  quit(status = 2)
}
require_pkg("jsonlite")
tool <- args[[1]]
params <- jsonlite::fromJSON(args[[2]], simplifyVector = TRUE)
result_path <- args[[3]]

started <- Sys.time()
captured_warnings <- character()
outcome <- tryCatch(
  withCallingHandlers(
    list(status = "ok", result = run_agent_tool(tool, params)),
    warning = function(w) {
      captured_warnings <<- c(captured_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ),
  error = function(e) list(status = "error", error = conditionMessage(e))
)
outcome$tool <- tool
outcome$elapsed_sec <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1)
if (length(captured_warnings)) outcome$r_warnings <- utils::head(unique(captured_warnings), 10)

jsonlite::write_json(outcome, result_path, auto_unbox = TRUE, null = "null", na = "null",
  digits = NA, dataframe = "rows", pretty = TRUE)
