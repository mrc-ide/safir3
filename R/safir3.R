#' safir3: squire and friends individual rewrite
#'
#' @description
#' safir3 is an individual-based simulation model of COVD-19 based on squire.
#'
#' @docType package
#' @name safir3
#'
#' @importFrom stats setNames
#' @importFrom utils getFromNamespace
#'
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib safir3
## usethis namespace: end
NULL

execute_any_process <- utils::getFromNamespace("execute_any_process", "individual")
execute_process <- utils::getFromNamespace("execute_process", "individual")
