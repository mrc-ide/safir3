#' safir: squire and friends individual rewrite
#'
#' @description
#' safir is an individual-based simulation model of COVD-19 based on squire.
#'
#' @docType package
#' @name safir
#'
#' @importFrom stats setNames
#' @importFrom utils getFromNamespace
#'
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib safir
## usethis namespace: end
NULL

execute_any_process <- utils::getFromNamespace("execute_any_process", "individual")
execute_process <- utils::getFromNamespace("execute_process", "individual")
