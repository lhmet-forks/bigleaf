########################
### helper functions ###
########################

#' checks columns in a data frame 
#' 
#' @description check if columns in a data.frame or matrix (or vector variables) 
#'              exist and if they are of the right type (numeric).
#' 
#' @param data         a data.frame or a matrix
#' @param column_name  the column name
#' 
#' 
#' 
check.columns <- function(data,column_name){
  
  if (is.character(column_name)){
    if (!missing(data)){
      if (length(column_name) == 1){
        if (column_name %in% colnames(data)){
          var <- data[,column_name]
          if (is.numeric(var)){
            return(unname(var))
          } else {
            stop("variable '",column_name,"' must be numeric",call.=FALSE)
          }
        } else {
          stop ("variable '",column_name,"' does not exist in the input matrix/data.frame",call.=FALSE)
        }
      } else {
        stop("name of variable '",deparse(substitute(column_name)),"' must have length 1",call.=FALSE)
      } 
    } else {
      if ("data" %in% names(formals(sys.function(which=-1)))){
        stop("variable '",deparse(substitute(column_name)),"' is of type character and interpreted as a column name, but no input matrix/data.frame is provided. Provide '",column_name,"' as a numeric vector, or an input matrix/data.frame with a column named '",column_name,"'",call.=FALSE)
      } else {
        stop("variable '",deparse(substitute(column_name)),"' must be numeric",call.=FALSE)
      }
    }
  } else {
    if (!missing(data)){
      if (is.numeric(column_name) & length(column_name) == nrow(data)){
        return(unname(column_name))
      } else if (is.numeric(column_name) & length(column_name) != nrow(data)) {
        stop("variable '",deparse(substitute(column_name)),"' must have the same number of observations as the input matrix/data.frame",call.=FALSE)
      } else if (!is.numeric(column_name)){
        stop("variable '",column_name,"' must be numeric",call.=FALSE)
      }
    } else {
      if (is.numeric(column_name)){
        return(unname(column_name))
      } else {
        stop("variable '",column_name,"' must be numeric",call.=FALSE)
      } 
    } 
  }
}


#' tests if variables have the same length
#' 
#' 
check.length <- function(...){
  
  args <- list(...)
  length.args <- sapply(args,length)
  
  if (length(unique(length.args)) > 1){
    stop("all input variables must have the same length",call.=FALSE)
  }
}

