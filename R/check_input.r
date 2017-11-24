#######################################
### helper functions to check input ###
#######################################

#' checks input for functions in the bigleaf package
#' 
#' @description Checks length and type of the provided input variables.
#' 
#' @param data   Input data.frame or matrix
#' @param ...    Input variables. Either a list or individual vectors
#' 
#' 
check.input <- function(data,...){

  check.length(list(...))
  
  if (missing(data)){
    data <- NULL
  }
  
  varlist  <- match.call()[-c(1,2)]
  varnames <- c(unlist(sapply(varlist,as.character)))
  varnames <- varnames[!varnames %in% c("c","list")]

  for (i in seq_along(varnames)){
    assign(varnames[i],check.columns(data,varnames[i]),envir=sys.frame(-1))
  } 
  
}


#' checks columns in a data frame 
#' 
#' @description check if columns in a data.frame or matrix (or vector variables)
#'              exist and if they are of the right type (numeric).
#' 
#' @param data     Input data.frame or a matrix
#' @param varname  Input variable. Must be a character 
#' 
#' @details The Input variable 'varname' must be of type character. The function
#'          searches for the object 'varname' in the parent environment and returns it.
#'          The function then ensures that the requirements of the input variable (E.g. does
#'          the column exist? Is the input numeric? etc.) is met.
#'          
check.columns <- function(data,varname){
  
  n <- -3
  
  var <- get0(varname,envir=sys.frame(n),ifnotfound="notfound")
  
  if (length(var) < 2){
    if (is.null(var)){
      return()
    } else if (is.na(var)){
      return()
    } else if (var == "notfound"){  # required for standalone check.columns() call
      var <- varname
    }
  }
  
  if (is.character(var)){
    if (!missing(data) & !is.null(data)){
      if (length(var) == 1){
        if (var %in% colnames(data)){
          var <- data[,var]
          if (is.numeric(var)){
            return(unname(var))
          } else {
            stop("variable '",var,"' must be numeric",call.=FALSE)
          }
        } else {
          stop ("variable '",var,"' does not exist in the input matrix/data.frame",call.=FALSE)
        }
      } else {
        stop("name of variable '",var,"' must have length 1",call.=FALSE)
      } 
    } else {
      if ("data" %in% names(formals(sys.function(which=n)))){
        if (var %in% as.character(unlist(match.call(definition=sys.function(n),call=sys.call(n))[-1]))){
          stop("variable '",var,"' is of type character and interpreted as a column name, but no input matrix/data.frame is provided. Provide '",var,"' as a numeric vector, or an input matrix/data.frame with a column named '",var,"'",call.=FALSE)
        } else {
          stop("variable '",var,"' is not provided",call.=FALSE)
        }
      } else {
        stop("variable '",varname,"' must be numeric",call.=FALSE)
      }
    }
  } else {
    if (!missing(data) & !is.null(data)){
      if (is.numeric(var) & length(var) == nrow(data)){
        return(unname(var))
      } else if (is.numeric(var) & length(var) != nrow(data)) {
        if (length(var) == 1){
          var <- rep(var,length=nrow(data))
          return(unname(var))
        } else {
          stop("variable '",varname,"' must have the same length as the input matrix/data.frame or length 1. Do NOT provide an input matrix/data.frame if none of its variables are used!",call.=FALSE)
        }
      } else if (!is.numeric(var)){
        stop("variable '",varname,"' must be numeric",call.=FALSE)
      }
    } else {
      if (is.numeric(var)){
        return(unname(var))
      } else {
        stop("variable '",varname,"' must be numeric",call.=FALSE)
      } 
    } 
  }
}


#' tests if variables have the same length
#' 
#' @param varlist List of variables for which the length has to be compared
#' 
#' @note This function only plays a role if no input data.frame or matrix are 
#'       provided. In this case it ensures that provided vectors have the same
#'       length to avoid trouble later up the function call.
#'       
check.length <- function(varlist){
  
  if (is.list(unlist(varlist,recursive=FALSE))){
    varlist <- unlist(varlist,recursive=FALSE)
  }
  
  length.vars <- sapply(varlist,length)
  length.vars <- length.vars[length.vars > 0]
  
  if (length(unique(length.vars)) >= 2){
    if (sort(unique(length.vars))[1] != 1 | length(unique(length.vars)) > 2){
      stop("All input variables must have the same length or a length of 1!",call.=FALSE)
    }
  }
}
