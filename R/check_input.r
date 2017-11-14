#######################################
### helper functions to check input ###
#######################################

#' checks input for functions in the bigleaf package
#' 
#' @description Checks length and type of the provided input variables.
#' 
#' @param data   Input data.frame or matrix
#' @param vars   Input variables. MUST be a list, provided as list(...)
#' 
#' 
check.input <- function(data,vars){
  
  check.length(vars)
  
  arg       <- deparse(substitute(vars))
  var_names <- trimws(unlist(strsplit(substr(arg,6,nchar(arg)-1),",")))
  
  for (i in seq_along(vars)){
    assign(var_names[i],check.columns(data,vars[[i]],var_names[i]),envir=sys.frame(-1))
  } 
  
}


#' checks columns in a data frame 
#' 
#' @description check if columns in a data.frame or matrix (or vector variables)
#'              exist and if they are of the right type (numeric).
#' 
#' @param data     Input data.frame or a matrix
#' @param var      Input variable 
#' @param var_name Name of the input variable
#' 
#' @details The Input variable can be of type character or numeric. If character,
#'          it is interpreted as the column name of the input data.frame / matrix.
#'          If numeric, it is interpreted as a direct input to the function.
#'          In both cases, the requirements of the input is checked (E.g. does
#'          the column exist? Is the input numeric? etc.)
check.columns <- function(data,var,var_name){
  
  if (missing(var_name)){
    var_name <- ""
  }
  
  if (is.character(var)){
    if (!missing(data)){
      if (length(var) == 1){
        if (var %in% colnames(data)){
          var <- data[,var]
          if (is.numeric(var)){
            return(unname(var))
          } else {
            stop("variable '",var_name,"' must be numeric",call.=FALSE)
          }
        } else {
          stop ("variable '",var_name,"' does not exist in the input matrix/data.frame",call.=FALSE)
        }
      } else {
        stop("name of variable '",var_name,"' must have length 1",call.=FALSE)
      } 
    } else {
      if ("data" %in% names(formals(sys.function(which=-2)))){
        stop("variable '",var_name,"' is of type character and interpreted as a column name, but no input matrix/data.frame is provided. Provide '",var_name,"' as a numeric vector, or an input matrix/data.frame with a column named '",var_name,"'",call.=FALSE)
      } else {
        stop("variable '",var_name,"' must be numeric",call.=FALSE)
      }
    }
  } else {
    if (!missing(data)){
      if (is.numeric(var) & length(var) == nrow(data)){
        return(unname(var))
      } else if (is.numeric(var) & length(var) != nrow(data)) {
        if (length(var) == 1){
          var <- rep(var,length=nrow(data))
          return(unname(var))
        } else {
          stop("variable '",var_name,"' must have the same number of observations as the input matrix/data.frame or length 1",call.=FALSE)
        }
      } else if (!is.numeric(var)){
        stop("variable '",var_name,"' must be numeric",call.=FALSE)
      }
    } else {
      if (is.numeric(var)){
        return(unname(var))
      } else {
        stop("variable '",var_name,"' must be numeric",call.=FALSE)
      } 
    } 
  }
}


#' tests if variables have the same length
#' 
#' @param input_list List of variables for which the length has to be compared
#' 
#' @note This function only plays a role if no input data.frame or matrix are 
#'       provided. In this case it ensures that provided vectors have the same
#'       length to avoid trouble later up the function call.
#'       
check.length <- function(input_list){
  
  length.vars <- sapply(input_list,length)
  
  if (length(unique(length.vars)) >= 2){
    if (unique(length.vars)[1] != 1 | length(unique(length.vars)) > 2){
      stop("All input variables must have the same length or have a length of 1!",call.=FALSE)
    }
  }
}

