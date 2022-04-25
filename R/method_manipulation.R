#' @importFrom dplyr %>%
NULL

#' Get Information of Methods
#'
#' @param method One or more method names. Default is all.
#'
#' @return A list.
#' @importFrom purrr map
#' @importFrom stringr str_split
#' @importFrom stats setNames
#' @importFrom utils lsf.str
#' @importFrom crayon red
#' @export
#'
#' @examples
#' get_method <- get_method(method = c("splat", "PROSSTT"))
get_method <-function(method = "all"){
  env <- asNamespace("simmethods")
  func_name <- as.character(lsf.str(env, pattern = "_definition"))
  right_name <- stringr::str_split(func_name, pattern = "_", simplify = T)[,1]
  if(method[1] != "all"){
    purrr::map(method, function(z){
      # Make sure the method name is valid
      if(!z %in% right_name) stop(crayon::red(paste(z, "is not wrapped in simmethods package. Please check your spelling or use () to get the right name")))
    })
    func_name <- method
  }
  # Get the information of methods
  methods_return <- purrr::map(func_name, .f = function(x){
    fun <- paste0(x,"()") %>% parse(text = .) %>% eval()
  }) %>% setNames(right_name)
  methods_return
}
