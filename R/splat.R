# splat_estimation <- function(){
#
#   requireNamespace("splatter")
#
#   # Set temporary files to record the error
#
#
#   # Read .h5 files to obtain the reference data
#
#
#   # Estimation
#   estimation_result <- splatter::splatEstimate(ref_data)
#
#   # Extract parameter
#   params_return <- get_S4_values(estimation_result) %>% setNames(slotNames(estimation_result))
#
#   # Save as .h5 file
#
#   # Return to user
#
# }
#
#
# get_S4_values <- function(S4Object){
#
#   purrr::map(slotNames(S4Object), function(x){
#
#     exp <- paste0(quote(S4Object),"@", x)
#
#     eval(parse(text = exp))
#
#   })
# }




