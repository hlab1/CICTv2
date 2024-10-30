createDataObj <- function(expression_matrix, ground_truth=NULL, association_matrix=NULL, rf_features=NULL, rf_outputs=NULL, grn=NULL, in_data_obj=NULL, out_existing_obj=FALSE, suppress_warnings=FALSE) {
  # user requests to create new output object
  if (!out_existing_obj) {
    out_data_obj <- list(expression_matrix = expression_matrix,
                         ground_truth = ground_truth,
                         association_matrix = association_matrix,
                         rf_features = rf_features,
                         rf_outputs = rf_outputs,
                         grn = grn)
    return(out_data_obj)
  }
  # steps for if user requests to use inputted object as output
  # if(is.null(cict_data_obj)) {
  #   print('Error: Tried to output existing cict_data_obj but one was not given')
  #   NO_CICT_OBJ()
  # }
  # try just checking if it's a list of class "cict_data_obj"
}
