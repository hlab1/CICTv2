# TODO: determine if needed

#' Ask user to confirm whether to proceed
#'
#' Asks user whether they wish to proceed and prompts user to enter `y` or `n`.
#' If the user enters `y`, return `TRUE`; if the user enters `n`, return
#' `FALSE`; if the user enters an invalid input, prompt again.
#'
#' @return A boolean value based on the user's input.
#'
#' @examples TODO
askUserProceed <- function() {
  input <- ""
  while(input != "y" && input != "n") {
    input = readline(prompt = "Do you want to proceed? (y/n)")
  }
  if(input == "y") {
    return(TRUE)
  }
  return(FALSE)
}
