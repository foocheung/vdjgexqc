#' Run the Shiny Application
#'
#' @param ... Arguments passed to [shiny::shinyApp()]
#'
#' @export
run_app <- function(...) {
  options(shiny.maxRequestSize = 5000 * 1024^2)

  shiny::shinyApp(
    ui = app_ui,
    server = app_server,
    ...
  )
}
