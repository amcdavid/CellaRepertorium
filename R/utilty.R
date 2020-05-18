##' Selectively muffle warnings based on output
##'
##' @param expr an expression
##' @param regexp a regexp to be matched (with str_detect)
##' @return the result of expr
##' @examples
##' CellaRepertorium:::hushWarning(warning('Beware the rabbit'), 'rabbit')
##' CellaRepertorium:::hushWarning(warning('Beware the rabbit'), 'hedgehog')
hushWarning <- function(expr, regexp){
    withCallingHandlers(expr, warning=function(w){
        if(grepl(regexp, conditionMessage(w))) invokeRestart("muffleWarning")
    })
}
