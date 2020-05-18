get_axis_labels = function(plt){
    if(!requireNamespace('ggplot2')) stop("Install ggplot2 >= 3.0.0")

    b = ggplot2::ggplot_build(plt)$layout$panel_params[[1]]
    if(packageVersion('ggplot2') >='3.3.0'){
        list(xlabels =  b$x$get_labels(), ylabels = b$y$get_labels())
    } else if(packageVersion('ggplot2') >= '3.0'){
        list(xlabels = b$x.label, ylabels = b$y.label)
    } else{
        stop("Install ggplot2 >= 3.0.0")
    }
}


if(requireNamespace('ggplot2')) {
    ScaleAxisColor = ggplot2::ggproto('ScaleAxisColor', ggplot2::ScaleDiscrete)
}

replace_empty = function(x, y) if(length(x) == 0) y else x

#' Color axis labels
#'
#' @param plt [ggplot2::ggplot()] object
#' @param label_data_x
#' @param label_data_y [data.frame()]s containing the mapping between x-axis/y-axis labels and
#' `aes_label`
#' @param aes_label character or bare symbol giving the column in `label_data` to be mapped
#' @param scale `ggplot2` discrete color
#' @return `plt` with axis text modified
#' @export
#'
#' @examples
#' require(ggplot2)
#' require(dplyr)
#' plt = ggplot(mpg, aes(x = manufacturer, y = drv)) + geom_jitter()
#' label_data = mpg %>% select(manufacturer) %>% unique() %>%
#' mutate(euro = manufacturer %in% c('audi', 'volkswagen'), drv = NA_character_)
#' map_axis_labels(plt, label_data, label_data, euro)
map_axis_labels = function(plt, label_data_x,  label_data_y, aes_label, scale = ggplot2::scale_color_hue(aesthetics = 'axis_color')){
    actual_labs = get_axis_labels(plt)
    aes_label_str = rlang::as_name(rlang::enquo(aes_label))
    label_data = bind_rows(label_data_x, label_data_y)
    data_x = hushWarning(left_join(tibble(X_ = actual_labs$xlabels), label_data_x, by = c(X_ = rlang::as_name(plt$mapping$x))), 'coercing into character')
    if(nrow(data_x) != length(actual_labs$xlabels)) stop("Bad join on xlabels!")
    data_y = hushWarning(left_join(tibble(Y_ := actual_labs$ylabels), label_data_y,  by = c(Y_ = rlang::as_name(plt$mapping$y))), 'coercing into character')
    if(nrow(data_y) != length(actual_labs$ylabels)) stop("Bad join on ylabels!")
    scale$train(label_data[[aes_label_str]])
    data_x$COLOR = replace_empty(scale$map(data_x[[aes_label_str]]),  'black')
    data_y$COLOR = replace_empty(scale$map(data_y[[aes_label_str]]), 'black')
    plt + hushWarning(ggplot2::theme(axis.text.x = ggplot2::element_text(color = data_x$COLOR), axis.text.y = ggplot2::element_text(color = data_y$COLOR)), 'Vectorized input')
}

