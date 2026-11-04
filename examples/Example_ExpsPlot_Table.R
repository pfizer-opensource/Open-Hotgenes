library(package = Hotgenes)

expression_data <- data.frame(
    Feature = rep(x = c("gene1", "gene2"), times = 2),
    Condition = rep(x = c("control", "treated"), each = 2),
    value = c(1, 2, 3, 4)
)
expression_plot <- ggplot2::ggplot(
    data = expression_data,
    mapping = ggplot2::aes(x = Condition, y = value)
)
ExpsPlot_Table(ggplot_obj = expression_plot, xVar = "Condition")
