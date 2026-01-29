# Helpers Functions for plotting Columbus results

## Layout Plotting Functions
plot_plate <-
  function(table,
           property,
           discrete_property = T,
           legend,
           title) {
    
    layout_plot <- ggplot(table,
                          aes(
                            x = column,
                            y = row,
                            fill = {{property}},
                          )) + 
      geom_tile(color = "grey30") +
      scale_y_reverse(breaks = 1:16, labels = LETTERS[1:16]) +
      scale_x_continuous(breaks = 0:24) +
      coord_equal() +
      xlab("Column") +
      ylab("Row") +
      ggtitle(title)
    
    if (discrete_property) {
      layout_plot + scale_fill_tableau(name = legend)
    }
    else{
      layout_plot + scale_fill_viridis_c(legend)
    }
  }

## Density Plotting Functions
plot_density <-
  function(table,
           property,
           property_legend,
           tint,
           tint_legend,
           title) {
    
    density_plot <- ggplot(table,
                           aes(
                             x = {{property}},
                             color = {{tint}},
                             fill =  {{tint}}
                           ))
    
    density_plot + geom_density(alpha = 0.3) +
      scale_x_continuous(property_legend, trans = "log2") +
      scale_color_viridis_d(name = tint_legend) +
      scale_fill_viridis_d(name = tint_legend) +
      facet_grid(vars(factor(ifn)), vars(probe_set)) +
      ggtitle(title)
    
  }

plot_crossbars <- function(table,
                           property,
                           property_legend,
                           title) {
  

    crossbar_plot <- ggplot(table, aes(x = .data$sgRNA,
                                       y = {{property}},
                                       color = .data$Dox))
    
    crossbar_plot + geom_jitter(alpha = 0.5,
                                width = 0.3) +
                    stat_summary(fun.data = 'mean_cl_boot',
                                 geom = "crossbar") +
                    scale_color_tableau() +
                    ylab(property_legend) +
                    xlab("sgRNA") +
                    ggtitle(title)
  }