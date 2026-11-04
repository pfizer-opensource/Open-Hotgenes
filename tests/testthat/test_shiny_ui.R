test_that(
  desc = "conditional plot panels control their own visibility",
  code = {
    panels <- list(
      volcano = DEstats_UI_main_support(id = "DEstats"),
      scatter = ExpsPlot_UI_main_support(id = "ExpsPlot"),
      expression = BoxPlot_UI_main_support(id = "BoxPlot")
    )

    panel_html <- lapply(X = panels, FUN = function(panel = NULL) {
      htmltools::renderTags(x = panel)$html
    })

    for (html in panel_html) {
      expect_false(
        object = grepl(
          pattern = 'style="display: none;"',
          x = html,
          fixed = TRUE
        )
      )
    }

    expect_match(
      object = panel_html$volcano,
      regexp = "input.DE_Contrasts",
      fixed = TRUE
    )
    expect_match(
      object = panel_html$scatter,
      regexp = "input.yVar",
      fixed = TRUE
    )
    expect_match(
      object = panel_html$expression,
      regexp = "input.NormSlot",
      fixed = TRUE
    )

    pca_html <- htmltools::renderTags(
      x = PCA_UI_main_support(id = "PCA")
    )$html
    expect_true(
      object = grepl(
        pattern = 'id="PCA-label_size"[^>]*data-from="3"',
        x = pca_html,
        perl = TRUE
      )
    )
  }
)

test_that(
  desc = "PCA label size has a readable default",
  code = {
    expect_identical(
      object = formals(fun = Shiny_Hotgenes_Server)$selected_label_size,
      expected = 3
    )
  }
)
