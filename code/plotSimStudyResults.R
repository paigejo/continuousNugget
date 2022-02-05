# function for making fancy colord table of simulation study results. 
# Refer to 
# https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_pdf.pdf
# for information about kableExtra
makeFancyTable = function(meanScoresDF) {
  # meanScoresDF contains the following variables:
  # beta, rho, nClustFac, nEAsFac, rmse, bias, 
  # intScore80, intScore90, intScore95, 
  # cvg89, cvg90, cvg95, 
  # width80, width90, with95
  
  require(kableExtra)
  scoreVars = c("rmse", "bias", 
                "intScore80", "intScore90", "intScore95", 
                "cvg89", "cvg90", "cvg95", 
                "width80", "width90", "with95")
  
  # make a fancy table for each scoreVar
  for(i in 1:length(scoreVars)) {
    thisScore = scoreVars[i]
    thisTab = meanScoresTab[[c("beta", "rho", "nClustFac", "nEAsFac", thisScore)]]
    
    # p. 15: nice color use illustration:
    # vs_dt <- iris[1:10, ]
    # vs_dt[1:4] <- lapply(vs_dt[1:4], function(x) {
    #   cell_spec(x, bold = T,
    #             color = spec_color(x, end = 0.9),
    #             font_size = spec_font_size(x))
    # })
    # vs_dt[5] <- cell_spec(vs_dt[[5]], color = "white", bold = T,
    #                       background = spec_color(1:10, end = 0.9, option = "A", direction = -1))
    # kbl(vs_dt, booktabs = T, escape = F, align = "c") %>%
    #   kable_classic("striped", full_width = F)
    
    #  p. 17: 3 headers with nice hlines
    # kbl(dt, booktabs = T) %>%
    #   kable_styling(latex_options = "striped") %>%
    #   add_header_above(c(" ", "Group 1" = 2, "Group 2" = 2, "Group 3" = 2)) %>%
    #   add_header_above(c(" ", "Group 4" = 4, "Group 5" = 2)) %>%
    #   add_header_above(c(" ", "Group 6" = 6), bold = T, italic = T
    
    #  p. 20: collapse_rows with correct hlines
    # collapse_rows_dt <- data.frame(C1 = c(rep("a", 10), rep("b", 5)),
    #                                C2 = c(rep("c", 7), rep("d", 3), rep("c", 2), rep("d", 3)),
    #                                C3 = 1:15,
    #                                C4 = sample(c(0,1), 15, replace = TRUE))
    # kbl(collapse_rows_dt, booktabs = T, align = "c") %>%
    #   column_spec(1, bold=T) %>%
    #   collapse_rows(columns = 1:2, latex_hline = "major", valign = "middle")
    
    
    
  }
}