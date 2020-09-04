#' Plot Amino Acid Change Table
#'
#' This function takes in the output data frame of makeAAChangeTable and plots it in bar graph form.
#' By default, only plots amino acid changes with frequency in Stage I patients greater than 1.
#' This can be changed by setting the minfreq parameter.
#'
#'@param table Data frame produced as output from makeAAChangeTable
#'@param minfreq Minimum frequency of amino acid changes from all Stage I patients to plot, defaults to 2.
#'
#'@export

plotAAChangeTable <- function(table, minfreq = 2) {
  t1 <- head(table, -1)
  t2 <- subset(t1, t1$LUAD_StageI >= minfreq)

  if(nrow(t2) == 0) {
    print(cat('There are no amino acid changes with frequency of', minfreq, 'or higher.\n'))
  }
  fin <- pivot_longer(data = t2, cols = c('HighRisk', 'IntRisk', 'LowRisk'),
                      names_to = 'Cohort', values_to = 'Frequency')

  p <- ggplot(fin, aes(x = AAChange, y = Frequency, fill = Cohort)) + geom_bar(stat = 'identity',
                                                                               position = position_dodge())
  p + scale_fill_brewer(palette = 'Blues')
}
