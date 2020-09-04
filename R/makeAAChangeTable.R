#' Make Amino Acid Change Table
#'
#' This function takes in the MAF data frame produced by GDCquery_Maf from TCGAbiolinks, a clinical information
#' table with the bspcacat risk classifer, and a gene name (Hugo symbol).
#' It return a data frame of the frequencies of each amino acid change in
#' high, intermediate, and low risk of recurrence patients, in addition to all Stage I LUAD patients and all LUAD
#' patients from the TCGA database.
#' Optionally, if csv is set to TRUE, it writes the data frame to a csv file in the folder 'AAChanges' in the current
#' directory.
#'
#' @param maf Data frame returned by GDCquery_Maf
#' @param riskcat Data frame containing the bspcacat risk classifier for each patient
#' @param gene Gene name (Hugo symbol) as a string
#' @param csv Logical that tells the function whether to write the result to a csv file or not. Defaults to FALSE
#'
#' @return Data frame containing computed frequencies of each amino acid change in high, intermediate, and low risk patients.
#' @export

makeAAChangeTable <- function(maf, riskcat, gene, csv = FALSE) {
  maf$bcr_patient_barcode <- gsub('.{16}$', '', maf$Tumor_Sample_Barcode)

  high <- subset(riskcat, riskcat$bspcacat == 'high')
  int <- subset(riskcat, riskcat$bspcacat == 'int')
  low <- subset(riskcat, riskcat$bspcacat == 'low')

  gene.all.maf <- subset(maf, maf$Hugo_Symbol == gene)
  gene.stageI.maf <- subset(maf, maf$bcr_patient_barcode %in% riskcat$bcr_patient_barcode & maf$Hugo_Symbol == gene)
  gene.high.maf <- subset(maf, maf$bcr_patient_barcode %in% high$bcr_patient_barcode & maf$Hugo_Symbol == gene)
  gene.int.maf <- subset(maf, maf$bcr_patient_barcode %in% int$bcr_patient_barcode & maf$Hugo_Symbol == gene)
  gene.low.maf <- subset(maf, maf$bcr_patient_barcode %in% low$bcr_patient_barcode & maf$Hugo_Symbol == gene)

  all.table <- as.data.frame(table(gene.all.maf$HGVSp_Short, dnn = 'AAChange'))
  stageI.table <- as.data.frame(table(gene.stageI.maf$HGVSp_Short, dnn = 'AAChange'))
  high.table <- as.data.frame(table(gene.high.maf$HGVSp_Short, dnn = 'AAChange'))
  int.table <- as.data.frame(table(gene.int.maf$HGVSp_Short, dnn = 'AAChange'))
  low.table <- as.data.frame(table(gene.low.maf$HGVSp_Short, dnn = 'AAChange'))

  all.table <- all.table %>% rename(LUAD_All = Freq)
  stageI.table <- stageI.table %>% rename(LUAD_StageI = Freq)
  high.table <- high.table %>% rename(HighRisk = Freq)
  int.table <- int.table %>% rename(IntRisk = Freq)
  low.table <- low.table %>% rename(LowRisk = Freq)

  t1 <- merge(high.table, int.table, by = 'AAChange', all = TRUE)
  t2 <- merge(t1, low.table, by = 'AAChange', all = TRUE)
  t3 <- merge(t2, stageI.table, by = 'AAChange', all = TRUE)
  combined <- merge(t3, all.table, by = 'AAChange', all = TRUE)

  combined[is.na(combined)] <- 0
  combined$AAChange <- str_replace(combined$AAChange, '^.{2}', '')
  combined <- combined %>% arrange(desc(LUAD_All))
  combined <- adorn_totals(combined, where = 'row', name = 'Total') # adorn_totals from janitor pkg

  if(csv) {
    dir.create('AAChange', showWarnings = FALSE)
    write.csv(combined, paste('AAChange/', gene, '.csv', sep = ''), row.names = FALSE)
  }

  return(combined)
}
