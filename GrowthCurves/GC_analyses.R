# Load growthCurveExperiment script
source("https://raw.githubusercontent.com/marcelo-nd/growthCurveExperiment/main/growthCurveExperiment.R")

source("C:/Users/marce/Documents/Github/growthCurveExperiment/growthCurveExperiment.R")


# Read Results Files
# BHI 1
gcbhi1 <- GrowthCurveExperiment(name = "BHI 1")

gcbhi1$create_gc_objects_from_table(gc_df_path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/GrowthCurves/(GC_BHI)(Inf)Cprop_Sepi_Slu__010224.xlsx",
                                     plate_reader_type = "Infinite",
                                     gc_range = "B44:CM173",
                                     strains_names = c("Blank1", "C. prop16", "C. prop17", "C. prop265", "S. epi28", "S. epi231",
                                                       "S. epi251","S. lug81", "S. lug115", "S. lug239", "Blank2"),
                                     strains_plate_cols = list("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"),
                                     strains_plate_rows = list("A", "B", "C", "D", "E", "F", "G", "H"),
                                     blank = TRUE, blank_col = 1, pr_correction = TRUE)

gcbhi1$strains_names

gcbhi1$plot_curves()

# BHI 2

gcbhi2 <- GrowthCurveExperiment(name = "BHI 2")

gcbhi2$create_gc_objects_from_table(gc_df_path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/GrowthCurves/(GC_BHI)(Syn)Cacc_Cpse_240207.xlsx",
                                     plate_reader_type = "Biotek",
                                     gc_range = "B220:BN365",
                                     strains_names = c("Blank1", "C. acc99", "C. acc157", "C. acc184", "C. pseDSM", "C. pse242",
                                                       "C. pse244", "Blank2"),
                                     strains_plate_cols = list("1", "2", "3", "4", "5", "6", "7", "8"),
                                     strains_plate_rows = list("A", "B", "C", "D", "E", "F", "G", "H"),
                                     blank = TRUE, blank_col = 1, pr_correction = TRUE)

gcbhi2$strains_names

gcbhi2$plot_curves(calculate_model = TRUE)


# BHI 3

gcbhi3 <- GrowthCurveExperiment(name = "BHI 3")

gcbhi3$create_gc_objects_from_table(gc_df_path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/GrowthCurves/(GC_BHI)(Syn)Dpig_Ctub_240208 .xlsx",
                                    plate_reader_type = "Biotek",
                                    gc_range = "B220:BF365",
                                    strains_names = c("D. pig21", "D. pig61", "D. pig245", "C. tubDSM", "C. tub102",
                                                      "C. tub223", "Blank1"),
                                    strains_plate_cols = list("1", "3", "4", "5", "6", "7", "8"),
                                    strains_plate_rows = list("A", "B", "C", "D", "E", "F", "G", "H"),
                                    blank = TRUE, blank_col = 8, pr_correction = TRUE)

gcbhi3$strains_names

gcbhi3$plot_curves()


#BHI C. acc
gcbhi_Cacc <- GrowthCurveExperiment(name = "BHI C. accolens")

gcbhi_Cacc$add_gco(gcbhi2$growthCurveObjects[2])
gcbhi_Cacc$add_gco(gcbhi2$growthCurveObjects[3])
gcbhi_Cacc$add_gco(gcbhi2$growthCurveObjects[4])

gcbhi_Cacc$strains_names

gcbhi_Cacc$plot_curves(calculate_model = TRUE)

#BHI C. pro

gcbhi_Cpro <- GrowthCurveExperiment(name = "BHI C. propinquum")

gcbhi_Cpro$add_gco(gcbhi1$growthCurveObjects[2])
gcbhi_Cpro$add_gco(gcbhi1$growthCurveObjects[3])
gcbhi_Cpro$add_gco(gcbhi1$growthCurveObjects[4])

gcbhi_Cpro$strains_names

gcbhi_Cpro$plot_curves()

#BHI C. pse

gcbhi_Cpse <- GrowthCurveExperiment(name = "BHI C. pseudodiphtericum")

gcbhi_Cpse$add_gco(gcbhi2$growthCurveObjects[5])
gcbhi_Cpse$add_gco(gcbhi2$growthCurveObjects[6])
gcbhi_Cpse$add_gco(gcbhi2$growthCurveObjects[7])

gcbhi_Cpse$strains_names

gcbhi_Cpse$plot_curves()

#BHI C. tub

gcbhi_Ctub <- GrowthCurveExperiment(name = "BHI C. tuberculostearicum")

gcbhi_Ctub$add_gco(gcbhi3$growthCurveObjects[4])
gcbhi_Ctub$add_gco(gcbhi3$growthCurveObjects[5])
gcbhi_Ctub$add_gco(gcbhi3$growthCurveObjects[6])

gcbhi_Ctub$strains_names

gcbhi_Ctub$plot_curves()

# BHI S. epi

gcbhi_Sepi <- GrowthCurveExperiment(name = "BHI S. epidermidis")

gcbhi_Sepi$add_gco(gcbhi1$growthCurveObjects[5])
gcbhi_Sepi$add_gco(gcbhi1$growthCurveObjects[6])
gcbhi_Sepi$add_gco(gcbhi1$growthCurveObjects[7])

gcbhi_Sepi$strains_names

gcbhi_Sepi$plot_curves()

# BHI S. lugdunensis

gcbhi_Slug <- GrowthCurveExperiment(name = "BHI S. lugdunensis")

gcbhi_Slug$add_gco(gcbhi1$growthCurveObjects[8])
gcbhi_Slug$add_gco(gcbhi1$growthCurveObjects[9])
gcbhi_Slug$add_gco(gcbhi1$growthCurveObjects[10])

gcbhi_Slug$strains_names

gcbhi_Slug$plot_curves()

# BHI D. pigrum

gcbhi_Dpig <- GrowthCurveExperiment(name = "BHI D. pigrum")

gcbhi_Dpig$add_gco(gcbhi3$growthCurveObjects[1])
gcbhi_Dpig$add_gco(gcbhi3$growthCurveObjects[2])
gcbhi_Dpig$add_gco(gcbhi3$growthCurveObjects[3])

gcbhi_Dpig$strains_names

gcbhi_Dpig$plot_curves()

