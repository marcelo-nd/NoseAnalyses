# Load microbiome graph helper
source("https://raw.githubusercontent.com/marcelo-nd/growthCurveExperiment/main/growthCurveExperiment.R")


# Read Results Files
# BHI 1
gcbhi1 <- GrowthCurveExperiment(name = "BHI 1")

gcbhi1$create_gc_objects_from_table(gc_df_path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/GrowthCurves/(GC_BHI)(Inf)Cprop_Sepi_Slu__010224.xlsx",
                                     plate_reader_type = "Infinite",
                                     gc_range = "B44:CM173",
                                     strains_names = c("Blank1", "C. prop1", "C. prop2", "C. prop3", "S. epi1", "S. epi2",
                                                       "S. epi3","S. lug1", "S. lug2", "S. lug3", "Blank2"),
                                     strains_plate_cols = list("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"),
                                     strains_plate_rows = list("A", "B", "C", "D", "E", "F", "G", "H"),
                                     blank = TRUE)

gcbhi1$strains_names

gcbhi1$plot_curves()


# BHI 2

gcbhi2 <- GrowthCurveExperiment(name = "BHI 2")

gcbhi2$create_gc_objects_from_table(gc_df_path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/GrowthCurves/(GC_BHI)(Syn)Cacc_Cpse_240207.xlsx",
                                     plate_reader_type = "Biotek",
                                     gc_range = "B220:BN365",
                                     strains_names = c("Blank1", "C. acc1", "C. acc2", "C. acc3", "C. pse1", "C. pse2",
                                                       "C. pse3", "Blank2"),
                                     strains_plate_cols = list("1", "2", "3", "4", "5", "6", "7", "8"),
                                     strains_plate_rows = list("A", "B", "C", "D", "E", "F", "G", "H"),
                                     blank = TRUE)

gcbhi2$strains_names

gcbhi2$plot_curves()


# BHI 3

gcbhi3 <- GrowthCurveExperiment(name = "BHI 3")

gcbhi3$create_gc_objects_from_table(gc_df_path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/GrowthCurves/(GC_BHI)(Syn)Dpig_Ctub_240208 .xlsx",
                                    plate_reader_type = "Biotek",
                                    gc_range = "B220:BF365",
                                    strains_names = c("D. pig1", "D. pig2", "D. pig3", "C. tub1", "C. tub2",
                                                      "C. tub3", "Blank1"),
                                    strains_plate_cols = list("1", "3", "4", "5", "6", "7", "8"),
                                    strains_plate_rows = list("A", "B", "C", "D", "E", "F", "G", "H"),
                                    blank = TRUE)

gcbhi3$strains_names

gcbhi3$plot_curves()


#BHI C. acc
gcbhi_Cacc <- GrowthCurveExperiment(name = "BHI C. accolens")

gcbhi_Cacc$add_gco(gcbhi2$growthCurveObjects[2])
gcbhi_Cacc$add_gco(gcbhi2$growthCurveObjects[3])
gcbhi_Cacc$add_gco(gcbhi2$growthCurveObjects[4])

gcbhi_Cacc$strains_names

gcbhi_Cacc$plot_curves()

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


########### Tests and debugging

