# Analyses of Growth Curves for Anaerobic Organisms in BHI and SNM3 + Mucins

# Load growthCurveExperiment script
source("https://raw.githubusercontent.com/marcelo-nd/growthCurveExperiment/main/growthCurveExperiment.R")

# Anaerobic Growth Curves

# First try, not used.
bhiacc <- GrowthCurveExperiment(name = "BHI ACC")

bhiacc$create_gc_objects_from_table(gc_df_path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/GrowthCurves/GC_ANA(GC_BHI)(Syn)Aoct_Cavidum_Cacnes_040324.xlsx",
                                      plate_reader_type = "Biotek",
                                      gc_range = "B220:CD365",
                                      strains_names = c("Blank", "A. oct133", "A. oct211", "A. oct239", "C. avidum32", "C. avidum181",
                                                        "C. avidum208","C. acnes33", "C. acnes86", "C. acnes149"),
                                      strains_plate_cols = list("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
                                      strains_plate_rows = list("A", "B", "C", "D", "E", "F", "G", "H"),
                                      blank = TRUE, blank_col = 1, pr_correction = TRUE)

bhiacc$strains_names

bhiacc$plot_curves(calculate_model = FALSE)


# BHI C. avidum 6.3.24

bhi_Cavi <- GrowthCurveExperiment(name = "BHI C. avidum 6.3.24")

bhi_Cavi$create_gc_objects_from_table(gc_df_path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/GrowthCurves/GC_ANA(GC_BHI)_Cavidum_060324.xlsx",
                                        plate_reader_type = "Biotek",
                                        gc_range = "B220:AH305",
                                        strains_names = c("Blank", "C. avi32", "C. avi181", "C. avi208"),
                                        strains_plate_cols = list("1", "2", "3", "4"),
                                        strains_plate_rows = list("A", "B", "C", "D", "E", "F", "G", "H"),
                                        blank = TRUE, blank_col = 1, pr_correction = TRUE)

bhi_Cavi$strains_names

bhi_Cavi$plot_curves(calculate_model = FALSE)


# BHI C. acnes 14.3.24 BHI

bhi_Cacn <- GrowthCurveExperiment(name = "BHI C. acnes")

bhi_Cacn$create_gc_objects_from_table(gc_df_path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/GrowthCurves/GC_ANA(GC_BHI_SNM)(Syn)C.acnes_140324.xlsx",
                                         plate_reader_type = "Biotek",
                                         gc_range = "B220:BN365",#how do I choose the BHI columns only?
                                         strains_names = c("BlankBHI", "C. acn33_BHI", "C. acn86_BHI", "C. acn149_BHI", "BlankSNM", "C. acn33_SNM", "C. acn86_SNM", "C. acnes149_SNM"),
                                         strains_plate_cols = list("1", "2", "3", "4", "5", "6", "7", "8"),
                                         strains_plate_rows = list("A", "B", "C", "D", "E", "F", "G", "H"),
                                         blank = TRUE, blank_col = 1, pr_correction = TRUE)

bhi_Cacn$strains_names

bhi_Cacn$plot_curves(calculate_model = FALSE)


# C. acnes 14.3.24 SNM

snm_Cacn <- GrowthCurveExperiment(name = "SNM C. acnes")

snm_Cacn$create_gc_objects_from_table(gc_df_path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/GrowthCurves/GC_ANA(GC_BHI_SNM)(Syn)C.acnes_140324.xlsx",
                                         plate_reader_type = "Biotek",
                                         gc_range = "B220:BN365",#how do I choose the BHI columns only?
                                         strains_names = c("BlankBHI", "C. acn33_BHI", "C. acn86_BHI", "C. acn149_BHI", "BlankSNM", "C. acn33_SNM", "C. acn86_SNM", "C. acnes149_SNM"),
                                         strains_plate_cols = list("1", "2", "3", "4", "5", "6", "7", "8"),
                                         strains_plate_rows = list("A", "B", "C", "D", "E", "F", "G", "H"),
                                         blank = TRUE, blank_col = 5, pr_correction = TRUE)

snm_Cacn$strains_names

snm_Cacn$plot_curves(calculate_model = FALSE)



# A. octavius BHI

bhi_Aoct <- GrowthCurveExperiment(name = "BHI A.octavius")

bhi_Aoct$create_gc_objects_from_table(gc_df_path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/GrowthCurves/GC_ANA(GC_BHI_SNM)(Syn)A.oct_C.avi.xlsx",
                                         plate_reader_type = "Biotek",
                                         gc_range = "B220:CL327",
                                         strains_names = c("BlankBHI", "A. oct133BHI", "A. oct211BHI", "A. oct259BHI", "BlankSNM", "A. oct133SNM", "A. oct211SNM", "A. oct259SNM", "C. avi32SNM", "C. avi181SNM", "C. avi208SNM"),
                                         strains_plate_cols = list("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"),
                                         strains_plate_rows = list("A", "B", "C", "D", "E", "F", "G", "H"),
                                         blank = TRUE, blank_col = 1, pr_correction = TRUE)

bhi_Aoct$strains_names

bhi_Aoct$plot_curves(calculate_model = FALSE)



# SNM for A. octavius and C. avidum

snm_Aoct_Cavi <- GrowthCurveExperiment(name = "SNM A.octavius C. avid")

snm_Aoct_Cavi$create_gc_objects_from_table(gc_df_path = "C:/Users/marce/OneDrive - UT Cloud/1_NoseSynCom Project/Experiments/GrowthCurves/GC_ANA(GC_BHI_SNM)(Syn)A.oct_C.avi.xlsx",
                                         plate_reader_type = "Biotek",
                                         gc_range = "B220:CL327",
                                         strains_names = c("BlankBHI", "A. oct133BHI", "A. oct211BHI", "A. oct259BHI", "BlankSNM", "A. oct133SNM", "A. oct211SNM", "A. oct259SNM", "C. avi32SNM", "C. avi181SNM", "C. avi208SNM"),
                                         strains_plate_cols = list("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"),
                                         strains_plate_rows = list("A", "B", "C", "D", "E", "F", "G", "H"),
                                         blank = TRUE, blank_col = 5, pr_correction = TRUE)

snm_Aoct_Cavi$strains_names

snm_Aoct_Cavi$plot_curves(calculate_model = FALSE)


#######################
# A. octavius BHI

gcbhi_Aoct <- GrowthCurveExperiment(name = "BHI A. octavius")

gcbhi_Aoct$add_gco(bhi_Aoct$growthCurveObjects[2])
gcbhi_Aoct$add_gco(bhi_Aoct$growthCurveObjects[3])
gcbhi_Aoct$add_gco(bhi_Aoct$growthCurveObjects[4])

gcbhi_Aoct$strains_names

gcbhi_Aoct$plot_curves(calculate_model = FALSE, yScalemin = 0, yScalemax = 3)


# A. octavius SNM

gcsnm_Aoct <- GrowthCurveExperiment(name = "SNM A. octavius")

gcsnm_Aoct$add_gco(snm_Aoct_Cavi$growthCurveObjects[6])
gcsnm_Aoct$add_gco(snm_Aoct_Cavi$growthCurveObjects[7])
gcsnm_Aoct$add_gco(snm_Aoct_Cavi$growthCurveObjects[8])

gcsnm_Aoct$strains_names

gcsnm_Aoct$plot_curves(calculate_model = FALSE, yScalemin = 0, yScalemax = 0.5)


# C. acnes BHI

gcbhi_Cacn <- GrowthCurveExperiment(name = "BHI C. acnes")

gcbhi_Cacn$add_gco(bhi_Cacn$growthCurveObjects[2])
gcbhi_Cacn$add_gco(bhi_Cacn$growthCurveObjects[3])
gcbhi_Cacn$add_gco(bhi_Cacn$growthCurveObjects[4])

gcbhi_Cacn$strains_names

gcbhi_Cacn$plot_curves(calculate_model = FALSE, yScalemin = 0, yScalemax = 3)

# C. acnes SNM

gcsnm_Cacn <- GrowthCurveExperiment(name = "SNM C. acnes")

gcsnm_Cacn$add_gco(snm_Cacn$growthCurveObjects[6])
gcsnm_Cacn$add_gco(snm_Cacn$growthCurveObjects[7])
gcsnm_Cacn$add_gco(snm_Cacn$growthCurveObjects[8])

gcsnm_Cacn$strains_names

gcsnm_Cacn$plot_curves(calculate_model = FALSE, yScalemin = 0, yScalemax = 0.5)


# BHI C. avidum

gcbhi_Cavi <- GrowthCurveExperiment(name = "BHI C. avidum")

gcbhi_Cavi$add_gco(bhi_Cavi$growthCurveObjects[2])
gcbhi_Cavi$add_gco(bhi_Cavi$growthCurveObjects[3])
gcbhi_Cavi$add_gco(bhi_Cavi$growthCurveObjects[4])

gcbhi_Cavi$strains_names

gcbhi_Cavi$plot_curves(calculate_model = FALSE, yScalemin = 0, yScalemax = 6)


# SNM C. avidum

gcbhi_Cavi_snm <- GrowthCurveExperiment(name = "SNM C. avidum")

gcbhi_Cavi_snm$add_gco(snm_Aoct_Cavi$growthCurveObjects[9])
gcbhi_Cavi_snm$add_gco(snm_Aoct_Cavi$growthCurveObjects[10])
gcbhi_Cavi_snm$add_gco(snm_Aoct_Cavi$growthCurveObjects[11])

gcbhi_Cavi_snm$strains_names

gcbhi_Cavi_snm$plot_curves(calculate_model = FALSE, yScalemin = 0, yScalemax = 1)

