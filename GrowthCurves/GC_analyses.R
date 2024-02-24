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
                                                       "S. pse3", "Blank2"),
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
gcbhi_Cacc <- GrowthCurveExperiment(name = "BHI C. Accolens")

gcbhi_Cacc$add_gco(gcbhi2$growthCurveObjects[2])
gcbhi_Cacc$add_gco(gcbhi2$growthCurveObjects[3])
gcbhi_Cacc$add_gco(gcbhi2$growthCurveObjects[4])

gcbhi_Cacc$strains_names

gcbhi_Cacc$plot_curves()

#BHI C. pro

#BHI C. pse

#BHI C. tub




gcbhi2$growthCurveObjects[[1]]$data




# Helper functions
parse_time_in_hours <- function(time_string){
  time_split <- strsplit(time_string, split = ":")
  hours <- as.numeric(time_split[[1]][1])
  mins <- as.numeric(time_split[[1]][2]) + (as.numeric(time_split[[1]][3])/60)
  return(hours + (mins/60))
}

parse_time_in_hours <- function(time_string){
  # Parse Time from Biotek-style results
  time_split <- strsplit(time_string, split = ":")
  hours <- as.numeric(time_split[[1]][1])
  mins <- as.numeric(time_split[[1]][2]) + (as.numeric(time_split[[1]][3])/60)
  return(hours + (mins/60))
}

parse_time_in_hours2 <- function(time_string){
  return(time_string * 24)
}

rowSDs <- function(x, cols_init, cols_end, na.rm=F) {
  # Vectorised version of variance filter
  return(apply(x[,cols_init:cols_end], MARGIN = 1, FUN = sd))
}

calculate_replicates_stats_per_strain <- function(plate_df, blank_col = 1, col_start = 2, n_samples = 11){
  # Rows are replicates, columns are strains.
  rows <- toupper(paste("", letters[1:8], sep=""))
  columns <- rep(col_start:(col_start+n_samples-1))
  # Create
  blank_wells <- paste(rows, as.character(blank_col), sep = "")
  #print(blank_wells)
  od_means_per_strain <- dplyr::select(plate_df, "Time")
  od_sds_per_strain <- dplyr::select(plate_df, "Time")
  #print(od_means_per_strain)
  col_names <- c("Time")
  for (strain_col in columns) {
    strain_wells <- paste(rows, as.character(strain_col), sep = "")
    strain_df <- dplyr::select(plate_df, any_of(strain_wells))
    strain_mean_df <- data.frame(rowMeans(strain_df))
    strain_sd_df <- data.frame(rowSDs(strain_df, 1, ncol(strain_df)))
    #print(head(strain_df_sd))
    #
    od_means_per_strain <- cbind(od_means_per_strain, strain_mean_df)
    od_sds_per_strain <- cbind(od_sds_per_strain, strain_sd_df)
    col_names <- c(col_names, strain_wells[1])
  }
  #print(dim(od_means_per_strain))
  #print(head(od_means_per_strain))
  #print(dim(od_sds_per_strain))
  #print(head(od_sds_per_strain))
  colnames(od_means_per_strain) <- col_names
  colnames(od_sds_per_strain) <- col_names
  return(list(od_means_per_strain, od_sds_per_strain))
}

calculate_growth_curve_models <- function(growth.curve.points){
  # Based on https://rpubs.com/angelov/growthcurver.
  summG <- function(x) {growthcurver::SummarizeGrowth(growth.curve.points$Time,x)}
  models.all <- lapply(growth.curve.points[2:ncol(growth.curve.points)], summG)
  
  df.predicted <- data.frame(time = growth.curve.points$Time)
  
  #print(head(df.predicted.plate))
  for (i in names(growth.curve.points[2:ncol(growth.curve.points)])) 
  {
    df.predicted[[i]] <- predict(models.all[[i]]$model)
  }
  colnames(df.predicted) <- c("Time", colnames(df.predicted)[2:ncol(df.predicted)])
  
  return(df.predicted)
}

create_growth_curve_plot_stats.df <- function(growth_curve_df, calculate_model = TRUE, col_start = 2, n_samples = 11, blank_col = 1, species_names = NULL){
  od_stats <- calculate_replicates_stats_per_strain(plate_df = growth_curve_df, n_samples = n_samples, blank_col = blank_col, col_start = col_start)
  
  od_means <- transform(od_stats[1], Time = as.numeric(Time))
  
  colnames(od_means) <- c("Time", species_names)
  
  od_sds <- transform(od_stats[2], Time = as.numeric(Time))
  
  colnames(od_sds) <- c("Time", species_names)
  
  print(head(od_means))
  
  print(head(od_sds))
  
  # .
  od_g <- tidyr::gather(od_means, key = "Species", value = "OD", 2:(n_samples + 1))
  od_sds_g <- tidyr::gather(od_sds, key = "Species", value = "SD", 2:(n_samples + 1))
  
  if (calculate_model) {
    #
    df.predicted.plate <- calculate_growth_curve_models(od_means)
    
    pred_g <- tidyr::gather(df.predicted.plate, key = "Species", value = "pred.od", 2:(n_samples + 1))
    
    od_means_sds_preds <- cbind(od_g, od_sds_g$SD, pred_g$pred.od)
    
    colnames(od_means_sds_preds) <- c("Time", "Species", "od", "sd", "pred.od")
    
    return(od_means_sds_preds)
    
  }else {
    od_means_sds_preds <- cbind(od_g, od_sds_g$SD)
    
    colnames(od_means_sds_preds) <- c("Time", "Species", "od", "sd")
    
    return(od_means_sds_preds)
  }
}

plot_curves_together <- function(gcurve.stats_df, plot_type = "curves_plot", color_scale = c("#AD7BE9", "#3E54AC", "#658864", "#6D9886", "#E96479", "#E69F00", "#FC7300",  "#183A1D", "#635985", "#F99417", "#FEA1BF")){
  # Create basic curves plot
  curves_plot <- gcurve.stats_df %>%
    ggplot(aes(x = Time, y = od, group = Species, color = Species)) +
    geom_errorbar(aes(ymin = od - sd, ymax = od + sd), width= 0.1) +
    geom_point(alpha=0.7) +
    scale_color_manual(values=color_scale)
  
  if (plot_type == "curves_plot2") {
    # Add model lines
    curves_plot +
      geom_line(aes(x = Time, y=pred.od), linewidth = 0.5)
  }
  else if (plot_type == "curves_plot3") {
    # Add line at 0.02 y intercept.
    curves_plot +
      geom_line(aes(x = Time, y = pred.od), linewidth = 0.5) +
      geom_hline(yintercept=0.02, linetype="dashed", color = "#0081B4")
  }else{
    curves_plot
  }
}





gc1 <- read.csv2("C:/Users/marce/Desktop/Mappe1.csv", sep = ";", header = TRUE)

summary(gc1)

strains_names = c("Spcs1", "Spcs2", "Spcs3", "Spcs4", "Spcs5")

growth.curve.stats_g1 <- create_growth_curve_plot_stats.df(gc1, calculate_model = TRUE, blank_col = 2, col_start = 3, n_samples = 5, species_names = strains_names)

plot1 <- plot_curves_together(gcurve.stats_df = growth.curve.stats_g1, plot_type = "curves_plot3")

plot1




###### tests with growthCurveExperiment

gcexp1 <- GrowthCurveExperiment(name = "test1")

gcexp1$create_gc_objects_from_table1(gc_df_path = "C:/Users/marce/Desktop/BestestGrowthCurves/GC 24hrs_231226_090916_0001.xlsx",
                                     plate_reader_type = "Biotek",
                                     gc_range = "B220:BN365",
                                     strains_names = c("Blank1", "Spcs1", "Spcs2", "Spcs3", "Spcs4",
                                                       "Spcs5","Spcs6", "Blank2"),
                                     strains_plate_cols = list("1", "2", "3", "4", "5", "6", "7", "8"),
                                     strains_plate_rows = list("A", "B", "C", "D", "E", "F", "G", "H"),
                                     blank = TRUE)


gcexp1$strains_names

gcexp1$plot_curves()

gcexp1$remove_gco("Blank2")

gcexp1$remove_gco("Blank1")


