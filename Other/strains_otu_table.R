df_strain <- df_joined_filtered

df_strain <- df_strain %>%
  mutate(Combined = paste(Species, Strain, sep = "_"))

df_strain2 <- tidyr::spread(data = df_strain, key = Sample, value = Abundance, fill = NA, convert = FALSE, drop = TRUE, sep = NULL)


df_strain2 <- df_strain2[, 3:50]


df_strain3 <- df_strain %>%
  mutate(Combined = paste(Species, Strain, sep = "_")) %>%
  spread(key = Sample, value = Abundance, fill = NA, convert = FALSE, drop = TRUE, sep = NULL) %>%
  select(3:50) %>%
  mutate_all(~ replace(., is.na(.), 0))


write.csv(x = df_strain3, file = "C:/Users/marce/Desktop/df_strain.csv",quote = FALSE, row.names = FALSE)


###############

chmpns_t0 <- select(strain_data, c("Species", "SC4", "SC7", "SC11", "SC25", "SC43"))

chmpns_t0_2 <- select(chmpns_t0, "SC4", "SC7", "SC11", "SC25", "SC43")

rownames(chmpns_t0_2) <- chmpns_t0$Species


barplot_from_feature_table(chmpns_t0_2)
