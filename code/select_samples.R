# Select samples for differential expression analysis

library(tidyverse)
options(tibble.print_max=Inf)
options(tibble.width=Inf)
metadata1 <- read_tsv(here("data", "atlas_metadata_stress_samples.tsv"))[,-10]
metadata2 <- read_tsv(here("data", "newsamples_metadata_stress_samples.tsv"))

#----Select only stress samples----
stress_metadata <- bind_rows(metadata1, metadata2) %>%
    filter(!is.na(Stress_info))

#----Create different objects for biotic and abiotic stress----
abiotic <- stress_metadata %>%
    filter(Stress_class == "abiotic")

biotic <- stress_metadata %>%
    filter(Stress_class == "biotic", !is.na(Sample_description))

#----Display BioProjects focusing on abiotic stress and their sample size
abiotic %>%
    group_by(BioProject, Stress_info) %>%
    summarise(N=n()) %>%
    arrange(Stress_info) %>%
    print(n=Inf)

#----Split tibbles into list of tibbles----
abiotic_list <- abiotic %>%
    select(c(2,4:7)) %>%
    arrange(BioProject, Sample_description) %>%
    group_by(Stress_info) %>%
    group_split()

# Rename 'control' class to avoid redundancy with other stress classes
abiotic_list[[12]]$Sample_description <- gsub("control", "control_root",
                                              abiotic_list[[12]]$Sample_description)
abiotic_list[[17]]$Sample_description <- paste(abiotic_list[[17]]$Sample_description,
                                               abiotic_list[[17]]$Tissue,sep="_")
abiotic_list[[10]]$Sample_description <- paste(abiotic_list[[10]]$Sample_description,
                                               abiotic_list[[10]]$Tissue,sep="_")
abiotic_list[[16]]$Sample_description <- paste(abiotic_list[[16]]$Sample_description,
                                               abiotic_list[[16]]$Tissue,sep="_")
abiotic_list[[15]]$Sample_description <- paste(abiotic_list[[15]]$Sample_description,
                                               abiotic_list[[15]]$Tissue,sep="_")
abiotic_list[[15]]$Sample_description <- gsub("block[0-9]_", "", 
                                              abiotic_list[[15]]$Sample_description)
abiotic_list[[18]]$Sample_description <- paste(abiotic_list[[18]]$Sample_description,
                                               abiotic_list[[18]]$Tissue,sep="_")


#----Select samples and create table of contrasts for abiotic stress----
abiotic_list

make_contrast_table <- function(control, case) {
    df <- data.frame(control=control,
                     case=case)
    return(df)
}

## Acidity ##
acidity <- bind_rows(
    abiotic_list[[1]][7:12, ],
    abiotic_list[[12]][c(9:11, 30:35, 41), ]
)
contrast_acidity <- make_contrast_table(c("control_9days_acidity", "control_root"),
                                        c("stress_9days_acidity", "stress_acidity"))

## Cadmium ##
cadmium <- abiotic_list[[2]]
contrast_cadmium <- make_contrast_table("control", "stress")

## Cold ##
cold <- bind_rows(
    abiotic_list[[3]],
    abiotic_list[[12]][c(9:11, 33:35, 41, 6:9, 38, 42), ]
)
contrast_cold <- make_contrast_table(c("control", "control"),
                                     c("stress_1h", "stress_24h"))

## Drought ##
drought <- bind_rows(
    abiotic_list[[4]][c(30:65, 70:72, 76:84, 88:93), ],
    abiotic_list[[5]][1:4, ],
    abiotic_list[[6]][c(1:3, 12:14), ],
    abiotic_list[[7]][c(1:6, 10:12, 16:18), ],
    abiotic_list[[12]][c(9:11, 33:35, 41, 25:27), ]
)
contrast_drought <- make_contrast_table(
    c("control_4h", "control_8h", 
      "control_12h", "control_16h",
      "control_20h", "control_24h",
      "control_mild", "control_verymild",
      "control_severe", "controlCO2_control",
      "control", "control_0h",
      "control_0h", "control_0h", 
      "control_root"),
    c("stress_4h", "stress_8h", 
      "stress_12h", "stress_16h",
      "stress_20h", "stress_24h",
      "stress_mild", "stress_verymild",
      "stress_severe", "controlCO2_50min",
      "stress_8h_drought", "stress_1h_drought",
      "stress_6h_drought", "stress_12h_drought",
      "stress_drought"))

## CO2 ##
co2 <- bind_rows(
    abiotic_list[[5]][c(3,4,7,8), ],
    abiotic_list[[12]][c(1:3, 9:11, 33:35, 41), ]
)
contrast_co2 <- make_contrast_table(c("control_root", "controlCO2_control"),
                                    c("stress_CO2", "stressCO2_control"))

## Heat ##
heat <- bind_rows(
    abiotic_list[[6]][c(1:3, 8:11, 17:19), ],
    abiotic_list[[8]][, ],
    abiotic_list[[12]][c(4:5, 9:11, 33:37, 41), ]
)
contrast_heat <- make_contrast_table(
    c("control", "control", 
      "control_6HAI_PI587982", "control_6HAI_PI636696",
      "control_germinated_PI587982", "control_germinated_PI636696",
      "control_mature_PI587982", "control_mature_PI636696",
      "control_root"),
    c("stress_8h_heat", "stress_24h_heat",
      "stress_6HAI_PI587982", "stress_germinated_PI636696",
      "stress_germinated_PI587982", "stress_germinated_PI636696",
      "stress_mature_PI587982", "stress_mature_PI636696",
      "stress_heat")
    )

## Salt ##
salt <- bind_rows(
    abiotic_list[[7]][c(1:3, 7:9, 13:15, 19:21), ],
    abiotic_list[[12]][c(9:11, 22:24, 33:35, 41), ],
    abiotic_list[[17]][c(1:8, 17:58), ]
)
contrast_salt <- make_contrast_table(
    c("control_0h", "control_0h", "control_0h",
      "control_root", "control_4h_leaves", "control_12h_leaves",
      "control_leaves", "control_leaves", "control_leaves",
      "control_leaves", "control_leaves",
      "control_root", "control_root", "control_root",
      "control_root", "control_root"),
    c("stress_1h_salt", "stress_6h_salt", "stress_12h_salt", 
      "stress_salt", "stress_4h_leaves", "stress_12h_leaves",
      "stress_1h_leaves", "stress_2h_leaves", "stress_4h_leaves",
      "stress_24h_leaves", "stress_48h_leaves",
      "stress_1h_root", "stress_2h_root", "stress_4h_root",
      "stress_24h_root", "stress_48h_root")
    )

## Herbicide ##
herbicide <- abiotic_list[[9]]
contrast_herbicide <- make_contrast_table(
    c("control_single", "control_stacked"),
    c("stress_single", "stress_stacked")
)

## Iron ##
iron <- bind_rows(
    abiotic_list[[10]][c(1:19, 23:25), ],
    abiotic_list[[11]][c(1:6, 20:23), ],
    abiotic_list[[12]][c(9:11, 33:35, 41, 20:21, 39), ]
)
iron$Sample_description[17:19] <- rep("control_root_iron", 3) 
contrast_iron <- make_contrast_table(
    c("control_1h_leaves", "control_1h_root", 
      "control_6h_leaves", "control_6h_root",
      "control_root_iron", "control_control_sens",
      "control_control_tol", "control_root"),
    c("stress_1h_leaves", "stress_1h_root",
      "stress_6h_leaves", "stress_6h_root",
      "stress_root", "stress_control_sens",
      "stress_control_tol", "stress_iron")
)

## Alkalinity ##
alkalinity <- bind_rows(
    abiotic_list[[11]][c(1:19), ],
    abiotic_list[[12]][c(9:11, 33:35, 41, 28:29, 40), ]
)
contrast_alkalinity <- make_contrast_table(
    c("control_control_sens", "control_control_sens",
      "control_control_tol", "control_control_tol"),
    c("control_stress2.5_sens", "control_stress5_sens",
      "control_stress2.5_tol", "control_stress5_tol")
)


## Phosphorus ##
phosphorus <- bind_rows(
    abiotic_list[[12]][c(9:11, 33:35, 41, 14:16), ],
    abiotic_list[[16]][c(1:22), ]
)
contrast_phosphorus <- make_contrast_table(
    c("control_root", "control_leaves", "control_nodule",
      "control_Bogao_root", "control_Nannong_root"),
    c("stress_phosphorus", "stress_leaves", "stress_nodule",
      "stress_Bogao_root", "stress_Nannong_root")
)


## Nitrogen ##
nitrogen <- abiotic_list[[12]][c(9:11, 33:35, 41, 12:13, 43), ]
contrast_nitrogen <- make_contrast_table("control_root", "stress_nitrogen")


## Potassium ##
potassium <- abiotic_list[[12]][c(9:11, 33:35, 41, 17:19, 44), ]
contrast_potassium <- make_contrast_table("control_root", "stress_potassium")


## Ozone ##
ozone <- bind_rows(
    abiotic_list[[13]][c(1:11, 15, 17, 20, 24), ],
    abiotic_list[[15]][c(1:22), ]
)
ozone$Sample_description <- gsub("[0-9]", "", ozone$Sample_description)
ozone$Sample_description <- gsub("_Rstage", "", ozone$Sample_description)

contrast_ozone <- make_contrast_table(
    c("control_block", "control_leaves", 
      "control_flower", "control_pod"),
    c("stress_block_ozone", "stress_leaves",
      "stress_flower", "stress_pod")
)

## Zinc ##
zinc <- abiotic_list[[18]]
contrast_zinc <- make_contrast_table(
    c("control_leaves", "control_root"),
    c("stress_leaves", "stress_root")
    )

#----Combine tables of abiotic stress-related samples and contrasts to be used----
abiotic_contrast_list <- list(
  acidity = contrast_acidity,
  alkalinity = contrast_alkalinity,
  cadmium = contrast_cadmium,
  co2 = contrast_co2,
  cold = contrast_cold,
  drought = contrast_drought,
  heat = contrast_heat,
  herbicide = contrast_herbicide,
  iron = contrast_iron,
  nitrogen = contrast_nitrogen,
  ozone = contrast_ozone,
  phosphorus = contrast_phosphorus,
  potassium = contrast_potassium,
  salt = contrast_salt,
  zinc = contrast_zinc
)

abiotic_samplelist <- list(
  acidity = acidity,
  alkalinity = alkalinity,
  cadmium = cadmium, 
  co2 = co2,
  cold = cold,
  drought = drought,
  heat = heat,
  herbicide = herbicide,
  iron = iron,
  nitrogen = nitrogen,
  ozone = ozone,
  phosphorus = phosphorus,
  potassium = potassium,
  salt = salt,
  zinc = zinc
)
abiotic_samplelist <- lapply(abiotic_samplelist, function(x) {
  return(x[!duplicated(x$BioSample), ])
})


##################################
#----Show number of samples involving biotic stress----
# By stress class
biotic %>%
    group_by(Stress_info) %>%
    summarise(N=n()) %>%
    arrange(Stress_info) %>%
    print(n=Inf)

# By pathogen
biotic %>%
    group_by(Pathogen) %>%
    summarise(N=n())


#----Split tibbles into list of tibbles----
biotic_list <- biotic %>%
    select(c(3,4:7)) %>%
    arrange(BioProject, Sample_description) %>%
    group_by(Pathogen) %>%
    group_split()


#----Select samples and create table of contrasts for biotic stress----
biotic_list

## Aphis glycines ##
aglycines <- bind_rows(
  biotic_list[[1]],
  biotic_list[[16]]
)
contrast_aglycines <- make_contrast_table(
    c("control_res", "control_res", "control_res",
      "control_sus", "control_sus", "control_sus",
      "control_24HAI_antibiotic", "control_24HAI_antixenotic",
      "control_24HAI_susceptible", "control_48HAI_antibiotic",
      "control_48HAI_antixenotic", "control_48HAI_susceptible",
      "control_96HAI_antibiotic", "control_96HAI_antixenotic",
      "control_96HAI_susceptible",
      "res_0h", "res_4h", "res_8h", "res_24h", "res_48h"),
    c("stress_res_6h", "stress_res_12h", "stress_res_48h",
      "stress_sus_6h", "stress_sus_12h", "stress_sus_48h",
      "stress_24HAI_antibiotic", "stress_24HAI_antixenotic",
      "stress_24HAI_susceptible", "stress_48HAI_antibiotic",
      "stress_48HAI_antixenotic", "stress_48HAI_susceptible",
      "stress_96HAI_antibiotic", "stress_96HAI_antixenotic",
      "stress_96HAI_susceptible",
      "sus_0h", "sus_4h", "sus_8h", "sus_24h", "sus_48h")
)


## Fusarium graminearum ##
fgraminearum <- biotic_list[[3]]
contrast_fgraminearum <- make_contrast_table(
  c("control_PI567301B", "control_wyandot"),
  c("stress_PI567301B", "stress_wyandot")
)

## Fusarium oxysporum ##
foxysporum <- biotic_list[[4]]
foxysporum$Sample_description[19:30] <- paste(foxysporum$Sample_description[19:30], 
                                       foxysporum$BioProject[19:30], 
                                       sep="_")
foxysporum$Sample_description <- gsub("PRJNA", "", foxysporum$Sample_description)
contrast_foxysporum <- make_contrast_table(
  c("control_72h", "control_72h",
    "control_96h", "control_96h",
    "control_393826", "control_599037"),
  c("stress_72h_FO36", "stress_72h_FO40",
    "stress_96h_FO36", "stress_96h_FO40",
    "stress_393826", "stress_599037")
)

## Fusarium virguliforme ##
fvirguliforme <- biotic_list[[6]]
contrast_fvirguliforme <- make_contrast_table(
  c("control_10-24days", "control_3-5days",
    "control_0dpi", "control_2dpi", "control_4dpi",
    "control_7dpi", "control_10dpi", "control_14dpi"),
  c("stress_10-24days", "stress_3-5days",
    "stress_0dpi", "stress_2dpi", "stress_4dpi",
    "stress_7dpi", "stress_10dpi", "stress_14dpi")
)

## Heterodera glycines ##
hglycines <- biotic_list[[7]][c(1:12, 25:36), ]
contrast_hglycines <- make_contrast_table(
  c("control_10days", "control_10days", "control_10days",
    "control_res", "control_sus"),
  c("stress_10days", "stress_15days", "stress_20days",
    "stress_res", "stress_sus")
)

## MAMP ##
mamp <- biotic_list[[8]]
contrast_mamp <- make_contrast_table(
  c("control_LD", "control_LDX", 
    "control_RIL-11268", "control_RIL-11272"),
  c("stress_LD", "stress_LDX",
    "stress_RIL-11268", "stress_RIL-11272")
)

## Macrophomina phaseolina ##
mphaseolina <- biotic_list[[9]]
contrast_mphaseolina <- make_contrast_table(
  c("control_res", "control_sus"),
  c("stress_res", "stress_sus")
)

## Phialophora gregata ##
pgregata <- biotic_list[[10]]
contrast_pgregata <- make_contrast_table(
  c("control_res", "control_sus", 
    "control_12h", "control_24h", "control_36h"),
  c("stress_res", "stress_sus",
    "stress_12h", "stress_24h", "stress_36h")
)

## Phakopsora pachyrhizi
ppachyrhizi <- biotic_list[[11]]
contrast_ppachyrhizi <- make_contrast_table(
  c("control_res_24h", "control_sus_24h",
    "control", "control", "control"), 
  c("stress_res_24h", "stress_sus_24h",
    "stress_7hai", "stress_10dai", "stress_48hai")
)

## Phytophthora sojae ##
psojae <- biotic_list[[12]][-c(28:38), ]
psojae$Sample_description[1:11] <- "control_rps"
psojae$Sample_description[12:22] <- "stress_rps"
psojae$Sample_description[23] <- "control_0hpi"
psojae$Sample_description <- gsub("_c.*", "", psojae$Sample_description)
psojae$Sample_description[123:125] <- "stress"
psojae$Sample_description <- gsub("_sloan", "", psojae$Sample_description)

contrast_psojae <- make_contrast_table(
  c("control_rps", 
    "control_0hpi", "control_0hpi", 
    "control_0hpi", "control_0hpi",
    "control"),
  c("stress_rps", 
    "stress_0.5hpi", "stress_3hpi", 
    "stress_6hpi", "stress_12hpi",
    "stress")
)

## Psojae glucan elicitor
psojae_gluc_elic <- biotic_list[[13]]
contrast_psojae_gluc_elic <- make_contrast_table(
  c("control_24h", "control_48h"),
  c("stress_24h", "stress_48h")
)

## Rotylenchulus reniformis ##
rreniformis <- biotic_list[[14]]
contrast_rreniformis <- make_contrast_table(
  c("control_3days", "control_6days", 
    "control_9days", "control_12days"),
  c("stress_3days", "stress_6days",
    "stress_9days", "stress_12days")
)

## Rhizoctonia solani ##
rsolani <- biotic_list[[15]]
contrast_rsolani <- make_contrast_table(
  c("control_24h"),
  c("stress_24h")
)

## Spodoptera litura ##
slitura <- biotic_list[[17]]
slitura$Sample_description <- gsub("_.*", "", slitura$Sample_description)
contrast_slitura <- make_contrast_table(
  "control",
  "stress"
)

## Soybean mosaic virus ##
smv <- biotic_list[[18]]
contrast_smv <- make_contrast_table(
  c("control", "control",
    "stress_0hpi_res", "stress_2hpi_res", "stress_4hpi_res",
    "stress_6hpi_res", "stress_8hpi_res"),
  c("stress_SMV-G7", "stress_SMV-L",
    "stress_0hpi_sus", "stress_2hpi_sus", "stress_4hpi_sus",
    "stress_6hpi_sus", "stress_8hpi_sus")
)

#----Combine tables of biotic stress-related samples to be used----
biotic_contrast_list <- list(
  aglycines = contrast_aglycines,
  fgraminearum = contrast_fgraminearum,
  foxysporum = contrast_foxysporum,
  fvirguliforme = contrast_fvirguliforme,
  hglycines = contrast_hglycines,
  mamp = contrast_mamp,
  mphaseolina = contrast_mphaseolina,
  pgregata = contrast_pgregata,
  ppachyrhizi = contrast_ppachyrhizi,
  psojae = contrast_psojae,
  psojae_gluc_elic = contrast_psojae_gluc_elic,
  rreniformis = contrast_rreniformis,
  rsolani = contrast_rsolani,
  slitura = contrast_slitura,
  smv = contrast_smv
)

biotic_samplelist <- list(
  aglycines = aglycines,
  fgraminearum = fgraminearum,
  foxysporum = foxysporum,
  fvirguliforme = fvirguliforme,
  hglycines = hglycines,
  mamp = mamp,
  mphaseolina = mphaseolina,
  pgregata = pgregata,
  ppachyrhizi = ppachyrhizi,
  psojae = psojae,
  psojae_gluc_elic = psojae_gluc_elic,
  rreniformis = rreniformis, 
  rsolani = rsolani,
  slitura = slitura,
  smv = smv
)
biotic_samplelist <- lapply(biotic_samplelist, function(x) {
  return(x[!duplicated(x$BioSample), ])
})


#----Export lists as .RData files----
save(
  abiotic_contrast_list, abiotic_samplelist,
  biotic_contrast_list, biotic_samplelist, 
  file = here("data", "stress_samples_for_DGE_analysis.RData"),
  compress="xz"
  )

