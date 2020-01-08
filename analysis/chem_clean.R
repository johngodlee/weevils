# Prepare a clean .csv of biochemical fingerprint data 
# John Godlee (johngodlee@gmail.com)
# 2019-01-08

# Preamble ----
# Remove old crap
rm(list=ls())
#dev.off()

# Set working directory to the location of the source file
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages
library(readxl)
library(dplyr)

# Import data ----
raw <- readxl::read_xlsx("~/data/Common Garden metabolomic NIR data_extract_for_JG.xlsx",
  sheet = 2)

# Clean data ----

tidy <- raw %>%
  dplyr::select(tree = Tree, 9:53, 58, 65, 68, 72,
    a_pinene_m = `alpha pinene m`,
    a_pinene_p = `alpha pinene p`,
    camphene_m = `camphene m`,
    camphene_p = `camphene p`,
    b_pinene_m = `beta pinene m`,
    b_pinene_p = `beta pinene p`,
    carene_p_3 = `3 carene p`,
    limonene_m = `limonene m`,
    limonene_p = `limonene p`,
    p_cymene = `p cymene`,
    terpinolene = `Terpinolene`,
    bornyl_acet_m = `bornyl acetate m`,
    b_elemene = `beta elemene`,
    t_caryophyllene = `trans caryophyllene`,
    a_humulene = `alpha humulene`,
    germacrene_d = `germacrene D`,
    result = `Reults Type`,                                  
    N_perw = `N (%w)`,                                       
    C_perw = `C (%W)`,                                       
    CT = `CT (Pine Equiv)`,                              
    phen_mgg = `Total Phenolics in Plant material (mg/g)`,     
    ox_phen_mgg = `Oxidisable Phenolics in Plant material (mg/g)`,
    glucose_mgg = `Glucose (mg/g DW)`,                            
    fructose_mgg = `Fructose (mg/g DW)`,                           
    sugar_mgg = `Total Sugars (mg/g DW)`,                       
    ADF_per = `ADF %`,                                        
    dry_mass_per = `dry mass %`,                                   
    mois_per = `moisture %`
    )

P_all <- tidy %>%
  dplyr::select(starts_with("P", ignore.case = FALSE)) %>%
  rowSums()

# Terpenes are named where a standard was found and just "TU" where no standard
TU_all <- tidy %>%
  dplyr::select(starts_with("TU", ignore.case = FALSE), 
    "a_pinene_m",
    "a_pinene_p",
    "camphene_m",
    "camphene_p",
    "b_pinene_p",
    "b_pinene_m",
    "carene_p_3",
    "limonene_m",
    "p_cymene",
    "limonene_p",
    "terpinolene",
    "bornyl_acet_m",
    "b_elemene",
    "t_caryophyllene",
    "a_humulene",
    "germacrene_d") %>%
  rowSums()

tidy$P_all <- P_all
tidy$TU_all <- TU_all

write.csv(tidy, "data/chem.csv", row.names = FALSE)
