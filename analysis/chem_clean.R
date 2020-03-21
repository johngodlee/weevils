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
raw <- readxl::read_xlsx("data/Common Garden metabolomic NIR data_extract_for_JG.xlsx",
  sheet = 2)

# Clean data ----

tidy <- raw %>%
  dplyr::select(tree = Tree, 9:53, 58, 65, 68, 72,
    terp_a_pinene_m = `alpha pinene m`,
    terp_a_pinene_p = `alpha pinene p`,
    terp_camphene_m = `camphene m`,
    terp_camphene_p = `camphene p`,
    terp_b_pinene_m = `beta pinene m`,
    terp_b_pinene_p = `beta pinene p`,
    terp_carene_p_3 = `3 carene p`,
    terp_limonene_m = `limonene m`,
    terp_limonene_p = `limonene p`,
    terp_p_cymene = `p cymene`,
    terp_terpinolene = `Terpinolene`,
    terp_bornyl_acet_m = `bornyl acetate m`,
    terp_b_elemene = `beta elemene`,
    terp_t_caryophyllene = `trans caryophyllene`,
    terp_a_humulene = `alpha humulene`,
    terp_germacrene_d = `germacrene D`,
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

# Phenolics, named "P*"
P_all <- tidy %>%
  mutate_at(vars(starts_with("P", ignore.case = FALSE)),
    .funs = list(std = ~(as.vector(scale(.,))))) %>%
  dplyr::select(matches("^P.*_std$")) %>%
  rowSums(na.rm = FALSE)

# Terpenes are named where a standard was found and "TU*" where no standard
TU_all <- tidy %>%
  dplyr::select(starts_with("TU", ignore.case = FALSE), 
    starts_with("terp_", ignore.case = FALSE)) %>%
  mutate_all(
    .funs = list(std = ~(as.vector(scale(.))))) %>%
  dplyr::select(ends_with("_std", ignore.case = FALSE)) %>%
  rowSums(na.rm = FALSE)

tidy$P_all <- P_all
tidy$TU_all <- TU_all

write.csv(tidy, "data/chem.csv", row.names = FALSE)
