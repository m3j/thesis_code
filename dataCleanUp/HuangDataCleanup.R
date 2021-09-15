library(tidyverse)

#include for housing free-range/cage
metadata <- Huang.Tabel.1 %>% mutate(Housing  = ifelse (location_name.farm. == "China:Guandong" | (location_name.farm. == "China:Hunan" & Feed == "not collected") , "free-range", "cage"))

#include breed based on farm and housing
metadata <- metadata %>% 
  mutate(Breed = case_when(location_name.farm. == "China:Hunan" & Housing == "cage" ~ "Local yellow-feather chickens",
                           location_name.farm. == "China:Hunan" & Housing == "free-range" ~ "Guangxi local chicken ",
                           location_name.farm. == "China:Shandong" ~ "Cobb 500",
                           location_name.farm. == "China:Shanxi" ~"Ross 308",
                           location_name.farm. == "China:Guandong" ~ "Yellow dwarf chicken",
                           location_name.farm. == "China:Henan" ~ "Hy-Line Variety Brown",
                           location_name.farm. == "China:Beijing" ~ "Arbor Acres broiler",
                            ))
#include Type layer/broiler
metadata <- metadata %>%  
  mutate(type = ifelse(Breed == "Hy-Line Variety Brown" | Breed == "Local yellow-feather chickens", "layer", "broiler"))

#clean days to only be the number (still as string)
metadata <- metadata %>%  
  mutate(age.days. = str_remove(age.days., "-day$"))


metadata <- metadata %>% 
   filter(!grepl('^pooled', Sample.replicate))

#save table.
write.csv(metadata,'/Users/Maja/Documents/DTU/thesis/Andet/Metadata/HuangMetadataRen.csv')
