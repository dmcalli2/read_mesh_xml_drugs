library(XML)
library(tidyverse)
library(stringr)

## There are 540 pharmacological actions
## All of these are stored in the substance (pharmacological action) file
## 422 are also in the descriptor file
## All of the descriptorui's from the pharma table are in the subst_code table
## The converse is not true.


## File names
# Record Type	Full File	Size	ZIP file	Size	GZ file	Size
# Descriptors	desc2017.xml	286 MB	desc2017.zip	15 MB	desc2017.gz	15 MB
# Supplementary Records	supp2017.xml	571 MB	supp2017.zip	39 MB	supp2017.gz	39 MB
# Qualifiers	qual2017.xml	274 KB
# Substances with Pharmacologic Action	pa2017.xml

## See https://www.nlm.nih.gov/mesh/xmlmesh.html for discussion of format

# Structure goes as follows. Each descriptor can have many concepts
# and each conept can have many terms

# Descriptor
#    Concept
#       Term

## Treenumbers and pharmacological action goe with the descriptor
## Concept ids, names, terms and relationships go with concepts

# Read and compress files then delete text files to save on space
# file_list <- c("desc2017.xml", "pa2017.xml", "qual2017.xml", "supp2017.xml")
# map(file_list, function(x) {
#   a <- read_lines(paste0("Data/", x))
#   x_out <- str_replace(x, "xml", "rds")
#   saveRDS(a, file = paste0("Data/", x_out))
# })

## Read in data, having compressed into R format
desc <- readRDS("data/desc2017.rds")
substance <- readRDS("data/pa2017.rds")
supp <- readRDS("data/supp2017.rds")

## Write small sample of code 
write_lines(desc[1:1000], path = "scratch_data/desc.xml")
write_lines(substance[1:1000], path = "scratch_data/subs.xml")
write_lines(supp[1:1000], path = "Scratch_data/supp.xml")

## Functions
ExtractXML <- function (my_nodeset, search_string) {
  ## Take list of nodes (a nodeset) and extact tree numbers
  a <- lapply(my_nodeset, function(my_node) getNodeSet(my_node, search_string))
  lapply(a, function(x_item) unique(sapply(x_item, xmlValue)))
}

## Parse tree into each record
x <- xmlParse(desc)
x <- xmlRoot(x)
x <- getNodeSet(x, "DescriptorRecord")
rm(desc)

## Extract treenumber, pharmacological action, etc
ui <- ExtractXML(x, search_string = "DescriptorUI")
ui <- unlist(ui)
tree_num <- ExtractXML(x, search_string = "TreeNumberList/TreeNumber")
pharmui <- ExtractXML(x, search_string = "PharmacologicalActionList/PharmacologicalAction/DescriptorReferredTo/DescriptorUI")
pharmstring <- ExtractXML(x, search_string = "PharmacologicalActionList/PharmacologicalAction/DescriptorReferredTo/DescriptorName")
rm(x)
names(pharmstring) <- ui

## Create dataframe
all_lists <- list(tree_num, pharmui, pharmstring)
all_lists <- map(all_lists, set_names, nm = ui)
names(all_lists) <- c("tree", "pharmui", "pharmstring")
rm(tree_num, pharmui, pharmstring)

all_dfs_pharm <- map(all_lists, stack)
all_dfs_pharm <- map2(all_dfs_pharm, names(all_lists), ~ set_names(.x, c(.y, "descriptorui")))
pharm <- reduce(all_dfs_pharm[c("pharmui", "pharmstring")], inner_join) %>%
  distinct()
pharm <- pharm %>% 
  distinct()
tree <- all_dfs_pharm$tree

## 422 strings pharmstrings of different actions
## 3470 pharmui, ranging from 1 to 22 per action, some Uis link to drug lists
## Must be more synonyms
## the descriptorui is for each specific agent
## No overlap between pharmui and descriptorUI
## pharmUI does not have a descriptor record, makes sense as not on the tree
pharm_unq <- pharm %>% 
  group_by(pharmstring) %>% 
  count()
rm(all_lists, all_dfs_pharm)

## Lookup some supplement IDs in pharm table
## Looks like this is the file for mapping the IDs from other database to get pharm IDs and tree numbers
ids <- c(SupplementalRecordUI = "C060863",
                        HeadingMappedTo1 = "D001388", # note my numbering
                        HeadingMappedTo2 = "D007294",
                        HeadingMappedTo3 = "Q000031",
                        HeadingMappedTo4 = "D050112",
                        IndexingInformation1 = "D001617", # note my numbering
                        IndexingInformation2 = "Q000037")
ids <- data.frame(name_id = names(ids), descriptorui = ids)
ids %>% 
  left_join(tree) %>% 
  left_join(pharm)
# Maps to tree for all D codes
# Maps to pharmacological action for only one of these


## Extract codes from supplemental concept records
x <- xmlParse(supp)
x <- xmlRoot(x)
x <- getNodeSet(x, "SupplementalRecord[@SCRClass = '1']")

supplementalrecordui <- ExtractXML(x, search_string = "SupplementalRecordUI")
supplementalrecordui <- unlist(supplementalrecordui)
headingmappedto <- ExtractXML(x,
                              "HeadingMappedToList/HeadingMappedTo/DescriptorReferredTo/DescriptorUI")
names(headingmappedto) <- supplementalrecordui
scr <- stack(headingmappedto)
names(scr) <- c("descriptorui", "supplementalrecordui")
scr <- scr %>%
  mutate(asterisk = str_sub(descriptorui, 1, 1) == "*",
         descriptorui = str_replace(descriptorui, ("^\\*"), ""))
rm(x, supplementalrecordui, headingmappedto)

## Extract codes from substance list for pharmacological classes, there are 540 of these
## Each has a list of descriptorus for different drugs, allowing you to
## relate drugs within a class, call this table substance_description
x <- xmlParse(substance)
x <- xmlRoot(x)
x <- getNodeSet(x, "PharmacologicalAction") # 540 pharmacological actions

subst_dui <- ExtractXML(x, search_string = "DescriptorReferredTo/DescriptorUI")
subst_name <- ExtractXML(x, search_string = "DescriptorReferredTo/DescriptorName")
subst_dui <- unlist(subst_dui)
subst_name <- unlist(subst_name)
substance_description <- data.frame(descriptorui_pharm_action540 = subst_dui, descriptorname = subst_name)
rm(subst_name)

subst_code <- ExtractXML(x, search_string = "PharmacologicalActionSubstanceList/Substance/RecordUI")
names(subst_code) <- subst_dui
subst_code <- stack(subst_code)
names(subst_code) <- c("descriptorui", "descriptorui_pharm_action540")

## Clean-up
rm(subst_combine)
rm(p1, s1, ids, smrms)
rm(a_tree, a_pharm, subst_dui, supp, ui, x)
rm(substance)

## Load terms required from rxnorm and eutils searching
my_terms <- readRDS("../../Fellowship/Trial_identify/clinical_trials_august_2017/Data/Automatic_atc_rxnorm_mesh.Rds")
srch <- my_terms %>% 
  select(MSH) %>% 
  na.omit() %>% 
  distinct() # 719 records

## Search in main database, supplemental records and substances with pharmacological action database
inpharm <- intersect(srch$MSH, pharm$descriptorui) # 397
inscr1 <- intersect(srch$MSH,  scr$supplementalrecordui) # 288
inscr2 <- intersect(srch$MSH, scr$headingmappedto) # 288
inscr_both <- union(inscr1, inscr2) ## 576
inscr_inpharm <- union(inscr_both, inpharm) ## 700
## Identify number of codes with one of 540 pharmacolgical action classes
insub1 <- intersect(srch$MSH, substance$descriptorui) # 3
insub2 <- intersect(srch$MSH, substance$recordui) # 630
insub_both <- union(insub1, insub2) # 633

## Identify number of pharmacological action codes which are in pharm
inpharmsub1 <- intersect(substance$descriptorui, pharm$descriptorui) # 2
inpharmsub2 <- intersect(substance$recordui,     pharm$descriptorui) # 2605
inpharmsub3 <- intersect(substance$descriptorui, pharm$pharmui) # 422
inpharmsub4 <- intersect(substance$recordui,     pharm$pharmui) # 2605
inpharmsub <- reduce(list(inpharmsub1, inpharmsub2, inpharmsub3, inpharmsub4), union) # 3027
check <- union(inpharmsub2, inpharmsub4) # 2605


## Link scr data to srch
scr_noasterisk <- scr %>% 
  mutate(headingmappedto = str_replace_all(headingmappedto, ("^\\*"), ""))

## Note all are distinct
srch2 <- srch %>% 
  left_join(scr_noasterisk, by = c("MSH" = "headingmappedto")) %>% 
  left_join(scr_noasterisk, by = c("MSH" = "supplementalrecordui")) %>% 
  gather(key = "scr", value = "scr_id", - MSH, na.rm = TRUE) 

## Link substance data to srch
substance2 <- substance %>% 
  select(substance_descriptorui = descriptorui,
         substance_recordui = recordui)

srch3 <- srch2 %>% 
  left_join(substance2, by = c("MSH" = "substance_recordui")) %>% 
  left_join(substance2, by = c("scr_id" = "substance_recordui")) %>% 
  gather(key = "substance_type", value = "substance_descriptorui", substance_descriptorui.x,
         substance_descriptorui.y, na.rm = TRUE) %>% 
  select(-substance_type) %>% 
  distinct()

# gather srch3 to make merge easier
srch3b <- srch3 %>%
  select(-scr) %>% 
  gather(key = "code_type", "ui", -MSH, na.rm = TRUE) %>% 
  distinct(MSH, ui, .keep_all = TRUE)

## Link pharm data to search
# Note zero overlap between pharmui and descriptorui
intersect(pharm$pharmui, pharm$descriptorui)

pharm_codes <- pharm %>% 
  select(descriptorui, pharmui, tree_num )
nrowreturn <- function(x) {
  print(nrow(x))
  x}
pharmcodeonly <- pharm %>% 
  select(pharmui, pharmstring) %>% 
  distinct()
treecodeonly <- pharm %>% 
  select(descriptorui, tree) %>% 
  distinct()

srch4 <- srch3b %>%
  left_join(pharm_codes, by = c("MSH" = "descriptorui")) %>% 
  nrowreturn() %>% 
  left_join(pharm_codes, by = c("MSH" = "pharmui")) %>% 
  nrowreturn() %>% 
  left_join(pharm_codes, by = c("ui" = "descriptorui")) %>% 
  nrowreturn() %>% 
  left_join(pharm_codes, by = c("ui" = "pharmui")) %>% 
  nrowreturn()
map(srch4, ~ all(is.na(.x)))

pharm_lkp <- srch4 %>% 
  select(MSH, pharmui.x, pharmui.y) %>% 
  gather(key = "pharm_codes", value = "pharmui", -MSH, na.rm = TRUE) %>% 
  select(-pharm_codes) %>% 
  distinct()

tree_lkp <- srch4 %>% 
  select(MSH,  descriptorui = descriptorui.y) %>% 
  distinct()

## Identifies all ids, can now link this id
## to any of the ids values in the other dataset
all_ids <- srch4 %>%
  select(-code_type) %>% 
  gather(key = "code_type", value = "ui", - MSH, na.rm = TRUE) %>% 
  select(-code_type) %>% 
  distinct()
all_ids <- all_ids %>% 
  add_row(MSH = unique(all_ids$MSH), ui = unique(all_ids$MSH)) %>% 
  distinct()

rm(srch, srch2, srch3, srch3b, srch4)
rm(pharm_codes, pharm_lkp, pharm2, pharmcodeonly, scr, scr_noasterisk)


## Trees
tree <- all_ids %>% 
  inner_join(pharm %>% 
               select(descriptorui, tree) %>% 
               distinct(),
             by = c("ui" = "descriptorui")) %>% 
  select(-ui) %>% 
  distinct()

