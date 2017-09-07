library(XML)
library(tidyverse)
library(stringr)

## There are 540 pharmacological actions
## All of these are stored in the substance (pharmacological action) file
## 422 are also in the descriptor file
## All of the descriptorui's from the pharma table are in the subst_code table
## The converse is not true.

## functions
nrowreturn <- function(x) {
  print(nrow(x))
  x}

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

## Save tables
saveRDS(scr, file= "Scratch_data/supplementaryconceptrecord.Rds")
saveRDS(tree, file= "Scratch_data/mesh_tree_numbers.Rds")
saveRDS(subst_code, file = "Scratch_data/substances_codes.Rds")
saveRDS(substance_description, file = "Scratch_data/substances_text.Rds")


## Load terms required from rxnorm and eutils searching
my_terms <- readRDS("../../Fellowship/Trial_identify/clinical_trials_august_2017/Data/Automatic_atc_rxnorm_mesh.Rds")
my_terms <- my_terms %>%
  filter(!is.na(MSH))

srch <- my_terms %>% 
  select(MSH) %>% 
  filter(!is.na(MSH)) %>% 
  mutate(rec_type = str_sub(MSH, 1, 1)) %>% 
  distinct() # 719 records

## Identify tree number for specific MESH IDs
## Direct join for D-codes (all drug codes)
my_tree_d <- srch %>% 
  filter(rec_type == "D") %>% 
  nrowreturn() %>% 
  inner_join(tree, by = c("MSH" = "descriptorui"))

## ConceptID related join
my_tree_c <- srch %>% 
  filter(rec_type == "C") %>% 
  nrowreturn() %>% 
  inner_join(scr,  by = c("MSH" = "supplementalrecordui")) %>% 
  nrowreturn() %>% 
  inner_join(tree, by = "descriptorui") %>% 
  nrowreturn()

## Check if have a tree for all MSH codes
my_tree <- bind_rows(my_tree_c, my_tree_d)
setdiff(my_terms$MSH, my_tree$MSH) # one missing only C063008, not found on eutils either


#### Find pharmacological action for all drugs
my_pa_direct1 <-  srch %>% 
  filter(MSH %in% subst_code$descriptorui_pharm_action540) %>% 
  mutate(descriptorui_pharm_action540 = MSH)

my_pa_direct2 <-  srch %>% 
  inner_join(subst_code, by = c("MSH" = "descriptorui")) 
my_pa_direct <- bind_rows(my_pa_direct1, my_pa_direct2)

setdiff(my_terms$MSH, my_pa_direct$MSH)

my_pa_indirect_c <- srch %>% 
  anti_join(my_pa_direct, by = "MSH") %>% 
  nrowreturn() %>% 
  filter(rec_type == "C") %>% 
  nrowreturn() %>% 
  inner_join(scr, by = c("MSH" = "supplementalrecordui")) %>%
  nrowreturn() %>% 
  inner_join(subst_code, by = "descriptorui") %>% 
  nrowreturn() 

my_pa_indirect_d <- srch %>% 
  anti_join(my_pa_direct, by = "MSH") %>% 
  nrowreturn()%>% 
  filter(rec_type == "D") %>% 
  nrowreturn() %>% 
  inner_join(scr, by = c("MSH" = "descriptorui")) %>% 
  inner_join(scr, by = "supplementalrecordui") %>% 
  inner_join(subst_code, by = "descriptorui")

my_pa_all <- bind_rows(my_pa_direct, my_pa_indirect_c, my_pa_indirect_d)
no_pa <- my_terms %>% 
  anti_join(my_pa_all) %>% 
  distinct(MSH, .keep_all = TRUE)

## Return strings for no PAs
# 41 MESH terms with no PAs, try using tree to identify these
my_tree_no_pa <- my_tree %>%
  filter(MSH %in% no_pa$MSH) %>% 
  select(MSH, tree) %>% 
  distinct()
sum(!duplicated(my_tree_no_pa$MSH)) # all but 1 has tree, as expected

## Detect tree lengths
tree_pa_link <- tree %>%
  inner_join(subst_code) %>% 
  select(-descriptorui) %>% 
  distinct()

## Identify where any of the trees match, bu trucating the trees in the
## MESH databse to the trees I ahve in my dataset
tree_lengths <- sort(nchar(my_tree_no_pa$tree) %>%  unique())
tree_lengths_res <- map(tree_lengths, function (tree_length) {
  tree_pa_link %>%
    mutate(tree_link = substr(tree, 1, tree_length)) %>% 
    filter(tree_link %in% my_tree_no_pa$tree) %>% 
    distinct()
})
names(tree_lengths_res) <- paste0("len", tree_lengths)
tree_lengths_res <- bind_rows(tree_lengths_res, .id = "tree_len")
tree_lengths_res <- distinct(tree_lengths_res)

my_pa_indirect_tree <- tree_lengths_res %>% 
  inner_join(my_tree_no_pa, by = c("tree_link" = "tree")) %>% 
  select(MSH, descriptorui_pharm_action540) %>% 
  unique()
sum(!duplicated(my_pa_indirect_tree$MSH)) # got an additional 25 with pharmacological actions

my_pa_all <- bind_rows(my_pa_all, my_pa_indirect_tree)

setdiff(my_terms$MSH, my_pa_all$MSH)
no_pa <- my_terms %>% 
  filter(!MSH %in% my_pa_all$MSH) %>% 
  group_by(MSH) %>% 
  summarise(str_srchd = paste(unique(str_to_lower(na.omit(str_srchd))), collapse = " | "))
write_csv(no_pa, path = "Data/assign_manual_pharacologic_action.csv")

## 16 items, manually assigned 7 pharmacological actions manually
rvd_no_pa <- read_csv(file = "Data/assignED_manual_pharacologic_action.csv")
rvd_no_pa <- rvd_no_pa %>% 
  filter(!is.na(descriptorui)) %>% 
  select(MSH, descriptorui) 
rvd_no_pa2 <- str_split(rvd_no_pa$descriptorui, pattern = fixed("|"))
names(rvd_no_pa2) <- rvd_no_pa$MSH 
 
my_pa_manual <- rvd_no_pa2 %>% 
  stack() %>% 
  set_names(c("descriptorui", "MSH")) %>% 
  inner_join(subst_code, by = "descriptorui") %>% 
  select(-descriptorui) %>% 
  distinct()
sum(!duplicated(my_pa_manual$MSH))

## Combine all mesh pharmacalogical actions
my_pa_all <- bind_rows(direct = my_pa_direct,
                       indirect_c = my_pa_indirect_c,
                       indirect_d = my_pa_indirect_d,
                       indirect_tree = my_pa_indirect_tree,
                       manual = my_pa_manual,
                       .id = "link_to_pa") %>% 
  select(MSH, link_to_pa, descriptorui_pharm_action540) %>% 
  distinct(MSH, descriptorui_pharm_action540, .keep_all = TRUE)
setdiff(my_terms$MSH, my_pa_all$MSH) # only 9 without MESh pharmacological action

my_pa_all <- my_pa_all %>%
  inner_join(substance_description)

saveRDS(my_tree, file = "Data/MSH_to_tree.Rds")
saveRDS(my_pa_all, file = "Data/MSH_to_pharmacological_action.Rds")
rm(list = ls())
my_tree <- readRDS(file = "Data/MSH_to_tree.Rds")
my_pa   <- readRDS(file = "Data/MSH_to_pharmacological_action.Rds")
