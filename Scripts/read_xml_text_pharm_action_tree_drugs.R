library(XML)
library(tidyverse)
library(stringr)


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

## Write small sample of code 
write_lines(desc[1:1000], path = "temp.txt")


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
concept <- ExtractXML(x, "ConceptList/Concept/ConceptUI")
children <-  ExtractXML(x,
            "ConceptList/Concept/ConceptRelationList/ConceptRelation[@RelationName='NRW']/*")
parents <-  ExtractXML(x,
            "ConceptList/Concept/ConceptRelationList/ConceptRelation[@RelationName='BRD']/*")
synonyms <-  ExtractXML(x,
            "ConceptList/Concept/ConceptRelationList/ConceptRelation[@RelationName='REL']/*")
rm(x)

## Select only variables with a drug action
all_lists <- list(tree_num, pharmui, pharmstring, concept, children, parents, synonyms)
all_lists <- map(all_lists, set_names, nm = ui)
rm(tree_num, pharmui, pharmstring, concept, children, parents, synonyms)
names(all_lists) <- c("tree", "pharmui", "pharmstring", "concept", "children", "parents",
                           "synonyms")
pharm_not_null<- map_lgl(all_lists$pharmui, ~ length(.x) >= 1)
mean(pharm_not_null)

## Select only pharm action data
all_lists_pharm <- map(all_lists, ~ .x[pharm_not_null])
all_dfs_pharm <- map(all_lists, stack)
all_dfs_pharm <- map2(all_dfs_pharm, names(all_lists), ~ set_names(.x, c(.y, "uid")))

## Select tree numbers where branch from pharm action ones
tree_pharm <- all_lists_pharm$tree %>% 
  stack() %>%  
  set_names(c("tree", "descriptor_id"))

tree <- all_lists$tree %>% 
  stack() %>%  
  set_names(c("tree", "descriptor_id"))
  
max(nchar(tree$tree))
trim_branches <- -seq(5, 48, 4)
tree_prune <- map(trim_branches, ~ str_sub(tree$tree, 1, .x))
tree_prune <- do.call(c, tree_prune)
names(tree_prune) <- tree$tree 
tree_prune <- data.frame(tree_prune = tree_prune)
tree_prune$tree_original <- tree$tree

tree_prune <- tree_prune %>% 
  filter(tree_prune != "") %>% 
  distinct(tree_prune, .keep_all = TRUE)
tree_prune2 <- tree_prune %>% 
  inner_join(tree, by = c("tree_prune" = "tree"))


map(all_lists, length) %>%  unique() # should be single value

## Merge into dataset where have tree and pharm class for all ids



pharm <- reduce(all_dfs, inner_join)
pharm <- pharm %>%
  rename(descriptorui = uid) %>% 
  distinct()

pharm_lng <- pharm %>%
  gather(key = "relation_type", value = "conceptid", concept, children, parents, synonyms) %>% 
  select(tree, descriptorui, pharmstring, conceptid) %>% 
  distinct()

my_terms <- readRDS("../../Fellowship/Trial_identify/clinical_trials_august_2017/Data/Automatic_atc_rxnorm_mesh.Rds")
srch <- my_terms %>% 
  select(MSH) %>% 
  na.omit() %>% 
  distinct()

a <- setdiff(srch$MSH, pharm_lng$descriptorui)
my_tree
my_tree <- stack(my_tree)

substance <- readRDS("data/pa2017.rds")
x <- xmlParse(substance)
x <- xmlRoot(x)
x <- getNodeSet(x, "PharmacologicalAction")

subst_code <- ExtractXML(x, search_string = "PharmacologicalActionSubstanceList/Substance/RecordUI")
subst_code <- unlist(subst_code)
b <- setdiff(srch$MSH, subst_code)
p <- setdiff(srch$MSH, names(pharm_check))
pharm_check <- all_lists$pharmstring[subst_code]
pharm_check <- pharm_check[!map_lgl(pharm_check, is.null)]

write_lines(substance[1:500], "temp2.txt")
