# add the lybrary
library(tidyverse)
library(readxl)
library(XLConnect)
library(enrichR)
require("gplots")
library(VennDiagram)

# filter the top disregulatd no matter whether they are shared or not
# this will be run over the list of lncRNA
# from the listDEG_test
# test<-listDEG_test[[1]]
# test%>%
#   group_by(Regulation)%>%
#   mutate(r=rank(desc(`Fold Change`),ties.method = "first"))%>%
#   select(log2FC,r)%>%
#   filter(r<11)%>%
#   arrange(desc(log2FC))

sheet_test <- excel_sheets("../Data Analysis Folder/All Comparison mRNAs.xlsx")
#####################################################################################################
# filter only the comparison of interest (exclude HRMECS)
# only some comparison are interesting
sheet2_test <- sheet_test[c(1:8)]

# save the table of comparison
listDEG_test <- lapply(sheet2_test,function(x){
  read_excel("../Data Analysis Folder/All Comparison mRNAs.xlsx",sheet = x,skip = 17)%>%
    # save the first 18 column
    .[,1:18]%>%
    # make the probename unique
    distinct()%>%
    # add the variable log2FC
    mutate(log2FC=ifelse(Regulation=="up",yes = log2(`Fold Change`),no = -log2(`Fold Change`)))
})

# add names
names(listDEG_test) <- sheet2_test

# sort in two separated list up and down
listDEG_up_test<- listDEG_test[str_detect(names(listDEG_test),pattern = "up")]
listDEG_down_test<- listDEG_test[str_detect(names(listDEG_test),pattern = "down")]

# put together up and down regualted in a single table
listDEG_tot_test <- pmap(list(listDEG_up_test,listDEG_down_test),function(x,y){
  bind_rows(x,y)
})

# fix the naming further by removing the up label
col_name2 <- names(listDEG_tot_test)%>%
  str_sub(start = 4,end = -1) 

names(listDEG_tot_test) <- col_name2

# this approah of the common one do not works


list_disregulated <- lapply(listDEG_tot_test,function(x){
  df<-x
  df%>%
    arrange(desc(abs(log2FC)))%>%
    .[1:300,]
})
#####################################################################################################

#####################################################################################################
# save the list of gene ordered by fold change from most disregulated to less disregulated
# create the xlsx files
# create the file
my_book <- loadWorkbook("test.xlsx",create = T)
name_v <- names(list_disregulated)
for(i in seq_along(list_disregulated)){
  page <- list_disregulated[[i]]
  name <- name_v[i]
  XLConnect::createSheet(object = my_book,name =name)
  writeWorksheet(object = my_book,data = page,sheet = name)
}
  
saveWorkbook(object = my_book,file = "mRNA_top_300.xlsx")

#####################################################################################################

#####################################################################################################
# make the enrichR qeury using the gene symbol

dbs <- c("TRANSFAC_and_JASPAR_PWMs","ARCHS4_TFs_Coexp")

common_top300 <- lapply(list_disregulated,function(x){
  unique(x$GeneSymbol)
})
# show the common genes in the top 300
#
#####################################################################################################

###################################################################################################
# make the query per comparison
list_enriched <- lapply(common_top300,function(x){
  enrichr(x, dbs)
})

# each element of the list will contain another list with data form both the databases.

# # from the list extract the TRANSFACT/JASPAR data and filter only for humand TF
# test <- list_enriched[[1]][[1]]
# 
# test%>%
#   tbl_df()%>%
#   filter(str_detect(Term,pattern = "human"))

# is IRF7 anywhere
lapply(list_enriched,function(x){
  x[[1]]%>%
    tbl_df()%>%
    filter(str_detect(Term,pattern = "human"))%>%
    # in each table split the column Overlap into present and total
    separate(Overlap,into = c("present","total"),sep = "/")%>%
    filter(str_detect(Term,pattern = "IRF"))
})
#####################################################################################################

#####################################################################################################
# try Reinhold approach of counting the total number of genes regulated by the transcription factors

tot_TF <- lapply(list_enriched,function(x){
  x[[1]]%>%
    tbl_df()%>%
    filter(str_detect(Term,pattern = "human"))%>%
    # in each table split the column Overlap into present and total
    separate(Overlap,into = c("present","total"),sep = "/")
})%>%
  bind_rows(.id = "contrast")%>%
  group_by(Term)%>%
  summarise(tot = sum(as.numeric(present)))%>%
  arrange(desc(tot))
# confirm whether the most aboundant TF are also the one that regulates the most genes

# to tot_TF add the present (how many gene a specific TF regulates)
tot_TF2 <- lapply(list_enriched,function(x){
  x[[1]]%>%
    tbl_df()%>%
    filter(str_detect(Term,pattern = "human"))%>%
    # in each table split the column Overlap into present and total
    separate(Overlap,into = c("present","total"),sep = "/")
})%>%
  bind_rows(.id = "contrast")%>%
  select(Term,total)%>%
  group_by(Term)%>%
  summarise(total=mean(as.numeric(total)))



# merge the two dataset to make a scatter plot
left_join(tot_TF,tot_TF2,"Term")%>%
  dplyr::rename(sum_present=tot,total_regulated=total)%>%
  ggplot(aes(total_regulated,sum_present))+geom_point()+xlab("total gene regulated by TF")+ylab("total by gene count in the top 300")

left_join(tot_TF,tot_TF2,"Term")%>%
  dplyr::rename(sum_present=tot,total_regulated=total)%>%
  write_csv("TF_rank_Reinhold_method.csv")

# this approach crlearly show that we are just picking the transcriptioon factor with the most aboundant number of regulated genes
#####################################################################################################

#####################################################################################################
# in my opinion a better approach is to select the enriched transcription factor
lapply(list_enriched,function(x){
  x[[1]]%>%
    tbl_df()
})
# all of the list contains ~200 transcription factor I will pick the top 20 of each list to see the common one
top20_enriched <- lapply(list_enriched,function(x){
  x[[1]]%>%
    tbl_df()%>%
    filter(str_detect(Term,pattern = "human"))%>%
    arrange(P.value)%>%
    .[1:20,]%>%
    select(Term)%>%
    .[[1]]
})
# see how they are organized

# enriched common is a list with the genes per group
VENN.LIST_tot <- top20_enriched

# choose a coclor palette
color_red = colorRampPalette(c("red","orange","yellow"),space="rgb")
# the data contain 4 groups, therefore I define 4 colors (evenly distributed)
red<-color_red(4)

# save the name of the groups
new_name <- names(VENN.LIST_tot)

# prodce the object to then plot the diagram
venn.plot_tot <- venn.diagram(VENN.LIST_tot,NULL ,alpha=0.7, cex = 2, cat.fontface=4, category.names=new_name, main="TFs",fill=red)


# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram

grid.draw(venn.plot_tot)

# To get the list of gene present in each Venn compartment we can use the gplots package
a_tot <- venn(VENN.LIST_tot, show.plot=FALSE)
# By inspecting the structure of the a object created, 
# you notice two attributes: 1) dimnames 2) intersections
# We can store the intersections in a new object named inters
inters_tot <- attr(a_tot,"intersections")



# if I use all of them is not that useful.
enriched_common <- lapply(list_enriched,function(x){
  x[[1]]%>%
    tbl_df()%>%
    filter(str_detect(Term,pattern = "human"))%>%
    .[[1]]
})


###################################################################################################
# see the common intersections between transcription factor
# start the venn diagram

# enriched common is a list with the genes per group
VENN.LIST_tot <- enriched_common

# choose a coclor palette
color_red = colorRampPalette(c("red","orange","yellow"),space="rgb")
# the data contain 4 groups, therefore I define 4 colors (evenly distributed)
red<-color_red(4)

# save the name of the groups
new_name <- names(VENN.LIST_tot)

# prodce the object to then plot the diagram
venn.plot_tot <- venn.diagram(VENN.LIST_tot,NULL ,alpha=0.7, cex = 2, cat.fontface=4, category.names=new_name, main="TFs",fill=red)


# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram

grid.draw(venn.plot_tot)

# To get the list of gene present in each Venn compartment we can use the gplots package
a_tot <- venn(VENN.LIST_tot, show.plot=FALSE)
# By inspecting the structure of the a object created, 
# you notice two attributes: 1) dimnames 2) intersections
# We can store the intersections in a new object named inters
inters_tot <- attr(a_tot,"intersections")