library(dplyr)
library(tidyverse)
library(ggplot2)
library(readr)
library(readxl)
library(ggrepel)
library(scales)
library(purrr)
library(writexl)
library(janitor)
library(lubridate)

setwd("~/Desktop/")
#upload data of all samples ran
df_wells <- read.csv("~/Desktop/df_total.csv")

#upload data of all repeated samples
repeats<- read.csv("~/Desktop/repeat_total.csv")

#make list of no repeats with only final repeats
#to be used for final stats
no_repeat_last<-df_wells %>%
  group_by(SampleBarcode) %>%
  slice(n()) %>%
  ungroup() 

#make list of no repeats with the first result to be used for plate mapping
no_repeat_first<-df_wells %>%
  group_by(SampleBarcode) %>%
  slice(1) %>%
  ungroup() 

#generates autoscoring
auto_score <- df_wells %>%
  mutate(Autoscore= ifelse(CtN1 > 40 & CtN2 >40 & CtRP < 38,"Not Detected","QNS")) %>%
  mutate(Autoscore= ifelse(CtN1 <= 40 | CtN2 <= 40,"Positive 2019-nCoV",Autoscore)) %>%
  mutate(Agree = if_else((aecAutomaticResult== Autoscore),1,0)) %>%
  write_csv("~/Desktop/results_well_agreement.csv")


#Generates extraction batch RVL of individual wells
#set the order of the extractions
so_1<- c("H01", "G01","F01","E01","D01","C01","B01","A01","H02","G02","F02","E02","D02","C02","B02","A02")
so_2<- c("H03", "G03","F03","E03","D03","C03","B03","A03","H04","G04","F04","E04","D04","C04","B04","A04")
so_3<- c("H05", "G05","F05","E05","D05","C05","B05","A05","H06","G06","F06","E06","D06","C06","B06","A06")
so_4<- c("H07", "G07","F07","E07","D07","C07","B07","A07","H08","G08","F08","E08","D08","C08","B08","A08")
so_5<- c("H09", "G09","F09","E09","D09","C09","B09","A09","H10","G10","F10","E10","D10","C10","B10","A10")
so_6<- c("H11", "G11","F11","E11","D11","C11","B11","A11","H12","G12","F12","E12","D12","C12","B12","A12")
b1<-"01|02"
b2<-"03|04"
b3<-"05|06"
b4<-"07|08"
b5<-"09|10"
b6<-"11|12"

output_rvl<- function(set,bounds,da){
  
  da<- da %>%
    dplyr::filter(grepl(bounds, RunWell)) %>%
    arrange(match(RunWell, set))
  
  da <- da %>%
    replace_na(list(RVL=0)) %>%
    mutate(HOT_EXT = sum(RVL))
  
  return(da)
}

t1 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_1,b1,.)) %>%
  bind_rows()
t2 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_2,b2,.)) %>%
  bind_rows()
t3 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_3,b3,.)) %>%
  bind_rows()
t4 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_4,b4,.)) %>%
  bind_rows()
t5 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_5,b5,.)) %>%
  bind_rows()
t6 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_6,b6,.)) %>%
  bind_rows()

output<-bind_rows(t1,t2,t3,t4,t5,t6) %>%
  dplyr::select(SampleBarcode,HOT_EXT) %>%
  write_csv("~/Desktop/BatchRLV_per_ext.csv")


#Geneartes Tecan RVL for individual wells
so_1<- "A01|A02|A03|B01|B02|B03"
so_2<- "C01|C02|C03|D01|D02|D03"
so_3<- "E01|E02|E03|F01|F02|F03"
so_4<- "G01|G02|G03|H01|H02|H03"

so_5<- "A04|A05|A06|B04|B05|B06"
so_6<- "C04|C05|C06|D04|D05|D06"
so_7<- "E04|E05|E06|F04|F05|F06"
so_8<- "G04|G05|G06|H04|H05|H06"

so_9<- "A07|A08|A09|B07|B08|B09"
so_10<- "C07|C08|C09|D07|D08|D09"
so_11<- "E07|E08|E09|F07|F08|F09"
so_12<- "G07|G08|G09|H07|H08|H09"

so_13<- "A10|A11|A12|B10|B11|B12"
so_14<- "C10|C11|C12|D10|D11|D12"
so_15<- "E10|E11|E12|F10|F11|F12"
so_16<- "G10|G11|G12|H10|H11|H12"

output_rvl<- function(set,da){
  
  da<- da %>%
    dplyr::filter(grepl(set, RunWell)) %>%
    replace_na(list(RVL=0)) %>%
    mutate(Tecan = sum(RVL))
  
  return(da)
}

t1 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_1,.)) %>%
  bind_rows()
t2 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_2,.)) %>%
  bind_rows()
t3 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_3,.)) %>%
  bind_rows()
t4 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_4,.)) %>%
  bind_rows()
t5 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_5,.)) %>%
  bind_rows()
t6 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_6,.)) %>%
  bind_rows()
t7 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_7,.)) %>%
  bind_rows()
t8 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_8,.)) %>%
  bind_rows()
t9 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_9,.)) %>%
  bind_rows()
t10 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_10,.)) %>%
  bind_rows()
t11 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_11,.)) %>%
  bind_rows()
t12 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_12,.)) %>%
  bind_rows()
t13 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_13,.)) %>%
  bind_rows()
t14 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_14,.)) %>%
  bind_rows()
t15 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_15,.)) %>%
  bind_rows()
t16 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_16,.)) %>%
  bind_rows()


output<-bind_rows(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16) %>%
  dplyr::select(SampleBarcode,Tecan) %>%
  write_csv("~/Desktop/batch_tecan.csv")

#Generates horizontal RVL batch
so_1<- "A01|A02|A03|A04|A05|A06|A07|A08|A09|A10|A11|A12"
so_2<- "B01|B02|B03|B04|B05|B06|B07|B08|B09|B10|B11|B12"
so_3<- "C01|C02|C03|C04|C05|C06|C07|C08|C09|C10|C11|C12"
so_4<- "D01|D02|D03|D04|D05|D06|D07|D08|D09|D10|D11|D12"
so_5<- "E01|E02|E03|E04|E05|E06|E07|E08|E09|E10|E11|E12"
so_6<- "F01|F02|F03|F04|F05|F06|F07|F08|F09|F10|F11|F12"
so_7<- "G01|G02|G03|G04|G05|G06|G07|G08|G09|G10|G11|G12"
so_8<- "H01|H02|H03|H04|H05|H06|H07|H08|H09|H10|H11|H12"


output_rvl<- function(set,da){
  
  da<- da %>%
    dplyr::filter(grepl(set, RunWell)) %>%
    replace_na(list(RVL=0)) %>%
    mutate(horizontal = sum(RVL))
  
  return(da)
}

t1 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_1,.)) %>%
  bind_rows()
t2 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_2,.)) %>%
  bind_rows()
t3 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_3,.)) %>%
  bind_rows()
t4 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_4,.)) %>%
  bind_rows()
t5 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_5,.)) %>%
  bind_rows()
t6 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_6,.)) %>%
  bind_rows()
t7 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_7,.)) %>%
  bind_rows()
t8 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_8,.)) %>%
  bind_rows()


output<-bind_rows(t1,t2,t3,t4,t5,t6,t7,t8) %>%
  dplyr::select(SampleBarcode,horizontal) %>%
  write_csv("~/Desktop/horizontal_tecan.csv")

#Generates vertical RVL batch
so_1<- "H01|G01|F01|E01|D01|C01|B01|A01"
so_2<- "H02|G02|F02|E02|D02|C02|B02|A02"
so_3<- "H03|G03|F03|E03|D03|C03|B03|A03"
so_4<- "H04|G04|F04|E04|D04|C04|B04|A04"
so_5<- "H05|G05|F05|E05|D05|C05|B05|A05"
so_6<- "H06|G06|F06|E06|D06|C06|B06|A06"
so_7<- "H07|G07|F07|E07|D07|C07|B07|A07"
so_8<- "H08|G08|F08|E08|D08|C08|B08|A08"
so_9<- "H09|G09|F09|E09|D09|C09|B09|A09"
so_10<- "H10|G10|F10|E10|D10|C10|B10|A10"
so_12<- "H11|G11|F11|E11|D11|C11|B11|A11"
so_12<- "H12|G12|F12|E12|D12|C12|B12|A12"


output_rvl<- function(set,da){
  
  da<- da %>%
    dplyr::filter(grepl(set, RunWell)) %>%
    replace_na(list(RVL=0)) %>%
    mutate(vertical = sum(RVL))
  
  return(da)
}

t1 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_1,.)) %>%
  bind_rows()
t2 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_2,.)) %>%
  bind_rows()
t3 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_3,.)) %>%
  bind_rows()
t4 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_4,.)) %>%
  bind_rows()
t5 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_5,.)) %>%
  bind_rows()
t6 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_6,.)) %>%
  bind_rows()
t7 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_7,.)) %>%
  bind_rows()
t8 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_8,.)) %>%
  bind_rows()
t9 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_9,.)) %>%
  bind_rows()
t10 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_10,.)) %>%
  bind_rows()
t11 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_11,.)) %>%
  bind_rows()
t12 <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~output_rvl(so_12,.)) %>%
  bind_rows()

output<-bind_rows(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12) %>%
  dplyr::select(SampleBarcode,vertical) %>%
  write_csv("~/Desktop/vertical_tecan.csv")

#calculates the surrounding RVL and FC variables for every well
library(platetools)
pm = matrix(c(
  "A01","B01","C01","D01","E01","F01","G01","H01",
  "A02","B02","C02","D02","E02","F02","G02","H02",
  "A03","B03","C03","D03","E03","F03","G03","H03",
  "A04","B04","C04","D04","E04","F04","G04","H04",
  "A05","B05","C05","D05","E05","F05","G05","H05",
  "A06","B06","C06","D06","E06","F06","G06","H06",
  "A07","B07","C07","D07","E07","F07","G07","H07",
  "A08","B08","C08","D08","E08","F08","G08","H08",
  "A09","B09","C09","D09","E09","F09","G09","H09",
  "A10","B10","C10","D10","E10","F10","G10","H10",
  "A11","B11","C11","D11","E11","F11","G11","H11",
  "A12","B12","C12","D12","E12","F12","G12","H12"
),nrow=8,ncol=12)

#finding the surrounding wells RVL and fc 
find_surrounding <- function(df)  {
  fcs=c()
  w_sur = c()
  ewt = c()
  nst = c()
  dt = c()
  d_f <- df %>%
    dplyr::select(RunWell,RVL) %>%
    arrange(RunWell)
  
  #arrange data into matrix with padding
  final <-fill_plate(d_f, "RunWell") %>%
    group_by(RunWell) %>%
    slice(1)
  
  final <-final %>% replace(is.na(.), 0)
  plate <- plate_matrix(data = final$RVL, well = final$RunWell)
  mat.pad = rbind(NA, cbind(NA, plate, NA), NA)
  
  
  #iterate through matrix and caluculate surrounding well RVL and fold change from well of interest
  for(col in 2:13){
    for(row in 2:9) {
      surrounding <- c(mat.pad[(row-1):(row+1),(col+1)],mat.pad[(row-1):(row+1),(col-1)],mat.pad[(row+1),(col)],mat.pad[(row-1),(col)])
      surrounding <-surrounding[!is.na(surrounding)]
      w_sur <- c(w_sur,sum(surrounding))
      well <- mat.pad[row,col]
      e_w <-c(mat.pad[row,(col+1)],mat.pad[row,(col-1)])
      n_s <-c(mat.pad[(row+1),(col)],mat.pad[(row-1),(col)])
      diagonals <-c(mat.pad[(row-1),(col-1)],
                    mat.pad[(row-1),(col+1)],
                    mat.pad[(row+1),(col-1)],
                    mat.pad[(row+1),(col+1)])
      e_w <-e_w[!is.na(e_w)]
      n_s <-n_s[!is.na(n_s)]
      diagonals<-diagonals[!is.na(diagonals)]
      ewt <- c(ewt,sum(e_w))
      nst <- c(nst,sum(n_s))
      dt <- c(dt,sum(diagonals))
      
      if(well>0) {
        fc=sum(surrounding)/well
        fcs <- c(fcs,fc)
      }
      else{
        fc=0
        fcs <- c(fcs,fc)
      }
    }
  }
  
  plate_output <- tibble(c(pm),c(plate),c(w_sur),c(fcs),c(ewt),c(nst),c(dt)) %>%
    rename(RunWell = `c(pm)`,Well_RVL = `c(plate)`,Surrounding_RVL = `c(w_sur)`,Fold_change = `c(fcs)`,
           East_west = `c(ewt)`,North_south = `c(nst)`,Diagonals = `c(dt)`)
  total <- left_join(df,plate_output,by = "RunWell")
}

#run function
surrounding <- no_repeat_first %>%
  split(.$RunName) %>%
  map(~find_surrounding(.))

df_surrounding <- bind_rows(surrounding) %>%
  write_csv("~/Desktop/surrounding_df.csv")

#identifies samples that have been repeated
# this scores the repeats for incorporation
repeats <- repeats %>%
  mutate(Autoscore= ifelse(CtN1 > 40 & CtN2 >40 & CtRP < 38,"Not Detected","Invalid")) %>%
  mutate(Autoscore= ifelse(CtN1 <= 40 | CtN2 <= 40,"Positive 2019-nCoV",Autoscore)) %>%
  dplyr::select(RunName,SampleBarcode,CtN1,CtN2,CtRP,ExtractionReplicate,Autoscore,rack,Result) %>%
  mutate(ratio = (CtN1-CtN2))

repeats$RunDate = substr(repeats$RunName,1,nchar(repeats$RunName)-3)

#bring in false positive results based on extraction order
proximity <- df_surrounding %>%
  dplyr::select(SampleBarcode,Surrounding_RVL,Fold_change,East_west,North_south,Diagonals,RVL)

#separate out based on number of repeats
repeatsx2 <- repeats %>%
  group_by(SampleBarcode) %>%
  filter(n()==2)
repeatsx3 <- repeats %>%
  group_by(SampleBarcode) %>%
  filter(n()==3) %>%
  slice(c(1, n())) 
repeatsx4 <- repeats %>%
  group_by(SampleBarcode) %>%
  filter(n()==4) %>%
  slice(c(1, n()))

repeatsx234<- bind_rows(repeatsx2,repeatsx3,repeatsx4,)

#make tables for manual and automated based on results, not autoscore
repeats_results_lims <- repeatsx234 %>%
  mutate(changes = paste(Result,lead(Result))) %>%
  mutate(Extraction = lead(ExtractionReplicate)) %>%
  filter(!grepl('NA',changes))

#manual re-extraciton
repeats_results_lims %>%
  filter(Extraction >= 2) %>%
  filter(rack>=90000) %>%
  ungroup()%>%
  add_count(changes) %>%
  group_by(changes) %>%
  slice(n()) %>%
  summarise(n) %>%
  arrange(desc(n)) %>%
  adorn_totals() %>%
  write_csv("~/Desktop/summarzied_results_repeats_tubee_with_incon.csv")
#automated re-extraction
repeats_results_lims  %>%
  ungroup() %>%
  filter(Extraction >= 2) %>%
  filter(rack < 90000) %>%
  add_count(changes) %>%
  group_by(changes) %>%
  slice(n()) %>%
  summarise(n) %>%
  arrange(desc(n)) %>%
  adorn_totals() %>%
  write_csv("~/Desktop/summarzied_results_repeats_racke_with_incon.csv")


#find out how many samples were changed
repeatsx2_o <- repeatsx234 %>%
  mutate(changes = paste(Autoscore,lead(Autoscore))) %>%
  mutate(repeatN1 = lead(CtN1)) %>%
  mutate(repeatN2 = lead(CtN2)) %>%
  mutate(Extraction = lead(ExtractionReplicate)) %>%
  filter(!grepl('NA',changes)) %>%
  mutate(N1 = if_else(CtN1 <46 & repeatN1==46,"N1 lost",
                      if_else(CtN1==46 & repeatN1==46,"N1 never det",
                              if_else(CtN1 < 46 & repeatN1 < 46,"N1 repeat+","N1 gained"))))%>% 
  mutate(N2 = if_else(CtN2 <46 & repeatN2==46,"N2 lost",
                      if_else(CtN2==46 & repeatN2==46,"N2 never det",
                              if_else(CtN2 < 46 & repeatN2 < 46,"N2 repeat+","N2 gained")))) 

rvl_batch <- read.csv("~/Desktop/BatchRLV_per_ext.csv")
vertical<- read.csv("~/Desktop/vertical_tecan.csv")
horizontal<- read.csv("~/Desktop/horizontal_tecan.csv")
tecan_batch<- read.csv("~/Desktop/batch_tecan.csv")

rvl_batch<- left_join(rvl_batch,vertical,by="SampleBarcode")%>%
  left_join(horizontal,by="SampleBarcode") %>%
  left_join(tecan_batch,by="SampleBarcode")


#write file of linked repeats samples and FP based on order
output_e <- left_join(repeatsx2_o,proximity,by="SampleBarcode") %>%
  filter(Extraction >=2) %>%
  left_join(rvl_batch, by="SampleBarcode") %>%
  dplyr::select(SampleBarcode,CtN1,CtN2,CtRP,changes,repeatN1,repeatN2,RVL,Surrounding_RVL,Fold_change,HOT_EXT,East_west,North_south,Diagonals,ratio,rack,vertical,horizontal,Tecan,RunDate) %>%
  write_csv("~/Desktop/redo_extract_FP.csv")


#make large file for total prediction
no_repeat_predictions<-no_repeat_first %>%
  left_join(proximity,by="SampleBarcode") %>%
  left_join(rvl_batch,by="SampleBarcode") %>%
  dplyr::select(SampleBarcode,CtN1,CtN2,CtRP,RVL.x,Surrounding_RVL,Fold_change,HOT_EXT,East_west,North_south,Diagonals,Result,rack,vertical,horizontal,Tecan) %>%
  rename(RVL=RVL.x) %>%
  mutate(ratio = (CtN1-CtN2)) %>%
  replace_na(list(RVL=0,Surrounding_RVL=0,Fold_change=0,Ext_batch_RVL=0,East_west=0,North_south=0,Diagonals=0))

##calculate the totals for the paper
#summary for singles
results <-no_repeat_last %>%
  add_count(Result)%>%
  group_by(Result) %>%
  slice(n()) %>%
  summarise(n) %>%
  arrange(desc(n))%>%
  adorn_totals() %>%
  write_csv("~/Desktop/summarzied_results_final.csv")
#summary for singles
results <-no_repeat_first %>%
  add_count(Result)%>%
  group_by(Result) %>%
  slice(n()) %>%
  summarise(n) %>%
  arrange(desc(n))%>%
  adorn_totals() %>%
  write_csv("~/Desktop/summarzied_results_first_results.csv")

#summary for extraction repeats
#automated
repeatsx2_o %>%
  ungroup() %>%
  filter(Extraction >= 2) %>%
  filter(rack < 90000) %>%
  add_count(changes) %>%
  group_by(changes) %>%
  slice(n()) %>%
  summarise(n) %>%
  arrange(desc(n)) %>%
  adorn_totals() %>%
  write_csv("~/Desktop/summarzied_results_repeats_racke.csv")

#manual
repeatsx2_o %>%
  filter(Extraction >= 2) %>%
  filter(rack>=90000) %>%
  ungroup()%>%
  add_count(changes) %>%
  group_by(changes) %>%
  slice(n()) %>%
  summarise(n) %>%
  arrange(desc(n)) %>%
  adorn_totals() %>%
  write_csv("~/Desktop/summarzied_results_repeats_tubee.csv")


#log transforms data and generates tables for GBM modeling
#program that makes 4 groups of rack-pcr, rack-extract, tube-pcr, tube-extract
#values are log transformed
extract <- read_csv("~/Desktop/redo_extract_FP.csv") %>%
  filter(changes == "Positive 2019-nCoV Positive 2019-nCoV"|changes == "Positive 2019-nCoV Not Detected") %>%
  mutate(True_positives = if_else(changes == "Positive 2019-nCoV Positive 2019-nCoV",1,0)) %>%
  mutate(type="Extraction repeat") %>%
  replace_na(list(RVL=0)) %>%
  replace_na(list(vertical=0))


extract_df <- extract %>%
  mutate(RVL=log(RVL+min(extract[extract$RVL>0,"RVL"]))) %>%
  mutate(Surrounding_RVL=log(Surrounding_RVL+min(extract[extract$Surrounding_RVL>0,"Surrounding_RVL"]))) %>%
  mutate(HOT_EXT=log(HOT_EXT+min(extract[extract$HOT_EXT>0,"HOT_EXT"]))) %>%
  mutate(East_west=log(East_west+min(extract[extract$East_west>0,"East_west"]))) %>%
  mutate(North_south=log(North_south+min(extract[extract$North_south>0,"North_south"]))) %>%
  mutate(Diagonals=log(Diagonals+min(extract[extract$Diagonals>0,"Diagonals"]))) %>%
  mutate(Fold_change=log(Fold_change+min(extract[extract$Fold_change>0,"Fold_change"]))) %>%
  mutate(horizontal=log(horizontal+min(extract[extract$horizontal>0,"horizontal"]))) %>%
  mutate(vertical=log(vertical+min(extract[extract$vertical>0,"vertical"]))) %>%
  mutate(Tecan=log(Tecan+min(extract[extract$Tecan>0,"Tecan"]))) %>%
  dplyr::select(-repeatN1) %>%
  dplyr::select(-repeatN2) %>%
  dplyr::select(-SampleBarcode)

#make TP a factor (may need be saved seperately for some modeling)
total <- extract_df
total %>%
  filter(True_positives==1) %>%
  write_csv("~/Desktop/total_repeats_calc_TP_extraction.csv")
total %>%
  filter(True_positives==0)%>%
  write_csv("~/Desktop/total_repeats_calc_FP_extraction.csv")

#seperate out rack and tubes and remove extraction batch from rack
#kingfisher
rack_pre<- total %>%
  filter((rack>=80000 & rack < 90000)) %>%
  dplyr::select(-HOT_EXT) %>%
  dplyr::select(-rack) %>%
  mutate(date=ymd(RunDate)) %>%
  filter(date < as.Date("2020-07-01")) %>%
  mutate(type= if_else(type == "PCR repeat",1,0))
#Tecan
rack_post<- total %>%
  filter((rack>=80000 & rack < 90000)) %>%
  dplyr::select(-HOT_EXT) %>%
  dplyr::select(-rack) %>%
  mutate(date=ymd(RunDate)) %>%
  filter(date >= as.Date("2020-07-01")) %>%
  mutate(type= if_else(type == "PCR repeat",1,0))
#manual
tube <-total %>%
  filter(rack>=90000) %>%
  dplyr::select(-Tecan) %>%
  dplyr::select(-vertical) %>%
  dplyr::select(-horizontal) %>%
  dplyr::select(-rack) %>%
  mutate(type = if_else(type == "PCR repeat",1,0))

#put into the 3 groups
tube_extract <- tube %>%
  filter(type==0) %>%
  dplyr::select(-type) %>%
  write_csv("~/Desktop/tube_extract.csv")

rack_pre_extract <- rack_pre %>%
  filter(type==0) %>%
  dplyr::select(-type)%>%
  write_csv("~/Desktop/rack_pre_extract.csv")

rack_post_extract <- rack_post %>%
  filter(type==0) %>%
  dplyr::select(-type)%>%
  write_csv("~/Desktop/rack_post_extract.csv")

#total lists
tube_total <- read_csv("~/Desktop/surrounding_df.csv") %>%
  filter(Result != "Not Detected") %>%
  filter((CtN1 <= 40 | CtN2 <=40) & CtRP < 38) %>%
  filter(rack>=90000) %>%
  replace_na(list(RVL=0)) %>%
  left_join(rvl_batch, by="SampleBarcode")

#calculate out log transformation for manual
tube_total <-tube_total %>%
  mutate(RVL=log(RVL+min(extract[extract$RVL>0,"RVL"]))) %>%
  mutate(Surrounding_RVL=log(Surrounding_RVL+min(extract[extract$Surrounding_RVL>0,"Surrounding_RVL"]))) %>%
  mutate(HOT_EXT=log(HOT_EXT+min(extract[extract$HOT_EXT>0,"HOT_EXT"]))) %>%
  mutate(East_west=log(East_west+min(extract[extract$East_west>0,"East_west"]))) %>%
  mutate(North_south=log(North_south+min(extract[extract$North_south>0,"North_south"]))) %>%
  mutate(Diagonals=log(Diagonals+min(extract[extract$Diagonals>0,"Diagonals"]))) %>%
  mutate(Fold_change=log(Fold_change+min(extract[extract$Fold_change>0,"Fold_change"]))) %>%
  mutate(horizontal=log(horizontal+min(extract[extract$horizontal>0,"horizontal"]))) %>%
  mutate(vertical=log(vertical+min(extract[extract$vertical>0,"vertical"]))) %>%
  mutate(Tecan=log(Tecan+min(extract[extract$Tecan>0,"Tecan"]))) %>%
  dplyr::select(-aecAutomaticResult) %>%
  dplyr::select(-RunType)%>%
  dplyr::select(-RunName) %>%
  write_csv("~/Desktop/tube_positive.csv")

#automated data with surrounding RVLs
rack_total <- read_csv("~/Desktop/surrounding_df.csv") %>%
  filter(Result != "Not Detected") %>%
  filter((CtN1 <= 40 | CtN2 <=40) & CtRP < 38) %>%
  filter((rack>=80000 & rack < 90000)) %>%
  replace_na(list(RVL=0)) %>%
  left_join(rvl_batch, by="SampleBarcode")

rack_total$RunDate = substr(rack_total$RunName,1,nchar(rack_total$RunName)-3)

rack_total_pre<- rack_total %>%
  mutate(date=ymd(RunDate)) %>%
  filter(date < as.Date("2020-07-01"))

rack_total_post<- rack_total %>%
  mutate(date=ymd(RunDate)) %>%
  filter(date >= as.Date("2020-07-01"))

#calculate log transformation for KingFisher
rack_total_pre <-rack_total_pre %>%
  mutate(RVL=log(RVL+min(extract[extract$RVL>0,"RVL"]))) %>%
  mutate(Surrounding_RVL=log(Surrounding_RVL+min(extract[extract$Surrounding_RVL>0,"Surrounding_RVL"]))) %>%
  mutate(East_west=log(East_west+min(extract[extract$East_west>0,"East_west"]))) %>%
  mutate(North_south=log(North_south+min(extract[extract$North_south>0,"North_south"]))) %>%
  mutate(Diagonals=log(Diagonals+min(extract[extract$Diagonals>0,"Diagonals"]))) %>%
  mutate(Fold_change=log(Fold_change+min(extract[extract$Fold_change>0,"Fold_change"]))) %>%
  mutate(horizontal=log(horizontal+min(extract[extract$horizontal>0,"horizontal"]))) %>%
  mutate(vertical=log(vertical+min(extract[extract$vertical>0,"vertical"]))) %>%
  mutate(Tecan=log(Tecan+min(extract[extract$Tecan>0,"Tecan"]))) %>%
  dplyr::select(-aecAutomaticResult) %>%
  dplyr::select(-RunType) %>%
  dplyr::select(-RunName) %>%
  write_csv("~/Desktop/rack_positive_pre.csv")

#calculate log transformation for Tecan
rack_total_post <-rack_total_post %>%
  mutate(RVL=log(RVL+min(extract[extract$RVL>0,"RVL"]))) %>%
  mutate(Surrounding_RVL=log(Surrounding_RVL+min(extract[extract$Surrounding_RVL>0,"Surrounding_RVL"]))) %>%
  mutate(East_west=log(East_west+min(extract[extract$East_west>0,"East_west"]))) %>%
  mutate(North_south=log(North_south+min(extract[extract$North_south>0,"North_south"]))) %>%
  mutate(Diagonals=log(Diagonals+min(extract[extract$Diagonals>0,"Diagonals"]))) %>%
  mutate(Fold_change=log(Fold_change+min(extract[extract$Fold_change>0,"Fold_change"]))) %>%
  mutate(horizontal=log(horizontal+min(extract[extract$horizontal>0,"horizontal"]))) %>%
  mutate(vertical=log(vertical+min(extract[extract$vertical>0,"vertical"]))) %>%
  mutate(Tecan=log(Tecan+min(extract[extract$Tecan>0,"Tecan"]))) %>%
  dplyr::select(-aecAutomaticResult) %>%
  dplyr::select(-RunType) %>%
  dplyr::select(-RunName) %>%
  write_csv("~/Desktop/rack_positive_post.csv")


#Modeling for paper
library(caret)

#function to generate GBM input data
regression <-function(name) {
  #make  data set for logistic regression
  df_main<-read_csv(paste("~/Desktop/",name,".csv", sep=""))
  
  total <- df_main %>%
    dplyr::select(-changes) %>%
    dplyr::select(-ratio) %>%
    dplyr::select(-RunDate)
  
  if(name=="rack_post_extract") {
    total <- total %>% dplyr::select(-date)
  } 
  
  if(name=="rack_pre_extract"){
    total <- total %>% dplyr::select(-date)
  }
  
  total$True_positives<-as.factor(total$True_positives)
  
  list(total,name)
}

list[total,name]<-regression("rack_post_extract")
#For each model, replace with the following: manual(tube_extract),Tecan(rack_post_extract),KingFisher(rack_pre_extract)

set.seed(90)
control <- rfeControl(functions = rfFuncs,
                      method = "repeatedcv",
                      repeats = 10,
                      verbose = FALSE)

#make all predictors
outcome <- "True_positives"
predictors<-names(total)[!names(total) %in% outcome]
#find which predictors are important
set.seed(4)
tp_pred<- rfe(total[,predictors], total$True_positives,
              rfeControl = control)

print(tp_pred)
predictors(tp_pred)

if(name=="rack_pre_extract") {
  predictors <- c("CtN2","CtN1","RVL","East_west","North_south","Diagonals","Fold_change")
} 
if(name=="rack_post_extract") {
  predictors <- c("CtN2","CtN1","RVL","East_west","North_south","Diagonals", "Fold_change","Tecan","horizontal","vertical")
} 

if(name=="tube_extract"){
  predictors <- c("CtN2","CtN1","RVL","East_west","North_south","Diagonals", "Fold_change","HOT_EXT")
}

#model
set.seed(198)
train<-trainControl(method = "repeatedcv",
                    number = 10,
                    repeats=10)
model_glm<-train(total[,predictors],total$True_positives,method='gbm',preProcess = c("center","scale"),trControl=train)
print(varImp(model_glm))
print(summary(model_glm))
predictions<-predict.train(object=model_glm,total[,predictors],type="raw")
table(predictions)
print(confusionMatrix(predictions,total$True_positives))

#generates partial dependency plots
library(pdp)           # for partial dependence plots
library(vip)           # for variable importance plots

set.seed(101)
rwb <- colorRampPalette(c("red", "white", "blue"))
if(name=="rack_post_extract") {
  n1n2 <- pdp::partial(model_glm, pred.var = c("CtN1","CtN2"),prob=TRUE)
  pdp1 <- plotPartial(n1n2, contour = TRUE, col.regions = rwb)
  pdp1
  
  n2ew <- pdp::partial(model_glm, pred.var = c("RVL","horizontal"),prob=TRUE)
  pdp2 <- plotPartial(n2ew, contour = TRUE, col.regions = rwb)
  pdp2
  
  n1ew <- pdp::partial(model_glm, pred.var = c("horizontal","Tecan"),prob=TRUE)
  pdp3 <- plotPartial(n1ew, contour = TRUE, col.regions = rwb)
  pdp3
  grid.arrange(pdp1,pdp2,ncol=2)
  grid.arrange(pdp1,pdp3,ncol=2)
  
} 

if(name=="rack_pre_extract") {
  n1n2 <- pdp::partial(model_glm, pred.var = c("CtN1","CtN2"),prob=TRUE)
  pdp1 <- plotPartial(n1n2, contour = TRUE, col.regions = rwb)
  pdp1
  
  n2ew <- pdp::partial(model_glm, pred.var = c("RVL","CtN2"),prob=TRUE)
  pdp2 <- plotPartial(n2ew, contour = TRUE, col.regions = rwb)
  pdp2
  
  n2ew <- pdp::partial(model_glm, pred.var = c("Surrounding_RVL","CtN1"),prob=TRUE)
  pdp3 <- plotPartial(n2ew, contour = TRUE, col.regions = rwb)
  pdp3
  
  grid.arrange(pdp1,pdp2,pdp3,ncol=2)
  
} 
if(name=="tube_extract") {
  
  n1n2 <- pdp::partial(model_glm, pred.var = c("CtN1","CtN2"),prob=TRUE)
  pdp1 <- plotPartial(n1n2, contour = TRUE, col.regions = rwb)
  pdp1
  
  hen2 <- pdp::partial(model_glm, pred.var = c("HOT_EXT","CtN2"),prob=TRUE)
  pdp2 <- plotPartial(hen2, contour = TRUE, col.regions = rwb)
  pdp2
  
  fcn2 <- pdp::partial(model_glm, pred.var = c("Fold_change","CtN2"),prob=TRUE)
  pdp3 <- plotPartial(fcn2, contour = TRUE, col.regions = rwb)
  pdp3
  
  ewn2 <- pdp::partial(model_glm, pred.var = c("CtN1","RVL"),prob=TRUE)
  pdp4 <- plotPartial(ewn2, contour = TRUE, col.regions = rwb)
  pdp4
  
  grid.arrange(pdp1,pdp2,ncol=2)
  grid.arrange(pdp3,pdp4,ncol=2)
  ew<-pdp::partial(model_glm, pred.var = "East_west",prob=TRUE)
  n1<-pdp::partial(model_glm, pred.var = "CtN1",prob=TRUE)
  n2<-pdp::partial(model_glm, pred.var = "CtN2",prob=TRUE)
  fc<-pdp::partial(model_glm, pred.var = "Fold_change",prob=TRUE)
  hot<-pdp::partial(model_glm, pred.var = "HOT_EXT",prob=TRUE)
  rvl<-pdp::partial(model_glm, pred.var = "RVL",prob=TRUE)
}


#bring in data generated from no repeats first-> surrounding data frame that has been log transformed for pos/inconclusive results
rack_post_total<-read_csv("~/Desktop/rack_positive_post.csv")
rack_pre_total<-read_csv("~/Desktop/rack_positive_pre.csv")
tube_total<-read_csv("~/Desktop/tube_positive.csv")
if(name=="tube_extract"){x=0.15
}
if(name=="rack_pre_extract"){x=0.1
}
if(name=="rack_post_extract"){x=0.25
}

#make new 2x2 with new cut-off
test<-predict.train(model_glm,type="prob")
train_table<-table(total$True_positives, test$`1` > x)
print(train_table)
test %>%
  bind_cols(total)%>%
  write_csv("~/Desktop/Covid/test.csv")

#takes either tube first repeat data and make predictions, save to total_tube/total_rack
if(name=="tube_extract"){
  predi_ves<-predict(model_glm,type="prob", newdata = tube_total)%>%
    bind_cols(tube_total)%>%
    write_csv("~/Desktop/Covid/total_tube.csv")
}

if(name=="rack_pre_extract"){
  predi_ves<-predict(model_glm,type="prob", newdata = rack_pre_total)%>%
    bind_cols(rack_pre_total)%>%
    write_csv("~/Desktop/Covid/total_rack_pre.csv")
}

if(name=="rack_post_extract"){
  predi_ves<-predict(model_glm,type="prob", newdata = rack_post_total)%>%
    bind_cols(rack_post_total)%>%
    write_csv("~/Desktop/Covid/total_rack_post.csv")
}

predi_ves1<-predi_ves %>%
  filter(predi_ves$`1`<= x)


#make list of all that were repeated, not just reextracted
df_wells <- read.csv("~/Desktop/df_total.csv")
repeats_all<-df_wells %>%
  group_by(SampleBarcode) %>%
  filter(n()>1) %>%
  ungroup()

#take difference from ones we have already repeated
out<-as_tibble(setdiff(predi_ves1$SampleBarcode,repeats_all$SampleBarcode))

# of samples extra that we predict to be FP
print(length(out$value))

