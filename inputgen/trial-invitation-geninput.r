########################################################################
## R Script to generate input files for invitation program
##
## CONSTANTS
##
## _/_/ Input files _/_/
## 1. csv list of chosen GPs with preferences weights
##  - CONST_FILE_GPLIST 
## 2. National opt out data by age/sex 
##    - CONST_FILE_DROPOUT_AGESEX 
## 3. National opt out data by GP
##    - CONST_FILE_DROPOUT_GP 
## 4. Invitations already sent, number already booked, csv file 
##    - CONST_FILE_ALREADYINVITED 
## 5. Invitations already requested by algorithm 
##    - CONST_FILE_ALREADYINVITED_ALGORITHM 
## _/_/ Input parameters _/_/
##  1. Date to use for national opt out (depends on latest file above)
##    - CONST_DROPDATE 
##  2. First (FALSE) or subsequent (TRUE) round - ie. needs updating for invitations already sent
##    - CONST_NOTFIRSTROUND 
##
## /_/_ Output /_/_
## 1. Number of patients in each GP practice by age / sex, for first run (input.txt)
## 2. Utility of patients by on age / sex / GP 
##
## _/_/ Other input files, not expected to change _/_/
## 1. expected-event-rate.csv : expected event rate 
##
## Info:
## Author: Adam Brentnall
## Last update: 20th July 2022
#########################################################################

myargs <- commandArgs(trailingOnly=TRUE)
if (length(myargs)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if(length(myargs)==1) {
    CONST_FILE_GPLIST <- myargs[1]
    CONST_ROUND <-1 
} else if(length(myargs)==2) {
    CONST_FILE_GPLIST <- myargs[1]
    CONST_FILE_ALREADYINVITED <- myargs[2] 
    CONST_ROUND <-2
}

## Libraries
library("tidyverse")
library("lubridate")


###############################
## Set up - input constants

## filenames for this site(s)
## current national drop out data
##https://digital.nhs.uk/data-and-information/publications/statistical/national-data-opt-out/
CONST_FILE_DROPOUT_AGESEX <- "configfiles/NDOP_age_gen_Sep_2021.csv"
CONST_FILE_DROPOUT_GP <- "configfiles/NDOP_reg_geo_Sep_2021.csv"
CONST_DROPDATE <- "01/09/2021" ## date to drop from
CONST_CHOOSEALLSITE <- 1 ## code might be used for multiple site (not intended for use in this way anymore though)

## uptake assumption file
CONST_FILE_UPTAKE <- "uptakefiles/uptake-gp-age-sex.csv"

##########################################
## Program
##########################################

##1. Setup

###########################################
##loop thru all sites
for(id_site in CONST_CHOOSEALLSITE){

    ## Load GP list for these sites
    testgp <- read_csv(CONST_FILE_GPLIST) %>%
        mutate(Rank = Overall_Ranked_Quintile) %>%
        select(ORG_CODE, ORG_NAME, F50to54, F55to59, F60to64, F65to69, F70to74, F75to77, M50to54, M55to59, M60to64, M65to69, M70to74, M75to77, IMD1, IMD2, IMD3, IMD4, IMD5, IMD_TOTAL, BoCaCoverage, LOCATION, POSTCODE.y, distance, Rank)

    testgp_locations <- unique(testgp$LOCATION)

    print(paste("Using location", testgp_locations[id_site]))

    testgp <- testgp %>%
        filter(LOCATION == testgp_locations[id_site])
    
   ## header names
    agesexnames <- c("F50to54", "F55to59", "F60to64", "F65to69", "F70to74", "F75to77", "M50to54", "M55to59", "M60to64", "M65to69", "M70to74", "M75to77")

    ## adjust capacity by expected drop out rate

    mydropout_agesex <- read_csv(CONST_FILE_DROPOUT_AGESEX) %>%
        filter(ACH_DATE==CONST_DROPDATE)

    my_overall <- mydropout_agesex %>%
        filter(GENDER=="All")%>%
        filter(AGE_BAND =="All") %>%
        mutate(CAPACITY = RATE/100) %>%
        select(CAPACITY) %>%
        unlist()  

    mymale_age <- mydropout_agesex %>%
        filter(GENDER=="Male")%>%
        filter(AGE_BAND %in% c("50-59", "60-69", "70-79")) %>%
        mutate(CAPACITY = (RATE/100) / my_overall) %>%
        select(CAPACITY) %>%
        unlist %>%
        rep(each=2)

    myfemale_age <- mydropout_agesex %>%
        filter(GENDER=="Female")%>%
        filter(AGE_BAND %in% c("50-59", "60-69", "70-79")) %>%
        mutate(CAPACITY = (RATE/100) / my_overall) %>%
        select(CAPACITY) %>%
        unlist %>%
        rep(each=2)

    ## gp opt out rates
    mydropout_gp <- read_csv(CONST_FILE_DROPOUT_GP) %>%
        filter(ACH_DATE==CONST_DROPDATE) %>%
        filter(LIST_SIZE !="Unallocated") %>%
        mutate(LIST_SIZE=as.numeric(LIST_SIZE))%>%
        mutate(OPT_OUT=as.numeric(OPT_OUT))%>%
        mutate(CAPACITY = OPT_OUT / LIST_SIZE) %>%
        select(GP_PRACTICE, CAPACITY)

    ##expected drop out by GP
    myexpdropout <- 1-cbind(t(c(myfemale_age) %*% t( mydropout_gp$CAPACITY)),
                                       t(c(mymale_age) %*% t( mydropout_gp$CAPACITY)) )

    myexpdropout_pos <- pmax(matrix(0.0000100001, ncol=12, nrow=nrow(myexpdropout)), myexpdropout)


    myexpcapacity<- data.frame(mydropout_gp,
                               myexpdropout_pos)
    
    colnames(myexpcapacity) <- c("ORG_CODE", "DROP_RATE", paste0("CAPACITY_", agesexnames))


    testgp2 <- left_join(testgp, myexpcapacity, by="ORG_CODE")

    newcapacity <- floor(
    (testgp2 %>% select(starts_with("CAPACITY"))) *
    (testgp2 %>% select(starts_with("F") | starts_with("M"))))

    ### make sure positive and greater than 0 - (latter needed for a technical reason in algorithm - generates a matrix for constrain in LP to ensure balance by age/sex/gp. If zero the constraint is that all zeroes in matrix. Could remove from matrix in algorithm, but more simple work around is to add a single patient, which won't have much effect on output in grand scheme) 
    newcapacity[,1:12] <- pmax(matrix(1, nrow=nrow(newcapacity), ncol=12), as.matrix(newcapacity))
    
    ## add absolute capacity projection
    testgp3 <- data.frame(testgp2 %>%
                          select(!((starts_with("F") | starts_with("M")))),
                          newcapacity)
    ##rename
    colnames(testgp3)[27:38] <- agesexnames

    testgp0 <- testgp
    testgp <- testgp3

#######################################################################
    ## date run for file names
    todayfiledate<- substr(str_replace_all(today(), "-", ""), 3,8)

############=############################
    ## 2. Input file 

    ## sites
    listsite <- unique(testgp$LOCATION)

    this_site_txt <- str_replace_all(as.character(listsite), ",", "") %>%
        str_replace_all(" ", "-") %>%
        str_replace_all("'","") %>%
        str_replace_all("\"", "") %>%
        str_replace_all("’","")

    ## if more than one site at this point (don't expect this as currently intended, but left over from previous defn)
    myinputs <- testgp %>%
        filter(LOCATION == listsite) %>%
        select( ORG_CODE, Rank, IMD_TOTAL, BoCaCoverage, M50to54, M55to59, M60to64, M65to69, M70to74, M75to77, F50to54, F55to59, F60to64, F65to69, F70to74, F75to77) %>%
        arrange(ORG_CODE)

    gps_bysite <- myinputs$ORG_CODE

    ## add info on what already happened - if second round+
    ## at moment only for one site (ie. myinputs length 1 - above to this point can work >1 site)

    if(CONST_ROUND>1){

        ##actual data, from central mailout
        mysite1_invited <- read_csv(CONST_FILE_ALREADYINVITED) %>%
            arrange(ORG_CODE) %>%
            filter(MAILOUT_TYPE=="CENTRAL") %>%
            filter(!is.na(GP_code))

        gp_invited <- mysite1_invited$ORG_CODE


        gp_book <- mysite1_invited %>%
            select(ORG_CODE, starts_with("Book")) %>%
            select(!Book_Total)

        gp_remain <- mysite1_invited %>%
            select(starts_with("REM")) 

        gp_total <- myinputs %>%
            arrange(ORG_CODE) %>%
            filter(ORG_CODE %in% gp_invited) %>%
            select(starts_with("M") | starts_with("F"))
        

        gp_remain[gp_remain=="NULL"]<-NA

        gp_remain_mtx<-matrix(as.numeric(unlist(data.frame(gp_remain))), ncol=ncol(gp_remain))

        gp_remain_prop <- gp_remain_mtx / gp_total

        ##theoretically can be more registered than we know about - set to max capacity 1 (could also correct the number registered data)
        gp_remain_prop[gp_remain_prop>1 & !is.na(gp_remain_prop)] <-1

        gp_invited_perc <- 1 - gp_remain_prop

        gp_invited_perc[is.na(gp_invited_perc)]<-0

        colnames(gp_invited_perc) <- paste0("PINV_", colnames(gp_invited_perc))
        
        gp_invited_perc <- gp_invited_perc %>%
            mutate(ORG_CODE=gp_invited)      
        
        gp_invited_summary <- left_join(gp_invited_perc, gp_book)

        myinput_file <- left_join(myinputs, gp_invited_summary, by="ORG_CODE")

        myinput_file[is.na(myinput_file)] <- 0.0        


        
    } else {
        
        myinput_file <- cbind(myinputs, matrix(0, ncol=24, nrow=nrow(myinputs)))

    }

    ##write the input file
    write_csv( path=paste0("genfiles/", todayfiledate, "_", this_site_txt, "_input.csv"), myinput_file)
    
    ## also need to output multiplication factor for algorithm- to scale number requested by to allow for opt outs
    testgp_multfactor <- 1/
        (testgp %>%
         select(starts_with("CAPACITY")))

    testgp_mfr <- data.frame(testgp %>%
                             select(ORG_CODE, LOCATION), testgp_multfactor)

    colnames(testgp_mfr)[3:14]<- agesexnames

    ## sites
    myinputs_mfr <- testgp_mfr %>%
            filter(LOCATION == listsite) %>%
            select( ORG_CODE, M50to54, M55to59, M60to64, M65to69, M70to74, M75to77, F50to54, F55to59, F60to64, F65to69, F70to74, F75to77) %>%        
            arrange(ORG_CODE)

    ##write the upscale files
    write_csv( path=paste0("genfiles/", todayfiledate, "_",
                           this_site_txt, "_upfactor.csv"), myinputs_mfr)

    

#############################################################################
    ## 3. Cost / utility by age/sex/gp
    ## expected event rate 
    myadcanrate <- read_csv("configfiles/expected-event-rate.csv") 

    imdgp<- testgp %>%
        select(ORG_CODE, Rank, starts_with("IMD")) %>%
        mutate(allIMD = IMD1+IMD2+IMD3+IMD4+IMD5) %>%
        mutate(IMDp1 = IMD1/allIMD) %>%
        mutate(IMDp2 = IMD2/allIMD) %>%
        mutate(IMDp3 = IMD3/allIMD) %>%
        mutate(IMDp4 = IMD4/allIMD) %>%
        mutate(IMDp5 = IMD5/allIMD) 

    imdp<- imdgp %>%
        select(ORG_CODE, starts_with("IMDp"))

    agegrps <- unique(myadcanrate$Age.group)

    sexgrp <- c(1,2)

    ngp<-nrow(imdgp)

    gps <- imdp %>% select(ORG_CODE) %>%unlist()  %>% as.character()

    ## Expected advanced cancer rate, assuming deprivation distribution constant by GP by age/sex
    fn.expt.inc <- function(inage, insex, ingp){
        
        thisrate <- myadcanrate %>%
            filter(Age.group==inage, Sex==insex) %>%
            select(rate) %>%
            unlist() %>%
            rev() #most deprived should be first
        

        thisexpt <- (imdp %>% select(starts_with("IMDp")) %>% as.matrix()) %*% matrix(thisrate)

        data.frame(ORG_CODE=ingp, age=rep(inage, length(thisexpt)), sex=rep(insex, length(thisexpt)), rate=thisexpt)

    }
     ##1=male,2=female

    counter<-0
    for(inage in agegrps){
        for(insex in sexgrp){
            counter<- counter+1
            if(counter>1){
                thisdta<- bind_rows(thisdta, fn.expt.inc(inage, insex, gps))
            }
            else thisdta <- fn.expt.inc(inage, insex, gps)
        }
    }

    thisdta2 <- left_join(thisdta, imdgp %>% select(ORG_CODE, Rank))

    thisdta2 <- thisdta2 %>%
        mutate(cost = 0.006/rate) %>%
        mutate(costrank = (cost) + seq(0,1000,by=10)[Rank])

    oututility <- thisdta2 %>%
        arrange(ORG_CODE, sex, age) %>%
        select(ORG_CODE, sex, age, costrank, rate)

    oututility_site <-  oututility %>%
            filter(ORG_CODE %in% gps_bysite) 


    ##write the utility input files for algorithm
    write_csv( path=paste0("genfiles/", todayfiledate, "_",
                           this_site_txt, "_utility.csv"),  oututility_site)


}
##END loop thru all sites (not expected any more)
###########################################

## and write in the prog dir, assumes 1 site per input file
write_csv( path="../prog/input/input.csv", myinput_file)
write_csv( path="../prog/input/utility.csv", oututility_site)
write_csv( path="../prog/input/upfactor.csv", myinputs_mfr)

## add expected uptake info
mydta <- read_csv(CONST_FILE_UPTAKE)
colnames(mydta)[1] <- "ORG_CODE"

thisgp <- testgp %>%
    select("ORG_CODE", Rank, starts_with("F"), starts_with("M")) 

thisgp_long <- thisgp %>%
    pivot_longer(
        col=(starts_with("M") | starts_with("F")),
        names_to=c("Sex", "AgeDes"),
        names_pattern= "(.)(.*)",
        values_to="capacity")

thisgp_long$Age<- as.character(as.integer(as.factor(thisgp_long$AgeDes)))

mysite <- thisgp %>%
    select(ORG_CODE, Rank) %>%
    left_join(mydta, on="ORG_CODE")

## gp without uptake data
missuptake<-is.na(mysite$Deprivation) | is.na(mysite$BowelUptake)

if(sum(missuptake)>0){
    
    avgupt<-mysite %>%
    filter(!is.na(M1)) %>%
    group_by(Rank) %>%
    summarise(across(starts_with("M")| starts_with("F"), mean)) %>%
    arrange(Rank) %>%
    select(starts_with("M"), starts_with("F"))
    
    missuptake_idx<-(1:nrow(mysite))[missuptake]

    print("/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_")
    print("Warning: some sites do not have deprivation / bowel data to predict uptake")

    for(idx in missuptake_idx){
        mysite[idx,6:17] <- avgupt[mysite$Rank[idx],]
        print(paste("GP: ", mysite$ORG_CODE[idx], " (priority rank ", mysite$Rank[idx], ")", sep=""))
        }
    print("/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_")

    }


## output for prog
mysite_input <- mysite %>%
    select(ORG_CODE, M1, M2, M3, M4, M5, M6, F1, F2, F3, F4, F5, F6) %>%
    arrange(ORG_CODE) 

mysite_input<-mysite_input %>%
    mutate(mysite_input %>% select(starts_with("M") | starts_with("F"))/100)



write_csv( path=paste0("genfiles/", todayfiledate, "_",
                           this_site_txt, "_uptake.csv"),  mysite_input)
write_csv( path="../prog/input/uptake.csv", mysite_input)


## plots
mysite_long <-   mysite %>%
    pivot_longer(
        col=(starts_with("M")| starts_with("F")),
        names_to=c("Sex","Age"),
        names_pattern= "(.)(.)",
        values_to="uptake")
                     

mysite_pdta<- left_join(mysite_long, thisgp_long, by=c("ORG_CODE", "Sex", "Age")) %>%
    mutate(RankAge = Rank.x + (as.integer(Age)-3.5)/12)

myggplot1<- 
    mysite_pdta %>%
    arrange(desc(capacity)) %>%
    filter(Sex=="M") %>%
    ggplot(aes(x=RankAge, y=uptake, size=capacity, color=AgeDes)) +
    geom_point(alpha=0.5) +
    scale_size(range = c(.1, 44), name="Capacity") +
    ylim(c(2,25)) +
    labs(title = "Predicted uptake by algorithm rank, age, GP",
         subtitle = "Male",
         caption = todayfiledate)

myggplot2<- 
    mysite_pdta %>%
    arrange(desc(capacity)) %>%
    filter(Sex=="F") %>%
    ggplot(aes(x=RankAge, y=uptake, size=capacity, color=AgeDes)) +
    geom_point(alpha=0.5) +
    scale_size(range = c(.1, 44), name="Capacity") +
    ylim(c(2,25)) +
    labs(title = "Predicted uptake by algorithm rank, age, GP",
         subtitle = "Female",
         caption = todayfiledate)

mysite_pdta_both<-
    mysite_pdta %>%
    group_by(ORG_CODE, Age) %>%
    summarise(capacity=sum(capacity), deprivation=Deprivation[1], rank=Rank.x[1], uptake=mean(uptake), rankage=RankAge[1], AgeDes=AgeDes[1], boweluptake=BowelUptake[1])



myggplot1a<- 
    mysite_pdta_both %>%
    arrange(desc(capacity)) %>%
    ggplot(aes(x=rankage, y=uptake, size=capacity, color=AgeDes)) +
    geom_point(alpha=0.5) +
    scale_size(range = c(.1, 44), name="Capacity") +
    ylim(c(2,25)) +
    labs(title = "Predicted uptake by algorithm rank, age, GP",
         subtitle = "Male+Female (mean)",
         caption = todayfiledate)


mysite_pdta_sum<-
    mysite_pdta_both %>%
    group_by(ORG_CODE) %>%
    summarise(capacity=sum(capacity), deprivation=deprivation[1], rank=rank[1], uptake=mean(uptake), boweluptake=boweluptake[1])



myggplot1b<- 
    mysite_pdta_sum %>%
    arrange(desc(capacity)) %>%
    ggplot(aes(x=rank, y=deprivation, size=capacity, color=ORG_CODE)) +
    geom_point(alpha=0.5) +
    scale_size(range = c(.1, 50), name="Total capacity") +
    labs(title = "Deprivation score by GP and algorithm rank",
         subtitle = "",
         caption = todayfiledate)


myggplot1c<- 
    mysite_pdta_sum %>%
    arrange(desc(capacity)) %>%
    ggplot(aes(x=rank, y=boweluptake, size=capacity, color=ORG_CODE)) +
    geom_point(alpha=0.5) +
    ylim(c(30,90)) +
    scale_size(range = c(.1, 50), name="Total capacity") +
    labs(title = "Bowel uptake by GP and algorithm rank",
         subtitle = "",
         caption = todayfiledate)




pdf("genfiles/uptake-report.pdf", width=16, height=8)
print(myggplot1a)
print(myggplot1b)
print(myggplot1c)
print(myggplot1)
print(myggplot2)
dev.off()
