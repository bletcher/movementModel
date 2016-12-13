
# to rebuild library:
# install.packages("devtools")
# devtools::install_github('Conte-Ecology/westBrookData/getWBData')

library(getWBData)
library(dplyr)
library(lubridate)



cdWB <- createCoreData(sampleType = "electrofishing", #"stationaryAntenna","portableAntenna"), 
                       whichDrainage = "west",
                       columnsToAdd=c("sampleNumber","river","riverMeter","survey",'observedLength','observedWeight')) %>%  
  addTagProperties( columnsToAdd=c("cohort","species","dateEmigrated","sex","species")) %>%
  dplyr::filter( !is.na(tag), area %in% c("trib","inside","below","above") ) %>% 
  createCmrData( maxAgeInSamples=20, inside=F) %>% 
  addSampleProperties() %>%
  addEnvironmental() %>%
  addKnownZ() %>%
  fillSizeLocation()

# get sites table
sitesIn <- data.frame(tbl(conDplyr,"data_sites") )
sites <- sitesIn %>% filter(is.na(quarter) & !is.na(quarter_length) & drainage == 'west') %>% select(-quarter)
sites$section <- as.numeric(sites$section)
############# 2_prepare data

library(dplyr)
library(lubridate)


# install dev version to fix NA problem with lag()
#if (packageVersion("devtools") < 1.6) {
#  install.packages("devtools")
#}
#devtools::install_github("hadley/lazyeval")
#devtools::install_github("hadley/dplyr")

# sample 2.5 = fyke net. not sure if we should keep it - prob with obs model
#cdWB <-  filter(cdWB, sampleNumber != 2.5 & sampleNumber != 10.1)

# some formatting fixes
cdWB$sectionOriginal <- cdWB$section
cdWB$section <- as.numeric( cdWB$section )
cdWB$inside <- ifelse( cdWB$section %in% 1:47 | cdWB$survey == "stationaryAntenna", T, F ) 

cdWB$year <- year(cdWB$detectionDate)
cdWB$yday <- yday(cdWB$detectionDate)

cdWB$riverOrdered <- factor(cdWB$river,levels=c('west brook', 'wb jimmy', 'wb mitchell',"wb obear"),labels=c("west brook","wb jimmy","wb mitchell","wb obear"), ordered=T)

cdWB <- cdWB %>%
  group_by(tag) %>%
  # arrange(tag,sampleNumber) %>%
  mutate( lagSection = lead(section),
          distMoved = section - lagSection,
          minSample = min(sampleNumber),
          maxSample = max(sampleNumber)) %>%
  ungroup()

cdWB$moveDir <- ifelse( cdWB$section == cdWB$lagSection, 0, ifelse( cdWB$section > cdWB$lagSection, 1,-1 ) )

cdWB$drainage <- "west"

# merge in riverMeter for sections
cdWB <- left_join(cdWB, sites, by = c("river","section","area")) 
cdWB$riverMeter <- ifelse( cdWB$survey == "shock" | cdWB$survey == "portableAntenna", cdWB$river_meter, cdWB$riverMeter )


save(cdWB, file='/home/ben/dataForModels/cdForMovementModelWB.RData')
cdWB2 <- cdWB%>%filter(sampleNumber==20)
save(cdWB2, file='/home/ben/movementModel/moveModelWithData/WB/temp.RData')







######
load(file='/home/ben/movementModel/moveModelWithData/cdForMovementModelWB.RData')

dat <- cdWB %>% filter(species=="bkt",distMoved<48,distMoved>0,enc==1)
counts <- hist(dat$distMoved, binwidth = 100, plot=F)
eFit <- lm( counts$counts+1 ~ lambda*exp(-lambda*counts$mids) ) #log(counts$counts+1)~counts$mids)

e <- function(x) {exp(eFit$coefficients[1]+eFit$coefficients[2]*x)}

ggplot(dat , aes(distMoved)) +
  geom_histogram(bins=200) +
  stat_function(fun = e, colour = "red")


cd <- cdWB %>% filter(cohort == 2003)












#cdWBHold <- cdWB

#cdWB <- cdWBHold %>% 
# filter(tag %in% c(  '257c683e2f', '00088d1c80', '00088d2e1f')) %>% 
#  select(tag,riverMeter,detectionDate,survey,sizeForGraph,observedLength)

counts <- cdWB %>% group_by(tag) %>% summarize(n=n()) %>% arrange(desc(n))
cdWB <- left_join(cdWB,counts) %>%
  filter(n > 1)  # remove single observations 


save(cdWB, file='/home/ben/shinyApps/pitAntenna/WB/cdForShiny.RData')

counts <- counts %>% filter(n > 1)
save(counts, file='/home/ben/shinyApps/pitAntenna/WB/counts.RData')

sampDateRange <- cdWB %>% group_by(sampleNumber) %>% 
                   summarize(minDate = min(detectionDate), maxDate = max(detectionDate), medDate = median(detectionDate))
save(sampDateRange, file='/home/ben/shinyApps/pitAntenna/WB/sampDateRange.RData')

countsByDay <- cdWB %>%
  filter(survey %in% c("stationaryAntenna")) %>%
  group_by(date=as.POSIXct(as.Date(detectionDate)),tag,survey,riverMeter) %>%
  summarize( nForInd = n() ) %>%
  ungroup() %>%
  group_by(date,survey,riverMeter) %>%
  summarize( n = n() )

portableSampleDays <- cdWB %>% filter(survey=="portableAntenna") %>% distinct(date = as.Date(detectionDate))

save(countsByDay, portableSampleDays, file='/home/ben/shinyApps/pitAntenna/WB/countsByDay.RData')

# for testing
#cdShort <- cdWB %>% 
#             filter(tag %in% c(  '257c683e2f', '00088d1c80', '00088d2e1f')) %>% 
#             select(tag,riverMeter,detectionDate,survey,sizeForGraph)

#save(cdShort,file='/home/ben/git/movementModel/shinyApp/cdForShinyShort.RData')


########## 3_envdata
library(dplyr)
library(getWBData)

# run this to reestablish connection
# createCoreData()

flowData <- tbl(conDplyr, "data_flow_extension") %>% collect(n = Inf) 
tempData <- tbl(conDplyr, "data_daily_temperature") %>% collect(n = Inf)

envData <- left_join(tempData, flowData, by=c("river","date")) %>% filter(river=="west brook")
#  filter( section == 11 ) %>%
  #      mutate(date = as.POSIXct(as.Date(datetime))) %>% 
#  mutate(date = format(datetime, "%Y-%m-%d")) %>% 
#  group_by(date) %>% 
#  mutate(tempMean = mean(temperature), depthMean = mean(depth,na.rm=T)) %>%
#  rename(temp15Min = temperature, daily_mean_temp = tempMean,
#         depth15Min = depth, qPredicted = depthMean) %>%
#  select(-datetime,-temp15Min,-depth15Min,-salinity,-logger) %>%
#  unique()

######
# temporary fix for negative depths
# 
envData$qPredicted <- ifelse(envData$qPredicted < 0, 0, envData$qPredicted)



# get median sample dates
# sampleName is the original, sampleNumber is consecutive
msd <- tbl(conDplyr,"data_seasonal_sampling") %>%
  filter(drainage=="west", seasonal) %>%
  select(sample_number,median_date,season,seasonal) %>%
  distinct() %>%
  collect() %>%
  arrange(median_date) %>%
 # mutate(date = format(median_date, "%Y-%m-%d")) %>%
  filter(sample_number !=2.5 & sample_number !=10.1)

envDataWB <- envData %>%
  left_join(.,msd, by = c("date" = "median_date"))

envDataWB$date <- as.POSIXct(envDataWB$date) # for row binding in the next step

envDataWB$seasonFill <- NA; envDataWB$sampleNumberFill <- NA
envDataWB$season[2] <- msd$season[1]; envDataWB$sample_number[2] <- msd$sample_number[1]
for (i in 2:(nrow(envDataWB))){
  if ( !is.na(envDataWB$season[i]) ) { envDataWB$seasonFill[i] = envDataWB$season[i]; envDataWB$sampleNumberFill[i] = envDataWB$sample_number[i] }
  else                               { envDataWB$seasonFill[i] = (envDataWB$seasonFill[i-1]); envDataWB$sampleNumberFill[i] = (envDataWB$sampleNumberFill[i-1]) }
}

envDataWB$drainage <- 'west'

save(envDataWB,file='/home/ben/shinyApps/pitAntenna/WB/envDataWBForMM.RData')


