#'Add transect data (estimated section volume) and fill in missing years by interpolation to a data frame. Also saves 'transects.RData' which contains the following objects: transects, vol, meansBySection, meansByYear
#'@return A data.frame appended with section volume properties
#'@param data A data.frame created with createCoreData(), which must include river and section as a column 
#'@param basecolumns Logical: include default columns? river,sampleName,width, maximum_depth
#'@param columnsToAdd A character vector of columns to inlcude; can replace or add to baseColumns
#'@export


# temporary
load('cdWBForTransects.RData')


addTransectData <- function(data){
  
  require(tidyverse)
  require(lubridate)
  require(getWBData)
  require(zoo)
  
  #rm( list=setdiff(ls(), c("conDplyr","con")) )
  reconnect()
  
  wb <- read.csv('wbCSV.csv', header = T, stringsAsFactors = F); wb$section <- as.numeric(wb$section)
  ob <- read.csv('obCSV.csv', header = T, stringsAsFactors = F); #ob$section <- as.numeric(as.character(ob$section))
  m <- read.csv('mCSV.csv', header = T, stringsAsFactors = F); #m$section <- as.numeric(as.character(m$section))
  j <- read.csv('jCSV.csv', header = T, stringsAsFactors = F); j$quarter <- as.numeric(j$quarter) # there is a quarter "5 pool"
  
  transects <- bind_rows(wb,ob,m,j)
  transects$date <- as.Date(transects$date, "%m/%d/%Y")
  transects$year <- year(transects$date)
  
  # fix based on volumeDis data graphs
  transects <- transects %>% filter( !(river == 'wb jimmy' & section == 5 & year == 2007 ))
  transects$section <- ifelse(transects$river == 'wb jimmy' & transects$section < 5 & transects$year == 2007,
                              transects$section + 1, transects$section)
  
  transects <- filter(transects,!(river == 'west brook' & date == '2009-10-19')) # this section was done at the end of a trib day and then repeated the next day. delete
  
  transects$depth[transects$river == 'wb mitchell' & transects$depth > 50] <- NA
  
  
  #r calc volumes}
  
  #transects$wD <- transects$width * transects$depth
  
  vol1 <- transects %>%
    group_by(river,year,date,section,quarter) %>%
    summarize( sumDepths = sum(depth), n = n(), maxDepth = max(depth), width = max(width) )
  
  vol <- vol1 %>%
    group_by(river,year,date,section) %>%
    mutate( volumeQ = sumDepths * width ) %>%
    summarize( volume = sum(volumeQ), nVol = n(), maxDepth = max(maxDepth), width = max(width) ) 
  
  
  #Link in discharge for transect dates and scale either by 1) discharge or 2) z-score of discharge within a river/year
  #r scale to flow}
  
  flow <- data.frame(tbl(conDplyr, "data_daily_discharge")) %>%
    #filter(date %in% samples$date) %>%
    select(date, river, discharge)
  
  #scale by discharge
  vol <- left_join(vol,flow) %>% mutate( volumeDis = volume/(discharge + 1) )
  
  #####################################################################
  # scale by z-score
  meansByYear <- vol %>%
    group_by(river,year) %>%
    summarize( meanV = mean(volume), stdV = sd(volume) ) %>%
    ungroup()
  
  vol <- vol %>%
    left_join(.,meansByYear) %>%
    mutate( volumeZ = (volume - meanV)/stdV + 2 ) # + 2 so all values are > 0
  
  #r, warning=FALSE, width = 9, height = 6}
  
  
  #Get means by section
  
  meansBySection <- vol %>%
    group_by(river,section) %>%
    summarize( meanVolumeDis = mean(volumeDis),
               meanVolumeZ = mean(volumeZ),
               meanVolumeBySection = mean(volume),
               meanMaxDepthBySection = mean(maxDepth),
               minYear = min(year))
  
  #Linearly interpolate missing years
  
  vol <- vol %>% 
    group_by(river,section) %>%
    complete(nesting(river,section), year = 1997:2015 ) %>%
    mutate( volumeDis = na.approx(volumeDis, rule = 2),
            volumeZ =   na.approx(volumeZ, rule = 2),
            volumeNAApprox = na.approx(volume, rule = 2),
            maxDepthNAApprox = na.approx(maxDepth, rule = 2))  
  
  vol <- vol %>%
    left_join(.,meansBySection) 
  
  vol$volumeDis <- ifelse( vol$year < vol$minYear, vol$meanVolumeDis, vol$volumeDis )
  vol$volumeZ   <- ifelse( vol$year < vol$minYear, vol$meanVolumeZ, vol$volumeZ )
  vol$volumeNAApprox   <- ifelse( vol$year < vol$minYear, vol$meanVolumeBySection, vol$volumeNAApprox )
  vol$maxDepthNAApprox <- ifelse( vol$year < vol$minYear, vol$meanMaxDepthBySection, vol$maxDepthNAApprox )
  
  # means by river
  volMeansByRiver <- vol %>%
    group_by(river) %>%
    summarize( meanMaxDepthByRiver = mean(maxDepth, na.rm = T),
               meanVolByRiver = mean(volume, na.rm = T) )
  
  vol <- left_join( vol, volMeansByRiver ) %>% mutate( pool01 = 1*(volumeNAApprox > meanVolByRiver & maxDepthNAApprox > meanMaxDepthByRiver )  )
  
  # means by isYOY,river,year,season,section,species
  # this will be our 'raw' mean biomass and counts
  
  data <- left_join(data,vol, by = c("river", "section", "year"))
  
  byYRYSS <- data %>% 
    group_by(isYOY,river,year,season,section,species) %>%
    summarize( meanBiomassYRYSS = mean(observedWeight,na.rm = T),
               meanBiomassByVolZYRYSS = mean(observedWeight/(volumeZ), na.rm = T), 
               varBiomassByVolZYRYSS = var(observedWeight/(volumeZ), na.rm = T), 
               nYRYSS = n(),
               meanVolZ = mean(volumeZ, na.rm = T),
               sumBiomassYRYSS = sum(observedWeight,na.rm = T)) %>%
    select(meanBiomassYRYSS,nYRYSS,meanBiomassByVolZYRYSS,varBiomassByVolZYRYSS,meanVolZ,sumBiomassYRYSS,isYOY,river,season,year,section,species)
  
  data <- left_join(data,byYRYSS)
  
  save(transects, vol, meansBySection, meansByYear, file = 'transects.RData')
  
  return(data)

}

