################ 
### 
###   Create data training and validatinfg for clasification modeling
### 
###                   
###                   
###   Author:  Andr?s Aguilar
###     CIAT - DAPA - AEPS
#########


rm(list=ls())
######## Functions
setwd("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/")

source(paste0("_scripts/R_Main_Functions_GrowCropIndentification.R"))


source("_scripts/process_main_function.R")


##### Load Libraries

libs=c("snowfall","stringr","raster", "rgeos","rgdal","dplyr", "ggplot2")
PackageReading(libs)
##### Set Locality Name

tile = "col_t3"
### Set Criteriums

## Period
DateStart="2015-07-07"
DateEnd="2016-01-23"

## Read optical data
inventory = get_inventory( tile)

vi_images = get_vi_layers(inventory = inventory , veg_indexes = c("NDVI","LSWI"),
                          DateStart, DateEnd,PercMaxBadPixels = 51)

plot(vi_images$LE07_008057_20150707_col_t3_10m.tif)


### Read Spatial data 


training_points = read.csv("model_inputs/phen_identification/data_coordinates.csv", row.names = 1) 
table(training_points$class)

### extract vegetation indexes info from images

names_optical = names(vi_images)
SentinelimageReference = vi_images[[1]][[1]]
SentinelimageReference[] = NA


training_points = SpatialPointsDataFrame(training_points[,1:2],data = training_points,
                                         proj4string = crs(SentinelimageReference))

training_points = Extract_raster_data(raster_images = vi_images , spatial_points = training_points , 
                                      spectral_bands = names(vi_images[[1]]))
names(training_points) = sapply(names_optical, function(name_imag) IdentiyImageDate(name_imag))


rm(vi_images)

### Join capturing dates 
i = 1
training_points = do.call(rbind , lapply(1:length(training_points) , function(i){
  dataper_date=training_points[[i]]
  dataper_date$capturing_date = names(training_points)[i]
  dataper_date$field_name = factor(as.character(dataper_date$field_name))
  return(dataper_date)
}))


tail(training_points)

### remove duplicated rows

data_WithOutDuplicated = data.frame(training_points %>%
                                      group_by(x,  y ,field_name,  date , class , capturing_date) %>% filter(row_number(y) == 1))

idrows =  with(data_WithOutDuplicated, paste0(x,  y,field_name,  date , class, capturing_date ,  "_"))

idrows[idrows %in%"517985428955none20151228other20150731_"]

## the date forma is changed
data_WithOutDuplicated$capturing_date = as.Date(data_WithOutDuplicated$capturing_date , format = "%Y%m%d")

#phase_stages = c("late_vegetative"  ,"early_vegetative","harvested","reproductive","ripening")
phase_stages = c("vegetative","harvested","reproductive","ripening")
# phasecolours = c ('harvested'='firebrick3', "ripening"="goldenrod", 
#    "early_vegetative"="cyan4", "reproductive" = "darkgreen","late_vegetative" = "chartreuse3 ,"soil" = "saddlebrown""
# )
phasecolours = c ('harvested'='firebrick3', "ripening"="goldenrod", 
                  "reproductive" = "darkgreen","vegetative" = "darkolivegreen3", "soil" = "saddlebrown"
)

## plot initial state
data_plot = data_WithOutDuplicated[ data_WithOutDuplicated$class %in%  
                                      phase_stages, ]



data_plot$ID = with(data_plot,paste(x , y , field_name,  date , class , sep = "_"))

data_date = do.call(rbind,lapply(c("20151127","20151211","20151228","20160114","20160123"), function(date_campaign){
  data_plot[data_plot$date %in% date_campaign & 
              data_plot$capturing_date<=(CheckDateFormat(date_campaign)),]
}))


m=ggplot(data_date, aes(capturing_date ,NDVI , colour=field_name, group = ID))+
  geom_point()+geom_line(alpha = 0.8)+ facet_grid(rows = vars(class), cols = vars(date))

ggsave(plot = m , filename = paste0("temp/testuncleanNDVIveg.jpg"),
       width = 34, height = 25, units = "cm")
rice_field = "21C919"

phen_ident_data=read.csv("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/field_data/training_fields/growthphases_campaigns_v2.csv")

campaign_dates = c("20151127","20151211","20151228","20160114","20160123")

data_plot2 = data_WithOutDuplicated

rice_field = "21C919"

lapply(levels(data_plot2$field_name), function(rice_field){
  cat(rice_field,"\n")
  subset_field = data_plot2[data_plot2$field_name%in%rice_field, ]
  subset_field$temp = with(subset_field,paste0(x,y,field_name))
  subset_fields = unique(subset_field$temp)[sample(1:length(unique(subset_field$temp)),round(length(unique(subset_field$temp)))*0.5)]
  if(!rice_field%in%c("other","Water","Urban_zones")){
    dates_camp = phen_ident_data[phen_ident_data$Field %in%rice_field,]
    dates_camp = dates_camp[order(as.Date(as.character(dates_camp$Date), format = "%Y%m%d")),]
    colours_vlines = sapply(as.character(dates_camp$state), function(x)
      switch(x,'harvested'='firebrick3', "ripening"="goldenrod", 
             "early_vegetative"="cyan4", "reproductive" = "darkgreen",
             "late_vegetative" = "chartreuse3", "soil" = "saddlebrown"))
    
    vlinekwars = data.frame(cdates = campaign_dates[campaign_dates%in%as.character(dates_camp$Date)],
                            ccols = unname(colours_vlines))
    
  } else{
    vlinekwars =  data.frame(cdates = "20151228",
                             ccols = "black")
  }
  m=ggplot(subset_field[subset_field$temp%in%subset_fields,], aes(capturing_date ,NDVI , colour=field_name, group = temp))+ylim(c(0,1))+
    geom_line(alpha = 0.8, colour = "gray60")+geom_point(size=1.2, colour= "gray30")+theme_bw()+
    geom_vline(colour = "black",xintercept =  as.Date(campaign_dates,format = "%Y%m%d"))+
    geom_vline(colour = vlinekwars$ccols,xintercept = as.Date(vlinekwars$cdates,format = "%Y%m%d"))
  
  ggsave(plot = m , filename = paste0("temp/testuncleanNDVIv2",rice_field,".jpg"),
         width = 18, height = 10, units = "cm")
})


### second graph : Boxplot
data_date = do.call(rbind,lapply(c("20151127","20151211","20151228","20160114","20160123"),
                                 function(date_campaign){
                                   
                                   data_plot = data_plot[data_plot$date %in% date_campaign & 
              data_plot$capturing_date<=CheckDateFormat(date_campaign)&
              data_plot$capturing_date>=(CheckDateFormat(date_campaign)-130),]
                                   
                                   
}))

lines_plot = data.frame(group_by(data_date, class,date,capturing_date) %>%summarise(ndvi = median(LSWI, na.rm = T)))
data_date$temp = paste0(data_date$capturing_date, data_date$class)

data_date$date = as.factor(data_date$date)
levels(data_date$date)= c("27-11-2015", "11-12-2015", 
                                                   "28-12-2015","14-01-2015","23-01-2016")


data_date$class = as.factor(as.character(data_date$class))
data_date$class = factor(as.character(data_date$class), levels = c("vegetative","reproductive","ripening","harvested"))
m = ggplot(data_date, aes(as.factor(capturing_date), NDVI ,group = temp)) +
  
  theme(text=element_text(size=12),
        axis.title.y=element_text(size = rel(1.1),colour = "black", face="bold"),
        axis.title.x=element_text(size = rel(1.1),colour = "black", face="bold"),
        axis.text.x  = element_text(angle=30, hjust=1),
        panel.background=element_rect(fill="white",colour = "black"),
        panel.grid.major = element_line(colour = "gray"),
        legend.title = element_text(face = "bold"))+ 
  labs( fill = "Growth Phase", x = "Satellite acquisition date")+ 
  geom_boxplot(aes(fill = class), alpha = .60)+lims(y = c(0,1))+ 
  stat_summary(fun.y=median, geom="line", aes(colour = class,group = class), size = 1.4)+ 
  stat_summary(fun.y=median, geom="point", aes(colour = class,group = class), size = 1.8)+
  facet_grid(date ~.)


m = m+ guides(colour=FALSE)+scale_colour_manual(values = phasecolours)+
  scale_fill_manual(values = phasecolours)

ggsave(plot = m , filename = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/figures/ndvi_time_seriesunclean_boxplot_veg.jpg"),
       width = 30, height = 20, units = "cm")



###### lines 
CheckDateFormat("20151127")
data_date = do.call(rbind,lapply(c("20151127","20151211","20151228","20160114"), function(date_campaign){
  training_points_veg[training_points_veg$Date %in% date_campaign & 
                        training_points_veg$Image_date<CheckDateFormat(date_campaign),]
}))


lines_plot = data.frame(group_by(data_date, Typ_Stg,Date,Image_date) %>%summarise(ndvi = median(LSWI, na.rm = T)))

lines_plot$temp = paste0(lines_plot$Image_date, lines_plot$Typ_Stg)


m = ggplot(lines_plot, aes(as.factor(Image_date), ndvi ,group = Typ_Stg)) + 
  geom_boxplot( alpha = .60)+
  geom_line(aes(color = Typ_Stg), alpha = .9)+
  geom_point(aes(color = Typ_Stg), alpha = .9)+
  lims(y = c(0,1))+ theme_bw() + labs( color = "Etapa de Crecimiento",
                                       y = "promedio de NDVI",
                                       x = "Fecha de captura de las imágenes satelitales")+ 
  facet_grid(Date ~.)

m+scale_colour_manual(values = c ('harvested'='firebrick3', "ripening"="goldenrod", 
                                                        "early_vegetative"="cyan4", "reproductive" = "darkgreen","late_vegetative" = "chartreuse3"
))+
  scale_fill_manual(values = c ('harvested'='firebrick3', "ripening"="goldenrod", 
                                "early_vegetative"="cyan4", "reproductive" = "darkgreen","late_vegetative" = "chartreuse3"
  ))



m=ggplot(training_points_veg[training_points_veg$iteration%in%1,] , aes(Image_date ,NDVI , colour=plygn_n, group = ID))+
  geom_point()+geom_line(alpha = 0.8)+ facet_grid(rows = vars(Typ_Stg), cols = vars(Date))



ggsave(plot = m , filename = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/temp/testcleanNDVI.jpg"),
       width = 34, height = 25, units = "cm")


head(data_WithOutDuplicated)
############################### end cleaning dataset
#------------------------------------------------------------------------

#-- Start smooth process

######### ---> Kernel Smooth process

## Date of Interest
ListIndex = c("NDVI")
bool_cross = F
## Set date period
datestoEvaluate = as.character(unique(data_WithOutDuplicated$date))
DateInt="20151228"
dataSmoothed = lapply(datestoEvaluate, function(DateInt){
  cat(DateInt,"\n")
  rangeDate=c(as.Date(DateInt,format="%Y%m%d")-96,(as.Date(DateInt,format="%Y%m%d")+0))
  nDaysPoints=length(rangeDate[1]:rangeDate[2])/16

  VI = "NDVI"
  data_reshape = lapply(ListIndex , function(VI){
    ## Select data for the interested date
    
    data_WithOutDuplicated_date = data_WithOutDuplicated[data_WithOutDuplicated$date %in% DateInt,]
    
    ## filter by interest date
    data_WithOutDuplicated_date=
      data_WithOutDuplicated_date[data_WithOutDuplicated_date$capturing_date %in% 
                                    (as.Date(DateInt,format=" %Y%m%d")-130):as.Date(DateInt,format=" %Y%m%d"),]
    
    if(dim(data_WithOutDuplicated_date)[1]>0){
      data_WithOutDuplicated_date$capturing_date = paste0(VI,"_", as.character(data_WithOutDuplicated_date$capturing_date, format = "%Y%m%d"))
      ## reshape dataset
      data_WithOutDuplicated_date$ID =
        with(data_WithOutDuplicated_date,paste0(as.character(x), as.character(y),"_p_" ,
                                                as.character(field_name),"_s_",as.character(class),"_d" ,as.character(date)))

      data_reshape = reshape2::dcast( data_WithOutDuplicated_date[,c(VI,"capturing_date","ID")] , ID ~ capturing_date, value.var = VI ,fun.aggregate=mean)
      ## assign ID code
      return(data_reshape)
      
    }else{
      return(NULL)
    }
  })
  
  if(!is.null(data_reshape[[1]])){
    ###### join
    
    data_reshape = plyr::join_all(data_reshape, by = "ID")
    
    ### remove those pixels with high percentage of NA in the time series
    
    pixelsPercentageNA=data.frame(Pixel=data_reshape$ID,
                                  Percentage=apply(data_reshape[,-which(names(data_reshape)%in% c("ID"))],
                                                   1,function(x){(sum(is.na(x))/length(x)*100)}))
    
    pixelsHighPercentage=which(pixelsPercentageNA$Percentage>55)
    if(length(pixelsHighPercentage)>0){
      data_reshape=data_reshape[-pixelsHighPercentage,]
      
    }else{
      data_reshape=data_reshape
    }
    
    
    ###############
    
    #### Smooth process
    
    rm(pixelsPercentageNA)
    row.names(data_reshape) = data_reshape$ID
    

    KernelSmootValues=SmoothAllPixels(DateInt = DateInt,ImagesInfo = data_reshape,
                                      Veg_Indexes = ListIndex,
                                      rangeDate = rangeDate,
                                      nDaysPoints = nDaysPoints,
                                      parallelProcess =T,ncores = 5)

    # rowsp = grep(row.names(ImagesInfo), pattern = "9855425265_p_21C919_s")
    # j = 816

    # lapply(rowsp,function(j){
    #   ImagesInfo_PerPixel=ImagesInfo[j,]
    #   cat(paste("interation ", j))
    #   data_PerPixel_VI  = ImagesInfo_PerPixel[grepl(names(ImagesInfo_PerPixel), pattern =VI)]
    #   timeseriesmetrics(data_PerPixel_VI,DateInt,rangeDate,nDaysPoints)
    # })
    # # 
    
    
    #-- End smooth process
    ########
    #----------------------------------
    
    ## bandwidth cross correlation 
    getcrosscorr_scrore = do.call(rbind,lapply(1:length(KernelSmootValues), function(x)
      {
      crossvalues = KernelSmootValues[[x]][[1]]$cross_corr_kernel
      lengthbands = length(seq(20,70,1))
      if(length(crossvalues)!= lengthbands){
        crossvalues = rep(NA, lengthbands)
      }
      data.frame(pixel=rep(x,lengthbands),
                 bandwidth = seq(20,70,1),
                 rsquare = crossvalues)
    }
      ))
    # ggplot(getcrosscorr_scrore, aes(rsquare))+geom_density()
    # ggplot(getcrosscorr_scrore, aes(factor(bandwidth), rsquare))+geom_boxplot()
    summary_cross=data.frame(getcrosscorr_scrore%>% group_by(bandwidth)%>%
                               summarise(meansquare = mean(rsquare, na.rm =T),
                                        sdvar = sd(rsquare, na.rm =T),
                                                           nacount = length(which(is.na(rsquare)))))
    print(summary_cross)
    ########
    #-- Start export step
    
    ####get smooth values
    
    TableValuesToClassify=do.call(cbind,lapply(1:length(ListIndex),function(index){
      ValTableperIndex=data.frame(do.call(rbind,lapply(lapply(KernelSmootValues,function(x){x[[index]]}),
                                                       function(x){x$regression_values})))
      names(ValTableperIndex)=paste0(ListIndex[index],"_Date_",1:length(names(ValTableperIndex)))
      return(ValTableperIndex)
    }))
    
    row.names(TableValuesToClassify)=row.names(data_reshape)
    ## get derivative values
    
    derivativeValues=do.call(cbind,lapply(1:1,function(index){
      ValTableperIndex=data.frame(do.call(rbind,lapply(lapply(KernelSmootValues,function(x){x[[index]]}),
                                                       function(x){x$derivative_values})))
      names(ValTableperIndex)=paste0(ListIndex[index],"_derivative_",1:length(names(ValTableperIndex)))
      return(ValTableperIndex)
    }))
    row.names(derivativeValues)=row.names(data_reshape)
    
    ## get summaryvariables
    
    summaryValues=do.call(cbind,lapply(1:1,function(index){
      ValTableperIndex=data.frame(do.call(rbind,lapply(lapply(KernelSmootValues,function(x){x[[index]]}),
                                                       function(x){x$summary_ts})))
      names(ValTableperIndex)=paste0(ListIndex[index],"_",names(ValTableperIndex))
      return(ValTableperIndex)
    }))
    row.names(summaryValues)=row.names(data_reshape)
    
    TableValuesToClassify = cbind(TableValuesToClassify,derivativeValues,summaryValues)
    
    ####### Assign class label to data frame
    
    num_pos = regexpr(pattern = "_s_",row.names(TableValuesToClassify))
    num_pos_end = regexpr(pattern = "_d2",row.names(TableValuesToClassify))
    TableValuesToClassify$class =str_sub(row.names(TableValuesToClassify),num_pos+3, num_pos_end-1)
    
    ### export training data
    ## remove NA
    TableValuesToClassify = TableValuesToClassify[!apply(TableValuesToClassify , 1 , function(column)T%in%is.na(column)),]
    
    ### smooth graphic

    sg_values=data.frame(do.call(rbind,lapply(lapply(KernelSmootValues,function(x){x[[1]]}),
                                                       function(x){x$sg_smooth})))
    
    sg_values$id = paste0(sg_values$id, "_da",sg_values$x) 
  
    data_orig = reshape2::melt(data_reshape)

    names(data_orig)[3] = "y_real" 
    
    data_orig$x = as.numeric(as.Date(IdentiyImageDate(data_orig$variable),"%Y%m%d"))
    data_orig$id = paste0(data_orig$ID, "_da",data_orig$x) 
    
    # data_plot = merge(data_orig, sg_values, by = 'id')
    # data_plot = data_plot[!is.na(data_plot$y_real),]
    # set.seed(123)
    # subsampleid = data_plot$ID[sample(1:length(unique(data_plot$ID)), 50)]
    # data_plot= data_plot[data_plot$ID%in%subsampleid,]
    # 
    # num_pos = regexpr(pattern = "_p_",data_plot$ID)
    # num_pos_end = regexpr(pattern = "_s",data_plot$ID)
    # 
    # data_plot$class =str_sub(data_plot$ID,num_pos+3, num_pos_end-1)
    # 
    # m =ggplot(data_plot, aes(x.x, y_real, group= ID))+ geom_point(colour = "green")
    # m+geom_point(data = data_plot, aes(x.x, y, group= ID))
    return(list(TableValuesToClassify,summary_cross))
  }else{
    return(NULL)
  }
  
})

### cross correlation band width
getcrosscorr_scrore = do.call(rbind,lapply(1:length(dataSmoothed), function(x)
  {
    dataSmoothed[[x]][[2]]$date = x
    return(dataSmoothed[[x]][[2]])}
))

ggplot(getcrosscorr_scrore, aes(sdvar, meansquare)) + geom_point(aes(color =bandwidth))  +
  scale_colour_gradientn(colours = terrain.colors(10))


ggplot(getcrosscorr_scrore, aes(meansquare))+geom_density()
ggplot(getcrosscorr_scrore, aes(factor(date),meansquare))+geom_boxplot()
summary_cross=data.frame(getcrosscorr_scrore%>% group_by(bandwidth)%>%
                           summarise(meansquare = mean(meansquare, na.rm =T),
                                     sdvarmean = mean(sdvar, na.rm =T),
                                     nacount =sum(nacount)))

summary_cross[order(summary_cross$meansquare,decreasing = T),]

####

TableValuesToClassify = do.call(rbind,lapply(1:length(dataSmoothed), function(x)
  dataSmoothed[[x]][[1]]
))

graphic_smoothed = function(test, vi = "NDVI",  limitbottom = 0.1, limitup = 1){
  test$ID = row.names(test)
  test = (reshape2::melt(test))
  num_pos = regexpr(pattern = "_p_",test$ID)
  num_pos_end = regexpr(pattern = "_s_",test$ID)
  test$plygn_n = str_sub(test$ID,num_pos+3, num_pos_end-1)
  
  num_pos = regexpr(pattern = "_s_",test$ID)
  num_pos_end = regexpr(pattern = "_d2",test$ID)
  test$Typ_Stg= str_sub(test$ID,num_pos+3, num_pos_end-1)
  
  num_pos = regexpr(pattern = "_d2",test$ID)

  test$Date= str_sub(test$ID,num_pos+2, num_pos+10)
  
   m=ggplot(test , aes(variable ,value , colour=plygn_n, group = ID))+ylim(limitbottom,limitup)+
    geom_point()+geom_line(alpha = 0.8)+ facet_grid(rows = vars(Typ_Stg), cols = vars(Date))
  m
}


test = TableValuesToClassify[,grepl(names(TableValuesToClassify), pattern = "vi_day")]
m = graphic_smoothed (test, vi = "NDVI",  limitbottom = 50, limitup = 0)

ggsave(plot = m , filename = paste0("temp/Data_v8_1_NDVI_days_veg2_.png"),
       width = 34, height = 25, units = "cm")

m = graphic_smoothed (test, vi = "NDVI")

test = TableValuesToClassify[,grepl(names(TableValuesToClassify), pattern = "NDVI_Da")]

test = test[grepl(row.names(test),pattern = "21D823_s_vegetative_d2015121"),]
#test = test[grepl(row.names(test),pattern = "vegetative"),]
graphic_smoothed (test, vi = "NDVI")
#### add data to other

TableValuesToClassify = TableValuesToClassify[!(TableValuesToClassify$class == "harvested" & grepl(row.names(TableValuesToClassify),pattern = "20151230")),]

TableValuesToClassify[(TableValuesToClassify$class == "ripening" & grepl(row.names(TableValuesToClassify),pattern = "20151230")),"class"] = "reproductive"

TableValuesToClassify[(TableValuesToClassify$class == "reproductive" & grepl(row.names(TableValuesToClassify),pattern = "p_11B011A_s_reproductive_d20160123")),"class"] = "ripening"


bad_positions_late =  row.names(TableValuesToClassify[grepl(row.names(TableValuesToClassify),pattern = "11B011A_s_vegetative_d20151"),])
TableValuesToClassify[row.names(TableValuesToClassify)%in%bad_positions_late,"class"] = "reproductive"

bad_positions_rep = row.names(TableValuesToClassify[TableValuesToClassify$class %in%  "reproductive" & 
                                                      (TableValuesToClassify$NDVI_Date_7<0.68),])

other_data = c(bad_positions_late,bad_positions_rep[grepl(bad_positions_rep, pattern = "1C914_")])

TableValuesToClassify[row.names(TableValuesToClassify)%in%other_data,"class"] = "other"


## remove bad signals

# bad_positions_late =  row.names(TableValuesToClassify[TableValuesToClassify$class %in%  "late_vegetative" & 
#                                                    (TableValuesToClassify$NDVI_Date_7<=0.50),])
# 
# bad_positions_early = row.names(TableValuesToClassify[TableValuesToClassify$class %in%  "early_vegetative" &
#                                     ((TableValuesToClassify$NDVI_Date_7<0.38) |
#                                       (TableValuesToClassify$NDVI_Date_4>0.5)|
#                                        (TableValuesToClassify$NDVI_Date_7>0.7)),])

bad_positions_veg = row.names(TableValuesToClassify[TableValuesToClassify$class %in%  
                                                        "vegetative" &((TableValuesToClassify$NDVI_Date_7<0.37) |
                                                                                              
                                                                                               (TableValuesToClassify$NDVI_Date_7>0.8)),])
bad_positions_ripe = row.names(TableValuesToClassify[TableValuesToClassify$class %in%  "ripening" & 
                                                       ((TableValuesToClassify$NDVI_Date_4<0.55) |
                                                          (TableValuesToClassify$NDVI_Date_5<0.6))
                                                    ,])


bad_positions_rep = row.names(TableValuesToClassify[TableValuesToClassify$class %in%  "reproductive" & 
                                    (TableValuesToClassify$NDVI_Date_5<0.40),])


bad_positions_harves = row.names(TableValuesToClassify[TableValuesToClassify$class %in%  "harvested" &
                                    (TableValuesToClassify$NDVI_Date_7>0.625),])

bad_positions_soil = row.names(TableValuesToClassify[TableValuesToClassify$class %in%  "soil" &
                                                         (TableValuesToClassify$NDVI_Date_7>0.40),])



wrongdata = c(bad_positions_veg,bad_positions_ripe,
              bad_positions_rep,bad_positions_harves,bad_positions_soil)


training_data = TableValuesToClassify[!row.names(TableValuesToClassify)%in%wrongdata,]

test = training_data[,grepl(names(training_data), pattern = "NDVI_Da")]

graphic_s = graphic_smoothed (test, vi =  "NDVI", limitbottom = 0.025,limitup = 0.9)

ggsave(plot = graphic_s , filename = paste0("temp/cleanData_v8_1_NDVI_DacaVeg2.png"),
       width = 34, height = 25, units = "cm")

#### remove by derivatives
test = training_data[,grepl(names(training_data), pattern = "NDVI_de")]

ggplot(training_data, aes(NDVI_derivative_1,NDVI_derivative_2, color = class)) + 
  geom_point()#+xlim(-0.03,0.03)#+ylim(-0.03,0.02)


bad_positions_early =  row.names(TableValuesToClassify[TableValuesToClassify$class %in%  "vegetative" & 
                                                        (TableValuesToClassify$NDVI_derivative_2 <= (-0.003)),])

bad_positions_ripe = row.names(TableValuesToClassify[TableValuesToClassify$class %in%  "ripening" & 
                                                       ((TableValuesToClassify$NDVI_derivative_2> (0.001)))
                                                     ,])


wrongdata = c(bad_positions_early,bad_positions_ripe)

training_data = training_data[!row.names(training_data)%in%wrongdata,]


training_data$class
#### Export data
table(training_data$class)


write.csv(training_data ,
          paste0("model_inputs/phen_identification/optical_data_ndvi_date_deri_veg.csv"))

data_training = read.csv(paste0("model_inputs/phen_identification/optical_data_ndvi_date_deri_veg.csv"),row.names =1)

data_training$class = as.character(data_training$class)
data_training$class [data_training$class %in% c("Urban_Zones", "Water")] = "other"
unique(data_training$class)

###

ggplot(data_training, aes(NDVI_derivative_1,NDVI_derivative_2, color = class)) + 
  geom_point()+xlim(-0.03,0.04)+ylim(-0.045,0.04)+
  scale_colour_manual(values = c(phasecolours, other = "gray19")
  ) + theme_bw() + 
  labs (x = "initial_NDVI_derivative", y = "ending_NDVI_derivative", colour = "Classification")+
  labs(y = "first_derivative_tsending", x = "first_derivative_tsstarting ", colour = "Growth Phase")+
    theme(text=element_text(size=12),
          axis.title.y=element_text(size = rel(1.1),colour = "black", face="bold"),
          axis.title.x=element_text(size = rel(1.1),colour = "black", face="bold"),
          axis.text.x  = element_text(angle=0, hjust=0.5),
          panel.background=element_rect(fill="white",colour = "black"),
          panel.grid.major = element_line(colour = "gray"),
          legend.title = element_text(face = "bold"))

###


ggplot(data_training[data_training$class %in% c("vegetative"),], aes(NDVI_maximumvi_day,NDVI_minimumvi_day, color = class)) + 
  geom_point()+xlim(-0.10,60)+ylim(-0.10,60)+
  scale_colour_manual(values = c ('harvested'='firebrick3', "ripening"="goldenrod", "other" = "gray", 
                                  "vegetative"="lightgreen", "reproductive" = "darkgreen", "soil" = "saddlebrown"
  )) + theme_bw() + labs (x = "maximum NDVI day", y = "minimum NDVI day", colour = "Classification")

###

ggplot(data_training[!data_training$class %in% c("other"),], aes(NDVI_min_ts,NDVI_max_ts, color = class)) + 
  geom_point()+#xlim(-0.10,50)+ylim(-0.10,50)+
  scale_colour_manual(values = c ('harvested'='firebrick3', "ripening"="goldenrod", "other" = "gray", 
                                  "vegetative"="lightgreen", "reproductive" = "darkgreen", "soil" = "saddlebrown"
  )) + theme_bw() + labs (x = "min NDVI", y = "max NDVI", colour = "Classification")

ggplot(data_training[!data_training$class %in% c("other"),], aes(NDVI_sd_ts,NDVI_minimumvi_day, color = class)) + 
  geom_point()+#xlim(-0.10,50)+ylim(-0.10,50)+
  scale_colour_manual(values = c ('harvested'='firebrick3', "ripening"="goldenrod", "other" = "gray", 
                                  "vegetative"="lightgreen", "reproductive" = "darkgreen", "soil" = "saddlebrown"
  )) + theme_bw() + labs (x = "SD NDVI", y = "Average NDVI", colour = "Classification")



### soil and other class
data_other = data_training[data_training$class %in%c("soil","other","vegetative"),]

mean_val= apply(data_other[,grepl(names(data_other), pattern = "Date")],1,mean)

sd_val= apply(data_other[,grepl(names(data_other), pattern = "Date")],1,sd)
dataplot = data.frame(mean_val, sd_val, class = data_other$class)
ggplot(dataplot, aes(mean_val, sd_val, colour= class)) + geom_point() 

othertochange = row.names(dataplot)[(dataplot$class %in% c("other"))&
                                      (dataplot$sd_val< 0.15 & dataplot$sd_val > 0.02)&
                                      (dataplot$mean_val> 0.02 & dataplot$mean_val < 0.5)]

data_other = data_other[row.names(data_other)%in%othertochange,
                        grepl(names(data_other), pattern = "Date")]

othertochange = row.names(data_other)[data_other$NDVI_Date_7>0.35&
                                        data_other$NDVI_Date_3<0.4&
                                        data_other$NDVI_Date_3>0.2&
                                        data_other$NDVI_Date_1<0.35&
                                        data_other$NDVI_Date_1>0.2&
                                        data_other$NDVI_Date_7<0.8&
                                        data_other$NDVI_Date_4<0.5&
                                        data_other$NDVI_Date_4>0.2&
                                        data_other$NDVI_Date_5<0.5&
                                        data_other$NDVI_Date_5>0.2]

data_other$id = row.names(data_other)

dataplot = reshape2::melt(data_other, by = "id")

num_pos = regexpr(pattern = "_s_",dataplot$id)
num_pos_end = regexpr(pattern = "_d2",dataplot$id)
dataplot$class= str_sub(dataplot$id,num_pos+3, num_pos_end-1)

ggplot(dataplot[dataplot$id %in%othertochange,], aes(group = id, variable, value, colour = class)) +geom_line()


data_training = data_training[!row.names(data_training)%in%othertochange,]


data_other = data_training[data_training$class %in%c("soil","other","vegetative"),]

mean_val= apply(data_other[,grepl(names(data_other), pattern = "Date")],1,mean)

sd_val= apply(data_other[,grepl(names(data_other), pattern = "Date")],1,sd)
dataplot = data.frame(mean_val, sd_val, class = data_other$class)
ggplot(dataplot, aes(mean_val, sd_val, colour= class)) + geom_point() 


othertochange = row.names(dataplot)[(dataplot$class %in% c("other"))&
                                      (dataplot$sd_val< 0.15 & dataplot$sd_val > 0.02)&
                                      (dataplot$mean_val> 0.02 & dataplot$mean_val < 0.5)]

data_other = data_other[row.names(data_other)%in%othertochange,
                        grepl(names(data_other), pattern = "Date")]

othertochange = row.names(data_other)[data_other$NDVI_Date_7<0.5&data_other$NDVI_Date_3<0.4&
                                        data_other$NDVI_Date_7>0.2&data_other$NDVI_Date_2>0.25&
                                        data_other$NDVI_Date_5<0.4&
                                        data_other$NDVI_Date_1<0.4]

data_other$id = row.names(data_other)

dataplot = reshape2::melt(data_other, by = "id")

num_pos = regexpr(pattern = "_s_",dataplot$id)
num_pos_end = regexpr(pattern = "_d2",dataplot$id)
dataplot$class= str_sub(dataplot$id,num_pos+3, num_pos_end-1)

ggplot(dataplot[dataplot$id %in%othertochange,], aes(group = id, variable, value, colour = class)) +
  geom_line()+ylim(0.1,0.8)

dataplot[dataplot$id %in%othertochange & (dataplot$value>0.39),]

data_training = data_training[!((row.names(data_training)%in%othertochange) & (data_training$NDVI_Date_7>0.39)),] 
data_training$class[row.names(data_training)%in%othertochange] = "soil"

data_training$class = factor(data_training$class, levels = c("vegetative","reproductive",
                                                             "ripening","harvested","soil","other"))


#############



### reproductive soil
data_other = data_training[data_training$class %in%c("reproductive","other","ripening"),]

mean_val= apply(data_other[,grepl(names(data_other), pattern = "Date")],1,mean)

sd_val= apply(data_other[,grepl(names(data_other), pattern = "Date")],1,sd)
dataplot = data.frame(mean_val, sd_val, class = data_other$class)
ggplot(dataplot, aes(mean_val, sd_val, colour= class)) + geom_point() 

othertochange = row.names(dataplot)[(dataplot$class %in% c("reproductive","other","ripening"))&
                                      (dataplot$sd_val< 0.15 & dataplot$sd_val > 0.06)&
                                      (dataplot$mean_val> 0.5 & dataplot$mean_val < 0.65)]

data_other = data_other[row.names(data_other)%in%othertochange,
                        grepl(names(data_other), pattern = "Date")]

data_other = data_other[data_other$NDVI_Date_7>0.55&
                          data_other$NDVI_Date_4>0.65&
                          data_other$NDVI_Date_3>0.4&
                          data_other$NDVI_Date_3<0.61&
                          data_other$NDVI_Date_5>0.6&
                          data_other$NDVI_Date_1<0.45 &
                          data_other$NDVI_Date_7<0.75,]

data_other$id = row.names(data_other)

dataplot = reshape2::melt(data_other, by = "id")

num_pos = regexpr(pattern = "_s_",dataplot$id)
num_pos_end = regexpr(pattern = "_d2",dataplot$id)
dataplot$class= stringr::str_sub(dataplot$id,num_pos+3, num_pos_end-1)

ggplot(dataplot, aes(group = id, variable, value, colour = class)) +geom_line()

data_training = data_training[!row.names(data_training)%in%dataplot$id[grepl(dataplot$id, pattern = "other")],]


#############33
data_other = data_training[data_training$class %in%c("reproductive","other","ripening"),]

mean_val= apply(data_other[,grepl(names(data_other), pattern = "Date")],1,mean)

sd_val= apply(data_other[,grepl(names(data_other), pattern = "Date")],1,sd)
dataplot = data.frame(mean_val, sd_val, class = data_other$class)
ggplot(dataplot, aes(mean_val, sd_val, colour= class)) + geom_point() 


othertochange = row.names(dataplot)[(dataplot$class %in% c("reproductive","other","ripening"))&
                                      (dataplot$sd_val< 0.16 & dataplot$sd_val > 0.06)&
                                      (dataplot$mean_val> 0.5 & dataplot$mean_val < 0.7)]


data_other = data_other[row.names(data_other)%in%othertochange,
                        grepl(names(data_other), pattern = "Date")]


data_other = data_other[data_other$NDVI_Date_7>0.64&
                          data_other$NDVI_Date_4>0.65&
                          data_other$NDVI_Date_3>0.6&
                          data_other$NDVI_Date_5>0.63&
                          data_other$NDVI_Date_1<0.55 &
                          data_other$NDVI_Date_1>0.4 &
                          data_other$NDVI_Date_6>0.63 &
                          data_other$NDVI_Date_2<0.65 &
                          data_other$NDVI_Date_7<0.75,]

data_other$id = row.names(data_other)

dataplot = reshape2::melt(data_other, by = "id")

num_pos = regexpr(pattern = "_s_",dataplot$id)
num_pos_end = regexpr(pattern = "_d2",dataplot$id)
dataplot$class= stringr::str_sub(dataplot$id,num_pos+3, num_pos_end-1)

ggplot(dataplot[dataplot$class %in% c("other" , "reproductive","ripening"),], aes(group = id, variable, value, colour = class)) +
  geom_line()+geom_point()

data_training = data_training[!row.names(data_training)%in%dataplot$id[grepl(dataplot$id, pattern = "reproductive")],] 


datostochange = unique(dataplot[dataplot$class %in% c("other"),"id"])
dataaux = data_training[row.names(data_training)%in%datostochange,]

valsaux = apply(dataaux[dataaux$class %in% c("other"),c("NDVI_Date_4","NDVI_Date_5","NDVI_Date_6","NDVI_Date_7")],1,function(x){
  mean(diff(x))
})

valsaux = names(valsaux[valsaux<0.02 & valsaux>-0.01])
ggplot(dataplot[dataplot$id %in% valsaux,], aes(group = id, variable, value, colour = class)) +
  geom_line()+geom_point()

### delete and change
todelete = dataplot$id[grepl(dataplot$id, pattern = "other") & (!dataplot$id%in%valsaux)]
data_training = data_training[!row.names(data_training)%in%todelete,]

data_training[row.names(data_training)%in%valsaux,"class"] = "ripening"



table(data_training$class)
#############

get_featuresmin_maxvals = function(dataf){
  ## number of features
  num_features = ncol(dataf)
  
  #calculate minimun and maximun values
  minmax_list =lapply(1:num_features, function(x){
    
    return(c(min(dataf[,x]),max(dataf[,x])))
  })
  
  ### rename list
  names(minmax_list) = names(dataf)
  return(minmax_list)
}
# 
feature_scaling =  function(feature, limits = c(0,1), type = "min_max"){
  
  ## calculate min max feature scaling
  if(type == "min_max"){
    feature_scaled = (feature  - limits[1])/(limits[2]-limits[1])
  }
  return(feature_scaled)
}
numfeatures = 14
minmaxvalues = get_featuresmin_maxvals(data_training[!data_training$class%in%"other",1:numfeatures])
dates_pos = grepl("Date", names(minmaxvalues))
#
# ###
#
min_dates = min(unlist(minmaxvalues[dates_pos]))
max_dates = max(unlist(minmaxvalues[dates_pos]))
#
# ## calculate min max derivatives
#
deriva_pos = grepl("deriva", names(minmaxvalues))

min_der = min(unlist(minmaxvalues[deriva_pos]))
max_der = max(unlist(minmaxvalues[deriva_pos]))


dates_scaled = do.call(cbind,lapply(which(dates_pos),function(x){
  feature_scaling(data_training[,x], c(min_dates, max_dates))
}))

derivatives_scaled = do.call(cbind,lapply(which(deriva_pos),function(x){
  feature_scaling(data_training[,x], c(min_der, max_der))
}))

summarize_scaled = do.call(cbind,lapply((1:numfeatures)[-c(which(dates_pos),which(deriva_pos))],
                                        function(x){
                                          
  feature_scaling(data_training[,x], minmaxvalues[[x]])
}))
# feature_scaled = quantable::robustscale(data_training[,1:14])
# 
# dates_scaled_robust = do.call(cbind,lapply(1:14,function(x){
#   feature_scaling(data_training[,x], type = "robust")
# }))

#
# ##
training_iter_scaled = data.frame(dates_scaled, derivatives_scaled,summarize_scaled,
                                  class = data_training$class)

# ## reassign names

row.names(training_iter_scaled) = row.names(data_training)
names(training_iter_scaled)[1:numfeatures] =
  names(data_training)[1:numfeatures]

pca_transform = FactoMineR::PCA(feature_scaled$data[,1:9],ncp = 3)
data_toplot = data.frame(pca_transform$ind$contrib,class =data_training$class)


datasub = data_toplot[!data_toplot$class%in%c("other","vegetative","soil"),]
ggplot(datasub, aes(Dim.1,Dim.2, color = class)) + 
  geom_point()+xlim(-0.000,0.002)+ylim(0.000,0.0012)

datasub = data_toplot#[data_toplot$class%in%c("reproductive","vegetative","other"),]
ggplot(datasub, aes(Dim.1,Dim.2, color = class)) + 
  geom_point()#+xlim(-0.000,0.0025)+ylim(0.000,0.0015)

ggplot(data_toplot, aes(Dim.1,Dim.2, color = class)) + 
  geom_point()+xlim(-0.00,0.002)+ylim(0.000,0.0015)

ggplot(data_toplot[data_toplot$class%in%c("reproductive","vegetative","other","ripening"),], 
       aes(Dim.1,Dim.3, color = class)) + 
  geom_point()+xlim(0.000,0.001)+ylim(0.000,0.002)

ggplot(data_toplot, aes(Dim.2,Dim.3, color = class)) + 
  geom_point()#+ylim(0.00125,0.004)+xlim(0.0001,0.0015)

rowsto_delete = row.names(data_toplot)[(data_toplot$Dim.1>0.0005 & data_toplot$Dim.1<0.001) &
                                         (data_toplot$Dim.2>0.0000 & data_toplot$Dim.2<0.0004) &
                                         grepl(row.names(data_toplot),pattern = "other")]

table(data_training$class)

data_other = data_training



data_other = data_other[data_other$NDVI_Date_1>0.4 ,
                        grepl(names(data_other), pattern = "Date")]


data_other$id = row.names(data_other)

data_other = data_other[data_other$id %in%rowsto_delete,]

data_other = data_other[data_other$NDVI_Date_7>0.60&
                          data_other$NDVI_Date_4>0.65&
                          data_other$NDVI_Date_3>0.6&
                          data_other$NDVI_Date_5>0.63&
                          data_other$NDVI_Date_1<0.65 &
                          data_other$NDVI_Date_1>0.35 &
                          data_other$NDVI_Date_6>0.63 &
                          data_other$NDVI_Date_2<0.7&
                          data_other$NDVI_Date_7<0.80,]

dataplot = reshape2::melt(data_other, by = "id")

num_pos = regexpr(pattern = "_s_",dataplot$id)
num_pos_end = regexpr(pattern = "_d2",dataplot$id)
dataplot$class= stringr::str_sub(dataplot$id,num_pos+3, num_pos_end-1)

ggplot(dataplot, aes(group = id, variable, value, colour = class)) +
  geom_line()+geom_point()

# dim(data_training)
data_training = data_training[!((row.names(data_training)%in% dataplot$id) & ( data_training$class%in%c("other"))),]
table(data_training$class)



rowsto_delete = row.names(data_toplot)[(data_toplot$Dim.1>0.0004 & data_toplot$Dim.1<0.0011) &
                                         (data_toplot$Dim.2>0.0000 & data_toplot$Dim.2<0.0001) &
                                         grepl(row.names(data_toplot),pattern = "other")]


data_other = data_training



data_other = data_other[ data_other$NDVI_Date_1<0.4,
                         grepl(names(data_other), pattern = "Date")]



data_other$id = row.names(data_other)

data_other = data_other[data_other$id %in%rowsto_delete,]

dataplot = reshape2::melt(data_other, by = "id")

num_pos = regexpr(pattern = "_s_",dataplot$id)
num_pos_end = regexpr(pattern = "_d2",dataplot$id)
dataplot$class= stringr::str_sub(dataplot$id,num_pos+3, num_pos_end-1)

ggplot(dataplot, aes(group = id, variable, value, colour = class)) +
  geom_line()+geom_point()
table(data_training[(row.names(data_training)%in% dataplot$id) , "class"])
table(data_training$class)



write.csv(data_training ,
          paste0("model_inputs/phen_identification/optical_data_ndvi_date_deri_veg.csv"))


data_training = read.csv(paste0("model_inputs/phen_identification/optical_data_ndvi_date_deri_veg.csv"),row.names =1)


graphic_s = graphic_smoothed (data_training, vi =  "NDVI_Date_", limitbottom = 0,limitup =1)
head(graphic_s$data)
dataplot = graphic_s$data[graphic_s$data$Typ_Stg%in%c ('harvested', "ripening" ,
                                                       "vegetative", "reproductive" ,"soil","other"),]


dataplot = dataplot[dataplot$variable%in% paste0("NDVI_Date_",1:7),]
dataplot$variable = as.factor(as.character(dataplot$variable)) 
levels(dataplot$variable ) = paste0("NDVI_",1:7)
dataplot$Typ_Stg = factor(dataplot$Typ_Stg, levels = c("vegetative", "reproductive", 
                                                       "ripening", "harvested","soil","other"))
m = ggplot(dataplot, 
           aes(as.factor(variable), value )) + 
  geom_boxplot(aes(fill = Typ_Stg, colour = Typ_Stg), alpha = 0.2)+
  lims(y = c(-0.035,1))+ theme_bw() + 
  labs( fill = "Classification", y = "NDVI", x = "")+ 
  stat_summary(fun.y=median, geom="line", aes(colour = Typ_Stg,group = Typ_Stg), size = 1.4)+ 
  stat_summary(fun.y=median, geom="point", aes(colour = Typ_Stg,group = Typ_Stg), size = 1.8)      



m=m+ guides(colour=FALSE)+
  scale_colour_manual(values = c(phasecolours, other = "gray19"))+
  scale_fill_manual(values = c(phasecolours, other = "gray19"))
m+labs(x = "Features", y = "NDVI", fill = "Growth Phases")+
    theme(text=element_text(size=12),
          axis.title.y=element_text(size = rel(1.1),colour = "black", face="bold"),
          axis.title.x=element_text(size = rel(1.1),colour = "black", face="bold"),
          axis.text.x  = element_text(angle=0, hjust=0.5),
          panel.background=element_rect(fill="white",colour = "black"),
          panel.grid.major = element_line(colour = "gray"),
          legend.title = element_text(face = "bold"))


+ stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                   geom="crossbar", width=0.5)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red")




########################


library(plotly)

plot3d = training_data[!training_data$class%in%c("other","Water","Urban_Zones"),]
plot3d$class = as.factor(as.character(plot3d$class))

p <- plot_ly(plot3d, x = ~NDVI_derivative_1, y = ~NDVI_derivative_2, z = ~NDVI_derivative_3, color = ~class,size = 0.8) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'derivative_1'),
                      yaxis = list(title = 'derivative_2'),
                      zaxis = list(title = 'derivative_3')))
p

test = (training_data)
num_pos = regexpr(pattern = "_p_",row.names(training_data))
num_pos_end = regexpr(pattern = "_s_",row.names(training_data))
test$plygn_n = str_sub(row.names(training_data),num_pos+3, num_pos_end-1)

table(test$class,test$plygn_n)

head(graphic_s$data)

levels(graphic_s$data$variable) = paste0("NDVI_",1:7)

m = ggplot(graphic_s$data[graphic_s$data$Typ_Stg%in%c ('harvested', "ripening", "other" ,
                                                       "vegetative", "reproductive" ,"soil"),], 
           aes(as.factor(variable), value )) + 
  geom_boxplot(aes(fill = Typ_Stg), alpha = 0.4)+
  lims(y = c(-0.035,1))+ theme_bw() + 
  labs( fill = "Classification", y = "NDVI", x = "")+ 
  stat_summary(fun.y=median, geom="line", aes(colour = Typ_Stg,group = Typ_Stg), size = 1.4)+ 
  stat_summary(fun.y=median, geom="point", aes(colour = Typ_Stg,group = Typ_Stg), size = 1.8)      



m=m+ guides(colour=FALSE)+scale_colour_manual(values = c ('harvested'='firebrick3', "ripening"="goldenrod", "other" = "gray", 
                                                        "early_vegetative"="cyan4", "reproductive" = "darkgreen",
                                                        "late_vegetative" = "chartreuse3", "soil" = "saddlebrown"
))+
  scale_fill_manual(values =  c ('harvested'='firebrick3', "ripening"="goldenrod","other" = "gray", 
                                 "early_vegetative"="cyan4", "reproductive" = "darkgreen",
                                 "late_vegetative" = "chartreuse3","soil" = "saddlebrown"
  ))+ stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                      geom="crossbar", width=0.5)+
   stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="pointrange", color="red")

ggsave(plot = m , filename = paste0("graphics/ndvi_time_features_boxplot.jpg"),
       width = 20, height = 12, units = "cm")




test = training_data[,grepl(names(training_data), pattern = "NDVI_de")]

graphic_s = graphic_smoothed (test, vi =  "NDVI", limitbottom = -0.025,limitup =0.025)


m = ggplot(graphic_s$data[graphic_s$data$Typ_Stg%in%c ('harvested', "ripening" ,
                                                       "vegetative", "reproductive" ,"soil"),], 
           aes(as.factor(variable), value )) + 
  geom_boxplot(aes(fill = Typ_Stg), alpha = 0.4)+
  lims(y = c(-0.035,0.035))+ theme_bw() + 
  labs( fill = "Classification", y = "", x = "")+ 
  stat_summary(fun.y=median, geom="line", aes(colour = Typ_Stg,group = Typ_Stg), size = 1.4)+ 
  stat_summary(fun.y=median, geom="point", aes(colour = Typ_Stg,group = Typ_Stg), size = 1.8)      


m=m+ guides(colour=FALSE)+scale_colour_manual(values = phasecolours)+
  scale_fill_manual(values =  phasecolours)+ stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                   geom="crossbar", width=0.5)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="pointrange", color="red")

ggsave(plot = m , filename = paste0("graphics/ndvi_derivative_features_boxplot.jpg"),
       width = 15, height = 10, units = "cm")




test = TableValuesToClassify[,grepl(names(TableValuesToClassify), pattern = "NDVI_Da")]
m = graphic_smoothed (test, vi = "NDVI", iter = "3")


ggsave(plot = m , filename = paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/temp/allData_v8_smoothed_NDVI_date.png"),
       width = 34, height = 25, units = "cm")



#############################################3
crossval_iter = 2



typeData = "Validation"
training_points_radar = SpatialPointsDataFrame(training_points_[,1:2],data = training_points_,
                                               proj4string = crs(SentinelimageReference))

########## ----- > Create Parametrics for radar images 
training_points = training_points_radar

datestoEvaluate = c("20151127","20151211","20151228", "20160114")
DateInt="20151127"
radar_Bands=c("Sigma0_VV_db", 'Sigma0_VV_GLCMMean', 'Sigma0_VV_GLCMVariance')


#### metrics per month

library(doParallel)
library(foreach)
ncores = 2
cl =  makeCluster(ncores)
registerDoParallel(cl)
Sys.time()->start
radar_path = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/satellite_imagery/radar_data/"
calculate_rastermetrics = function(DateStart, DateEnd,radar_bands,tile,
         fpath = "D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/satellite_imagery/radar_data/",
         metrics = c("mean", "min", "max")){
  ## set capturing period
  Date_interval = DateStart : DateEnd
  # filter filenames
  fnames = list.files(fpath, pattern = paste0(tile , ".RData$"))
  
  ## Filter files by dates
  
  fnames = fnames[as.Date(IdentiyImageDate(fnames), format = "%Y%m%d")%in% Date_interval]
  fnames = fnames[order(as.Date(IdentiyImageDate(fnames), format = "%Y%m%d"))]
  
  ## read imagery
  radar_images = lapply(1:length(fnames), function(imag_index){
    file_name = fnames[imag_index]
    
    load(file = paste0(  fpath,file_name))
    
    bands_toSelect = names(sat_imag)[names(sat_imag)%in%radar_Bands]
    sat_imag = sat_imag[[bands_toSelect]]
    cat(str_sub(file_name ,1, -7) ," loaded\n")
    return(sat_imag)
  })
  
  ## calculate metrics
  metrics_features = lapply(radar_bands, function(band){
    band_info = stack(lapply(radar_images , function(x) x[[band]]))
    
    band_ifo = metrics_(band_info,metrics)
  })
  
  stack(metrics_features)
  
}


data_radar = foreach(DateInt = datestoEvaluate, 
                     .export = c("calculate_rastermetrics","radar_Bands", "training_points")) %dopar% {
                       
  library(raster)
  library(stringr)

  ## 30 days metrics
  DateEnd=as.Date(DateInt, format = "%Y%m%d")
  DateStart= DateEnd - 30

  metrics_30days = calculate_rastermetrics ( DateEnd - 30,DateEnd ,radar_bands , tile)

  ## 60  days metrics
  DateStart= (DateEnd ) - ( 60)

  metrics_60days = calculate_rastermetrics (DateStart, ( DateEnd - 31),radar_Bands, tile)
  
  ## 90  days metrics
  DateStart= (DateEnd ) - ( 90)

  metrics_90days = calculate_rastermetrics (DateStart, ( DateEnd - 61),radar_Bands, tile)
  
  ## stack metrics
  metrics= c("mean","min","max")
  radar_data = stack(metrics_30days, metrics_60days, metrics_90days)
  names(radar_data) = do.call(c,lapply(c("30days","60days","90days"),function(y)
    paste0(do.call(c,lapply(radar_Bands,function(x) 
      (paste0(x, "_", metrics)))),"_",y)))

  radar_data_extracted = Extract_raster_data(raster_images = radar_data, spatial_points = training_points [training_points$Date %in% DateInt,] , 
                                             spectral_bands = names(radar_data))
  
  # graph_values_limits = reshape2::melt(radar_data_extracted[[1]])
  # 
  # ggplot(graph_values_limits[grepl(graph_values_limits$variable, pattern = "Sigma0_VV_db_min"),] , aes(variable, value, fill = Typ_Stg)) + geom_boxplot()
  

  radar_data_extracted[[1]]
  
  
}
print(Sys.time()-start)
stopCluster(cl)

#### metrics quantile


Sys.time()->start

data_radar = foreach(DateInt = datestoEvaluate, 
                     .export = c("calculate_raster_quantiles","radar_Bands", "training_points")) %do% {
                       
    library(raster)
    library(stringr)
                       
                       
  DateEnd=as.Date(DateInt, format = "%Y%m%d")
  DateStart= DateEnd - 60
  Date_interval = DateStart : DateEnd
  # filter filenames
  fnames = list.files(radar_path, pattern = paste0(tile , ".RData$"))
  
  ## Filter files by dates
  
  fnames = fnames[as.Date(IdentiyImageDate(fnames), format = "%Y%m%d")%in% Date_interval]
  fnames = fnames[order(as.Date(IdentiyImageDate(fnames), format = "%Y%m%d"))]
  
  ## read imagery
  radar_images = lapply(1:length(fnames), function(imag_index){
    file_name = fnames[imag_index]
    
    load(file = paste0(  radar_path,file_name))
    
    bands_toSelect = names(sat_imag)[names(sat_imag)%in%radar_Bands]
    sat_imag = sat_imag[[bands_toSelect]]
    cat(str_sub(file_name ,1, -7) ," loaded\n")
    return(sat_imag)
  })
  
  ncores = 4
  
  quantile_imags = calculate_raster_quantiles(radar_images,radar_Bands,SentinelimageReference, division = 20, ncores = 4)
  
  quantile_imags = stack(quantile_imags)
  
  names(quantile_imags) = do.call(c,lapply(radar_Bands,function(x) 
      (paste0(x, c("_p05","_p25","_p50","_p75","_p95")))))
  
  
  
  radar_data_extracted = Extract_raster_data(raster_images = quantile_imags, 
                                             spatial_points = training_points [training_points$Date %in% DateInt,] , 
                                             spectral_bands = names(quantile_imags))[[1]]
  
  radar_data_extracted
  #radar_data_extracted$diff = radar_data_extracted$Sigma0_VV_db_p95 - radar_data_extracted$Sigma0_VV_db_p05
  
}

graph_values_limits = reshape2::melt(radar_data_extracted[[1]])
# 
ggplot(graph_values_limits[grepl(graph_values_limits$variable, pattern = "saldana_rasterref"),] , aes(variable, value, fill = Typ_Stg)) + geom_boxplot()

data_ras = data_radar[[1]]
data_ras = do.call(rbind,lapply(data_radar, function(data_ras){
  data_ras$Class = data_ras$Typ_Stg
  
  data_ras$ID = 
    with(data_ras,paste0(as.character(x), as.character(y),"_p_" ,as.character(plygn_n),"_s_",as.character(type),"_type_",as.character(Typ_Stg),"_d" ,as.character(Date), "_i_", as.character(iteration), "_t_",data_type))
  data_ras = data_ras[,-which(names(data_ras)%in% c("x", "y", "plygn_n", "type","Typ_Stg", "Date" , "iteration","data_type"))]
  duplicated_info = which(duplicated(data_ras$ID ))
  if(length(duplicated_info)>0)
    data_ras = data_ras[-duplicated_info,]
  row.names(data_ras) = data_ras$ID
  data_ras = data_ras[,-ncol(data_ras)]
  print(dim(data_ras))
  return(data_ras)
}))

head(data_ras) 
dataras_back = data_ras
data_ras = dataras_back
for(band in radar_Bands){
  data_ras= data_ras[,-ncol(data_ras)]
  data_ras[, paste0(band, "_diffmin_30_60days")] = data_ras[, paste0(band, "_min_30days")] - data_ras[, paste0(band, "_min_60days")]
  data_ras[, paste0(band, "_diffmax_30_60days")] =data_ras[, paste0(band, "_max_30days")] - data_ras[, paste0(band, "_max_60days")]
  
  data_ras[, paste0(band, "_diffmin_60_90days")] =data_ras[, paste0(band, "_min_60days")] - data_ras[, paste0(band, "_min_90days")]
  data_ras[, paste0(band, "_diffmax_60_90days")] =data_ras[, paste0(band, "_max_60days")] - data_ras[, paste0(band, "_max_90days")]
  
  data_ras[, paste0(band, "_diffmin_30_90days")] =data_ras[, paste0(band, "_min_30days")] - data_ras[, paste0(band, "_min_90days")]
  data_ras[, paste0(band, "_diffmax_30_90days")] =data_ras[, paste0(band, "_max_30days")] - data_ras[, paste0(band, "_max_90days")]
  data_ras$Class = dataras_back$Class
}

band = radar_Bands[1]
lapply(radar_Bands ,function(band){
  graph_ = data_ras[,c(grep(names(data_ras) , pattern = paste0(band,"_m")), ncol(data_ras))]
  
  graph_ = reshape2::melt(graph_ )
  t_print = unique(graph_$Class)[1]
  m=ggplot(graph_, aes(variable , value , fill = Class))+geom_boxplot()+
    theme(text=element_text(size=10),
          axis.title.y=element_text(size = rel(1.2),colour = "#999999"),
          axis.title.x=element_text(size = rel(1.2),colour = "#888888"),
          axis.text.x  = element_text(angle=60, hjust=1),
          panel.background=element_rect(fill="white",colour = "black"),
          panel.grid.major = element_line(colour = "gray"))
  ggsave(plot = m , filename = paste0("process/exploratory_images/radar_images/phen_identification_",typeData,"_iterationCV_", crossval_iter ,"_",band,".png"),width = 30, height = 14, units = "cm")
  
} )

str_length(row.names(TableValuesToClassify))
write.csv(data_ras ,paste0("D:/OneDrive - Universidad Nacional de Colombia/MScPhil/phen_identification/model_inputs/phen_identification/RADAR_conf8_radar.csv"))

table(str_sub(row.names(data_ras),str_length(row.names(data_ras))-2,str_length(row.names(data_ras))))
windows()



# 
# DateStart= (DateEnd ) - ( 90)
# Date_interval = DateStart :( DateEnd - 61)
# 
# list_radar_names = list.files(paste0(Main_Folder , "process/rdata/radar_info/") , pattern = paste0(Locality , ".RData$"))
# 
# ## Filter files by dates
# 
# list_radar_names = list_radar_names[as.Date(IdentiyImageDate(list_radar_names), format = "%Y%m%d")%in% Date_interval]
# list_radar_names = list_radar_names[order(as.Date(IdentiyImageDate(list_radar_names), format = "%Y%m%d"))]
# 
# 
# radar_images = read_radar_Rdata (list_radar_names, radar_Bands,
#                                  SentinelimageReference, paste0(Main_Folder , "process/rdata/radar_info/"))
# 
# 
# max_60 = max(stack(radar_images))
# min_60 = min(stack(radar_images))
# mean_60 = mean(stack(radar_images))
# plot(max_30 - max_60)
# radar_data_extracted = Extract_raster_data(raster_images = stack(abs(min_30- max_60),abs(min_60- max_90)), spatial_points = training_points , 
#                                            spectral_bands = names( stack((min_30- min_60),(min_60- min_90))))
# 
# graph_values_limits = reshape2::melt(radar_data_extracted[[1]])
# 
# ggplot(graph_values_limits[grepl(graph_values_limits$variable, pattern = "layer"),] , aes(variable, value, fill = Typ_Stg)) + geom_boxplot()
# 
# 
# dates_images = as.Date(IdentiyImageDate(list_radar_names), format = "%Y%m%d")
# 
# wavelengths =  sapply(dates_images , function(x) x - dates_images[1])
# names(radar_images) = IdentiyImageDate(list_radar_names)
# 
# Sys.time()->start
# data_to_integrated = do.call(cbind,lapply(radar_data_extracted, function(x) x[,"Sigma0_VV_db"] ))
# data_integral = apply(data_to_integrated, 1, function(valraster)calc_integral(wavelengths,valraster))
# 
# head()[,1:6]
# datatoplot = data.frame(value = data_integral, typecov = radar_data_extracted[[1]]$Typ_Stg)
# ggplot(datatoplot, aes(typecov, value))+geom_boxplot()
# 
# dataModel_inegral=lapply (radar_data_extracted,function(image_info){
#   image_info = image_info[,names(central_wavelengths)]/2^12
# })
# print(Sys.time()-start)
# image_info = optical_data_extracted[[1]]
# 
# dataModel_inegral=lapply (optical_data_extracted,function(image_info){
#   image_info = image_info[,names(central_wavelengths)]/2^12
#   data_integral = apply(image_info, 1, function(valraster)calc_integral(wavelengths,valraster))
# })
# 
# rm(quantile_imags)
# 
# rm(radar_images)
# quantile_imags = raster::stack(quantile_imags)
# radar_names_database = do.call(c,lapply(radar_Bands, function(band) paste0(band,"_pa_",1:length(c(0.05,.25, .50, .75,0.95)))))
# names(quantile_imags) = radar_names_database
# 
# radar_min =  list(min_30, min_60,min_90)
# meanValues = mean(stack(radar_min))
# SumbandDiff=list()
# for (j in 1:length(radar_min)){
#   SumbandDiff[[j]]=(radar_min[[j]]-meanValues)*(radar_min[[j]]-meanValues)
# }
# SumbandDiff=sum(stack(SumbandDiff),na.rm=T)
# 
# stdValues=((SumbandDiff*(1/(length(radar_min)-1)))^(1/2))
# 
# 
# rm(raster_temp)
# 
# dataModel_Param_radar = cbind(dataModel_Param_radar, dataModel_Param_radar_conf2)

k <- floor(7/2)
Fm <- matrix(0., 7, 7)

for (row  in  1:(k+1)) {

  Ce <- ( ((1:7)-3) %*% matrix(1, 1, 3+1) ) ^ ( matrix(1, 7) %*% (0:3) )
  
  A <- MASS::ginv(Ce, tol = .Machine$double.eps)
  Fm[row,] <- A[1+0,]
} 

Fm[(k+2):7,] <- (-1)^0 * Fm[k:1,7:1]
x = df_to_fit$y
n = 7
k = floor(n/2)
z = filter(Fm[k+1,n:1], 1, x)
c(Fm[1:k,] %*% x[1:n], z[n:len], Fm[(k+2):n,] %*% x[(len-n+1):len])

