##################
library("dplyr")

library("Kendall")

library("corrplot")

library("ggplot2")
library("ggpubr")
library("ggcorrplot")

##################

create_folderpath = function(pathsnippet1, markername, pathsnippet2, subfolder)
{
  return(paste(pathsnippet1, markername, pathsnippet2, sep = ""))
}

empty_dataframe = function()
{
  data_all = data.frame(X = integer(),
                        Area = numeric(),
                        Mean1 = numeric(),
                        Min1 = integer(),
                        Max1 = integer(),
                        Mean2 = numeric(),
                        Min2 = integer(),
                        Max2 = integer(),
                        Mean3 = numeric(),
                        Min3 = integer(),
                        Max3 = integer(),
                        Origin = character())
}

empty_dataframe_nucl = function()
{
  data_all = data.frame(X = integer(),
                        Area = numeric(),
                        Mean1 = numeric(),
                        Min1 = integer(),
                        Max1 = integer(),
                        Mean2 = numeric(),
                        Min2 = integer(),
                        Max2 = integer(),
                        Mean3 = numeric(),
                        Min3 = integer(),
                        Max3 = integer(),
                        Origin = character(),
                        nucleusArea = integer(),
                        nucleusMean = numeric(),
                        nucleusMin = integer(),
                        nucleusMax = integer(),
                        nucleusROI = numeric(),
                        nucleusIndex = integer())
}

label_nucToROI = function(folderpath, dataset, max_length)
{
  data_nuc = read.csv(paste(paste(folderpath, "nucleus\\", sep = ""), "nuclei_stardist_", dataset, sep = ""), header = TRUE, stringsAsFactors = FALSE)
  data_roi = empty_dataframe()
  associatedROI = data.frame(nucleus = integer()
                             ,roi = integer())
  
  fils = list.files(paste(folderpath, "roi\\", sep = ""))
  
  print("appending ROI datasets")
  for (ROIset in fils)
  {
    data_0=read.csv(paste(paste(folderpath, "roi\\", sep = ""), ROIset, sep = ""), header = TRUE, stringsAsFactors = FALSE)
    data_roi = rbind(data_roi, data_0)
    rm(data_0)
    print(ROIset)
  }
  
  print("labeling nuclei")
  for (row in 1:nrow(data_nuc))
  {
    roiInput = data_roi$iROI[which(data_roi$ROIpixelsX == round(data_nuc$X[row]) & data_roi$ROIpixelsY == round(data_nuc$Y[row]))]
    if (length(roiInput) == 0)
    {
      roiInput = -1
    }
    associatedROI = rbind(associatedROI, list(data_nuc$X.1[row], roiInput))
    print(paste(associatedROI[row,1], associatedROI[row,2]))
  }
  data_nuc$associatedROIcheck = associatedROI[,1]
  data_nuc$associatedROI = associatedROI[,2]
  rm(associatedROI)
  for (check in 1:max_length)
  {
    if (length(which(data_nuc$associatedROI == check)) == 0)
    {
      data_nuc = rbind(data_nuc, list(0, 0, 0, 0, 0, 0, 0, 0, check))
    } else if (length(which(data_nuc$associatedROI == check)) > 1)
    {
      liste = which(data_nuc$associatedROI == check)
      data_nuc = rbind(data_nuc, list(-2, sum(data_nuc$Area[liste]), sum(data_nuc$Mean[liste]*data_nuc$Area[liste])/sum(data_nuc$Area[liste]), min(data_nuc$Min[liste]), max(data_nuc$Max[liste]), mean(data_nuc$X[liste]), mean(data_nuc$Y[liste]), -2, check))
      for (element in liste)
      {
        
        data_nuc$associatedROI[element] = -1
      }
      print(length(data_nuc))
      print(duplicated(data_nuc$associatedROI))
      rm(liste)
    }
  }
  data_nuc = subset(data_nuc,associatedROI != -1)
  
  data_nuc = data_nuc[order(data_nuc$associatedROI),]
  
  print("got nuc data")
  
  rm(data_roi, roiInput)
  return(data_nuc)
}

fill_dataframe = function(folderpath, markerlist, chosenmarker)
{
  
  fils = list.files(paste(folderpath , "protein\\", sep = ""))
  num_files <- length(fils)
  
  data_all = empty_dataframe()
  
  for (dataset in fils)
  {
    data_0=read.csv(paste(paste(folderpath , "protein\\", sep = ""), dataset, sep = ""), header = TRUE, stringsAsFactors = FALSE)
    max_length = nrow(data_0)/3
    data_1 = data_0[which(data_0$Ch == 1),]
    data_2 = data_0[which(data_0$Ch == 2),]
    data_3 = data_0[which(data_0$Ch == 3),]
    X = data_1$X
    Area = data_1$Area
    Mean1 = data_1$Mean
    Min1 = data_1$Min
    Max1 = data_1$Max
    Mean2 = data_2$Mean
    Min2 = data_2$Min
    Max2 = data_2$Max
    Mean3 = data_3$Mean
    Min3 = data_3$Min
    Max3 = data_3$Max
    Origin = rep(dataset, max_length)
    
    print(summary(data_1))
    print(summary(data_2))
    print(summary(data_3))
    
    print(paste(X, Area, Mean1, Min1, Max1, Mean2, Min2, Max2, Mean3, Min3, Max3))
    print(Origin)
    
    print("got all data")
    print(max_length)
    
    data_new = data.frame(X, Area, Mean1, Min1, Max1, Mean2, Min2, Max2, Mean3, Min3, Max3, Origin)
    rm(data_0, data_1, data_2, data_3, X, Area, Mean1, Min1, Max1, Mean2, Min2, Max2, Mean3, Min3, Max3, Origin)
    data_all = rbind(data_all, data_new)
    rm(data_new)
  }
  #  data_all$adjustedMean2 = 0
  comment(data_all) = markerlist[[chosenmarker]]
  return(data_all)
}

fill_dataframe_nucl = function(folderpath, markerlist, chosenmarker)
{
  
  fils = list.files(paste(folderpath , "protein\\", sep = ""))
  num_files <- length(fils)
  
  data_all = empty_dataframe_nucl()
  
  for (dataset in fils)
  {
    data_0=read.csv(paste(paste(folderpath , "protein\\", sep = ""), dataset, sep = ""), header = TRUE, stringsAsFactors = FALSE)
    max_length = nrow(data_0)/3
    data_1 = data_0[which(data_0$Ch == 1),]
    data_2 = data_0[which(data_0$Ch == 2),]
    data_3 = data_0[which(data_0$Ch == 3),]
    X = data_1$X
    Area = data_1$Area
    Mean1 = data_1$Mean
    Min1 = data_1$Min
    Max1 = data_1$Max
    Mean2 = data_2$Mean
    Min2 = data_2$Min
    Max2 = data_2$Max
    Mean3 = data_3$Mean
    Min3 = data_3$Min
    Max3 = data_3$Max
    Origin = rep(dataset, max_length)
    
    print(summary(data_1))
    print(summary(data_2))
    print(summary(data_3))
    
    print(paste(X, Area, Mean1, Min1, Max1, Mean2, Min2, Max2, Mean3, Min3, Max3))
    print(Origin)
    
    print("got all data")
    print(max_length)
    
    data_new = data.frame(X, Area, Mean1, Min1, Max1, Mean2, Min2, Max2, Mean3, Min3, Max3, Origin)
    rm(data_0, data_1, data_2, data_3, X, Area, Mean1, Min1, Max1, Mean2, Min2, Max2, Mean3, Min3, Max3, Origin)
    data_all = rbind(data_all, data_new)
    rm(data_new)
  }
#  data_all$adjustedMean2 = 0
  comment(data_all) = markerlist[[chosenmarker]]
  return(data_all)
}

test_plot = function(input_data)
{
  ggplot(data=input_data, aes(x=Mean1, y=Mean2)) +
    labs(
      title = comment(input_data),
      x = "Calbindin",
      y = "Marker",
      caption = "average 8-bit cell-staining intensity",
      colour = "image origin"
    )+
    theme_light()+
    geom_point(aes(color=Origin))+
    geom_quantile()+
    geom_rug()+
    stat_cor(method = "kendall")
}

batch_import_data = function(base = "", middle, ending = "")
{
  for(middle_entry in middle)
  {
    folder_path = create_folderpath(base, middle_entry, ending)
    assign(middle_entry, fill_dataframe(folder_path))
    test_plot(as.name(middle_entry))
  }
}

##################

basepath = "D:\\Documents\\git\\Forschungsmodul\\img\\input\\"

marker = list("bassoon","betatub3","chat","eaat1","gad65_67","map2","nestin","psd95","synaptophysin","th","vglut","test","testvglut")

pathend = "\\data\\"


#folder = 
test_folder = list()

##################
pick = 13
#assign(marker, empty_dataframe())
#assign(marker, fill_dataframe(create_folderpath(marker, basepath, pathend)))
assign(marker[[pick]], fill_dataframe(create_folderpath(basepath, marker[[pick]], pathend), marker, pick))
summary(eval(as.name(marker[[pick]])))
test_plot(eval(as.name(marker[[pick]])))

#create corrplot https://stackoverflow.com/questions/67608684/how-to-use-corrplot-with-correlation-matrix-created-by-hand-of-type-list
summary(eval(as.name(marker[[pick]])))
meanNucleusArea = mean(eval(as.name(marker[[pick]]))$nucleusArea)
sddevNucleusArea = sd(eval(as.name(marker[[pick]]))$nucleusArea)
upperlimit = meanNucleusArea + sddevNucleusArea
lowerlimit = meanNucleusArea - sddevNucleusArea
test_plot(eval(as.name(marker[[pick]]))[which(eval(as.name(marker[[pick]]))$nucleusArea < upperlimit & eval(as.name(marker[[pick]]))$nucleusArea > lowerlimit & eval(as.name(marker[[pick]]))$nucleusIndex > 0),])

plot(eval(as.name(marker[[pick]]))$Mean1, eval(as.name(marker[[pick]]))$Mean2)

plot(eval(as.name(marker))$Mean1,eval(as.name(marker))$Mean2)

ggcorrplot(cor(select(eval(as.name(marker[[pick]]))[which(eval(as.name(marker[[pick]]))$Mean1 > 5),], Area, Mean1, Mean2, Mean3), method = "kendall"))

corrplot(cor(select(eval(as.name(marker[[pick]])), Area, Mean1, Mean2, Mean3), method = "kendall"),
         diag = FALSE,
         method = 'color',
         order = 'alphabet',
         col = colorRampPalette(c("red", "white", "blue"))(20),
         title = comment(eval(as.name(marker[[pick]]))),
         mar=c(0,0,2,0))

batch_import_data(basepath, marker, pathend)
batch_import_data("", test_folder, "")

#linear regression model
lmMarkerCalb = lm(adjustedMean2 ~ Mean1, data = eval(as.name(marker[[pick]])))
summary(lmMarkerCalb)




#cor.test(data1$Mean, data2$Mean[148:294], method="spearman")
#cor.test(data_all$Mean1, data_all$Mean2, method="kendall")
Kendall(eval(as.name(marker[[pick]]))$Mean1,eval(as.name(marker[[pick]]))$Mean2)
print(marker[[pick]])

ggplot(data=eval(as.name(marker[[pick]]))[which(eval(as.name(marker[[pick]]))$Mean1 > 5),], aes(x=Mean1, y=Mean2)) +
  labs(
    title = comment(eval(as.name(marker[[pick]]))),
    x = "Calbindin",
    y = "Marker",
    caption = "average 8-bit cell-staining intensity",
    colour = "image origin"
  )+
  theme_light()+
  geom_point(aes(color=Origin))+
#  geom_point(aes())+
#  geom_point()+
#  geom_line(aes(y=Mean2, color = "red"))
  geom_quantile()+
  geom_rug()+
#  ggscatterhist(data=eval(as.name(marker[[pick]])), x="nucleusArea", y="Mean1",
#                margin.params = list(fill = "red"))
#                margin.plot = "boxplot"
#              )
  stat_cor(method = "kendall")
