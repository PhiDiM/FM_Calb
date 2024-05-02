library("dplyr")

library("Kendall")

library("corrplot")

library("ggplot2")
library("ggpubr")
library("ggcorrplot")

create_folderpath = function(pathsnippet1, markername, pathsnippet2)
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


fill_dataframe = function(folderpath, markerlist, chosenmarker)
{
  
  fils = list.files(folderpath)
  num_files <- length(fils)
  
  data_all = empty_dataframe()
  
  for (dataset in fils)
  {
    data_0=read.csv(paste(folderpath, dataset, sep = ""), header = TRUE, stringsAsFactors = FALSE)
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
    Origin = dataset
    data_new = data.frame(X, Area, Mean1, Min1, Max1, Mean2, Min2, Max2, Mean3, Min3, Max3, Origin)
    rm(data_0, data_1, data_2, data_3, X, Area, Mean1, Min1, Max1, Mean2, Min2, Max2, Mean3, Min3, Max3)
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

basepath = "basepath"

marker = list("bassoon","betatub3","chat","eaat1","gad65_67","map2","nestin","psd95","synaptophysin","th","vglut","test")

pathend = "\\data\\"

#folder = 
test_folder = list()

##################

pick = 7#5,10,11,7
#assign(marker, empty_dataframe())
#assign(marker, fill_dataframe(create_folderpath(marker, basepath, pathend)))
assign(marker[[pick]], fill_dataframe(test_folder[[1]], marker, pick))
summary(eval(as.name(marker[[pick]])))
test_plot(eval(as.name(marker[[pick]])))#[eval(as.name(marker[[12]]))$Area > 195,])

#create corrplot https://stackoverflow.com/questions/67608684/how-to-use-corrplot-with-correlation-matrix-created-by-hand-of-type-list



plot(eval(as.name(marker[[pick]]))$Mean1, eval(as.name(marker[[pick]]))$adjustedMean2)

plot(eval(as.name(marker))$Mean1,eval(as.name(marker))$Mean2)

ggcorrplot(cor(select(eval(as.name(marker[[pick]])), Area, Mean1, Mean2, Mean3, adjustedMean2), method = "kendall"))

corrplot(cor(select(eval(as.name(marker[[pick]])), Area, Mean1, Mean2, Mean3, adjustedMean2), method = "kendall"),
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



#data2=read.csv("C:\\Users\\Phili\\Documents\\git\\Forschungsmodul\\img\\PL_LSM_24.04.2024\\PL LSM 24.04.2024\\testResultsTHlow1.csv",header=TRUE,stringsAsFactors = FALSE)

#cor.test(data1$Mean, data2$Mean[148:294], method="spearman")
#cor.test(data_all$Mean1, data_all$Mean2, method="kendall")
Kendall(eval(as.name(marker[[pick]]))$Mean1,eval(as.name(marker[[pick]]))$adjustedMean2)
print(marker[[pick]])

ggplot(data=eval(eval(as.name(marker[[pick]]))), aes(x=Mean1, y=adjustedMean2)) +
  labs(
    title = deparse(substitute(marker[[pick]])),
    x = "Calbindin",
    y = "Marker",
    caption = "average 8-bit cell-staining intensity",
    colour = "image origin"
  )+
  theme_light()+
  geom_point(aes(color=Origin))+
  geom_quantile()+
  geom_rug()+
  ggscatterhist(data=eval(as.name(marker[[pick]])), x="Mean1", y="adjustedMean2",
                margin.params = list(fill = "red"))
#                margin.plot = "boxplot"
              )
  stat_cor(method = "kendall")
