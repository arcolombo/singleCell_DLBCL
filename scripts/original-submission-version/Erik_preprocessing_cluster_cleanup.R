
path="/Volumes/G-DRIVE mobile SSD R-Series/USC DLBCL exp31/custom_gates_0/"


do.call_rbind_fread <- function(path, pattern = "*.csv") {
  files = list.files(path, pattern, full.names = TRUE)
  do.call(rbind, lapply(files, function(x) fread(x, stringsAsFactors = FALSE)[,c(1:50)]))
  #only include the first 50 columns since the number of neighbours(1-10 columns) varies
}


listan=do.call_rbind_fread(path)



names=unique(listan$ImageId)
files = list.files(path, "*.csv", full.names = TRUE)
imageID=data.frame(names,files)
imageID$files=sub(path,"",imageID$files)
imageID$files = substr(imageID$file,2,6)
imageID$files=sub("_","",imageID$files)






listan$ROIID=imageID$files[match(unlist(data.frame(listan$ImageId)), imageID$names)]

b=listan$ROIID
#exclude ROI
c=b!="4058" & b!="4054" & b!="4051"
#table.all.cur=table.all[c,]

listan.exl= listan[c,]

basic.channels=c("Cell_CD20" ,"Cell_CD3" ,"Cell_CD31" ,"Cell_CD4","Cell_CD45RA" , "Cell_CD45RO","Cell_CD68", "Cell_CD8","Cell_FOXP3")

listan.exl$sum=apply(listan.exl[,c("Cell_CD20" ,"Cell_CD3" ,"Cell_CD31" ,"Cell_CD4","Cell_CD45RA" , "Cell_CD45RO","Cell_CD68", "Cell_CD8","Cell_FOXP3")], 1, sum)

#remove cells with no sign signal
listan.exl_tresh50000=listan.exl[listan.exl$sum>50000,] 

set.seed(42)
Rphenograph.listan.exl_tresh50000.45 <- Rphenograph(listan.exl_tresh50000[,c("Cell_CD20" ,"Cell_CD3" ,"Cell_CD31" ,"Cell_CD4","Cell_CD45RA" , "Cell_CD45RO","Cell_CD68", "Cell_CD8","Cell_FOXP3")], k = 45)


phenograph_cluster <- factor(membership(Rphenograph.listan.exl_tresh50000.45[[2]]))
final.list=cbind(listan.exl_tresh50000, phenograph_cluster)


t=final.list$phenograph_cluster

cluster.name=NA
cluster.name[t=="1"] = "CD8"
cluster.name[t=="2"]= "Tumor"
cluster.name[t=="3"]= "Dendritic"
cluster.name[t=="4"]= "Macrophages"
cluster.name[t=="5"]= "Endothelial"
cluster.name[t=="6"]="Tumor"
cluster.name[t=="7"]="Stroma"
cluster.name[t=="8"]="Tumor"
cluster.name[t=="9"]="CD8"
cluster.name[t=="10"]="Treg"
cluster.name[t=="11"]="Tumor"
cluster.name[t=="12"]="Tumor"
cluster.name[t=="13"]="CD4"
cluster.name[t=="14"]="Tumor"
cluster.name[t=="15"]="Tumor"
cluster.name[t=="16"]="Tumor"
cluster.name[t=="17"]="Tumor"
cluster.name[t=="18"]="CD4"
cluster.name[t=="19"]="Tumor"
cluster.name[t=="20"]="CD4"
cluster.name[t=="21"]="Tumor"
cluster.name[t=="22"]="Tumor"
cluster.name[t=="23"]="CD8"
cluster.name=as.factor(cluster.name)

final.list$cell.type=cluster.name

final.list$global.id=c(1:732689)





channels= c("Cell_BCL2"     ,     "Cell_BCL6"  ,
            "Cell_C"         ,    "Cell_CCR4"       ,   "Cell_CD134"    ,     "Cell_CD20"  ,
            "Cell_CD206"    ,     "Cell_CD3"          , "Cell_CD30"    ,      "Cell_CD31"  ,
            "Cell_CD4"         ,  "Cell_CD45RA"     ,   "Cell_CD45RO"    ,    "Cell_CD68"   ,
            "Cell_CD8"      ,     "Cell_CXCR3"         ,
            "Cell_EphrinB2"    ,  "Cell_FOXP3"       ,  "Cell_Granzym"   ,    "Cell_HLADR" ,
            "Cell_ICOS"   ,       "Cell_Ki67"       ,   "Cell_Lag3"      ,
            "Cell_PD1"   ,        "Cell_PDL1"     ,     "Cell_PDL2"  ,
            "Cell_Tbet"     ,     "Cell_Tim3"      ,    "Cell_Vimentin" ,     "Cell_Vista"  ,
            "Cell_p"           )



##

#run phenograph on one sample
#
final.list$phenograph_ind=NA
basic.channels=c("Cell_CD20" ,"Cell_CD3" ,"Cell_CD31" ,"Cell_CD4","Cell_CD45RA" , "Cell_CD45RO","Cell_CD68", "Cell_CD8","Cell_FOXP3")
ROI.ID=unique(final.list$ROIID)
k=15


##add ind clusters column to final list 
for (i in ROI.ID){
  ROI=data.frame(subset(final.list, ROIID==i))
  set.seed(42)
  #CD4.channels=c("Cell_CCR4", "Cell_CXCR3" ,"Cell_ICOS","Cell_Lag3",  "Cell_PD1", "Cell_Tbet","Cell_Tim3","Cell_Vista")
  Rphenograph<- Rphenograph(ROI[,basic.channels], k = k)
  phenograph_cluster<- factor(membership(Rphenograph[[2]]))
  final.list$phenograph_ind[intersect(final.list$global.id,ROI$global.id)]=phenograph_cluster
}

meta=final.list %>% group_by(ROIID, phenograph_ind) %>% summarise(count = n(),Cell_CD20_centroid=median(Cell_CD20), 
                                                                  Cell_CD3_centroid=median(Cell_CD3), 
                                                                  Cell_CD31_centroid=median(Cell_CD31), 
                                                                  Cell_CD4_centroid=median(Cell_CD4), 
                                                                  Cell_CD45RA_centroid=median(Cell_CD45RA), 
                                                                  Cell_CD45RO_centroid=median(Cell_CD45RO),
                                                                  Cell_CD68_centroid=median(Cell_CD68),
                                                                  Cell_CD8_centroid=median(Cell_CD8),
                                                                  Cell_FOXP3_centroid=median(Cell_FOXP3)
)

meta=data.frame(meta)
Rphenograph<- Rphenograph(meta[,c(4:12)], k = k)
phenograph_cluster_meta<- factor(membership(Rphenograph[[2]]))
meta=cbind(meta, phenograph_cluster_meta)
set.seed(42)
tsne <- Rtsne(meta[,c(4:12)], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
list.tsne=as.data.frame(tsne$Y)
colnames(list.tsne)=c("tSNE1.2", "tSNE2.2")
meta=cbind(meta, list.tsne)

#ggplot(meta, aes(x=tSNE1.2, y=tSNE2.2, col=phenograph_cluster_meta, shape = phenograph_cluster_meta)) + geom_point(size = 1)+theme_bw() +scale_shape_manual(values = c(1:length(levels(meta$phenograph_cluster_meta))))





final.list.new=left_join(final.list, meta, by = c("ROIID", "phenograph_ind"))




##rename clusters
t=final.list.new$phenograph_cluster_meta
cluster.name=NA
cluster.name[t=="1"] = "Endothelial"
cluster.name[t=="2"]= "Tcyt"
cluster.name[t=="3"]= "Th"
cluster.name[t=="4"]= "Tumor"
cluster.name[t=="5"]= "Tcyt"
cluster.name[t=="6"]="Th"
cluster.name[t=="7"]="Treg"
cluster.name[t=="8"]="Macrophages"
cluster.name[t=="9"]="Tumor"
cluster.name[t=="10"]="Tumor"
cluster.name[t=="11"]="Th"
cluster.name[t=="12"]="Tumor"
cluster.name[t=="13"]="Tumor"
cluster.name[t=="14"]="Tumor"

cluster.name=as.factor(cluster.name)

final.list.new$cell.type.new=cluster.name








