#*********************************************************************
# R code for drawing map output variables of SEIB-DGVM ver. 2.60      
#                                                                     
# All rights are reserved by Dr. Hisashi SATO                         
# @Graduate School of Environmental Studies, Nagoya University        
#*********************************************************************

#____________________ Set environment ____________________ 
   #Working directory
   setwd('result_visualized/') 
   
   #Missing value
   Missing  =  0
   
   #Grid size (length on a side @ deg)
   #GridSize   = 0.5
   #GridSize   = 1.0
   GridSize   = 2.0
   
   #Pixel size of a grid
   #PixSize   = 1.0 #Suited for 0.5deg simualtion
   #PixSize   = 2.0 #Suited for 1.0deg simualtion
   PixSize   = 4.0 #Suited for 2.0deg simualtion
   
   #Visualizing area (designated by grid numbers)
   #(Global @ 2.0deg mesh)
   LatNoStart =   1 
   LatNoEnd   =  90
   LonNoStart =   1
   LonNoEnd   = 180
   
   #(African continent @ 2.0deg mesh)
   #LatNoStart =  27
   #LatNoEnd   =  63
   #LonNoStart =  80
   #LonNoEnd   = 116
    
   #(African continent @ 1.0deg mesh)
   # LatNoStart =  54
   # LatNoEnd   = 126
   # LonNoStart = 160
   # LonNoEnd   = 232
   
   #(African continent @ 0.5deg mesh)
   #LatNoStart = 108
   #LatNoEnd   = 252
   #LonNoStart = 320
   #LonNoEnd   = 464
   
#____________________ Common procedure ____________________ 

#Compute valiables for coodination and image size
   #Row and column numbers of this coordination system
   LatNoMax   = 180 / GridSize #grid numbers @ vertical   axis
   LonNoMax   = 360 / GridSize #grid numbers @ horizontal axis
   
   #Row and column numbers of drawing area
   Lat      = LatNoEnd - LatNoStart + 1
   Lon      = LonNoEnd - LonNoStart + 1
   
   #Vertical and horizontal picture size @ pixel
   width_size  = Lon * PixSize
   height_size = Lat * PixSize
   
   #Latitude and longitude at the center of grid of the most south-west grid cell
   Lat1_loc =   90.0 - (LatNoEnd  -0.5) * GridSize      #(North:+ , Sourth:-)
   Lon1_loc = -180.0 + (LonNoStart-0.5) * GridSize      #(West :- , East  :+)
   
   #Latitude coordinate
   y <- array(0.0, dim=c(Lat))                          #Prepare array y
   for(i in 1:Lat)   y[i] <- (i-1)*GridSize + Lat1_loc  #Input latitude from south to north
   
   #Longitude coordinate
   x <- array(0.0, dim=c(Lon))                          #Prepare attay x
   for(i in 1:Lon)   x[i] <- (i-1)*GridSize + Lon1_loc  #Input longitude from west to east
   
#Other preparation
   #Activate map library
   library(maps)
   
   #Prepare data aray
   z  <- array(Missing, dim=c(Lon,Lat))
   z1 <- array(Missing, dim=c(Lon,Lat))
   z2 <- array(Missing, dim=c(Lon,Lat))

#____________________ Subroutines for color palette ____________________ 
set_color_topo <- function(num_col) {
   col <-  c(num_col)                #Prepare color variable (this must be division number - 1)
   col <- c(topo.colors(num_col))    #Give color from a palette
   col <- col[length(col):1]         #Turn upside down for color palette
   col[1]  = "white"                 #Change minimum class color
   return(col)                       }
   
set_color_heat <- function(num_col) {
   col <-  c(num_col)                #Prepare color variable (this must be division number - 1)
   col <- c(heat.colors(num_col))    #Give color from a palette
   col <- col[length(col):1]         #Turn upside down for color palette
   col[1]  = "white"                 #Change minimum class color
   return(col)                       }
   
set_color_terrain <- function(num_col) {
   col <-  c(num_col)                #Prepare color variable (this must be division number - 1)
   col <- c(terrain.colors(num_col)) #Give color from a palette
   col <- col[length(col):1]         #Turn upside down for color palette
   col[1]  = "antiquewhite1"         #Change minimum class color
   return(col)                       }

#____________________ Subroutines for data reading and coordination conversion ____________________ 
read_data <- function(fname) {
   d <- read.csv(fname, header=F)
   for (i in 1:Lat) {
   for (j in 1:Lon) {
      #Tuen upside down for latitude, and replace row and column
      z[j,i] <- d[Lat-i+1,j] 
   }
   }
   return(z)
}

#____________________ Subroutine for drawing color pannel ____________________
# LabelName :Label on the top of color pannel
# DivedNum  :Cell number of color pannel
# LabelNum  :Number of value below the color pannel
# IncreStep :Increment of color pannel
draw_panel <- function(LabelName, DivedNum, LabelNum, IncreStep) {
   x_start <-  -15
   x_width <-    5
   y_start <-  -47
   y_width <-    5
   
   #Write label on top of the color pannel
   text(x_start, y_start+2, pos=4, LabelName)
   
   #Draw white belt under color pannel
   x_end <-  x_start+x_width*LabelNum
   polygon( c(x_start, x_end, x_end, x_start), c(y_start, y_start, y_start-y_width, y_start-y_width), col='white') 
   
   #hDraw color pannel
   if (DivedNum==0 || LabelNum==0) {return}
   i <- floor(DivedNum/ LabelNum)
   for (j in 1:DivedNum) {
      polygon( c(x_start, x_start+x_width, x_start+x_width, x_start), c(y_start, y_start, y_start-y_width, y_start-y_width), col=col[j])
      if (floor(j/i) == j/i) {text(x_start+0.5*x_width, y_start-y_width, pos=1, j*IncreStep)}
      x_start <- x_start+x_width
   }
   }

#____________________ Subroutine for drawing a distribution map 1 ____________________
draw_dist <- function(DivedNum,PannelStep,col,x,y,z) {
   PannelMax  <- DivedNum * PannelStep                #カラーパネルの最大値を算出
   br  <- c(seq(from=0, to=PannelMax, by=PannelStep)) #値の分割点の設定１(内部の仕切り)
   br[DivedNum+1]  = PannelMax*10                     #値の分割点の設定２(天井)
   par(mar = c(0,0,0,0) )                             #下・左・上・右の順で内部マージンを設定
   frame()                                            #画面のクリア
   image(x, y, z, breaks=br, col=col, xlab='', ylab='', axes = FALSE) #マッピング
   map(add=T, interior=FALSE)                                         #地図を重ね書き
   }

#____________________ Subroutine for drawing a distribution map 2 ____________________
draw_dist_LAI <- function(DivedNum,PannelStep,col,x,y,z) {
   PannelMax  <- DivedNum * PannelStep                #カラーパネルの最大値を算出
   br  <- c(seq(from=0, to=PannelMax, by=PannelStep)) #値の分割点の設定１(内部の仕切り)
   br[DivedNum+1]  = PannelMax*10                     #値の分割点の設定２(天井)
   br[1]           = 0.01                             #値の分割点の設定３(始点)
   par(mar = c(0,0,0,0) )                             #下・左・上・右の順で内部マージンを設定
   frame()                                            #画面のクリア
   image(x, y, z, breaks=br, col=col, xlab='', ylab='', axes = FALSE) #マッピング
   map(add=T, interior=FALSE)                                         #地図を重ね書き
   }
   
#____________________ Biome ____________________ 
   DivedNum   <- 17 #カラーパネルの分割数
   PannelStep <- 1  #カラーパネルの増分
   
   z <- read_data('out_biome.txt') #データ読みだしと、整形
   
   #色の設定を行う
   col <- set_color_topo(DivedNum) 
   col[ 0] <- "white"         # 0: - water -
   col[ 1] <- "white"         # 1: Polar desert
   col[ 2] <- "violet"        # 2: Arctic/Alpine-tundra
   col[ 3] <- "blue"          # 3: tropical evergreen forest (wet in any month)
   col[ 4] <- "red"           # 4: tropical deciduous forest
   col[ 5] <- "darkseagreen1" # 5: temperate conifer forest
   col[ 6] <- "darkseagreen3" # 6: temperate broad-leaved evergreen forest
   col[ 7] <- "darkseagreen4" # 7: temperate deciduous forest
   col[ 8] <- "skyblue4"      # 8: boreal evergreen forest / woodland
   col[ 9] <- "skyblue1"      # 9: boreal deciduous forest / woodland
   col[10] <- "orange"        #15: xeric woodland / scrub
   col[11] <- "pink"          #16: Grassland / Savanna/ Arid shrubland / Steppe
   col[12] <- "gray"          #17: Desert
   
#   col <- set_color_topo(DivedNum) 
#   col[ 0] <- "white"         # 0: - water -
#   col[ 1] <- "white"         # 1: Polar desert
#   col[ 2] <- "gray"          # 2: Arctic/Alpine-tundra
#   col[ 3] <- "blue"          # 3: tropical rain forest (wet in any month)
#   col[ 4] <- "blue"          # 4: tropical rain forest (seasonaly cycle of water situation)
#   col[ 5] <- "red"           # 5: tropical deciduous forest
#   col[ 6] <- "blue"          # 6: temperate conifer forest
#   col[ 7] <- "blue"          # 7: temperate broad-leaved evergreen forest
#   col[ 8] <- "red"           # 8: temperate deciduous forest
#   col[ 9] <- "blue"          # 9: boreal evergreen forest / woodland
#   col[10] <- "red"           #10: boreal deciduous forest / woodland
#   col[11] <- "green"         #11: short grass land
#   col[12] <- "green"         #12: tall grass land
#   col[13] <- "yellow"        #13: moist savannas
#   col[14] <- "yellow"        #14: dry savannas
#   col[15] <- "yellow"        #15: xeric woodland / scrub
#   col[16] <- "gray"          #16: Arid shrubland / steppe
#   col[17] <- "gray"          #17: Desert
   
   png('out_biome.png', width=width_size, height=height_size) #デバイスドライバ開く
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #図本体の描画
   dev.off()                                                  #デバイスドライバ閉じる
   
#____________________ Fire number ____________________ 
   DivedNum   <- 20     #カラーパネルの分割数
   PannelStep <- 0.025  #カラーパネルの増分
   
   z   <- read_data('out_fire.txt')                               #データ読みだしと、整形
   col <- set_color_topo(DivedNum)                                #色の設定を行う
   
   png('out_fire.png', width=width_size, height=height_size)      #デバイスドライバ開く
   draw_dist (DivedNum, PannelStep, col, x, y, z)                 #図本体の描画
   draw_panel('Fire intensity (n/year)', DivedNum, 4, PannelStep) #カラーパネルを重ね書き
   dev.off()                                                      #デバイスドライバ閉じる
   
#____________________ Biomass ____________________ 
   DivedNum   <- 20   #カラーパネルの分割数
   PannelStep <- 0.75 #カラーパネルの増分
   
   z <- read_data('out_wbiomass.txt')                            #データ読みだしと、整形
   col <- set_color_heat(DivedNum)                               #色の設定を行う
   
   png('out_wbiomass.png', width=width_size, height=height_size) #デバイスドライバ開く
   draw_dist (DivedNum, PannelStep, col, x, y, z)                #図本体の描画
   draw_panel('Biomass (KgC/m2)', DivedNum, 4, PannelStep)       #カラーパネルを重ね書き
   dev.off()                                                     #デバイスドライバ閉じる
   
#____________________ LAI ____________________ 
   DivedNum   <- 14                    #カラーパネルの分割数
   PannelStep <- 0.5                   #カラーパネルの増分
   
   col <- set_color_terrain(DivedNum)  #分割点と色の設定
   
   #描画a1
   z <- read_data('out_lai_amean.txt')                             #データ読みだしと、整形
   png('out_lai_amean.png', width=width_size, height=height_size)  #デバイスドライバ開く
   draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)              #図本体の描画
   draw_panel('Annual Mean LAI (m2/m2)', DivedNum, 7, PannelStep)  #カラーパネルを重ね書き
   dev.off()                                                       #デバイスドライバ閉じる
   
   #描画a2
   z <- read_data('out_lai_amean_t.txt')                                        #データ読みだしと、整形
   png('out_lai_amean_t.png', width=width_size, height=height_size)             #デバイスドライバ開く
   draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)                           #図本体の描画
   draw_panel('Annual Mean LAI of Woody PFTs(m2/m2)', DivedNum, 7, PannelStep)  #カラーパネルを重ね書き
   dev.off()                                                                    #デバイスドライバ閉じる
   
   #描画a3
   z <- read_data('out_lai_amean_g.txt')                                        #データ読みだしと、整形
   png('out_lai_amean_g.png', width=width_size, height=height_size)             #デバイスドライバ開く
   draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)                           #図本体の描画
   draw_panel('Annual Mean LAI of Grass PFT (m2/m2)', DivedNum, 7, PannelStep)  #カラーパネルを重ね書き
   dev.off()                                                                    #デバイスドライバ閉じる
   
#   #描画a4
#   z <- read_data('out_lai_amean_5.txt')                                             #データ読みだしと、整形
#   png('out_lai_amean_5.png', width=width_size, height=height_size)                  #デバイスドライバ開く
#   draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)                                #図本体の描画
#   draw_panel('Annual Mean LAI of PFT5 (m2/m2)', DivedNum, 7, PannelStep)            #カラーパネルを重ね書き
#   dev.off()                                                                         #デバイスドライバ閉じる
#   
#   #描画a5
#   z <- read_data('out_lai_amean_6.txt')                                             #データ読みだしと、整形
#   png('out_lai_amean_6.png', width=width_size, height=height_size)                  #デバイスドライバ開く
#   draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)                                #図本体の描画
#   draw_panel('Annual Mean LAI of PFT6 (m2/m2)', DivedNum, 7, PannelStep)            #カラーパネルを重ね書き
#   dev.off()                                                                         #デバイスドライバ閉じる
   
   #描画b1
   z <- read_data('out_lai_max.txt')                               #データ読みだしと、整形
   png('out_lai_max.png', width=width_size, height=height_size)    #デバイスドライバ開く
   draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)              #図本体の描画
   draw_panel('LAI (m2/m2)', DivedNum, 7, PannelStep)              #カラーパネルを重ね書き
   dev.off()                                                       #デバイスドライバ閉じる
   
   #描画b2
   z <- read_data('out_lai_max_t.txt')                             #データ読みだしと、整形
   png('out_lai_max_t.png', width=width_size, height=height_size)  #デバイスドライバ開く
   draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)              #図本体の描画
   draw_panel('LAI of tree PFTs (m2/m2)', DivedNum, 7, PannelStep) #カラーパネルを重ね書き
   dev.off()                                                       #デバイスドライバ閉じる
   
   #描画b3
   z  <- read_data('out_lai_max_g.txt')
   png('out_lai_max_g.png', width=width_size, height=height_size)   #デバイスドライバ開く
   draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)               #図本体の描画
   draw_panel('LAI of grass PFTs (m2/m2)', DivedNum, 7, PannelStep) #カラーパネルを重ね書き
   dev.off()                                                        #デバイスドライバ閉じる
   
#   #描画b4
#   z <- read_data('out_lai_max_5.txt')                             #データ読みだしと、整形
#   png('out_lai_max_5.png', width=width_size, height=height_size)  #デバイスドライバ開く
#   draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)              #図本体の描画
#   draw_panel('LAI of PFT5 (m2/m2)', DivedNum, 7, PannelStep)      #カラーパネルを重ね書き
#   dev.off()                                                       #デバイスドライバ閉じる
#   
#   #描画b5
#   z <- read_data('out_lai_max_6.txt')                             #データ読みだしと、整形
#   png('out_lai_max_6.png', width=width_size, height=height_size)  #デバイスドライバ開く
#   draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)              #図本体の描画
#   draw_panel('LAI of PFT6 (m2/m2)', DivedNum, 7, PannelStep)      #カラーパネルを重ね書き
#   dev.off()                                                       #デバイスドライバ閉じる
   
#____________________ NPP ____________________ 
   DivedNum   <- 12                #カラーパネルの分割数
   PannelStep <- 0.1               #カラーパネルの増分
   
   z <- read_data('out_npp.txt')   #データ読みだしと、整形
   col <- set_color_topo(DivedNum) #色の設定を行う
   
   png('out_npp.png', width=width_size, height=height_size)   #デバイスドライバ開く
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #図本体の描画
   draw_panel('NPP (gC/m2/yr)', DivedNum, 6, PannelStep*1000) #カラーパネルを重ね書き
   dev.off()                                                  #デバイスドライバ閉じる
   
#____________________ GPP ____________________ 
   DivedNum   <- 12                    #カラーパネルの分割数
   PannelStep <- 0.3                   #カラーパネルの増分
   
   z <- read_data('out_gpp.txt')       #データ読みだしと、整形
   col <- set_color_topo(DivedNum)     #色の設定を行う
   
   png('out_gpp.png', width=width_size, height=height_size)   #デバイスドライバ開く
   draw_dist (DivedNum, PannelStep, col, x, y, z)             #図本体の描画
   draw_panel('GPP (gC/m2/yr)', DivedNum, 6, PannelStep*1000) #カラーパネルを重ね書き
   dev.off()                                                  #デバイスドライバ閉じる
   
#____________________ Soil water content in soil layers 1 to 5 ____________________ 
   DivedNum   <- 20                 #カラーパネルの分割数
   PannelStep <- 0.05               #カラーパネルの増分
   
   z <- read_data('out_water1.txt') #データ読みだしと、整形
   col <- set_color_topo(DivedNum)  #色の設定を行う
   
   png('out_water1.png', width=width_size, height=height_size)                            #デバイスドライバ開く
   draw_dist (DivedNum, PannelStep, col, x, y, z)                                         #図本体の描画
   draw_panel('Soil water saturation @ 0-50cm depth (fraction)', DivedNum, 5, PannelStep) #カラーパネルを重ね書き
   dev.off()                                                                              #デバイスドライバ閉じる
   
#____________________ LAI_every_month ____________________ 
   DivedNum   <- 14                    #カラーパネルの分割数
   PannelStep <- 0.5                   #カラーパネルの増分
   col <- set_color_terrain(DivedNum)  #分割点と色の設定
   
   for (month in 1:12) {
      if (month <10) {z <- read_data(paste('out_lai_month_0', month, '.txt', sep=""))}
      else           {z <- read_data(paste('out_lai_month_' , month, '.txt', sep=""))}
      
      png( paste('out_lai_month',month,'.png'), width=width_size, height=height_size)
         draw_dist_LAI (DivedNum, PannelStep, col, x, y, z)
      dev.off()
   }

#____________________ 降水量と草本LAI_every_month ____________________ 
#  年降水量と年NPPとのプロット
#   Higgins et al. (2000)によると、アフリカのサバナ地帯において
#   草本の地上部生産量(KgDM/ha/yr) = 3.37 × 降水量(mm/yr)
#   
#   d <- read.csv('out_analysis2.txt', header=F)

