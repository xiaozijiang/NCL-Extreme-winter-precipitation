;####文件说明：寻找冬季暴雨的研究区域；对研究区域EOF；第一模态时间序列进行天气特征分析（8.3）。

;（1）输出图：中国2144站站点图：1-station；

;（2）计算：冬季极端降水阈值（95百分位），station_95

;（3）输出图：中国极端降水降水频次（2.1-sta_fq）、标准差（2.2-sta_stdv）以及95百分位阈值（2.3-sta_ep95）；

;（4）输出图：中国东南地区EOF（3.1-eof），ts（3.2-eof_ts），REOF(3.3-reof)；
 
;（5）输出柱状图和数据：面积平均的每年极端降水站次（EPI）及其柱状图（4-fq_series）以及EOF6_ts（2-eof_ts6.txt）和EPI指数（2-EPI.txt）。


load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"                                             ;绘制南海和长江黄河

begin

startime            =    get_cpu_time()
;########################START##############一.工作站、变量定义###########################################################################################################################################    

;   datadir         =   "/home/jupiter/jiangsh/"
   datadir         =   "/home/tom/"   
   wks_type         =   "eps"

    minlat_c        =   15                                                                          ; 中国地图边界范围
    maxlat_c        =   55
    minlon_c        =   72
    maxlon_c        =   136

    minlat          =   20                                                                          ; 东南地区地图边界范围（不包括海南岛）
    maxlat          =   33
    minlon          =   106
    maxlon          =   123          
    
    syear           =   1961                                                                        ; 实际开始年份
    eyear           =   2020                                                                        ; 实际结束年份
    yer             =   ispan(syear,eyear,1)    
    sey             =   dimsizes(yer) 
    
    m               =   0     
    n               =   0
;    rscan           =   (/2.0,1.0,0.5/)                                                            ; 连续的有效半径大小，最大为5，依次递减，插值参数       
;    rscan           =   (/10,7,4,2/)                                                               ; 连续的有效半径大小，最大为5，依次递减，插值参数，0.25度的分辨率       

    stdy            =   new(sey,float)                                                              ; 标准差1.0
    stdy!0          =   "yer"
    stdy&yer        =   yer  
    stdy            =   1.0 
    stdya           =   -stdy                                                                       ; 标准差-1.0                        

;#######################START############### .数据读取########################################################################    
  
;1.0; mask land and ocean----------------------------------------------------------------------------

    lsm_f           =   addfile(datadir+"data/5-sfc-land/land-sea-mask.nc","r")                       ; ≤0.5是海洋，＞0.5是陆地；分辨率：0.25度*0.25度 ; short型
    lsm_land        =   short2flt(lsm_f->lsm)
    land            =   lsm_land(0,{minlat:maxlat},{minlon:maxlon})
printVarSummary(land)    
;1.1; read  station-lat-lon  data 2144站点经纬度数据---------分割符为逗号,,----------------------------------------------------------------------------
    name_station  =   asciiread(datadir+"data/1-Interannual/1-f_station2144.txt",-1,"string")   
    delim         =   ","    
    name_num6     =   stringtoint(str_get_field(name_station,1,delim))
    name_lat6     =   stringtofloat(str_get_field(name_station,2,delim))
    name_lon6     =   stringtofloat(str_get_field(name_station,3,delim))
    name_id6      =   stringtoint(str_get_field(name_station,4,delim))                              ; 各省的ID
    name_wgt6     =   stringtofloat(str_get_field(name_station,5,delim))                            ; 各省的面积权重


    stn_2144        =   dimsizes(name_num6)                                                         ; 全国所有站点数（90%无缺测），用于计算：全国暴雨频次  
    stsn_2144       =   ispan(1,stn_2144,1)  
    print(stn_2144)           
;1.2; read pricipation station  data ---1951-2016-------冬季（0-2月）-----------分隔符为逗号，，，，，-------------------
    f_1           =   asciiread(datadir+"data/1-Interannual/1-f_winterain2144_61-20.txt",-1,"string")                                                                   
    delim         =   ","    
    stanum        =   stringtoint(str_get_field(f_1,1,delim))  
    stalat        =   stringtofloat(str_get_field(f_1,2,delim)) 
    stalon        =   stringtofloat(str_get_field(f_1,3,delim))        
    stayear       =   stringtoint(str_get_field(f_1,4,delim))  
    stamon        =   stringtoint(str_get_field(f_1,5,delim)) 
    staday        =   stringtoint(str_get_field(f_1,6,delim))    
    starain       =   stringtofloat(str_get_field(f_1,7,delim)) 
      
    delete(f_1)    
         
;--1.3 ----START------绘制全国站点图:2144站--------------------------------------------------------------------------------------------
   wksds             =    gsn_open_wks(wks_type,datadir+"pic/1-Interannual/2-spatial-temporal/1-stations")                         ; 2144站站点图    
   gsn_define_colormap(wksds,"BlueDarkRed18")        
    
       resg                             =   True
       
       resg@gsnFrame                    =   False
       resg@gsnDraw                     =   False    
       resg@gsnMaximize                 =   False                                                    ; --在工作台中尽可能大绘图（最大化绘图）            
           
       resg@mpGridAndLimbOn             =   False                                                   ; --绘制网格线
       resg@mpGridLineColor             =   "gray"                                                  ; --网格线颜色
       resg@mpGridLatSpacingF           =   5                                                       ; --纬度间隔
       resg@mpGridLonSpacingF           =   5                                                       ; --经度间隔
       resg@tmXTOn                                 = False                                                     ; 去除X轴上部的刻度
       resg@tmYROn                                = False                                                     ; 去除Y轴右侧的刻度
       resg@tiXAxisString                       = ""  ; 去除 x 轴标题
       resg@tiYAxisString                       = ""  ; 去除 y 轴标题                  
       
       resg@mpDataSetName               =   "Earth..4"                                              ; This new database contains
       resg@mpDataBaseVersion           =   "Mediumres"                                             ; Medium resolution database                      
       resg@pmTickMarkDisplayMode       =   "Always"                                                ; --使坐标轴更美观   
              
       resg@mpOutlineOn                 =   True 
       resg@mpOutlineSpecifiers         =   (/"China","Taiwan"/)                                    ; 中国边界有问题，缺藏南、台湾
;       resg@mpOutlineBoundarySets     =   "NoBoundaries"
       resg@mpGeophysicalLineColor      =   "black"                                                 ; 地图边界线的颜色
       
       resg@mpAreaMaskingOn             =   True                                                    ; 使能填充覆盖
       resg@mpMaskAreaSpecifiers        =   (/"China:states","Taiwan","Disputed area between India and China","India:Arunachal Pradesh"/) ; 除括号内的区域其他被抹去（不绘制）  

       resg@mpFillColors                =   (/"transparent","transparent","gray","transparent"/)

       resg@mpLimitMode                 =   "LatLon"
       resg@mpMinLatF                   =   minlat_c                                                ; Asia limits
       resg@mpMaxLatF                   =   maxlat_c
       resg@mpMinLonF                   =   minlon_c
       resg@mpMaxLonF                   =   maxlon_c
       
       resg@tmXBLabelFontHeightF        = 0.025     ;set the size of x-Axis words
       resg@tmYLLabelFontHeightF        = 0.025    ;set the size
; 你还可以使用这两行代码来加粗边界线。

       resg@mpGeophysicalLineThicknessF =   1.0                                                      ; double the thickness of geophysical boundaries
       resg@mpNationalLineThicknessF    =   1.0                                                      ; double the thickness of national boundaries

; 默认的底图投影方式是等经纬度投影，画出来的中国地图比较扁，我们常看到的中国地图，投影方式是兰伯特投影，因此需要对投影方式进行修改

;       resg@mpProjection                =   "LambertConformal"                                     ; 兰伯特投影
;       resg@mpLambertMeridianF          =   110.0
;       resg@mpLimitMode                 =   "LatLon"
;       resg@mpLambertParallel1F         =   0.001                                                  ; Default: .001 ;可以自己改一改，看看投影有什么不同，挺有趣的
;       resg@mpLambertParallel2F         =   89.999                                                 ; Default: 89.999
; 最后将填充区域设定在中国行政区域图之内，如果使用默认效果，等值线会对整个矩形区域填充颜色，因此需要去掉中国边境范围外的填充颜色

       resg@mpOceanFillColor            =   "skyblue"                                                       ; 用白色填充海洋  0是colormap的索引值
       resg@mpInlandWaterFillColor      =   0                                                       ; 用白色填充内陆湖水
       resg@mpLandFillColor             =   "gray88"       

       res_fq                           =   resg                                                    ; 定义暴雨频次的res                     
       map_c                            =   gsn_csm_map(wksds,resg)        
              
;; 中国地图河流等，使用了chinamap数据
       cnres                            =   True
       cnres@china                      =   True                                                    ; draw china map or not
       cnres@river                      =   False                                                    ; draw changjiang&huanghe or not
       cnres@province                   =   False                                                   ; draw province boundary or not
       cnres@nanhai                     =   True                                                    ; draw nanhai or not 
       cnres@diqu                       =   False                                                   ; draw diqujie or not

       chinamap                         =   add_china_map(wksds,map_c,cnres)                        ; use shapefile   
            
; 全国站点图，使用

       mkres                            =   True
       mkres@gsMarkerColor              =   "black"
       mkres@gsMarkerSizeF              =   0.002
       mkres@gsMarkerIndex              =   7                            ; 7 空心三角，16实心圆点
       station_plot                     =   gsn_add_polymarker(wksds,map_c,name_lon6,name_lat6,mkres)
 
       draw(map_c)
       frame(wksds)

;--1.3 ----OVER------绘制全国站点图:2144站------------------------------------

;########################START############## 二.计算：全国极端降水阈值##################################################################################  
 
;2.0------START------1961-2020年冬季全国极端降水阈值（95百分位数）----------------
    
    epi           =   0.95

    station_ep95  =   new(stn_2144,float,-9999)
    station_ep95!0     =   "stsn"
    station_ep95&stsn  =   stsn_2144    
    
    do  i         =    0   ,   stn_2144-1                                                           ; 站点数              
                            
        if(any(stanum.eq.name_num6(i))) then
                
           ind_n :=  ind(stanum.eq.name_num6(i))                                                    ; 某个站
                
           if (any(starain(ind_n).ge.0.1)) then                                                     ; 判断雨量≥0.1
                                             
               ind_n01 := ind(starain(ind_n).ge.0.1)
                          
               starain_q:=starain(ind_n(ind_n01))
              
               qsort(starain_q)                                                                     ; 排序 
                
               epn     :=  dimsizes(starain_q)  
               
               ep95     =  round(epn*epi,3)                                                        ; 95百分位数，3指定ep95为整型              
               
               rain_ep95= starain_q(ep95-1)                                           
               
               station_ep95(i) = rain_ep95 
                                                                                             
           end if      
                                            
        end if
        
        print(i)
                                      
    end do

    printMinMax(station_ep95,True)            
    print("staion_ep95")
    
    
;2.0------END------1961-2020年冬季全国极端降水阈值（95百分位数）----------------

;########################START############## 三.绘图：极端降水频次、标准差和极端降水阈值（95百分位）##################################################################################  

;3.1------START------计算：冬季全国极端降水平均频次、标准差----------------

;3.1.1---计算冬季极端降水频次
    sta_00        =   new((/sey,stn_2144/),float,-9999)                                          ; 每个站每年的站点数量            
    sta_00!0      =   "yer"
    sta_00&yer    =   yer     
    sta_00!1      =   "s_2144"
    sta_00&s_2144 =   stsn_2144        
     
    sta_00         =       0
     
   sta_95         =   sta_00                                                                     ; 每个站每年极端降水的数量
    
    do  i         =    0   ,   eyear-syear                                                       ; 1961-2020冬季

        if (any(stayear.eq.(i+syear)))  then                                                     ; 判断此年是否存在
            
            ind_y:=   ind(stayear.eq.(i+syear))                                                  ; 某一年               
            
            do  j =   0   ,    stn_2144-1                                                        ; 某个站
                            
                if(any(stanum(ind_y).eq.name_num6(j))) then
                
                   ind_yn :=   ind(stanum(ind_y).eq.name_num6(j))                                ; 某个站
                
                   if (any(starain(ind_y(ind_yn)).ge.0.0)) then                                  ; 判断此年是否有站点，只有存在雨量≥0，说明有站
                
                       flag= 1                                                                   ; 此站存在就+1
                 
                       if (any(starain(ind_y(ind_yn)).ge.station_ep95(j))) then                  ; 判断是否有极端降水站    
                            
                           ind_yn95 := ind(starain(ind_y(ind_yn)).ge.station_ep95(j))
                           
                           sta_95(i,j)= num(.not.(ismissing(ind_yn95)))                          ; 每年每个站的极端降水的站数
                           
                       else
                       
                           sta_95(i,j)= 0
                       
                       end if      
                         
                   else
                   
                       flag=0
                   
                   end if
                  
                  
                  sta_00(i,j) = flag                                                             ; 每年每个站是否存在
                  
              
                end if
                
            end do

        end if
        
        print(tostring(1961+i))
                    
    end do
         
    sta_95s       =    dim_sum_n_Wrap(sta_95,0)                                                  ; 暴雨站数，时间维求和
    sta_00s       =    dim_sum_n_Wrap(sta_00,0)                                                  ; 站点年数，时间维求和
    sta_00s       =    where(sta_00s.eq.0,sta_00s@_FillValue,sta_00s)                            ; 如果有站不存在，变成缺失值
    station_fq    =    sta_00s
    station_fq    =    sta_95s/sta_00s                                                           ; 每年的极端降水频次
    
    sta_00        =    where(sta_00.eq.0,sta_00@_FillValue,sta_00)                               ; 如果有站不存在，变成缺失值    
    sta_fq        =    sta_95/sta_00
       
    printVarSummary(station_fq)
    printMinMax(station_fq,True)                     ; max=3.08
    
    print(sta_00(0:20,0:100))
           
;3.1.2---计算冬季暴雨标准差

    sta_95sd      =    dim_stddev_n_Wrap(sta_95,0)    

    printMinMax(sta_95sd,True)                                                                   ; max=0.988       
    
;3.2------START------绘图：1951-2016年冬季全国暴雨平均频次、暴雨标准差绘图station_fq-------2078站数据绘图----------------

;;;;3.2.2.---------作图中国地图上站点图（按照数据大小划分站点）----------------------------- 
 wks_fq        =   gsn_open_wks(wks_type,datadir+"pic/1-Interannual/2-spatial-temporal/2.1-station_fq")                            ; 全国暴雨频次站点图      
 wks_stdv      =   gsn_open_wks(wks_type,datadir+"pic/1-Interannual/2-spatial-temporal/2.2-station_stdv")                          ; 全国暴雨标准差站点图  
    
 gsn_define_colormap(wks_fq,  "BlueDarkRed18")      
 gsn_define_colormap(wks_stdv,"BlueDarkRed18")   
 
;;;;3.2.2.1---------极端降水频次站点图---------------------------------------------------------- 
   
 procedure draw_legend(wks_fq,lat[*][*]:numeric,lon[*][*]:numeric, arr[*]:numeric,colors)
 local gsres, txres, xleg, xtxt, yleg, ytxt, i, labels, nitems

begin

  narr   = dimsizes(arr)
  nmarkers = narr+0
  labels = new(nmarkers,string)

;---Generate the labels for each marker.
  do i = 0, nmarkers-1
  
     if (i.eq.0) then
      labels(i) = arr(0) + " < x < " + arr(1)
     end if
            
     if (i.eq.nmarkers-1) then
      labels(i-1) = "x >= " + max(arr)
     end if
     
     if (i.gt.0.and.i.lt.nmarkers-1) then       
      labels(i) = arr(i) + " <= x < " + arr(i+1)
     end if
     
  end do

;  Create logical variables to hold the marker and text resources.
;  These markers are different than the XY markers, because they are not
;  associated with an XY plot. You can put these markers on any plot.
;
  gsres               = True
  gsres@gsMarkerIndex = 7          ; Use filled dots for markers.

  txres               = True
  txres@txFontHeightF = 0.015

;
; Loop through each grouping of markers, and draw them one set at
; a time, assigning the proper color and size with gsn_marker.
;
; At the same time, draw a legend showing the meaning of the markers.


  xleg = (/0.07,0.07,0.32,0.32,0.57,0.57,0.82,0.82,1.07/)   ; Location of
  xtxt = (/0.16,0.16,0.41,0.41,0.66,0.66,0.91,0.91,1.16/)   ; legend markers
;  yleg = (/0.22,0.17,0.22,0.17,0.22,0.17/)   ; and text
;  ytxt = (/0.22,0.17,0.22,0.17,0.22,0.17/)   ; strings.

  yleg = (/0.12,0.18,0.12,0.18,0.12,0.18,0.12,0.18,0.12/)   ; and text
  ytxt = (/0.12,0.18,0.12,0.18,0.12,0.18,0.12,0.18,0.12/)   ; strings.

  do i = 0, dimsizes(lat(:,0))-1
  
    if (.not.ismissing(lat(i,0)))
      
      if(i.eq.0) then
      
        gsres@gsMarkerIndex    = 7
        
      else
      
        gsres@gsMarkerIndex    = 7
        
      end if        
      
      gsres@gsMarkerColor      = colors(i)
      gsres@gsMarkerThicknessF = 0.8*(i+1)
      gsres@gsMarkerSizeF      = i+2
;---Add marker and text for the legend.
      gsn_polymarker_ndc(wks_fq, xleg(i),yleg(i),gsres)
      gsn_text_ndc      (wks_fq,labels(i),xtxt(i),ytxt(i),txres)
      
    end if
    
  end do

end
;----------------------------------------------------------------------
; This function attaches a labelbar to the given plot.
;----------------------------------------------------------------------

 function attach_labelbar(wks_fq,map,arr[*]:numeric,colors[*])
 local lbres, vph, vpw, nboxes
 
 begin
 
  getvalues map
    "vpHeightF" : vph
    "vpWidthF"  : vpw
  end getvalues

  nboxes = dimsizes(colors)
   
  lbres                    = True          ; labelbar only resources
  lbres@lbAutoManage       = False          ; Necessary to control sizes
  lbres@lbFillColors       = colors
  lbres@vpWidthF           = 1.1 * vpw     ; labelbar width
  lbres@vpHeightF          = 0.3 * vph     ; labelbar height
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.025          ; font height. default is small
  lbres@lbOrientation      = "horizontal"
  lbres@lbPerimOn          = False
  lbres@lbLabelAlignment   = "ExternalEdges"
;--- labelbar
    lbres@lbLabelBarOn                 = True    
   lbres@pmLabelBarOrthogonalPosF     = 0.02    ; ＞0为朝下移动
    lbres@pmLabelBarParallelPosF       = 0.50    ; ＞0为朝右移动    default=0.5，表示图中间
    lbres@pmLabelBarWidthF             = 0.6     ; 设置色标宽度
    lbres@pmLabelBarHeightF            = 0.18    ; 设置色标高度   
  
  lbid = gsn_create_labelbar(wks_fq,nboxes,""+arr,lbres)

  amres                  = True
;  amres@amJust           = "TopCenter"
;  amres@amParallelPosF   =  0.0    ; Center
  amres@amOrthogonalPosF =  0.8    ; Bottom
  annoid = gsn_add_annotation(map,lbid,amres)
  return(annoid)
  
 end
;----------------------------------------------------------------------
; Main code
;----------------------------------------------------------------------
;---Set some needed arrays
  arr    = (/0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5/)   ; bin settings (bin0 = < 0., ;  arr    = (/0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5/)   ; bin settings (bin0 = < 0., 
                                           ; bin1 = 0.:4.999, etc.)
  narr   = dimsizes(arr)
  labels = new(narr,string)                ; Labels for legend.

  npts   = dimsizes(name_lat6)             ; Number of points.
  lat    = name_lat6                       ; Create some dummy latitude
  lon    = name_lon6                       ; and longitude data that
                                           ; will contain the position of
                                           ; our markers.

  R   =  station_fq                        ; This is dummy data for determining
                                           ; how to color the markers.

;
; Create X and Y arrays to hold the points for each range and initialize
; them to missing values.  We want to use num_distinct_markers
; different colors, so we need num_distinct_markers sets of X and
; Y points.
;
  num_distinct_markers = dimsizes(arr)+0        ; number of distinct markers
  lat_new = new((/num_distinct_markers,dimsizes(R)/),float,-999)
  lon_new = new((/num_distinct_markers,dimsizes(R)/),float,-999)

;---Group the points according to which range they fall in.
  do i = 0, num_distinct_markers-1     
  
     if (i.eq.0) then
      indexes = ind(R.gt.arr(0).and.R.lt.arr(1))      
     end if
     
      if (i.eq.num_distinct_markers-1) then
         indexes = ind(R.ge.max(arr))
      end if
      
      if (i.gt.0.and.i.lt.num_distinct_markers-1) then 
         indexes = ind(R.ge.arr(i).and.R.lt.arr(i+1))      
      end if
      
; Now that we have the set of indexes whose values fall within 
; the given range, take the corresponding lat/lon values and store
; them, so later we can color this set of markers with the appropriate
; color.

    if (.not.any(ismissing(indexes))) then 
      npts_range = dimsizes(indexes)   ; # of points in this range.
      
      lat_new(i,0:npts_range-1) = lat(indexes)
      lon_new(i,0:npts_range-1) = lon(indexes)
      
    end if
    delete(indexes)            ; Necessary b/c "indexes" may be a different
                               ; size next time.
  end do

;----------------------------------------------------------------------
; Begin plotting section.
;----------------------------------------------------------------------

   map                               =   gsn_csm_map(wks_fq,resg)                      
; 中国地图河流等，使用了chinamap数据
   chinamap                          =   add_china_map(wks_fq,map,cnres)                        ; use shapefile     
 
;添加边框

  boxlat=(/maxlat,maxlat,minlat,minlat,maxlat/) 
  boxlon=(/minlon,maxlon,maxlon,minlon,minlon/)
  
  lnres= True
  
   lnres@gsLineColor="red"
   lnres@gsLineThicknessF=2.0
   lnres@gsLineDashPattern=0
  
  dum=gsn_add_polyline(wks_fq,map,boxlon,boxlat,lnres)  
     
  gsres               = True
  gsres@gsMarkerIndex = 7                 ; Use filled dots for markers.

;---Get nice spacing through color map for marker colors
  getvalues wks_fq
    "wkColorMapLen" : clen                 ; number of colors in color map
  end getvalues

;  nstep = (clen-2)/narr
;  colors = ispan(2,clen-4,nstep)
  colors = (/6,4,2,13,14,16,18/)
;---Loop through each "bin" and attach the markers to the map.

  do i = 0, num_distinct_markers-1
    if (.not.ismissing(lat_new(i,0)))
    
      if(i.eq.0)
         gsres@gsMarkerIndex = 7          ; Use filled dots for markers.
      else
         gsres@gsMarkerIndex = 7          ; Use filled dots for markers.
      
      end if  
      gsres@gsMarkerColor      = colors(i)
      gsres@gsMarkerSizeF      = i+2
      dumstr = unique_string("marker")
      map@$dumstr$ = gsn_add_polymarker(wks_fq,map,lon_new(i,:),lat_new(i,:),gsres)
    end if
  end do

;----------------------------------------------------------------------
; First version of this plot has a legend at the bottom
;----------------------------------------------------------------------
;  draw(map)                  ; Drawing the map will draw the markers
;  draw_legend(wks_fq,lat_new,lon_new,arr,colors) ; Draw a legend at the bottom

;  frame(wks_fq)    ; Advance the frame 

;----------------------------------------------------------------------
; Second version of this plot has a labelbar added.
;----------------------------------------------------------------------
  lbid = attach_labelbar(wks_fq,map,arr,colors)                                                     ; Attach a labelbar
  
  chinamap                          =   add_china_map(wks_fq,map,cnres)                             ; use shapefile       

  draw(map)                                                                                         ; Drawing the map will draw everything
  frame(wks_fq)                                                                                     ; Advance the frame.

;;;;3.2.2.2----------画图-极端降水频次的标准差------------------------------------------------------------------------ 

procedure draw_legend2(wks_stdv,lat[*][*]:numeric,lon[*][*]:numeric, arr[*]:numeric,colors)
 local gsres, txres, xleg, xtxt, yleg, ytxt, i, labels, nitems

begin

  narr   = dimsizes(arr)
  nmarkers = narr+0
  labels = new(nmarkers,string)

;---Generate the labels for each marker.
  do i = 0, nmarkers-1
  
     if (i.eq.0) then
      labels(i) = arr(0) + " < x < " + arr(1)
     end if
            
     if (i.eq.nmarkers-1) then
      labels(i-1) = "x >= " + max(arr)
     end if
     
     if (i.gt.0.and.i.lt.nmarkers-1) then       
      labels(i) = arr(i) + " <= x < " + arr(i+1)
     end if
     
  end do
;  Create logical variables to hold the marker and text resources.
;  These markers are different than the XY markers, because they are not
;  associated with an XY plot. You can put these markers on any plot.
;
  gsres               = True
  gsres@gsMarkerIndex = 7          ; Use filled dots for markers.

  txres               = True
  txres@txFontHeightF = 0.015

;
; Loop through each grouping of markers, and draw them one set at
; a time, assigning the proper color and size with gsn_marker.
;
; At the same time, draw a legend showing the meaning of the markers.


  xleg = (/0.07,0.07,0.32,0.32,0.57,0.57,0.82,0.82,1.07/)   ; Location of
  xtxt = (/0.16,0.16,0.41,0.41,0.66,0.66,0.91,0.91,1.16/)   ; legend markers
;  yleg = (/0.22,0.17,0.22,0.17,0.22,0.17/)   ; and text
;  ytxt = (/0.22,0.17,0.22,0.17,0.22,0.17/)   ; strings.

  yleg = (/0.12,0.18,0.12,0.18,0.12,0.18,0.12,0.18,0.12/)   ; and text
  ytxt = (/0.12,0.18,0.12,0.18,0.12,0.18,0.12,0.18,0.18/)   ; strings.

  do i = 0, dimsizes(lat(:,0))-1
  
    if (.not.ismissing(lat(i,0)))
      
      if(i.eq.0) then
      
        gsres@gsMarkerIndex    = 7
        
      else
      
        gsres@gsMarkerIndex    = 7
        
      end if        
      
      gsres@gsMarkerColor      = colors(i)
      gsres@gsMarkerThicknessF = 0.8*(i+1)
      gsres@gsMarkerSizeF      = i+3
;---Add marker and text for the legend.
      gsn_polymarker_ndc(wks_stdv, xleg(i),yleg(i),gsres)
      gsn_text_ndc      (wks_stdv,labels(i),xtxt(i),ytxt(i),txres)
      
    end if
    
  end do

end
;----------------------------------------------------------------------
; This function attaches a labelbar to the given plot.
;----------------------------------------------------------------------

 function attach_labelbar2(wks_stdv,map,arr[*]:numeric,colors[*])
 local lbres, vph, vpw, nboxes
  begin
 
  getvalues map
    "vpHeightF" : vph
    "vpWidthF"  : vpw
  end getvalues

  nboxes = dimsizes(colors)
   
  lbres                    = True          ; labelbar only resources
  lbres@lbAutoManage       = False          ; Necessary to control sizes
  lbres@lbFillColors       = colors
  lbres@vpWidthF           = 1.1 * vpw     ; labelbar width
  lbres@vpHeightF          = 0.3 * vph     ; labelbar height
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.025          ; font height. default is small
  lbres@lbOrientation      = "horizontal"
  lbres@lbPerimOn          = False
  lbres@lbLabelAlignment   = "ExternalEdges"
;--- labelbar
    lbres@lbLabelBarOn                 = True    
   lbres@pmLabelBarOrthogonalPosF     = 0.10    ; ＞0为朝下移动
    lbres@pmLabelBarParallelPosF       = 0.50    ; ＞0为朝右移动    default=0.5，表示图中间
    lbres@pmLabelBarWidthF             = 0.6     ; 设置色标宽度
    lbres@pmLabelBarHeightF            = 0.18    ; 设置色标高度     
  
  lbid = gsn_create_labelbar(wks_stdv,nboxes,""+arr,lbres)
  amres                  = True
;  amres@amJust           = "TopCenter"
;  amres@amParallelPosF   =  0.0    ; Center
; amres@amOrthogonalPosF =  0.6    ; Bottom
  annoid = gsn_add_annotation(map,lbid,amres)
  return(annoid)
  
 end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;标注；；；；；；；；；；；；；；；；；；；；；；；；；；；；；；
;----------------------------------------------------------------------
; Main code
;----------------------------------------------------------------------
;---Set some needed arrays
  arr    = (/0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5/); bin settings (bin0 = < 0., 
                                           ; bin1 = 0.:4.999, etc.)
  narr   = dimsizes(arr)
  labels = new(narr,string)                ; Labels for legend.

  npts   = dimsizes(name_lat6)              ; Number of points.
  lat    = name_lat6                        ; Create some dummy latitude
  lon    = name_lon6                        ; and longitude data that
                                           ; will contain the position of
                                           ; our markers.

  R   =  sta_95sd                          ; This is dummy data for determining
                                           ; how to color the markers.

;
; Create X and Y arrays to hold the points for each range and initialize
; them to missing values.  We want to use num_distinct_markers
; different colors, so we need num_distinct_markers sets of X and
; Y points.
;
  num_distinct_markers = dimsizes(arr)        ; number of distinct markers
  lat_new = new((/num_distinct_markers,dimsizes(R)/),float,-999)
  lon_new = new((/num_distinct_markers,dimsizes(R)/),float,-999)

;---Group the points according to which range they fall in.
  do i = 0, num_distinct_markers-1     
  
     if (i.eq.0) then
      indexes = ind(R.gt.arr(0).and.R.lt.arr(1))      
     end if
     
      if (i.eq.num_distinct_markers-1) then
         indexes = ind(R.ge.max(arr))
      end if
      
      if (i.gt.0.and.i.lt.num_distinct_markers-1) then 
         indexes = ind(R.ge.arr(i).and.R.lt.arr(i+1))      
      end if
; Now that we have the set of indexes whose values fall within 
; the given range, take the corresponding lat/lon values and store
; them, so later we can color this set of markers with the appropriate
; color.
    if (.not.any(ismissing(indexes))) then 
      npts_range = dimsizes(indexes)   ; # of points in this range.
      
      lat_new(i,0:npts_range-1) = lat(indexes)
      lon_new(i,0:npts_range-1) = lon(indexes)
      
    end if
    delete(indexes)            ; Necessary b/c "indexes" may be a different
                               ; size next time.
  end do

;----------------------------------------------------------------------
; Begin plotting section.
;----------------------------------------------------------------------

;  map = gsn_csm_map(wks_stdv,mpres)
;        resg@gsnLeftString     = "(c) Standard deviation "
   map                        =   gsn_csm_map(wks_stdv,resg)                      
; 中国地图河流等，使用了chinamap数据
   chinamap                   =   add_china_map(wks_stdv,map,cnres)                        ; use shapefile       
    
;;添加边框
;
;  boxlat=(/maxlat,maxlat,minlat,minlat,maxlat/)
;  
;  boxlon=(/minlon,maxlon,maxlon,minlon,minlon/)
;  
;  lnres= True
;  
;   lnres@gsLineColor="black"
;   lnres@gsLineThicknessF=2.0
;   lnres@gsLineDashPattern=0
  
  dum=gsn_add_polyline(wks_stdv,map,boxlon,boxlat,lnres)  
     

  gsres               = True
  gsres@gsMarkerIndex = 7          ; Use filled dots for markers.

;---Get nice spacing through color map for marker colors
  getvalues wks_stdv
    "wkColorMapLen" : clen     ; number of colors in color map
  end getvalues

;  nstep = (clen-2)/narr
;  colors = ispan(2,clen-4,nstep)
  colors = (/6,4,2,13,14,16,18/)
;---Loop through each "bin" and attach the markers to the map.
  do i = 0, num_distinct_markers-1
    if (.not.ismissing(lat_new(i,0)))
    
      if(i.eq.0)
         gsres@gsMarkerIndex = 7          ; Use filled dots for markers.
      else
         gsres@gsMarkerIndex = 7          ; Use filled dots for markers.
      
      end if  
      gsres@gsMarkerColor      = colors(i)
      gsres@gsMarkerSizeF      = i+2
      dumstr = unique_string("marker")
      map@$dumstr$ = gsn_add_polymarker(wks_stdv,map,lon_new(i,:),lat_new(i,:),gsres)
    end if
  end do
   
  lbid = attach_labelbar(wks_stdv,map,arr,colors)   ; Attach a labelbar

  chinamap                          =   add_china_map(wks_stdv,map,cnres)                           ; use shapefile       
    
  draw(map)                                ; Drawing the map will draw the markers
                                           ;  draw_legend(wks_stdv,lat_new,lon_new,arr,colors) ; Draw a legend at the bottom   
  frame(wks_stdv)
;3.2.2------OVER------1958-2016年冬季全国暴雨频次、标准差绘图station_fq-------2078站------------


;;;;2.2.2.3----------画图-极端降水阈值95百分位数------------------------------------------------------------------------ 

 wks_ep95       =   gsn_open_wks(wks_type,datadir+"pic/1-Interannual/2-spatial-temporal/2.3-sta_ep95")                          ; 全国极端降水阈值95百分位数  
 gsn_define_colormap(wks_ep95,"BlueDarkRed18")   
 

procedure draw_legend3(wks_ep95,lat[*][*]:numeric,lon[*][*]:numeric, arr[*]:numeric,colors)
 local gsres, txres, xleg, xtxt, yleg, ytxt, i, labels, nitems

begin

  narr   = dimsizes(arr)
  nmarkers = narr+0
  labels = new(nmarkers,string)

;---Generate the labels for each marker.
  do i = 0, nmarkers-1
  
     if (i.eq.0) then
      labels(i) = arr(0) + " < x < " + arr(1)
     end if
            
     if (i.eq.nmarkers-1) then
      labels(i-1) = "x >= " + max(arr)
     end if
     
     if (i.gt.0.and.i.lt.nmarkers-1) then       
      labels(i) = arr(i) + " <= x < " + arr(i+1)
     end if
     
  end do
;  Create logical variables to hold the marker and text resources.
;  These markers are different than the XY markers, because they are not
;  associated with an XY plot. You can put these markers on any plot.
;
  gsres               = True
  gsres@gsMarkerIndex = 7          ; Use filled dots for markers.

  txres               = True
  txres@txFontHeightF = 0.015

;
; Loop through each grouping of markers, and draw them one set at
; a time, assigning the proper color and size with gsn_marker.
;
; At the same time, draw a legend showing the meaning of the markers.


  xleg = (/0.07,0.07,0.32,0.32,0.57,0.57,0.82,0.82,1.07/)   ; Location of
  xtxt = (/0.16,0.16,0.41,0.41,0.66,0.66,0.91,0.91,1.16/)   ; legend markers
;  yleg = (/0.22,0.17,0.22,0.17,0.22,0.17/)   ; and text
;  ytxt = (/0.22,0.17,0.22,0.17,0.22,0.17/)   ; strings.

  yleg = (/0.12,0.18,0.12,0.18,0.12,0.18,0.12,0.18,0.12/)   ; and text
  ytxt = (/0.12,0.18,0.12,0.18,0.12,0.18,0.12,0.18,0.18/)   ; strings.

  do i = 0, dimsizes(lat(:,0))-1
  
    if (.not.ismissing(lat(i,0)))
      
      if(i.eq.0) then
      
        gsres@gsMarkerIndex    = 7
        
      else
      
        gsres@gsMarkerIndex    = 7
        
      end if        
      
      gsres@gsMarkerColor      = colors(i)
      gsres@gsMarkerThicknessF = 0.8*(i+1)
      gsres@gsMarkerSizeF      = i+2
;---Add marker and text for the legend.
      gsn_polymarker_ndc(wks_ep95, xleg(i),yleg(i),gsres)
      gsn_text_ndc      (wks_ep95,labels(i),xtxt(i),ytxt(i),txres)
      
    end if
    
  end do

end
;----------------------------------------------------------------------
; This function attaches a labelbar to the given plot.
;----------------------------------------------------------------------

 function attach_labelbar3(wks_ep95,map,arr[*]:numeric,colors[*])
 local lbres, vph, vpw, nboxes
  begin
 
  getvalues map
    "vpHeightF" : vph
    "vpWidthF"  : vpw
  end getvalues

  nboxes = dimsizes(colors)
   
  lbres                    = True          ; labelbar only resources
  lbres@lbAutoManage       = False          ; Necessary to control sizes
  lbres@lbFillColors       = colors
  lbres@vpWidthF           = 1.1 * vpw     ; labelbar width
  lbres@vpHeightF          = 0.3 * vph     ; labelbar height
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.025          ; font height. default is small
  lbres@lbOrientation      = "horizontal"
  lbres@lbPerimOn          = False
  lbres@lbLabelAlignment   = "ExternalEdges"
;--- labelbar
    lbres@lbLabelBarOn                 = True    
   lbres@pmLabelBarOrthogonalPosF     = 0.10    ; ＞0为朝下移动
    lbres@pmLabelBarParallelPosF       = 0.50    ; ＞0为朝右移动    default=0.5，表示图中间
    lbres@pmLabelBarWidthF             = 0.6     ; 设置色标宽度
    lbres@pmLabelBarHeightF            = 0.18    ; 设置色标高度     
  
  lbid = gsn_create_labelbar(wks_ep95,nboxes,""+arr,lbres)
  amres                  = True
;  amres@amJust           = "TopCenter"
;  amres@amParallelPosF   =  0.0    ; Center
;  amres@amOrthogonalPosF =  0.6    ; Bottom
  annoid = gsn_add_annotation(map,lbid,amres)
  return(annoid)
  
 end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;标注；；；；；；；；；；；；；；；；；；；；；；；；；；；；；；
;----------------------------------------------------------------------
; Main code
;----------------------------------------------------------------------
;---Set some needed arrays
  arr   := (/0,5,10,15,20,25,30,35,40/)       ; bin settings (bin0 = < 0., 
                                           ; bin1 = 0.:4.999, etc.)
  narr   = dimsizes(arr)
  labels:= new(narr,string)                ; Labels for legend.

  npts   = dimsizes(name_lat6)              ; Number of points.
  lat    = name_lat6                        ; Create some dummy latitude
  lon    = name_lon6                        ; and longitude data that
                                           ; will contain the position of
                                           ; our markers.

  R   :=  station_ep95                      ; This is dummy data for determining
                                           ; how to color the markers.

;
; Create X and Y arrays to hold the points for each range and initialize
; them to missing values.  We want to use num_distinct_markers
; different colors, so we need num_distinct_markers sets of X and
; Y points.
;
  num_distinct_markers = dimsizes(arr)        ; number of distinct markers
  lat_new  := new((/num_distinct_markers,dimsizes(R)/),float,-999)
  lon_new  := new((/num_distinct_markers,dimsizes(R)/),float,-999)

;---Group the points according to which range they fall in.
  do i = 0, num_distinct_markers-1     
  
     if (i.eq.0) then
      indexes = ind(R.gt.arr(0).and.R.lt.arr(1))      
     end if
     
      if (i.eq.num_distinct_markers-1) then
         indexes = ind(R.ge.max(arr))
      end if
      
      if (i.gt.0.and.i.lt.num_distinct_markers-1) then 
         indexes = ind(R.ge.arr(i).and.R.lt.arr(i+1))      
      end if
; Now that we have the set of indexes whose values fall within 
; the given range, take the corresponding lat/lon values and store
; them, so later we can color this set of markers with the appropriate
; color.
    if (.not.any(ismissing(indexes))) then 
      npts_range = dimsizes(indexes)   ; # of points in this range.
      
      lat_new(i,0:npts_range-1) = lat(indexes)
      lon_new(i,0:npts_range-1) = lon(indexes)
      
    end if
    delete(indexes)            ; Necessary b/c "indexes" may be a different
                               ; size next time.
  end do

;----------------------------------------------------------------------
; Begin plotting section.
;----------------------------------------------------------------------

;  map = gsn_csm_map(wks_ep95,mpres)
;         resg@gsnLeftString     = "(a) 95 percentile "  
   map                        =   gsn_csm_map(wks_ep95,resg)                      
; 中国地图河流等，使用了chinamap数据
  chinamap                   =   add_china_map(wks_ep95,map,cnres)                        ; use shapefile   
  
;;添加边框
;
;  boxlat=(/maxlat,maxlat,minlat,minlat,maxlat/)
;  
;  boxlon=(/minlon,maxlon,maxlon,minlon,minlon/)
;  
;  lnres= True
;  
;   lnres@gsLineColor="black"
;   lnres@gsLineThicknessF=2.0
;   lnres@gsLineDashPattern=0
  
  dum=gsn_add_polyline(wks_ep95,map,boxlon,boxlat,lnres)  
     

  gsres               = True
  gsres@gsMarkerIndex = 7          ; Use filled dots for markers.

;---Get nice spacing through color map for marker colors
  getvalues wks_ep95
    "wkColorMapLen" : clen     ; number of colors in color map
  end getvalues

;  nstep = (clen-2)/narr
;  colors = ispan(2,clen-4,nstep)
  colors := (/8,6,4,2,13,14,16,18/)
;---Loop through each "bin" and attach the markers to the map.
  do i = 0, num_distinct_markers-1
    if (.not.ismissing(lat_new(i,0)))
    
      if(i.eq.0)
         gsres@gsMarkerIndex = 7          ; Use filled dots for markers.
      else
         gsres@gsMarkerIndex = 7          ; Use filled dots for markers.
      
      end if  
      gsres@gsMarkerColor      = colors(i)
      gsres@gsMarkerSizeF      = i+2
      dumstr = unique_string("marker")
      map@$dumstr$ = gsn_add_polymarker(wks_ep95,map,lon_new(i,:),lat_new(i,:),gsres)
    end if
  end do
   
  lbid = attach_labelbar(wks_ep95,map,arr,colors)   ; Attach a labelbar

  chinamap                          =   add_china_map(wks_ep95,map,cnres)                           ; use shapefile       
    
  draw(map)                                ; Drawing the map will draw the markers
                                           ;  draw_legend(wks_ep95,lat_new,lon_new,arr,colors) ; Draw a legend at the bottom   
  frame(wks_ep95)
;3.3------OVER------1961-2021年冬季全国极端降水阈值station_ep95-------2144站------------


;########################START############## 四.计算：研究区域：东南区域##################################################################################  

;4.1------START------ 筛选出东南区域：(20-33N,106,123E)----779站----------------------------------

   ind_d           =    ind((name_lat6.ge.minlat).and.(name_lat6.le.maxlat))
   
   ind_dn          =    ind((name_lon6(ind_d).ge.minlon).and.(name_lon6(ind_d).le.maxlon))
   
   dns             =    dimsizes(ind_dn)                                                            ; 东南站点的个数：650

;--4.1.1 经纬度数据，转换为650站数据，-------------------   
   name_num       :=    name_num6(ind_d(ind_dn))                                                    ; 转换为东南650站
   name_lat       :=    name_lat6(ind_d(ind_dn))
   name_lon       :=    name_lon6(ind_d(ind_dn))
   name_id        :=    name_id6(ind_d(ind_dn))
   name_wgt       :=    name_wgt6(ind_d(ind_dn))
   
   ep_95           =    station_ep95(ind_d(ind_dn))                                                 ; 东南地区的极端降水阈值   
   st_fq           =    sta_fq(:,ind_d(ind_dn))                                                     ; 东南地区的极端降水频次
   
   print(name_num)   
   printVarSummary(st_fq)     
   
   stn             =   dns                                                                              ; 东南区（20-33N，106-123E）站点数，用于计算：东南地区站点 
   stsn            =   ispan(1,dns,1) 
   
   s_fq            =     dim_sum_n_Wrap(st_fq,1)                                                  ; 站点年数，
print(s_fq)
;--4.1.2 雨量等数据，转换为650站数据，-------------------stnum,,,779站数据:st开头的变量
   ny              =   4205613
   
   stnum           =   new(ny,integer,-9999)                                                        ; 东南669站数据的行数  
   styear          =   new(ny,integer,-9999)         
   stmon           =   new(ny,integer,-9999)         
   stday           =   new(ny,integer,-9999)         
   strain          =   new(ny,float,-9999)         
   stlat           =   new(ny,float,-9999)         
   stlon           =   new(ny,float,-9999)  
              
   sx              =    0                                                                           ; 计数器  

   do    i         =    0  ,   stn-1
     
         if(any(stanum.eq.name_num(i))) then
         
           indx   :=    ind(stanum.eq.name_num(i))
           
           n       =    dimsizes(indx)
           
           stnum(sx:sx+n-1)  = stanum(indx)
           
           styear(sx:sx+n-1) = stayear(indx)
           
           stmon(sx:sx+n-1)  = stamon(indx)
           
           stday(sx:sx+n-1)  = staday(indx)
         
           strain(sx:sx+n-1) = starain(indx)
           
           stlat(sx:sx+n-1)  = stalat(indx)
           
           stlon(sx:sx+n-1)  = stalon(indx)
           
;           sta_95s(:,sx:sx+n-1) = sta_95
           
           sx      =    sx + n
                            
         end if
         
   end do  
   
           print(sx)     
;4.1------OVER------- 筛选出东南区域：(20-33N,106-123E)----已全部转化为st(stnum)开头（779站），以st开头的变量                                                          

;########################START############## 五.绘图：东南区域极端降水频次的EOF及ts##################################################################################  

;5.1------计算所有年份极端降水频次------------hr_frq 

    hr_frq          =    new((/sey,stn/),integer)
    hr_frq!0        =    "year"
    hr_frq!1        =    "station"
    hr_frq&year     =    yer
    hr_frq&station  =    name_num
    
    do  i           =    0   ,   eyear-syear
    
        if(any(styear.eq.(syear+i))) then                                                              ; 某年
           
           ind_i   :=    ind(styear.eq.(syear+i))
           
           do  j    =    0   ,   stn-1                                                              ; 站点数              
                
               n    =    0  
                            
               if(any(stnum(ind_i).eq.name_num(j))) then
                
                 ind_ij :=  ind(stnum(ind_i).eq.name_num(j))                                        ; 某个站                 
                                
                 if (any(strain(ind_i(ind_ij)).ge.ep_95(j))) then                       ; 判断雨量≥极端降水阈值
                                             
                    ind_ijh := ind(strain(ind_i(ind_ij)).ge.ep_95(j))
                      
                    n        = num(.not.ismissing(ind_ijh))
                                                                                    
                 end if
                                            
               end if

               hr_frq(i,j) = n 
                                                    
           end do
           
        end if
    
    end do



    frq_grid  = hr_frq  

;5.2.1------计算：极端降水频次的EOF和REOF---------------------------------------------------------------------------------------------------------    
    neof                      =    10
    
    eof_hr                    =    eofunc_n_Wrap(frq_grid,neof,False,0)
    eofts_hr                  =    eofunc_ts_n_Wrap(frq_grid,eof_hr,False,0)
     
    print(eofts_hr(0,:))
    print(max(eof_hr))
    
;;;;;EOF_ts的标准化和滑动平均
    eofts_sd                  =    dim_standardize_n(eofts_hr,1,1)                                  ; 时间序列标准化
    dim_eofts_sd              =    runave_n_Wrap(eofts_sd,9,0,1)                                    ; 时间序列九年平滑   

;;;;;EOF_ts扩大100倍，便于画图

    eof_hr1                   =    eof_hr                                                           ; 扩大100倍，便于画图
    eof_hr1                   =    eof_hr*100
    eof_hr                    =    eof_hr1                                                          ; 扩大100倍，便于画图
    

    dtrd_tshr     =    dtrend_msg_n(yer,eofts_hr,False,True,1)                                     ; 时间序列去趋势

;;;;;north检验
    
    sig_pcv                   =    eofunc_north(eof_hr@pcvar,sey,False)                             ; north检验,前四个模态显著

;;;;;REOF分析
    
    sig_neof                  =    num(sig_pcv)                                                     ; EOF显著的模态数量
    reof_hr                   =    eofunc_varimax_Wrap(eof_hr,1)                                    ; REOF分析
    
    eof_varimax_reorder (reof_hr)                                                                   ; reorder the information
    printVarSummary(reof_hr)

;5.2.2------绘图：EOF和REOF极端降水频次------------  

    wks_eof                   =   gsn_open_wks(wks_type,datadir+"pic/1-Interannual/2-spatial-temporal/3.1-eof")                      ; eof    
    gsn_define_colormap(wks_eof,"rainbow+gray")                                                     ; define a different colormap.
    
    wks_eofts                 =   gsn_open_wks(wks_type,datadir+"pic/1-Interannual/2-spatial-temporal/3.2-eof_ts")                   ; eof    

    wks_reof                  =   gsn_open_wks(wks_type,datadir+"pic/1-Interannual/2-spatial-temporal/3.3-reof")                     ; eof    
    gsn_define_colormap(wks_reof,"rainbow+gray")                                                    ; define a different colormap. 
    
    plot_eof                  = new(neof,"graphic")      
    plot_eofts                = new(neof,"graphic")    
    mean_plot                 = new(neof,"graphic")      
    
;;;;EOF绘图 差值为1度*1度进行绘图

;;5.2.1-----------站点格点化------0.25度*0.25度-----中国东南地区范围：(20-33N,106-123N)--------------------     

    ratio                     =   4                                                                 ; 分辨率0.25=1.0/4
      
    blon                      =   fspan(minlon,maxlon,(maxlon-minlon)*ratio+1)                      ; 106-123；(123-106)*1+1=121            
    blat                      =   fspan(minlat,maxlat,(maxlat-minlat)*ratio+1)                      ; 20-30  ； (30-20)*8+1=81    
    
    blon!0                    =   "lon"
    blon@units                =   "degrees-east"
    blon&lon                  =   blon

    blat!0                    =   "lat"
    blat@units                =   "degrees_north"
    blat&lat                  =   blat
;;         
;;;----每年的极端降水频次插值    
    rscan                     =   (/4.0,2.0,0.5/)                                                     ; 连续的有效半径大小，依次递减，插值参数，0.25度的分辨率
           
    eof_grid                  =   new((/neof,dimsizes(blat),dimsizes(blon)/),float,-9999)             ; 每年极端降水日数格点化 
    eof_grid!0                =   "nnum"
    eof_grid&nnum             =   ispan(1,10,1)
    eof_grid!1                =   "lat"
    eof_grid&lat              =    blat          
    eof_grid!2                =   "lon"
    eof_grid&lon              =    blon        
      
    do  i =   0 , neof-1                                                                              ; 避免站点对应不上，先每年站点暴雨日数插值
    
        grid                  = obj_anal_ic_Wrap(name_lon,name_lat,eof_hr(i,:),blon,blat,rscan, False)    ; Creanm插值 插值成1度*1度
               
        eof_grid(i,:,:)       = where(land.gt.0.5,grid,eof_grid@_FillValue)
        
              print(i)
               printVarSummary(eof_grid(i,:,:))
               printMinMax(eof_grid(i,:,:),0)
               print(num(.not.ismissing(eof_grid(i,:,:))))
        
    end do
    printVarSummary(eof_grid)    
 
    map_eof                   = new(neof,"graphic")  
    
;---Set up some map resources.   
    mpres                     =   True

    mpres@gsnFrame            = False                                                               ; Don't advance the frame
    mpres@gsnDraw             = False                                                               ; Don't advance the frame
    mpres@gsnMaximize         = False                                                               ; Maximize plot in frame.
    
    mpres@mpMinLatF           = minlat
    mpres@mpMaxLatF           = maxlat
    mpres@mpMinLonF           = minlon
    mpres@mpMaxLonF           = maxlon
  
    mpres@mpDataSetName               =   "Earth..4"                                           ; This new database contains
    mpres@mpDataBaseVersion           =   "Mediumres"                                          ; Medium resolution database                      
    mpres@pmTickMarkDisplayMode       =   "Always"                                             ; --使坐标轴更美观   
           
    mpres@mpGridAndLimbOn             =   False                                                ; --绘制网格线
    mpres@mpGridLineColor             =   "gray50"                                               ; --网格线颜色
    mpres@mpGridLatSpacingF           =   5                                                    ; --纬度间隔
    mpres@mpGridLonSpacingF           =   5                                                    ; --经度间隔
       
    mpres@mpOutlineOn                 =   True 
    mpres@mpOutlineSpecifiers         =   (/"China","Taiwan","Disputed area between India and China","India:Arunachal Pradesh"/)     ; 中国边界有问题，缺藏南、台湾
    mpres@mpGeophysicalLineColor      =   "gray50"  
  
    mpres@mpAreaMaskingOn             =   True                                                 ;使能填充覆盖
    mpres@mpMaskAreaSpecifiers        =   (/"China:states","Taiwan","Disputed area between India and China","India:Arunachal Pradesh"/) ; 除括号内的区域其他被抹去（不绘制）  

    mpres@mpFillColors                = (/"transparent","transparent","gray","transparent"/)
;-------------------------------
    res_eof                   = True
    res_eof@gsnMaximize       = False 
    res_eof@gsnFrame          = False                                                               ; Don't advance the frame
    res_eof@gsnDraw           = False                                                               ; Don't advance the frame  
;    
;    res_eof@sfXArray          = name_lon                                                           ; 站点资料时需要定义坐标！！！
;    res_eof@sfYArray          = name_lat                 

    res_eof@cnLinesOn            =   True   
    res_eof@cnSmoothingOn        =   True                                                           ; 平滑线
    res_eof@cnLineLabelsOn       =   True
    res_eof@cnLineThicknessF     =   0.7                                                            ; 等值线粗细  
    res_eof@cnInfoLabelOn        = False                                                            ; 打开图右下方的等值线信息标签
    res_eof@cnSmoothingDistanceF = 0.005                                                            ; 平滑曲线
    res_eof@cnSmoothingTensionF  = 0.005
    res_eof@cnLineLabelBackgroundColor     = -1
    res_eof@cnLineLabelPlacementMode       ="Constant"
    res_eof@cnLineLabelFontHeightF         = 0.01
    res_eof@cnLineLabelInterval = 1                                                                 ; 每条线上都有标签    

   sign = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)"/)
   
    do n=0,neof-1
    
       res_eof@gsnLeftString  = "(a) "+"EOF "+(n+1)
       res_eof@gsnRightString = sprintf("%5.1f", eof_hr@pcvar(n)) +"%"
       res_eof@cnFillOn       = False
       res_eof@cnLinesOn      = True                                                                ; 绘制等值线    
       res_eof@cnLevelSelectionMode = "ExplicitLevels"                                              ; AutomaticLevels  ;ExplicitLevels
       res_eof@cnLevels      := (/-7,-5,-3,-1,0,1,3,5,7/)                                           ; 颜色        
       plot_eof(n)            = gsn_csm_contour(wks_eof,eof_grid(n,:,:),res_eof)                      ; 
       plot_eof(n)            = ColorNegDashZeroPosContour(plot_eof(n),"blue","black","red")       
       map_eof(n)             = gsn_csm_map(wks_eof,mpres)
       
       overlay(map_eof(n),plot_eof(n))      
       
    end do
    
;;  panel plot only resources
    resP                      = True                                                                ; modify the panel plot
    resP@gsnMaximize          = True                                                                ; large format
;  resP@gsnPanelLabelBar    = True                                                                  ; add common colorbar
    resP@gsnPanelRowSpec      = True
  
print("1234")    
  gsn_panel(wks_eof,map_eof(0:sig_neof-1),(/sig_neof/),resP)                                        ; only plot the 1st four   

; ----------------time series (principal component) plot------------------------------------------------------------------------
   rts                        = True
   rts@gsnDraw                = False                                                               ; don't draw yet
   rts@gsnFrame               = False                                                               ;don't advance frame yet

   rts@vpHeightF              = 0.35                                                                ; Changes the aspect ratio
   rts@vpWidthF               = 0.50
   rts@vpXF                   = 0.10                                                                ; change start locations
   rts@vpYF                   = 0.75                                                                ; the plot

  rts@trXMinF	                = min(yer)-1                                                          ; 设定X轴坐标范围
  rts@trXMaxF	                = max(yer)+1

  res2                        = rts   
  
  rts@gsnXYBarChart           = True                                                                ; 直方柱bar的宽度
  rts@gsnXYBarChartBarWidth   = 0.5
  rts@tmXBOn                  = True                                                                ; X轴刻度
  rts@tmXBMode                = "Manual"
  rts@tmXBTickStartF          = syear
  rts@tmXBTickSpacingF        = 5
  rts@tmXBLabelFontHeightF    = 0.010
  rts@gsnYRefLine             = 0.                                                                  ; reference line   
  rts@gsnAboveYRefLineColor   = "black"                                                             ; above ref line fill red
  rts@gsnBelowYRefLineColor   = "gray"                                                              ; below ref line fill blue
; panel plot only resources
  rtsP                        = True                                                                ; modify the panel plot
  rtsP@gsnMaximize            = True                                                                ; large format
   
  do n    = 0, neof-1
  
     rts@gsnLeftString        = "PC"+(n+1)
     rts@gsnRightString       = sprintf("%5.1f", eof_hr@pcvar(n)) +"%"
     plot_eofts(n)            = gsn_csm_xy (wks_eofts,yer,eofts_sd(n,:),rts)                        ; 标准化后的时间序列
          
     res2@xyLineThicknessF    = 3.0
     
     mean_plot(n)             = gsn_csm_xy(wks_eofts,yer,dim_eofts_sd(n,:),res2)                    ; 标准化时间序列后---再9年滑动平均；
     
     overlay(plot_eofts(n),mean_plot(n))     
     
  end do
  
  rtsP@gsnPanelRowSpec        = True 
   
  gsn_panel(wks_eofts,plot_eofts(0:sig_neof-1),(/sig_neof/),rtsP)                                   ; draw all 'neof' as one plot  
  

;; EOF1-PC1--time series (principal component) plot----------绘制柱状图;------------------------------------
  wks_pc1                     =  gsn_open_wks(wks_type,datadir+"pic/1-Interannual/2-spatial-temporal/3.4-PC1")                       ; eof    

  resbar                      = True       
   resbar@gsnDraw             = False                                                               ; 暂不画
   resbar@gsnFrame            = False                                                               ; 暂不翻页
   resbar@tmXTOn              = False                                                               ; 去除X轴上部的刻度
   resbar@tmYROn              = False                                                               ; 去除Y轴右侧的刻度   
   
   res2                       = resbar          
       
   resbar@vpHeightF           = 0.35                                                                ; 改变viewport的大小
   resbar@vpWidthF 	          = 0.6                                                                 ; 
   resbar@trXMinF	            = min(yer)-1                                                          ; 设定X轴坐标范围
   resbar@trXMaxF	            = max(yer)+1
   resbar@trYMaxF             = 4.0                                                                 ; 设定Y轴坐标范围
   resbar@trYMinF             = -2.0 
   resbar@gsnXYBarChart       = True                                                                ; 直方柱bar的宽度
   resbar@gsnXYBarChartBarWidth= 0.5
   
   resbar@tmXBOn              = True                                                                ; X轴刻度
   resbar@tmXBMode            = "Manual"
   resbar@tmXBTickStartF      = syear
   resbar@tmXBTickSpacingF    = 5
   resbar@tmXBLabelFontHeightF= 0.010

   resbar@gsnYRefLine         = 0.                                                                  ; 设定参考值。  ;设定bar的颜色
   resbar@gsnAboveYRefLineColor= "red"                                                            ; 大于该参考值所用的颜色
   resbar@gsnBelowYRefLineColor= "blue"                                                             ; 小于该参考值所用的颜色

  ;resbar@gsnXYBarChartColors = (/"black"/)                                                         ; 如果不设定参考值,则可通过该resbarources来设定	颜色 
   
;       resbar@tiXAxisString	           = tostring(yer(indstd1))
;       resbar@tiYAxisString	           = tostring(yer(indstd2))                                   ; 添加X轴名称   
;       resbar@tiYAxisString	= "precipitation (units: mm)"                                         ; 添加Y轴名称 
  
   top_plot                   =  gsn_csm_xy(wks_pc1,yer,eofts_sd(0,:),resbar)                        ; 标准化的PC1
           
   res2@xyLineColor           =  "black"
   res2@xyLineThicknessF      =  2.0
       
   std_plot                   =  gsn_csm_xy(wks_pc1,yer,stdy,res2)                                  ; 0.6标准差
   overlay(top_plot,std_plot)      
   stda_plot                  = gsn_csm_xy(wks_pc1,yer,stdya,res2)
   overlay(top_plot,stda_plot)
       
   res2@xyLineThicknessF      = 3.0
   mean9_plot                 = gsn_csm_xy(wks_pc1,yer,dim_eofts_sd(0,:),res2)                      ; 9年滑动平均
   overlay(top_plot,mean9_plot)      
    
   draw(top_plot)                                                                                   ; 绘图
   frame(wks_pc1) 
 
;;    rEOF绘图----------------------------------------------------------------------------------------------------------
    map_reof = new(neof,"graphic")  
    
     res_reof            = res_eof
;    res_reof            = True     
;    res_reof@gsnMaximize= False 
;    res_reof@gsnFrame   = False            ; Don't advance the frame
;    res_reof@gsnDraw    = False            ; Don't advance the frame            
;    res_reof@sfXArray   = name_lon      
;    res_reof@sfYArray   = name_lat
;    res_reof@cnMinLevelValF       =   -0.1
;    res_reof@cnMaxLevelValF       =   0.5
;    res_reof@cnLevelSpacingF      =   0.1 
 
    plot_reof           = new(neof,"graphic")    
    
    do n=0,neof-1
    
     res_reof@gsnLeftString     = "EOF "+(n+1)
     
     res_reof@gsnRightString    = sprintf("%5.1f", reof_hr@pcvar_varimax(n)) +"%"
     
     plot_reof(n)=gsn_csm_contour(wks_reof,reof_hr(n,:),res_reof)
     
     plot_reof(n)=ColorNegDashZeroPosContour(plot_reof(n),"black","black","black")
     
     map_reof(n) = gsn_csm_map(wks_reof,mpres)
       
     overlay(map_reof(n),plot_reof(n))      
     
    end do
  gsn_panel(wks_reof,map_reof(0:sig_neof-1),(/sig_neof/),resP)     ; now draw as one plot

;5------OVER------1958-2016年冬季东南的暴雨频次EOF及REOF分区，绘图hr_reof----650站数据绘图------------
    
;##六###################START########(暴雨频次&面积平均)的年际变化##############################################################################    

;--6.1----------------计算每年每省极端降水频次的面积平均的时间序列
  
        s_wgt     =    new(sey,float,-9999)                                                     ; 每年的极端降水频次面积平均
        s_wgt!0        =   "yer"
        s_wgt&yer      =   yer

    do  i         =    0   ,   eyear-syear                                                      ; 1961-2020年冬季

        if (any(styear.eq.(i+syear)))  then                                                     ; 判断此年是否存在
            
            indwy:=   ind(styear.eq.(i+syear))                                                  ; 某年
             
            p_fq    =  new(31,float,-9999)                                                      ; 某一年东南区域站点频次
            p_st    =  new(31,float,-9999)                                                      ; 区域内省站点占省总站点比例
            p_wgt   =  new(31,float,-9999)                                                      ; 区域内省站点占省总站点比例
            
          ;6.1.1 start-省区域内极端降水频次
          
           do j  =    0   ,   30                                                                ; 判断在某个省的极端降水频次
            
               indid  :=  ind(name_id.eq.(j+1))                                                 ; 东南区域的省内站点数               
               indp   :=  ind(name_id6.eq.(j+1))                                                ; 省内站点数               
                              
               if(all(.not.(ismissing(indid))))then          

                 p_fq(j)  =  sum(st_fq(i,indid))                                                ; 某省的极端降水频次和
                 
                 ps_num   =  dimsizes(indid)                                                    ; 区域省站点数
                 p_num    =  dimsizes(indp)                                                     ; 省站点数
                 
                 ;6.1.2 start-区域省站站省总面积权重                                 
                 p_st(j)  =  ps_num/p_num                                                       ; 区域内省站点占省总站点比例

                 ;6.1.3 start-省面积权重                 
                 p_wgt(j) = name_wgt(indid(0))                                                  ; 某个省的面积权重                                                
                   
               end if
  
           end do 

           
           ; 综合指数：极端降水频次面积平均
           s_wgt(i) = sum(p_fq * p_st * p_wgt)         
           
        end if

    end do

;--6.2----------------极端降水频次标准化
    sd            =  1.0
     
    dim_wgt       =  dim_standardize_n_Wrap(s_wgt,0,0)                                          ; 数据标准化
    dim_mean      =  runave_Wrap(dim_wgt,9,0)                                                   ; 九年平滑
    
    dtrd_fq       =  dtrend_msg_n(s_wgt&yer,s_wgt,False,True,0)                                 ; 去趋势    
    dtst_fq       =  dim_standardize_n_Wrap(dtrd_fq,0,0)                                        ; 去趋势后标准化
    
    
    indstd1       =  ind(dim_wgt.gt.sd)                                                         ; 0.8个标准差
    indstd2       =  ind(dim_wgt.lt.-sd)                                                        ; -0.8个标准差
     
    print(yer(indstd1))
    print(dim_wgt(indstd1))
    
    print(yer(indstd2))
    print(dim_wgt(indstd2))  

    print(dim_wgt)
    
    print(dtst_fq)

;####6.3################################绘制柱状图;##############################################

    wksstd          =   gsn_open_wks(wks_type,datadir+"/pic/1-Interannual/2-spatial-temporal/4-fq_series")             ; 极端降水站点数标准化      
    gsn_define_colormap(wksstd,"BlueDarkRed18")                                                     ; define a different colormap.
    
       res_bar                           = True       

       res_bar@gsnDraw                   = False                                                  ; 暂不画
       res_bar@gsnFrame                  = False                                                  ; 暂不翻页
       res_bar@tmXTOn                    = False                                                  ; 去除X轴上部的刻度
       res_bar@tmYROn                    = False                                                  ; 去除Y轴右侧的刻度 
       
       res2                              = res_bar          
       
       res_bar@vpHeightF                 = 0.35                                                   ; 改变viewport的大小
       res_bar@vpWidthF 	               = 0.6                                                    ; 

       res_bar@trXMinF	                 = min(yer)-1                                             ; 设定X轴坐标范围
       res_bar@trXMaxF	                 = max(yer)+1
       res_bar@trYMaxF                   = 3.5                                                    ; 设定Y轴坐标范围
       res_bar@trYMinF                   = -2.0 
       res_bar@gsnXYBarChart             = True                                                   ; 直方柱bar的宽度
       res_bar@gsnXYBarChartBarWidth     = 0.5

       res_bar@tmXBOn                    = True                                                   ; X轴刻度
       res_bar@tmXBMode                  = "Manual"
       res_bar@tmXBTickStartF            = min(yer)-1
       res_bar@tmXBTickSpacingF          = 5
       res_bar@tmXBLabelFontHeightF      = 0.010

; 设定bar图的颜色
       res_bar@gsnYRefLine               = 0.                                                     ; 设定参考值。  ;设定bar的颜色
       res_bar@gsnAboveYRefLineColor     = "red"                                                  ; 大于该参考值所用的颜色 40的颜色比60深
       res_bar@gsnBelowYRefLineColor     = "blue"                                                 ; 小于该参考值所用的颜色

  ;res_bar@gsnXYBarChartColors = (/"black"/)                                                      ; 如果不设定参考值,则可通过该res_barources来设定	颜色 
  
;; 设定x轴和y轴的名称   
;       res_bar@tiXAxisString	           = "year"
;;       res_bar@tiYAxisString	= "precipitation (units: mm)"                                      ; 添加Y轴名称 

;EPI指数（标准化，未去趋势）   
       t_plot = gsn_csm_xy(wksstd,yer,dim_wgt,res_bar)                                           ; 极端降水站点数标准化绘图
     
;正负1倍标准差(线)   
    
       res2@xyLineColor                 = "red"
       res2@xyLineThicknessF            = 1.0
       
       std_plot = gsn_csm_xy(wksstd,yer,stdy,res2)
       overlay(t_plot,std_plot)   
       
       res2@xyLineColor                 = "blue"
       
       stda_plot= gsn_csm_xy(wksstd,yer,stdya,res2)
       overlay(t_plot,stda_plot)

;EPI9年滑动平均指数       
;       res2@xyLineThicknessF            = 3.0
;       mavg_plot:= gsn_csm_xy(wksstd,yer,dim_mean,res2) 
;       overlay(t_plot,mavg_plot)  

;eof_ts1(标准化指数)         
       res2@xyLineColor                 = "gray40"
       res2@xyLineThicknessF            = 2.0
       
       ts1_plot = gsn_csm_xy(wksstd,yer,eofts_sd(0,:),res2)
       overlay(t_plot,ts1_plot)    
    
       draw(t_plot)                                                                            ; 绘图
       frame(wksstd) 

;####6.4################################ EPI 绘制柱状图;##############################################

; 线性趋势图
    rc   =  regCoef(yer,dim_wgt)                         ; 线性相关系数
    
    yreg =  rc*yer + rc@yintercept                       ; 计算的y值
    
    copy_VarCoords(dim_wgt,yreg)
    print(rc)
    print(yreg)
    
 ; 线性趋势检验
 
    p_reg = student_t(rc@tval,rc@nptxy-2)
    print(p_reg)   

    wks_epi          =   gsn_open_wks(wks_type,datadir+"/pic/1-Interannual/2-spatial-temporal/5-epi")             ; 极端降水站点数标准化      
    gsn_define_colormap(wksstd,"BlueDarkRed18")                                                     ; define a different colormap.
    
;EPI指数（标准化，未去趋势）   

       t_plot = gsn_csm_xy(wks_epi,yer,dim_wgt,res_bar)                                           ; 极端降水站点数标准化绘图
     
;;正负1倍标准差(线)   
;    
;       res2@xyLineColor                 = "red"
;       res2@xyLineThicknessF            = 1.0
;       
;       std_plot = gsn_csm_xy(wks_epi,yer,stdy,res2)
;       overlay(t_plot,std_plot)   
;       
;       res2@xyLineColor                 = "blue"
;       
;       stda_plot= gsn_csm_xy(wks_epi,yer,stdya,res2)
;       overlay(t_plot,stda_plot)

; 线性趋势
       res_reg                           = True       

       res_reg@gsnDraw                   = False                                                  ; 暂不画
       res_reg@gsnFrame                  = False                                                  ; 暂不翻页
       res_reg@tmXTOn                    = False                                                  ; 去除X轴上部的刻度
       res_reg@tmYROn                    = False                                                  ; 去除Y轴右侧的刻度 
       res_reg@xyLineColor               = "black"       
       res_reg@xyLineThicknessF          = 2.0
       
        res_reg@gsnLeftString     = "(b)  Time series of frequency of extreme precipitation "

       reg_plot = gsn_csm_xy(wks_epi,yer,yreg,res_reg) 
       
       overlay(t_plot,reg_plot)  
    
       draw(t_plot)                                                                            ; 绘图
       frame(wks_epi) 

;##七###################相关计算(二者必须都去趋势)########EOF_ts1，与EPI（极端降水面积平均指数）################################################################    

;;;;;;;;;相关计算

   r_pc1EPI        =  escorc(dtrd_tshr(0,:),dtrd_fq)
   p_pc1EPII        =  rtest(r_pc1EPI,sey,0)  
   
   print(r_pc1EPI)
   print(p_pc1EPII)

;##八###################START########输出：EOF_ts1，以及极端降水面积平均指数：EPI##############################################################################    

  EPI =  s_wgt                                      ; 未经过任何处理的极端降水面积平均指数     
    
  ;;;;;输出原始EOF_ts的元数据
    alist                     =    [/eofts_hr(0,:),eofts_hr(1,:),eofts_hr(2,:),eofts_hr(3,:),eofts_hr(4,:),eofts_hr(5,:)/]     
    file_name                 =    datadir+"data/1-Interannual/2-eof_ts6.txt"                       ; 建立冬季降水的文件
    system("/bin/rm -f "+file_name)                                                                                        ; 移除存在的文件                                                          
    write_table(file_name,"w",alist,"%8.5f,%8.5f,%8.5f,%8.5f,%8.5f,%8.5f")  


   ;;;;;输出原始EPI数据
    alist                     =    [/EPI/]     
    file_name                 =    datadir+"data/1-Interannual/2-EPI.txt"                           ;建立冬季降水的文件
    system("/bin/rm -f "+file_name)                                                                 ;移除存在的文件                                                          
    write_table(file_name,"w",alist,"%8.5f")     

   ;;;;;输出标准化EPI数据;;;EPFI
    alist                     =    [/dim_wgt/]     
    file_name                 =    datadir+"data/1-Interannual/3-EPFI.txt"                           ;建立冬季降水的文件
    system("/bin/rm -f "+file_name)                                                                 ;移除存在的文件                                                          
    write_table(file_name,"w",alist,"%8.5f")     

   ;;;;;输出东南地区极端降水频次数据st_fq(年份，站点)    779站
    
    yer!0             =   "yer"
    yer@long_name     =   "years"
    yer@units         =   "y"    

    name_num!0        =   "name_num"
    name_num@long_name=   "779 station names"
    name_num@units    =   "s" 

    name_lat!0        =   "name_lat"
    name_lat@long_name=   "779 station lat"
    name_lat@units    =   "lat"     
    
    name_lon!0        =   "name_lon"
    name_lon@long_name=   "799 station lon"
    name_lon@units    =   "lon"     
    
    st_fq!0           =   "year"
    st_fq&year        =   yer     
    st_fq!1           =   "s_779"
    st_fq&s_779       =   name_num     
        
    
    st_fq@long_name   = "extreme precipitation frequency"
    st_fq@short_name  = "EPF"
    st_fq@units       = "frequency"
    st_fq@_FillValue  = -9999
    
    printVarSummary(st_fq)       
 
    file_name         =    datadir+"data/1-Interannual/3-EPF.nc"                                   ; 极端降水频次文件（年份，站点）
    system("/bin/rm -rf "+file_name)                                                                                   ; 移除存在的文件                
    
    fout              = addfile(file_name,"c")
    fout@title        = "NetCDF output EPF"
    fout@creation_date= systemfunc("date")
    filedimdef(fout,"year",-1,True)                                         
    fout->EPF         = st_fq
    fout->yr          = yer
    fout->n_num       = name_num
    fout->n_lat       = name_lat    
    fout->n_lon       = name_lon
    
  endtime            =    get_cpu_time()
  
  print("   runing time: " + (endtime-startime)/60.0+"minutes")    
  
  print("the number of dns is "+tostring(dns))

end