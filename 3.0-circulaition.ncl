;####1.合成图（极端降水多年-少年）

load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"                                             ;绘制南海和长江黄河
;<------------------------------------------------------------------------------------------------------------------------------->
begin

   startime            =    get_cpu_time()

   datadir                             =   "/home/jupiter/jiangsh/"
;   datadir                             =   "/home/tom/"

    syear                               =   1961                                                    ; 实际开始年份
    eyear                               =   2020                                                    ; 实际结束年份      

    it_ss                               =   196112                                                  ; 资料的起始月份
    it_ee                               =   202111
    yearse                              =   ispan(1961,2020,1)                                      ; 气候所有的年份
    sey                                 =   dimsizes(yearse)                        
        
    it_s                                =   1991                                                    ; 气候态起始年 ; 81/82 - 10/11气候态1981-2010
    it_e                                =   2020                                                    ; 气候态结束年

    levl                                =   (/1000,925,850,700,600,500,400,300,250,200,150,100,70/) ; 13层                    
    levlq                               =   (/1000,925,850,700,600,500,400,300/)                    ; 水汽层，所有的层次，共8层

    minlat                              =   -90                                                     ; 地图边界范围
    maxlat                              =   90
    minlon                              =   -180
    maxlon                              =   180
        
    ref8                                =   10                                                      ; 参考箭头所表示的风速大小：850hPa平均场
    ref8a                               =   5                                                       ; 参考箭头所表示的风速大小：850hPa异常场
    ref2a                               =   20                                                      ; 参考箭头所表示的风速大小：200hPa异常场
    reffx                               =   400                                                     ; 参考箭头所表示的风速大小：整层水汽通量参考箭头:平均场
    refplum                             =   15                                                      ; 波作用通量参数 
    
    wks_type                            =   "eps"
    test_sd                             =   0.1                                                     ; 90%的信度
;    wks_type@wkHeight=1080                                                                         ; 分辨率
;    wks_type@wkWidth=1920   

    sign                                =   95                                                      ; 显著性检验 95 significance

;0;plot 因为要循环画图，最好把不需要循环的变量（地图）提到循环外面来;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

;0.1****map****;;res地图设置;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    
  res                                   = True                                                      ; 调整地图和等值线, 每个子图均需该res  
;  res@gsnAddCyclic                      = False                                                    ; 添加循环点，否则会在0度经线左侧出现一根白条
  res@gsnMaximize                       = False                                                      ;--在工作台中尽可能大绘图（最大化绘图）       
  res@gsnDraw                           = False        
  res@gsnFrame                          = False        
  res@gsnLeftString                     = ""
  res@gsnRightString                    = ""
  res@tmXTOn                            = False                                                     ; 去除X轴上部的刻度
  res@tmYROn                            = False                                                     ; 去除Y轴右侧的刻度
  res@gsnLeftStringOrthogonalPosF       = 0.015
  res@gsnLeftStringFontHeightF          = 0.012                                                     ; 标题字体大小。后面不修改该值，则每幅图的主标题字体均是此大小
  res@gsnRightStringOrthogonalPosF      = 0.015
  res@gsnRightStringFontHeightF         = 0.012           
  res@tmBorderThicknessF                = 0.8  
  rest                                  = res                                                       ; 检验填色
  resd                                  = res                                                       ; 检验打点
  restc                                 = res                                                       ; 检验0.01线
  resv                                  = res                                                       ; 箭头
  resc                                  = res                                                       ; 等值线及填色
  resz                                  = res                                                       ; 青藏高原
  resf                                  = res                                                       ; 波作用通量
  resea                                 = res                                                       ; 海洋地图
  resea_nino                            = res
  res_p                                 = res
  res_xy                                = res
  res_eqt                               = res
  res@mpDataSetName                     = "Earth..4"                                                ; This new database contains  
  res@mpDataBaseVersion                 = "MediumRes"                                               ; Medium resolution database
  res@mpFillOn                          = False                                                     ; 不填色地图
  res@mpGeophysicalLineThicknessF       = 0.3                                                       ; 地图边界的粗细
  res@mpGeophysicalLineColor            = "grey60"                                                  ; 地图底图为浅色
  res@mpOutlineOn                       = True                                                      ; Turn on map outlines  
  res@mpOutlineSpecifiers               = (/"China","Taiwan"/) 
;  res@mpOutlineSpecifiers               = (/"China","Taiwan","Disputed area between India and China","India:Arunachal Pradesh"/)   
  res@mpNationalLineColor               = "grey60" 
  res@mpNationalLineThicknessF          = 0.3
  
  res@mpMinLatF                         = minlat
  res@mpMaxLatF                         = maxlat
  res@mpMinLonF                         = minlon
  res@mpMaxLonF                         = maxlon                       
;>----------------------------add China map----------------------------------------------------------------------------------------<
;  res@china                             = False                                                    ; draw china map or not
;  res@river                             = True                                                     ; draw changjiang&huanghe or not
;  res@province                          = False                                                    ; draw province boundary or not
;  res@nanhai                            = False                                                    ; draw nanhai or not
;  res@diqu                              = False                                                    ; draw diqujie or not
;>---------------------------------------------------------------------------------------------------------------------------------<  
; 海温地图坐标
;   resea@gsnAddCyclic                   = False                                                     ; 添加循环点，否则会在0度经线左侧出现一根白条
   resea@mpDataSetName                  = "Earth..4"                                                ; This new database contains  
   resea@mpDataBaseVersion              = "MediumRes"                                               ; Medium resolution database

   resea@mpCenterLonF                   = 160                                                       ; 地图的中心经度  
   resea@mpMinLatF                      = -20
   resea@mpMaxLatF                      = 70
   resea@mpMinLonF                      = 40
   resea@mpMaxLonF                      = 280  

   resea@tmYLMode                       = "Explicit"                                                ; 对Y轴的坐标进行设置,显示自己想要的值
   resea@tmYLValues                     = (/-20,0,20,40,60/)
   resea@tmYLLabels                     = (/"20~S~o~N~S","0","20~S~o~N~N","40~S~o~N~N","60~S~o~N~N"/)
   resea@tmXBMode                       = "Explicit"
   resea@tmXBValues                     = (/40,80,120,160,200,240,280/)
   resea@tmXBLabels                     = (/"40~S~o~N~E","80~S~o~N~E","120~S~o~N~E","160~S~o~N~E","160~S~o~N~W","120~S~o~N~W","80~S~o~N~W"/)   
      
;0.2****绘制青藏高原大于1500m地形****;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  resz@cnFillOn                         = True
  resz@cnLinesOn                        = False
  resz@cnLevelSelectionMode             = "ExplicitLevels"
  resz@cnLevels                         = (/2000,10000/)
  resz@cnLineColor                      = "red"                                                     ; 颜色
  resz@cnLineThicknessF                 = 1.0                                                       ; 粗细
  resz@cnLineLabelsOn                   = False                                                     ; 关闭标签
  resz@cnLineDashPattern                = 16                                                        ; 线型为16的虚线
  resz@cnInfoLabelOn                    = False                                                     ; 打开图右下方的等值线信息标签  
  resz@cnSmoothingOn                    = True                                                      ;--平滑线条
  resz@lbLabelBarOn                     = False                                                     ; 关闭labelbar
  resz@cnFillPalette                    = "GMT_gray"                                                ; 用GMT_gray 进行填色。即调用了第2种色板
  resz@cnFillColors                     = (/-1,11,-1/)                                              ; -1 为透明 

;   resz@gsLineColor                    = "grey55"
;   resz@cnLineThicknessF               = 8.0
  
;0.3.1****填色显著性0.01，0.05****;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;               
  rest@cnFillOn                         = True                                                      ; 填色等值线
  rest@cnFillMode                       = "RasterFill"                                              ; 平滑线
  rest@cnRasterSmoothingOn              = True
  rest@cnLinesOn                        = False                                                     ; 绘制等值线
  rest@cnLineLabelsOn                   = False                                                     ; 关闭标签 
  rest@cnLevelSelectionMode             = "ExplictLevels"
  rest@cnLevels                         = (/95,99/)
  rest@cnFillPalette                    = "GMT_gray"                                                ; 用GMT_gray 进行填色。即调用了第2种色板
  rest@cnFillColors                     = (/-1,7,5/)                                                ; -1 为透明 
  rest@cnInfoLabelOn                    = False                                                     ; 关闭图右下方的等值线信息标签
  rest@lbLabelBarOn                     = False                                                     ; 关闭labelbar
  print("1")
  
;0.3.2****打点显著性0.05****;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;               
  resd@cnLinesOn                        = False
  resd@cnLineLabelsOn                   = False                                                     ; do not use line labels             
  resd@cnInfoLabelOn                    = False                                                                                                        
  resd@cnFillOn                         = True                                                      ; dot fill                                 
  resd@lbLabelBarOn                     = False                                                                                                            
  resd@cnLevelSelectionMode             = "ExplicitLevels"                                                                                             
  resd@cnFillDrawOrder                  = "postDraw"                                                                                                     
  resd@cnLevels                         = (/95,99/)                                                                                                
  resd@cnMonoFillPattern                = False                                                     ; want multiple patterns                                                               
  resd@cnFillPatterns                   = (/-1,17,17/)                                                ; the patterns                                                                        
  resd@cnMonoFillScale                  = False                                                     ; want different densities                                                            
  resd@cnFillScales                     = (/0.7,0.7,0.7/)                                           ; change densities                                                                        
  resd@cnMonoFillColor                  = True                                                                                                     
  resd@cnFillDotSizeF                   = -0.1                                                                                                                                 

;0.3.3****画显著线0.01****绿色;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;               
  restc@cnLinesOn                       = True
  restc@cnFillOn                        = False                  
  restc@cnSmoothingOn                   = True                                                      ; 平滑线
  restc@cnLineLabelsOn                  = False                                                     ; do not use line labels 关闭标签            
  restc@cnInfoLabelOn                   = False                                                                                                        
  restc@cnLineThicknessF                = 1.0                         
  restc@lbLabelBarOn                    = False                                                                                                            
  restc@cnLevelSelectionMode            = "ExplicitLevels"                                                                                             
  restc@cnLineColor                     = "green"                                                                                               
  restc@cnLevels                        = sign                                                        ; 0.1显著性画线                                                                                                
  restc@cnFillDrawOrder                 = "postDraw"       
  print("2")
;0.4****合成暴雨异常场等高线****;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;              

  resc@cnFillOn                         = False                                                     ; 关闭等值线填色 
  resc@cnLineThicknessF                 = 0.7                                                       ; 等值线粗细  
  resc@gsnContourZeroLineThicknessF     = 0.                                                        ; 设置0值线粗细。0则不画
  resc@cnLineLabelsOn                   = False                                                     ; 关闭标签
  resc@cnInfoLabelOn                    = False                                                     ; 打开图右下方的等值线信息标签
  resc@lbLabelBarOn                     = True 
  resc@cnSmoothingOn                    = True                                                      ; 平滑线
  resc@cnSmoothingDistanceF             = 0.01 
  resc@cnSmoothingTensionF              = 0.1
  resc@cnLineLabelBackgroundColor       = -1
  resc@cnLineLabelPlacementMode         = "Constant"
  resc@cnLineLabelFontHeightF           = 0.008
  resc@cnLineLabelInterval              = 1
  resc@lbLabelFontHeightF               = 0.012                                                     ; 色标标签字体大小
  resc@pmLabelBarOrthogonalPosF         = 0.08
  resc@pmLabelBarHeightF                = 0.09                                                      ; labelbar高度
    
;  resc@cnInfoLabelOrthogonalPosF = 0.05                                                            ; 移动等值线信息标签的位置
  
;0.5.1****合成异常场箭头****;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
;  resv@vcPositionMode                  = "ArrowTail”"                                              ; 箭头尾部对应着格点的位置
  resv@vcGlyphStyle                     = "LineArrow"                                               ; 其余三种选项为“LineArrow”（默认值）、“WindBarb” 、“CurlyVector”
  resv@vcMinDistanceF                   = 0.03                                                      ; 箭头之间的最小距离（在单位平方中）
;--------------------------------------------------------------------------------------------------
  resv@vcFillArrowEdgeThicknessF        = 1                                                         ; 箭头边界粗细
  resv@vcFillArrowEdgeColor             = "white"                                                   ; 及颜色
  resv@vcFillArrowFillColor             = "darkslateblue"                                           ; 箭头内部填充颜色
  resv@vcFillArrowWidthF                = 0.07                                                      ; 箭头宽度
  resv@vcFillArrowHeadXF                = 0.4                                                       ; 请参考附录中Fillarrow箭头示意图
  resv@vcFillArrowHeadYF                = 0.10                                                      ;
  resv@vcFillArrowHeadInteriorXF        = 0.20                                                      ;            
  
  resv@vcFillArrowMinFracWidthF         = 1.0 
  resv@vcFillArrowHeadMinFracXF         = 1.0  
  resv@vcFillArrowHeadMinFracYF         = 1.0 
  
;0.5.2****设定参考箭头****
   resv@vcRefAnnoOn                     =  True  
   resv@vcRefLengthF                    =  0.040                                                    ; 标准长度箭头在单位平方中的大小
   resv@vcRefAnnoBackgroundColor        =  -1                                                       ; 背景颜色，透明     
   resv@vcRefAnnoPerimOn                =  False                                                    ; 关闭边框                                           
   resv@vcRefAnnoFontHeightF            =  0.012                                                    ; 参考箭头标签字体大小         
   resv@vcRefAnnoString1On              =  False                                                    ; 设定参考箭头上的字符        
   resv@vcRefAnnoString2On              =  True                                                     ; 这里仅设定其下方的字符           
   resv@vcRefAnnoSide                   =  "Top"                                                    ; 参考箭头放至图形上方
   resv@vcRefAnnoOrthogonalPosF         =  -0.040                                                   ; 调整箭头垂直位置(>0为向上)
   resv@vcRefAnnoParallelPosF           =  1.00                                                     ; 调整箭头水平位置
   resv@vcRefAnnoArrowSpaceF            =  1.3                                                      ; 调整箭头和文本间的距离
 
;########################################### 一.气象和气候数据分析###########################################################################################################################################            

;1.1 eof1_ts 为 eof 分解后第一模态时间序列的元数据--------------------------------------------------------------------------------------

;    f_eofts                             =   asciiread(datadir+"data/1-Interannual/1-eof_ts6.txt",-1,"string")   
;    delim                               =   ","    
;    eof1_ts                             =   stringtofloat(str_get_field(f_eofts,1,delim))           ; 提取EOF1的PC1
;    print(eof1_ts)

    f_epi                               =   asciiread(datadir+"data/1-Interannual/2-EPI.txt",-1,"string")   
    delim                               =   ","    
    eof1_ts                             =   stringtofloat(str_get_field(f_epi,1,delim))           ; 未经过任何处理的极端降水频次面积平均指数
    print(eof1_ts)
    
    dp_epi                              =   dim_rmvmean_n_Wrap(eof1_ts,0)     ; 去趋势

    print(dp_epi)
    

    st_pc1                              =    dim_standardize_n(eof1_ts,1,0)                         ; 标准化    
    st_pc1!0                            =   "year"    
    st_pc1&year                         =   yearse    

    std_flag                            =   1.0                                                     ; 正负1.0个标准差  
    ind_stdy1                           =   ind(st_pc1.ge.std_flag)                                 ; 1982,1984,1988,1989,1994,1997,2002,2015,2019
    ind_stdy2                           =   ind(st_pc1.lt.-std_flag)                                ; 1962,1967,1969,1976,1983,1985,1995,1998,2008 
    
    stdy1                               =  yearse(ind_stdy1)
    stdy2                               =  yearse(ind_stdy2)
    print(stdy1)
    print(stdy2)

;;;;;
    std_epf                             =   dim_standardize_n(eof1_ts,1,0)                         ;去趋势后，再标准化   
    ind_epf1                            =   ind(std_epf.ge.std_flag)                                 ; 1968,1982,1988,1991,1994,1997,2013,2015,2017,2019
    ind_epf2                            =   ind(std_epf.lt.-std_flag)                                ; 1959,1962,1965,1967,1969,1973,1976,1983,1986,1995,1998,2008 
    
    year_epf1                               =  yearse(ind_epf1)
    year_epf2                               =  yearse(ind_epf2)
    print(year_epf1)
    print(year_epf2)

;;;计算结果表明：很巧合，去趋势与不去趋势的1Std年份是一样的。
;1.2; read 青藏高原地形数据----------------------------------------------------------------------------------------------------------

    fz            =   addfile(datadir+"data/5-sfc-land/hgt.sfc.nc","r")
    z1            =   fz->hgt(0,:,:)
    z             =   lonFlip(z1)                                                                   ; 将0~360度经度转化为-180~180度经度
    delete(fz)
    delete(z1)
print("z")                  
       

;1.4;;;;;;; read climate winter ncepdata;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;数据长度：194801-202111
   f_hgt            =  addfile(datadir+"data/3-ncep/monthly/hgt.mon.mean.nc","r")                     ; float hgt ( time, level, lat, lon )
   f_uwd            =  addfile(datadir+"data/3-ncep/monthly/uwnd.mon.mean.nc","r") 
   f_vwd            =  addfile(datadir+"data/3-ncep/monthly/vwnd.mon.mean.nc","r")
   f_omega          =  addfile(datadir+"data/3-ncep/monthly/omega.mon.mean.nc","r")                   ; float omega ( time, level, lat, lon )
   f_q              =  addfile(datadir+"data/3-ncep/monthly/shum.mon.mean.nc","r")                    ; float shum ( time, level, lat, lon ) grams/kg level(8层，1000,925,850,700,600,500,400,300)
   f_pres           =  addfile(datadir+"data/3-ncep/monthly/pres.mon.mean.nc","r")                    ; float pres ( time, lat, lon )  units :millibars
   f_air            =  addfile(datadir+"data/3-ncep/monthly/air.mon.mean.nc","r")                     ; 多气压层的空气温度
   f_rh             =  addfile(datadir+"data/3-ncep/monthly/rhum.mon.mean.nc","r")                    ; 相对湿度

   time             =  f_hgt->time                                                                  ; 读取气候态日期
   YYYYMM           =  cd_calendar(time,-1)                                                         ; 转换成公历日期
   rec_s            =  ind(it_ss.eq.YYYYMM)
   rec_e            =  ind(it_ee.eq.YYYYMM)
   
   hgt_ca1          =  f_hgt->hgt(rec_s:rec_e,{levl},:,:)                                           ; 截取指定时间段196112-202111
   uwd_ca1          =  f_uwd->uwnd(rec_s:rec_e,{levl},:,:)
   vwd_ca1          =  f_vwd->vwnd(rec_s:rec_e,{levl},:,:)
   pres_ca1         =  f_pres->pres(rec_s:rec_e,:,:)                                                ; 地面气压   pres ( time, lat, lon )            
   air_ca1          =  f_air->air(rec_s:rec_e,{levl},:,:)                                           ; 大气温度   air(time,lvel,lat,lon) 
   omega_ca1        =  f_omega->omega(rec_s:rec_e,:,:,:)
   omega_ca1        =  omega_ca1*100
   
   rh_ca1           =  f_rh->rhum(rec_s:rec_e,{levlq},:,:)                                          ; 相对湿度 rhum(time,lvel,lat,lon)levl层和比湿相同               
   q_ca3            =  f_q->shum(rec_s:rec_e,{levlq},:,:)                                           ; 比湿 
   q_ca1            =  q_ca3/1000.0                                                                 ; 转化为kg/kg
   copy_VarMeta(q_ca3,q_ca1)
         
   hgt_ca1         :=  lonFlip(hgt_ca1)                                                             ; 将经度0~360度转换为-180~180度
   uwd_ca1         :=  lonFlip(uwd_ca1)
   vwd_ca1         :=  lonFlip(vwd_ca1)  
   omega_ca1       :=  lonFlip(omega_ca1)
   q_ca1           :=  lonFlip(q_ca1)
   pres_ca1        :=  lonFlip(pres_ca1)
   air_ca1         :=  lonFlip(air_ca1)   
   rh_ca1          :=  lonFlip(rh_ca1)      
   
   hgt_ca           =  hgt_ca1(:,:,{minlat:maxlat},{minlon:maxlon})
   uwd_ca           =  uwd_ca1(:,:,{minlat:maxlat},{minlon:maxlon})
   vwd_ca           =  vwd_ca1(:,:,{minlat:maxlat},{minlon:maxlon})
   omega_ca         =  omega_ca1(:,:,{minlat:maxlat},{minlon:maxlon})
   q_ca             =  q_ca1(:,:,{minlat:maxlat},{minlon:maxlon})
   pres_ca          =  pres_ca1(:,{minlat:maxlat},{minlon:maxlon})
   air_ca           =  air_ca1(:,:,{minlat:maxlat},{minlon:maxlon})
   air_ca           =  air_ca+273.15
   rh_ca            =  rh_ca1(:,:,{minlat:maxlat},{minlon:maxlon})
;--------------------------------------------------------------------------------------------------------------------------------------------------

;;1.4.1 流函数和势函数;;;;;;;; 
   sfvp_ca          =  uv2sfvpF(uwd_ca,vwd_ca)                                                      ; 计算流函数和势函数，sfvp比原来多一维数据。
   sf_ca            =  sfvp_ca(0,:,:,:,:)                                                           ; 第一维为流函数，没有坐标
   vp_ca            =  sfvp_ca(1,:,:,:,:)
   
   copy_VarCoords(uwd_ca,sf_ca)
   copy_VarCoords(uwd_ca,vp_ca)

;;1.4.2 散度场;;;;;;;;  
   div_ca           =  uv2dv_cfd(uwd_ca,vwd_ca,uwd_ca&lat,uwd_ca&lon,2) 
   div_ca           =  div_ca*1000000                                                                ; 扩大10的6次方倍                        
   copy_VarCoords(uwd_ca,div_ca)     

;;1.4.3 涡度场;;;;;;;; 
   vor_ca           =  uv2vr_cfd(uwd_ca,vwd_ca,uwd_ca&lat,uwd_ca&lon,2) 
   vor_ca           =  vor_ca*1000000                                                               ; 扩大10的6次方倍                        
   copy_VarCoords(uwd_ca,vor_ca)

;;1.4.4 water vapor flux;;;;;;;; 
  
   qu1              =  new(dimsizes(q_ca1),"float")
   qu1              =  q_ca1*uwd_ca1(:,{levlq},:,:)
   qv1              =  q_ca1*vwd_ca1(:,{levlq},:,:)
   copy_VarCoords(q_ca1,qu1)
   copy_VarCoords(q_ca1,qv1)
   printVarSummary(qu1)

   ;;水汽500hPa层数太少，导致积分错误，所以多设置几层。
   
   qu_new1   =  new(dimsizes(uwd_ca1),"float")
   qv_new1   =  new(dimsizes(uwd_ca1),"float")

   do      i = 0, dimsizes(levlq)-1
     qu_new1(:,i,:,:)= qu1(:,i,:,:)
     qv_new1(:,i,:,:)= qv1(:,i,:,:)
   end do
         
   do     j  = dimsizes(levlq),dimsizes(levl)-1
     qu_new1(:,j,:,:)= 0.0
     qv_new1(:,j,:,:)= 0.0
   end do 
    
   qu_new   = qu_new1(:,:,{minlat:maxlat},{minlon:maxlon})   
   qv_new   = qv_new1(:,:,{minlat:maxlat},{minlon:maxlon})
   
   p           :=  levl                                                                                         ; 气压的层
   ps           =  pres_ca                                                                                      ; 地面气压hPa（逐月资料为hPa，逐日为pa需要转化为hPa）
   linlog       =  1
   pbot         =  1100
   ptop         =  300      
   g            =  9.80665
   
   vint_qu_ca   =  vibeta(p,qu_new(time|:,lat|:,lon|:,level|:),linlog,ps,pbot,ptop)/g*100                       ; 整层水汽通量
   copy_VarCoords(uwd_ca(:,0,:,:),vint_qu_ca)                                                                   ; 将风速u场坐标复制给u场
   printVarSummary(vint_qu_ca)
   
   print("maxvint="+max(vint_qu_ca)) 
   
   vint_qv_ca   =  vibeta(p,qv_new(time|:,lat|:,lon|:,level|:),linlog,ps,pbot,ptop)/g*100
   copy_VarCoords(vwd_ca(:,0,:,:),vint_qv_ca)
   
;;1.4.5 整层水汽通量散度;;;;;;;;  
   div_quv_ca  =  uv2dv_cfd(vint_qu_ca,vint_qv_ca,uwd_ca&lat,uwd_ca&lon,2) 
   div_quv_ca  =  div_quv_ca*100000                                                                                   ; 扩大10的5次方倍                        
   copy_VarCoords(vint_qu_ca,div_quv_ca)     
      
;;1.4.6 theta_se;;;;;;;;;假相当位温    
   
   levv4d     :=  conform(air_ca(:,{levlq},:,:),levlq,1)
   copy_VarCoords(air_ca(:,{levlq},:,:),levv4d)
  
  es    := where((air_ca(:,{levlq},:,:)-273.15).lt.0,6.112*exp(17.67*(air_ca(:,{levlq},:,:)-273.15)/(air_ca(:,{levlq},:,:)-29.65)),6.11*10^(7.5*(air_ca(:,{levlq},:,:)-273.15)/( air_ca(:,{levlq},:,:)-35.85)))                                                            ; 注意t温度为绝对温度 
  qq    := rh_ca * (0.62197*es/( levv4d-0.378*es))/100.
  e     := levv4d * qq/(0.62197+qq)+10^(-10)
  tlcl  := 55.0+2840.0/(3.5*log(air_ca(:,{levlq},:,:))-log(e)-4.805)
  theta := air_ca(:,{levlq},:,:)*(1000. / levv4d )^(0.2854*(1.0-0.28*qq))
  eqt_ca:= theta*exp(((3376./tlcl)-2.54)*qq*(1.0+0.81*qq))                                           ; 假相当位温
  
  copy_VarCoords(air_ca(:,{levlq},:,:),eqt_ca)   
  
;------------------------------------------------------------------------------------------------------    
  delete([/hgt_ca1,uwd_ca1,vwd_ca1,omega_ca1,q_ca1,pres_ca1,air_ca1,rh_ca1,qu_new1,qv_new1/])
   
print("c")
; 1.5 -------------------------------------------climatic season 冬季气候平均和及异常计算---------------------------------------------------------------------< 
;;1.5.1air;;;;;;;;;;
   air_DJF          =  month_to_season(air_ca,"JFM")                  ;每年冬季;;DJF 平均;JFM季节平均，实际是12/1/2月三个月平均，因为从1981年12月开始截取;;;;;;;  
   copy_VarMeta(air_ca(0,:,:,:),air_DJF(0,:,:,:))
   air_DJF!0        =  "year"
   air_DJF&year     =  yearse
   air_c            =  dim_avg_n_Wrap(air_DJF({it_s:it_e},:,:,:),0)    ;30年冬季平均(1981-2010) float air_c (level, lat, lon )  [level | 3] x [lat | 33] x [lon | 144]
   printVarSummary(air_c)

;;1.5.2pres;;;;;;;;;;地面气压pres ( time, lat, lon ) 
   pres_DJF         =  month_to_season(pres_ca,"JFM")                 ; 每年冬季;;DJF 平均;JFM季节平均，实际是12/1/2月三个月平均，因为从1951年12月开始截取;;;;;;;  
   copy_VarMeta(pres_ca(0,:,:),pres_DJF(0,:,:))
   pres_DJF!0       =  "year"
   pres_DJF&year    =  yearse                                         ; 1951-2016年的冬季
   pres_c           =  dim_avg_n_Wrap(pres_DJF({it_s:it_e},:,:),0)    ; 30年冬季平均 float pres_c (lat, lon )  [lat | 33] x [lon | 144]
   printVarSummary(pres_c)
    
;;1.5.3hgt;;;;;;;;;;
   hgt_DJF          =  month_to_season(hgt_ca,"JFM")                  ;每年冬季;;DJF 平均;JFM季节平均，实际是12/1/2月三个月平均，因为从1981年12月开始截取;;;;;;;  
   copy_VarMeta(hgt_ca(0,:,:,:),hgt_DJF(0,:,:,:))
   hgt_DJF!0        =  "year"
   hgt_DJF&year     =  yearse
   hgt_c            =  dim_avg_n_Wrap(hgt_DJF({it_s:it_e},:,:,:),0)   ;30年冬季平均 float hgt_c (level, lat, lon )  [level | 3] x [lat | 33] x [lon | 144]
   printVarSummary(hgt_c)
   
   hgt_ano          =  hgt_DJF
   hgt_ano          =  hgt_DJF-conform(hgt_DJF,hgt_c,(/1,2,3/))
   hgt_ano_d        =  dtrend_msg_n(yearse,hgt_ano,True,True,0)                                      ; 去趋势
   copy_VarCoords(hgt_DJF,hgt_ano_d)

   printVarSummary(hgt_ano_d)
   
   hgt_ave1         =  dim_avg_n_Wrap(hgt_DJF({stdy1},:,:,:),0)                                     ; 多暴雨年的均值   
   hgt_ave2         =  dim_avg_n_Wrap(hgt_DJF({stdy2},:,:,:),0)                                     ; 少暴雨年的均值
   var_hgt1         =  dim_variance_n_Wrap(hgt_DJF({stdy1},:,:,:),0)                                ; 多暴雨年的方差
   var_hgt2         =  dim_variance_n_Wrap(hgt_DJF({stdy2},:,:,:),0)                                ; 少暴雨年的方差
   
   hgt_ano12        =  hgt_ave1
   hgt_ano12        =  hgt_ave1-hgt_ave2                                                            ; 多暴雨年均值-少暴雨年均值

   hgt_ano1         =  dim_avg_n_Wrap(hgt_ano({stdy1},:,:,:),0)                                     ; 多暴雨年的异常平均
   hgt_ano2         =  dim_avg_n_Wrap(hgt_ano({stdy2},:,:,:),0)                                     ; 少暴雨年的异常平均
   
;;1.5.4uwd;;;;;;;;
   uwd_DJF          =  month_to_season(uwd_ca,"JFM") 
   copy_VarMeta(uwd_ca(0,:,:,:),uwd_DJF(0,:,:,:))
   uwd_DJF!0        =  "year"
   uwd_DJF&year     =  yearse
   uwd_c            =  dim_avg_n_Wrap(uwd_DJF({it_s:it_e},:,:,:),0)                                 ; 30年冬季平均
   
   uwd_ano          =  uwd_DJF
   uwd_ano          =  uwd_DJF-conform(uwd_DJF,uwd_c,(/1,2,3/))
   uwd_ano_d        =  dtrend_msg_n(yearse,uwd_ano,True,True,0)                                     ; 去趋势
   copy_VarCoords(uwd_DJF,uwd_ano_d)
        
   uwd_ave1         =  dim_avg_n_Wrap(uwd_DJF({stdy1},:,:,:),0)                                     ; 多暴雨年的均值
   uwd_ave2         =  dim_avg_n_Wrap(uwd_DJF({stdy2},:,:,:),0)
   var_uwd1         =  dim_variance_n_Wrap(uwd_DJF({stdy1},:,:,:),0)                                ; 多暴雨年的方差
   var_uwd2         =  dim_variance_n_Wrap(uwd_DJF({stdy2},:,:,:),0)                                ; 少暴雨年的方差   
      
   uwd_ano12        =  uwd_ave1
   uwd_ano12        =  uwd_ave1-uwd_ave2                                                            ; 正负1个标准差相减
   printVarSummary(uwd_ano12)
   
   uwd_ano1         =  dim_avg_n_Wrap(uwd_ano({stdy1},:,:,:),0)                                     ; 多暴雨年的异常平均
   uwd_ano2         =  dim_avg_n_Wrap(uwd_ano({stdy2},:,:,:),0)                                     ; 少暴雨年的异常平均
   
;;1.5.4vwd;;;;;;;;
   vwd_DJF          =  month_to_season(vwd_ca,"JFM")                                                ; 每年冬季
   copy_VarMeta(vwd_ca(0,:,:,:),vwd_DJF(0,:,:,:))
   vwd_DJF!0        =  "year"
   vwd_DJF&year     =  yearse
   vwd_c            =  dim_avg_n_Wrap(vwd_DJF({it_s:it_e},:,:,:),0)                                 ; 30年冬季平均
   
   vwd_ano          =  vwd_DJF
   vwd_ano          =  vwd_DJF-conform(vwd_DJF,vwd_c,(/1,2,3/))
   vwd_ano_d        =  dtrend_msg_n(yearse,vwd_ano,True,True,0)                                     ; 去趋势
   copy_VarCoords(vwd_DJF,vwd_ano_d)  
   
   vwd_ave1         =  dim_avg_n_Wrap(vwd_DJF({stdy1},:,:,:),0)                                     ; 正1个标准差平均
   vwd_ave2         =  dim_avg_n_Wrap(vwd_DJF({stdy2},:,:,:),0)
   var_vwd1         =  dim_variance_n_Wrap(vwd_DJF({stdy1},:,:,:),0)                                ; 多暴雨年的方差
   var_vwd2         =  dim_variance_n_Wrap(vwd_DJF({stdy2},:,:,:),0)                                ; 少暴雨年的方差  
      
   vwd_ano12        =  vwd_ave1   
   vwd_ano12        =  vwd_ave1-vwd_ave2                                                            ; 正负1个标准差相减

   vwd_ano1         =  dim_avg_n_Wrap(vwd_ano({stdy1},:,:,:),0)                                     ; 多暴雨年的异常平均
   vwd_ano2         =  dim_avg_n_Wrap(vwd_ano({stdy2},:,:,:),0)                                     ; 少暴雨年的异常平均
;;1.5.5omega;;;;;;;;
   omega_DJF        =  month_to_season(omega_ca,"JFM")                        
   copy_VarMeta(omega_ca(0,:,:,:),omega_DJF(0,:,:,:))
   omega_DJF!0      =  "year"
   omega_DJF&year   =  yearse
   omega_c          =  dim_avg_n_Wrap(omega_DJF({it_s:it_e},:,:,:),0)        
   printVarSummary(omega_c)
                                         
   omega_ano        =  omega_DJF
   omega_ano        =  omega_DJF-conform(omega_DJF,omega_c,(/1,2,3/))
   omega_ano_d      =  dtrend_msg_n(yearse,omega_ano,True,True,0)                                   ; 去趋势
   copy_VarCoords(omega_DJF,omega_ano_d)

   
   omega_ave1       =  dim_avg_n_Wrap(omega_DJF({stdy1},:,:,:),0)                                   ; 正1个标准差平均
   omega_ave2       =  dim_avg_n_Wrap(omega_DJF({stdy2},:,:,:),0)
   var_omega1       =  dim_variance_n_Wrap(omega_DJF({stdy1},:,:,:),0)                              ; 多暴雨年的方差
   var_omega2       =  dim_variance_n_Wrap(omega_DJF({stdy2},:,:,:),0)                              ; 少暴雨年的方差       
   
   omega_ano12      =  omega_ave1
   omega_ano12      =  omega_ave1-omega_ave2                                                        ; 正负1个标准差相减

   omega_ano1       =  dim_avg_n_Wrap(omega_ano({stdy1},:,:,:),0)                                   ; 多暴雨年的异常平均
   omega_ano2       =  dim_avg_n_Wrap(omega_ano({stdy2},:,:,:),0)                                   ; 少暴雨年的异常平均
   
;;1.5.6theta_se;;;;;;;;eqt   

   eqt_DJF          =  month_to_season(eqt_ca,"JFM") 
   copy_VarMeta(eqt_ca(0,:,:,:),eqt_DJF(0,:,:,:))
   eqt_DJF!0        =  "year"
   eqt_DJF&year     =  yearse
   eqt_c            =  dim_avg_n_Wrap(eqt_DJF({it_s:it_e},:,:,:),0)
   
   eqt_ano          =  eqt_DJF
   eqt_ano          =  eqt_DJF-conform(eqt_DJF,eqt_c,(/1,2,3/))
   eqt_ano_d        =  dtrend_msg_n(yearse,eqt_ano,True,True,0)                                     ; 去趋势
   copy_VarCoords(eqt_DJF,eqt_ano_d)
        
   eqt_ave1         =  dim_avg_n_Wrap(eqt_DJF({stdy1},:,:,:),0)                                     ; 正1个标准差平均
   eqt_ave2         =  dim_avg_n_Wrap(eqt_DJF({stdy2},:,:,:),0)
   var_eqt1         =  dim_variance_n_Wrap(eqt_DJF({stdy1},:,:,:),0)                                ; 多暴雨年的方差
   var_eqt2         =  dim_variance_n_Wrap(eqt_DJF({stdy2},:,:,:),0)                                ; 少暴雨年的方差  
  
   eqt_ano12        =  eqt_ave1
   eqt_ano12        =  eqt_ave1-eqt_ave2                                                            ; 正负1个标准差相减
   printVarSummary(eqt_ano12)
   eqt_ano1         =  dim_avg_n_Wrap(eqt_ano({stdy1},:,:,:),0)                                     ; 多暴雨年的异常平均
   eqt_ano2         =  dim_avg_n_Wrap(eqt_ano({stdy2},:,:,:),0)                                     ; 少暴雨年的异常平均 

;;1.5.7water vapor flux整层水汽通量;;;;;;;; 
; qu
   vint_qu_DJF      =  month_to_season(vint_qu_ca,"JFM")                                            ; 每年冬季
   copy_VarMeta(vint_qu_ca(0,:,:),vint_qu_DJF(0,:,:))
   vint_qu_DJF!0    =  "year"
   vint_qu_DJF&year =  yearse
   vint_qu_c        =  dim_avg_n_Wrap(vint_qu_DJF({it_s:it_e},:,:),0)                               ; 30年冬季平均
   printVarSummary(vint_qu_c)
   
   vint_qu_ano      =  vint_qu_DJF
   vint_qu_ano      =  vint_qu_DJF-conform(vint_qu_DJF,vint_qu_c,(/1,2/)) 
   vint_qu_ano_d    =  dtrend_msg_n(yearse,vint_qu_ano,True,True,0)   
   copy_VarCoords(vint_qu_DJF,vint_qu_ano_d)  
; qv
   vint_qv_DJF      =  month_to_season(vint_qv_ca,"JFM")                                          
   copy_VarMeta(vint_qv_ca(0,:,:),vint_qv_DJF(0,:,:))
   vint_qv_DJF!0    =  "year"
   vint_qv_DJF&year =  yearse
   vint_qv_c        =  dim_avg_n_Wrap(vint_qv_DJF({it_s:it_e},:,:),0)                               ; 30年冬季平均
   printVarSummary(vint_qv_c)
   print("maxvintqv_c="+max(vint_qv_c)) 
   
   vint_qv_ano      =  vint_qv_DJF
   vint_qv_ano      =  vint_qv_DJF-conform(vint_qv_DJF,vint_qv_c,(/1,2/)) 
   vint_qv_ano_d    =  dtrend_msg_n(yearse,vint_qv_ano,True,True,0)                                 ; 去趋势
   copy_VarCoords(vint_qv_DJF,vint_qv_ano_d) 
      
;   wvflux_x     =  sqrt(vint_qu_c^2+vint_qv_c^2)
;   copy_VarMeta(vint_qu_c,wvflux_c)   

   vint_qu_ave1     =  dim_avg_n_Wrap(vint_qu_DJF({stdy1},:,:),0)                                   ; 多暴雨年标准差平均
   vint_qu_ave2     =  dim_avg_n_Wrap(vint_qu_DJF({stdy2},:,:),0)
   var_vint_qu1     =  dim_variance_n_Wrap(vint_qu_DJF({stdy1},:,:),0)                              ; 多暴雨年的方差
   var_vint_qu2     =  dim_variance_n_Wrap(vint_qu_DJF({stdy2},:,:),0)                              ; 少暴雨年的方差  
      
   vint_qu_ano12    =  vint_qu_ave1
   vint_qu_ano12    =  vint_qu_ave1-vint_qu_ave2       

   vint_qv_ave1     =  dim_avg_n_Wrap(vint_qv_DJF({stdy1},:,:),0)                                   ; 多暴雨年标准差平均
   vint_qv_ave2     =  dim_avg_n_Wrap(vint_qv_DJF({stdy2},:,:),0)
   var_vint_qv1     =  dim_variance_n_Wrap(vint_qv_DJF({stdy1},:,:),0)                              ; 多暴雨年的方差
   var_vint_qv2     =  dim_variance_n_Wrap(vint_qv_DJF({stdy2},:,:),0)                              ; 少暴雨年的方差  
   
   vint_qv_ano12    =  vint_qv_ave1
   vint_qv_ano12    =  vint_qv_ave1-vint_qv_ave2       

   vint_qu_ano1         =  dim_avg_n_Wrap(vint_qu_ano({stdy1},:,:),0)                               ; 多暴雨年的异常平均
   vint_qu_ano2         =  dim_avg_n_Wrap(vint_qu_ano({stdy2},:,:),0)                               ; 少暴雨年的异常平均
   vint_qv_ano1         =  dim_avg_n_Wrap(vint_qv_ano({stdy1},:,:),0)                               ; 多暴雨年的异常平均
   vint_qv_ano2         =  dim_avg_n_Wrap(vint_qv_ano({stdy2},:,:),0)                               ; 少暴雨年的异常平均   

;;1.5.8整层水汽通量散度

   div_quv_DJF      =  month_to_season(div_quv_ca,"JFM")                                            ; 每年冬季
   copy_VarMeta(div_quv_ca(0,:,:),div_quv_DJF(0,:,:))
   div_quv_DJF!0    =  "year"
   div_quv_DJF&year =  yearse
   div_quv_c        =  dim_avg_n_Wrap(div_quv_DJF({it_s:it_e},:,:),0)                               ; 30年冬季平均
   printVarSummary(div_quv_c)
   
   div_quv_ano      =  div_quv_DJF
   div_quv_ano      =  div_quv_DJF-conform(div_quv_DJF,div_quv_c,(/1,2/)) 
   div_quv_ano_d    =  dtrend_msg_n(yearse,div_quv_ano,True,True,0)   
   copy_VarCoords(div_quv_ano,div_quv_ano_d)  

   div_quv_ave1     =  dim_avg_n_Wrap(div_quv_DJF({stdy1},:,:),0)                                   ; 正0.6个标准差平均
   div_quv_ave2     =  dim_avg_n_Wrap(div_quv_DJF({stdy2},:,:),0)
   var_div_quv1     =  dim_variance_n_Wrap(vint_qv_DJF({stdy1},:,:),0)                              ; 多暴雨年的方差
   var_div_quv2     =  dim_variance_n_Wrap(vint_qv_DJF({stdy2},:,:),0)                              ; 少暴雨年的方差 

   
   div_quv_ano12    =  div_quv_ave1
   div_quv_ano12    =  div_quv_ave1-div_quv_ave2       
   
   div_quv_ano1         =  dim_avg_n_Wrap(div_quv_ano({stdy1},:,:),0)                               ; 多暴雨年的异常平均
   div_quv_ano2         =  dim_avg_n_Wrap(div_quv_ano({stdy2},:,:),0)                               ; 少暴雨年的异常平均

;;1.5.9div;;;;;;;;;;   
   div_DJF          =  month_to_season(div_ca,"JFM") 
   copy_VarMeta(div_ca(0,:,:,:),div_DJF(0,:,:,:))
   div_DJF!0        =  "year"
   div_DJF&year     =  yearse
   div_c            =  dim_avg_n_Wrap(div_DJF({it_s:it_e},:,:,:),0)               
   
   div_ano          =  div_DJF
   div_ano          =  div_DJF-conform(div_DJF,div_c,(/1,2,3/))
   div_ano_d        =  dtrend_msg_n(yearse,div_ano,True,True,0)                                     ; 去趋势
   copy_VarCoords(div_DJF,div_ano_d)
        
   div_ave1         =  dim_avg_n_Wrap(div_DJF({stdy1},:,:,:),0)                                     ; 正1个标准差平均
   div_ave2         =  dim_avg_n_Wrap(div_DJF({stdy2},:,:,:),0)
   var_div1         =  dim_variance_n_Wrap(div_DJF({stdy1},:,:,:),0)                                ; 多暴雨年的方差
   var_div2         =  dim_variance_n_Wrap(div_DJF({stdy2},:,:,:),0)                                ; 少暴雨年的方差 
      
   div_ano12        =  div_ave1
   div_ano12        =  div_ave1-div_ave2                                                            ; 正负1个标准差相减
   printVarSummary(div_ano12)
   
   div_ano1         =  dim_avg_n_Wrap(div_ano({stdy1},:,:,:),0)                                     ; 多暴雨年的异常平均
   div_ano2         =  dim_avg_n_Wrap(div_ano({stdy2},:,:,:),0)                                     ; 少暴雨年的异常平均
;;1.5.10vor;;;;;;;;;;
   vor_DJF          =  month_to_season(vor_ca,"JFM")                                                ; 每年冬季;;DJF 平均;JFM季节平均，实际是12/1/2月三个月平均，因为从1981年12月开始截取;;;;;;;  
   copy_VarMeta(vor_ca(0,:,:,:),vor_DJF(0,:,:,:))
   vor_DJF!0        =  "year"
   vor_DJF&year     =  yearse
   vor_c            =  dim_avg_n_Wrap(vor_DJF({it_s:it_e},:,:,:),0)           
   printVarSummary(vor_c)
   
   vor_ano          =  vor_DJF
   vor_ano          =  vor_DJF-conform(vor_DJF,vor_c,(/1,2,3/))
   vor_ano_d        =  dtrend_msg_n(yearse,vor_ano,True,True,0)                                     ; 去趋势
   copy_VarCoords(vor_DJF,vor_ano_d)

   printVarSummary(vor_ano_d)
   
   vor_ave1         =  dim_avg_n_Wrap(vor_DJF({stdy1},:,:,:),0)                                     ; 正1个标准差平均   
   vor_ave2         =  dim_avg_n_Wrap(vor_DJF({stdy2},:,:,:),0)                                     ; 负1个标准差平均
   var_vor1         =  dim_variance_n_Wrap(vor_DJF({stdy1},:,:,:),0)                                ; 多暴雨年的方差
   var_vor2         =  dim_variance_n_Wrap(vor_DJF({stdy2},:,:,:),0)                                ; 少暴雨年的方差 
   
   vor_ano12        =  vor_ave1
   vor_ano12        =  vor_ave1-vor_ave2                                                            ; 正负1个标准差相减  
   
   vor_ano1         =  dim_avg_n_Wrap(vor_ano({stdy1},:,:,:),0)                                     ; 多暴雨年的异常平均
   vor_ano2         =  dim_avg_n_Wrap(vor_ano({stdy2},:,:,:),0)                                     ; 少暴雨年的异常平均


   vor_anos         =  new((/3,dimsizes(vor_DJF&level),dimsizes(vor_DJF&lat),dimsizes(vor_DJF&lon)/),"float")  ; 便于后期画图
   vor_anos(0,:,:,:) = vor_ano1  
   vor_anos(1,:,:,:) = vor_ano2  
   vor_anos(2,:,:,:) = vor_ano12  
   
   vor_anos!0       =  "nub"
   vor_anos&nub     =  (/1,2,3/)
 
;;1.5.11流函数--;;;;;;;;;;

   sf_DJF           =  month_to_season(sf_ca,"JFM")                                                 ; 每年冬季
   copy_VarMeta(sf_ca(0,:,:,:),sf_DJF(0,:,:,:))
   sf_DJF!0         =  "year"
   sf_DJF&year      =  yearse
   sf_c             =  dim_avg_n_Wrap(sf_DJF({it_s:it_e},:,:,:),0)                                  ; 30年冬季平均 流函数
     
   sf_ano           =  sf_DJF
   sf_ano           =  sf_DJF-conform(sf_DJF,sf_c,(/1,2,3/))
   sf_ano_d         =  dtrend_msg_n(yearse,sf_ano,True,True,0)                                      ; 去趋势
   copy_VarCoords(sf_DJF,sf_ano_d)
   
   sf_ave1          =  dim_avg_n_Wrap(sf_DJF({stdy1},:,:,:),0)                                      ; 多暴雨年的均值   
   sf_ave2          =  dim_avg_n_Wrap(sf_DJF({stdy2},:,:,:),0)                                      ; 少暴雨年的均值
   var_sf1          =  dim_variance_n_Wrap(sf_DJF({stdy1},:,:,:),0)                                 ; 多暴雨年的方差
   var_sf2          =  dim_variance_n_Wrap(sf_DJF({stdy2},:,:,:),0)                                 ; 少暴雨年的方差
   
   sf_ano12         =  sf_ave1
   sf_ano12         =  sf_ave1-sf_ave2                                                              ; 正负1个标准差相减

   sf_ano1          =  dim_avg_n_Wrap(sf_ano({stdy1},:,:,:),0)                                      ; 多暴雨年的异常平均
   sf_ano2          =  dim_avg_n_Wrap(sf_ano({stdy2},:,:,:),0)                                      ; 少暴雨年的异常平均
   
;;1.5.12 散度风和速度势  

; 辐散风----------再将辐散风u,v分量分开

   divuv_DJF        =  uv2dvG_Wrap(uwd_DJF,vwd_DJF)           ; u,v ==> divergence  

   divu_DJF         = new ( dimsizes(uwd_DJF), typeof(uwd_DJF), uwd_DJF@_FillValue )
   divv_DJF         = new ( dimsizes(vwd_DJF), typeof(vwd_DJF), vwd_DJF@_FillValue )

   dv2uvg(divuv_DJF,divu_DJF,divv_DJF) ; div  ==> divergent  wind components
   copy_VarCoords(uwd_DJF, divu_DJF ) 
   copy_VarCoords(uwd_DJF, divv_DJF ) 
   divu_DJF@long_name= "Zonal Divergent Wind"
   divu_DJF@units    = uwd_DJF@units
   divv_DJF@long_name= "Meridional Divergent Wind"
   divv_DJF@units    = vwd_DJF@units
  
   divuwd_c         =  dim_avg_n_Wrap(divu_DJF({it_s:it_e},:,:,:),0)                                ; 30年冬季平均   
   divuwd_ano       =  divu_DJF
   divuwd_ano       =  divu_DJF-conform(divu_DJF,divuwd_c,(/1,2,3/))
   divuwd_ano_d     =  dtrend_msg_n(yearse,divuwd_ano,True,True,0)                                  ; 去趋势
   copy_VarCoords(divu_DJF,divuwd_ano_d)  
  

   divvwd_c         =  dim_avg_n_Wrap(divv_DJF({it_s:it_e},:,:,:),0)                                ; 30年冬季平均   
   divvwd_ano       =  divv_DJF
   divvwd_ano       =  divv_DJF-conform(divv_DJF,divvwd_c,(/1,2,3/))
   divvwd_ano_d     =  dtrend_msg_n(yearse,divvwd_ano,True,True,0)                                  ; 去趋势
   copy_VarCoords(divv_DJF,divvwd_ano_d)    
  
; calculate 速度势CHI
   chi_DJF          = ilapsG_Wrap ( divuv_DJF , 0)
   chi_DJF          = (/chi_DJF/1e4/)
   chi_DJF@long_name= "velocity potential"
   chi_DJF@units    = "m/s" 

   chi_c            =  dim_avg_n_Wrap(divv_DJF({it_s:it_e},:,:,:),0)                                ; 30年冬季平均   
   chi_ano          =  divv_DJF
   chi_ano          =  divv_DJF-conform(chi_DJF,chi_c,(/1,2,3/))
   chi_ano_d        =  dtrend_msg_n(yearse,chi_ano,True,True,0)                                     ; 去趋势
   copy_VarCoords(chi_DJF,chi_ano_d)    
  
;;T-N波作用通量-------------------------------------------------------------------
   lat              =  hgt_DJF&lat
   lon              =  hgt_DJF&lon

; make arrays for climatology that are suitable for calculation
; bzvar <- bzvarorg
; btvar <- btvarorg

czvar = new(dimsizes(hgt_DJF),float,hgt_DJF@_FillValue)
ctvar = new(dimsizes(hgt_DJF),float,hgt_DJF@_FillValue)
cuvar = new(dimsizes(hgt_DJF),float,hgt_DJF@_FillValue)
cvvar = new(dimsizes(hgt_DJF),float,hgt_DJF@_FillValue)

czvar = conform(hgt_DJF,hgt_c,(/1,2,3/))
ctvar = conform(air_DJF,air_c,(/1,2,3/))
cuvar = conform(uwd_DJF,uwd_c,(/1,2,3/))
cvvar = conform(vwd_DJF,vwd_c,(/1,2,3/))

;height anomaly(所有冬季-气候平均)
zavar = hgt_DJF - czvar

;  Gas constant
gc=290
;  Gravitational acceleration
ga=9.80665

;  Radius of the earth
re=6378388

; scale height
sclhgt=8000.

; pi
pi = atan(1.0)*4.

; Coriolis parameter
f =  2.*2.*pi/(60.*60.*24.)*sin(pi/180. * lat(:))
f!0 = "lat"
f&lat = lat
f@_FillValue = hgt_DJF@_FillValue

; missing for 10S - 10N
do ilat = 0, dimsizes(lat)-1
if (abs(lat(ilat) ).lt. 10. ) then
f(ilat)= f@_FillValue
end if
end do

; cosine
coslat = cos(lat(:)*pi/180.)

; 1-D -> 4-D
leveltmp = conform_dims(dimsizes(hgt_DJF),levl,1)
coslattmp = conform_dims(dimsizes(hgt_DJF),coslat,2)
ftmp = conform_dims(dimsizes(hgt_DJF),f,2)

; vertical gradient of potential temperature  (K/m)
dthetadz = center_finite_diff_n(ctvar*(1000./leveltmp)^0.286,-sclhgt*log(levl/1000),False,0,1)

; Brunt Vaisala frequency
NN = (gc*(leveltmp/1000.)^0.286)/sclhgt * dthetadz
NN@_FillValue = ctvar@_FillValue
NN = where(NN .gt. 0, NN, NN@_FillValue)

; magnitude of climatological wind
cumag = sqrt(cuvar^2 + cvvar^2)
cumag@_FillValue = cuvar@_FillValue
cumag = where(cumag .gt. 0, cumag, cumag@_FillValue)

; QG steam function for anomaly
psidev = zavar*ga /ftmp

;dpsidev/dlon
dpsidevdlon =  center_finite_diff_n(psidev,lon*pi/180.,True,0,3)

;ddpsidev/dlonlon
ddpsidevdlonlon =  center_finite_diff_n(dpsidevdlon,lon*pi/180.,True,0,3)

;dpsidev/dlat
dpsidevdlat = center_finite_diff_n(psidev, lat*pi/180., False,0,2)

;ddpsidev/dlonlat
ddpsidevdlonlat =  center_finite_diff_n(dpsidevdlon,lat*pi/180.,False,0,2)

;ddpsidev/dlatdlat
ddpsidevdlatlat = center_finite_diff_n(dpsidevdlat, lat*pi/180.,False,0,2)

;dpsidev/dz
dpsidevdz = center_finite_diff_n(psidev, -sclhgt*log(levl/1000),False,0,1)

;ddpsidev/dlondz
ddpsidevdlonz = center_finite_diff_n(dpsidevdlon, -sclhgt*log(levl/1000),False,0,1)

;ddpsidev/dlatdz
ddpsidevdlatz = center_finite_diff_n(dpsidevdlat, -sclhgt*log(levl/1000),False,0,1)

xuterm = (dpsidevdlon*dpsidevdlon - psidev*ddpsidevdlonlon)
xvterm = (dpsidevdlon*dpsidevdlat - psidev*ddpsidevdlonlat)
;yuterm = xvterm
yvterm = (dpsidevdlat*dpsidevdlat - psidev*ddpsidevdlatlat)

zuterm = (dpsidevdlon*dpsidevdz - psidev*ddpsidevdlonz)
zvterm = (dpsidevdlat*dpsidevdz - psidev*ddpsidevdlatz)

; Mask out where westerlies is small or negative (less than 5 m/s).
;  by using mask

;x-component of (38)
Fx = mask( leveltmp/1000./(2.*cumag*re*re)*( cuvar/coslattmp * xuterm + cvvar * xvterm), cuvar.lt.5,False)

;y-component 
Fy = mask(leveltmp/1000./(2.*cumag*re*re)*( cuvar*xvterm + coslattmp*cvvar*yvterm),  cuvar.lt.5,False)

;z-component 
Fz = mask(leveltmp/1000.*ftmp*ftmp/(2.*cumag*NN*re)*(cuvar*zuterm + coslattmp*cvvar*zvterm), cuvar.lt.5, False)

; for output
Fx!0 = "year"
Fx&year = yearse
Fx!1 = "level"
Fx&level = levl
Fx!2 = "lat"
Fx&lat = lat
Fx!3 = "lon"
Fx&lon = lon

Fy!0 = "year"
Fy&year = yearse
Fy!1 = "level"
Fy&level = levl
Fy!2 = "lat"
Fy&lat = lat
Fy!3 = "lon"
Fy&lon = lon

Fz!0 = "year"
Fz&year = yearse
Fz!1 = "level"
Fz&level = levl
Fz!2 = "lat"
Fz&lat = lat
Fz!3 = "lon"
Fz&lon = lon

psidev!0 = "year"
psidev&year = yearse
psidev!1 = "level"
psidev&level = levl
psidev!2 = "lat"
psidev&lat = lat
psidev!3 = "lon"
psidev&lon = lon

NN!0 = "year"
NN&year = yearse
NN!1 = "level"
NN&level = levl
NN!2 = "lat"
NN&lat = lat
NN!3 = "lon"
NN&lon = lon

Fx@units = "m^2/s^2"
Fx@units = "m^2/s^2"
Fz@units = "m^2/s^2"

psidev@units = "m^2/s"

NN@units = "1/s^2" 
 
;Fx 
   Fx_ave1          =  dim_avg_n_Wrap(Fx({stdy1},:,:,:),0)                                          ; 多暴雨年的均值
   Fx_ave2          =  dim_avg_n_Wrap(Fx({stdy2},:,:,:),0)                                          ; 少暴雨年的均值   
   var_Fx1          =  dim_variance_n_Wrap(Fx({stdy1},:,:,:),0)                                     ; 多暴雨年的方差
   var_Fx2          =  dim_variance_n_Wrap(Fx({stdy2},:,:,:),0)                                     ; 少暴雨年的方差        
   Fx_ano12         =  Fx_ave1
   Fx_ano12         =  Fx_ave1-Fx_ave2                                                              ; 正负标准差相减
;Fy
   Fy_ave1          =  dim_avg_n_Wrap(Fy({stdy1},:,:,:),0)                                          ; 多暴雨年的均值
   Fy_ave2          =  dim_avg_n_Wrap(Fy({stdy2},:,:,:),0)                                          ; 少暴雨年的均值 
   var_Fy1          =  dim_variance_n_Wrap(Fy({stdy1},:,:,:),0)                                     ; 多暴雨年的方差
   var_Fy2          =  dim_variance_n_Wrap(Fy({stdy2},:,:,:),0)                                     ; 少暴雨年的方差        
   Fy_ano12         =  Fy_ave1
   Fy_ano12         =  Fy_ave1-Fy_ave2                                                              ; 正负标准差相减
;Fz
   Fz_ave1          =  dim_avg_n_Wrap(Fz({stdy1},:,:,:),0)                                          ; 多暴雨年的均值
   Fz_ave2          =  dim_avg_n_Wrap(Fz({stdy2},:,:,:),0)                                          ; 少暴雨年的均值 
   var_Fz1          =  dim_variance_n_Wrap(Fz({stdy1},:,:,:),0)                                     ; 多暴雨年的方差
   var_Fz2          =  dim_variance_n_Wrap(Fz({stdy2},:,:,:),0)                                     ; 少暴雨年的方差     
   Fz_ano12         =  Fz_ave1
   Fz_ano12         =  Fz_ave1-Fz_ave2        
   
   Fx_ave           =  vor_anos                                                                     ; 便于后期画图
   Fx_ave(0,:,:,:)  =  Fx_ave1  
   Fx_ave(1,:,:,:)  =  Fx_ave2  
   Fx_ave(2,:,:,:)  =  Fx_ano12                                                                     ; 12为二者的差值场，与1，2的均值场不同

   Fy_ave           =  vor_anos                                                                     ; 
   Fy_ave(0,:,:,:)  =  Fy_ave1  
   Fy_ave(1,:,:,:)  =  Fy_ave2  
   Fy_ave(2,:,:,:)  =  Fy_ano12                                                                     ; 12为二者的差值场，与1，2的均值场不同
    
   Fz_ave           =  vor_anos                                                                     ; 
   Fz_ave(0,:,:,:)  =  Fz_ave1  
   Fz_ave(1,:,:,:)  =  Fz_ave2  
   Fz_ave(2,:,:,:)  =  Fz_ano12                                                                     ; 12为二者的差值场，与1，2的均值场不同
   
   hgt_anos          =  vor_anos                                                                     ; 便于后期画图
   hgt_anos(0,:,:,:) = hgt_ano1  
   hgt_anos(1,:,:,:) = hgt_ano2  
   hgt_anos(2,:,:,:) = hgt_ano12                                                                    ; 12为二者的差值场，与1，2的均值场不同

   sf_anos           =  vor_anos                                                                     ; 便于后期画图
   sf_anos(0,:,:,:) = sf_ano1  
   sf_anos(1,:,:,:) = sf_ano2  
   sf_anos(2,:,:,:) = sf_ano12  
   
   
print("TN-1219")   
   
;;2.--------------标准差的两组均值显著性分析-ttest检验;即std-6需要检验--------------------------------------------------------------------<
   iflag            = False
   x1               = dimsizes(stdy1)
   x2               = dimsizes(stdy2)

;;2.1--ttest-hgt位势高度---多暴雨和少暴雨年的差异检验    
   prob_hgt        := ttest(hgt_ave1, var_hgt1, x1, hgt_ave2, var_hgt2, x2, iflag, False)           ; 平均值，方差，样本数（维度同平均值）
   copy_VarCoords(hgt_ave1, prob_hgt)    
   prob_hgt         =   where(prob_hgt.eq.0,10,100*(1-prob_hgt))                                    ; 因为数据参与rtest运算，缺失值得到结果为0
;   prob_hgt         = 100*(1-prob_hgt)                                                                                      ; 因为数据参与rtest运算，缺失值得到结果为0    

;----------------------ttest;;;;;;hgt-positive
  hgt_std1          = dim_variance_n_Wrap(hgt_ano({stdy1},:,:,:),0)                                 ; 标准差
  hgt_std1          = sqrt(hgt_std1/x1)                                                             ; 标准差
  hgt_std1          = where(hgt_std1.eq.0,hgt_DJF@_FillValue,hgt_std1)  
  t_hgt1            = hgt_ano1/hgt_std1       
  prob_hgt1         = hgt_ano1
  prob_hgt1         = student_t(t_hgt1, x1-1)  
  prob_hgt1         = where(prob_hgt1.eq.0,10,100*(1-prob_hgt1))  
;----------------------ttest;;;;;;hgt-negtive
  hgt_std2          = dim_variance_n_Wrap(hgt_ano({stdy2},:,:,:),0)                                 ; 标准差
  hgt_std2          = sqrt(hgt_std2/x2)                                                             ; 标准差
  hgt_std2          = where(hgt_std2.eq.0,hgt_DJF@_FillValue,hgt_std2)  
  t_hgt2            = hgt_ano2/hgt_std2       
  prob_hgt2         = hgt_ano2
  prob_hgt2         = student_t(t_hgt2, x2-1)  
  prob_hgt2         = where(prob_hgt2.eq.0,10,100*(1-prob_hgt2))  
   
;;2.2--ttest-U水平风速 
   prob_uwd        := ttest(uwd_ave1, var_uwd1, x1, uwd_ave2, var_uwd2, x2, iflag, False)           ; 平均值，方差，样本数（维度同平均值）
   copy_VarCoords(uwd_ave1, prob_uwd) 
   
   prob_uwd         = where(prob_uwd.eq.0,10,100*(1-prob_uwd))                                      ; 因为数据参与rtest运算，缺失值得到结果为0
;   prob_uwd        = 100*(1-prob_uwd)                                                              ; 因为数据参与rtest运算，缺失值得到结果为0     
;----------------------uwd-positive
  uwd_std1          = dim_variance_n_Wrap(uwd_ano({stdy1},:,:,:),0)                                 ; 标准差
  uwd_std1          = sqrt(uwd_std1/x1)                                                             ; 标准差
  uwd_std1          = where(uwd_std1.eq.0,uwd_DJF@_FillValue,uwd_std1)  
  t_uwd1            = uwd_ano1/uwd_std1       
  prob_uwd1         = uwd_ano1
  prob_uwd1         = student_t(t_uwd1, x1-1)  
  prob_uwd1         = where(prob_uwd1.eq.0,10,100*(1-prob_uwd1))  
;----------------------uwd-negtive
  uwd_std2          = dim_variance_n_Wrap(uwd_ano({stdy2},:,:,:),0)                                 ; 标准差
  uwd_std2          = sqrt(uwd_std2/x2)                                                             ; 标准差
  uwd_std2          = where(uwd_std2.eq.0,uwd_DJF@_FillValue,uwd_std2)  
  t_uwd2            = uwd_ano2/uwd_std2       
  prob_uwd2         = uwd_ano2
  prob_uwd2         = student_t(t_uwd2, x2-1)  
  prob_uwd2         = where(prob_uwd2.eq.0,10,100*(1-prob_uwd2))  
  
;;2.3--ttest-V水平风速   
   prob_vwd        := ttest(vwd_ave1, var_vwd1, x1, vwd_ave2, var_vwd2, x2, iflag, False)           ; 平均值，方差，样本数（维度同平均值）
   copy_VarCoords(vwd_ave1, prob_vwd) 
   
   prob_vwd         = where(prob_vwd.eq.0,10,100*(1-prob_vwd))
;   prob_vwd         = 100*(1-prob_vwd)                                                             ; 因为数据参与rtest运算，缺失值得到结果为0    

   uwd_ano12_pb     = uwd_ano12
   vwd_ano12_pb     = vwd_ano12   
   uwd_ano12_pb     = where((prob_uwd.ge.sign).or.(prob_vwd.ge.sign),uwd_ano12,hgt_DJF@_FillValue)      ; 显著的风场(多暴雨年-少暴雨年)   
   vwd_ano12_pb     = where((prob_uwd.ge.sign).or.(prob_vwd.ge.sign),vwd_ano12,hgt_DJF@_FillValue)   

   uwd_ano12_pbn    = uwd_ano12
   vwd_ano12_pbn    = vwd_ano12
   uwd_ano12_pbn    = where((prob_uwd.lt.sign).and.(prob_vwd.lt.sign),uwd_ano12,hgt_DJF@_FillValue)     ; 不显著的风场(多暴雨年-少暴雨年)     
   vwd_ano12_pbn    = where((prob_uwd.lt.sign).and.(prob_vwd.lt.sign),vwd_ano12,hgt_DJF@_FillValue)

;----------------------vwd-positive
  vwd_std1          = dim_variance_n_Wrap(vwd_ano({stdy1},:,:,:),0)                                 ; 标准差
  vwd_std1          = sqrt(vwd_std1/x1)                                                             ; 标准差
  vwd_std1          = where(vwd_std1.eq.0,vwd_DJF@_FillValue,vwd_std1)  
  t_vwd1            = vwd_ano1/vwd_std1       
  prob_vwd1         = vwd_ano1
  prob_vwd1         = student_t(t_vwd1, x1-1)  
  prob_vwd1         = where(prob_vwd1.eq.0,10,100*(1-prob_vwd1))  
;----------------------vwd-negtive
  vwd_std2          = dim_variance_n_Wrap(vwd_ano({stdy2},:,:,:),0)                                 ; 标准差
  vwd_std2          = sqrt(vwd_std2/x2)                                                             ; 标准差                                 
  vwd_std2          = where(vwd_std2.eq.0,vwd_DJF@_FillValue,vwd_std2)  
  t_vwd2            = vwd_ano2/vwd_std2       
  prob_vwd2         = vwd_ano2
  prob_vwd2         = student_t(t_vwd2, x2-1)  
  prob_vwd2         = where(prob_vwd2.eq.0,10,100*(1-prob_vwd2))   
  
;----多暴雨年风场------------------------------------------------------------------------------------------------------------
   uwd_ano1_pb     = uwd_ano1
   vwd_ano1_pb     = vwd_ano1   
   uwd_ano1_pb     = where((prob_uwd1.ge.sign).or.(prob_vwd1.ge.sign),uwd_ano1,hgt_DJF@_FillValue)      ; 显著的风场(多暴雨年-气候场)     
   vwd_ano1_pb     = where((prob_uwd1.ge.sign).or.(prob_vwd1.ge.sign),vwd_ano1,hgt_DJF@_FillValue)      ; 显著的风场(多暴雨年-气候场)     

   uwd_ano1_pbn    = uwd_ano1
   vwd_ano1_pbn    = vwd_ano1
   uwd_ano1_pbn    = where((prob_uwd1.lt.sign).and.(prob_vwd1.lt.sign),uwd_ano1,hgt_DJF@_FillValue)     ; 不显著的风场 (多暴雨年-气候场)   
   vwd_ano1_pbn    = where((prob_uwd1.lt.sign).and.(prob_vwd1.lt.sign),vwd_ano1,hgt_DJF@_FillValue)     ; 不显著的风场 (多暴雨年-气候场)   
   
;----少暴雨年风场------------------------------------------------------------------------------------------------------------
   uwd_ano2_pb     = uwd_ano2
   vwd_ano2_pb     = vwd_ano2   
   uwd_ano2_pb     = where((prob_uwd2.ge.sign).or.(prob_vwd2.ge.sign),uwd_ano2,hgt_DJF@_FillValue)      ; 显著的风场(多暴雨年-气候场)     
   vwd_ano2_pb     = where((prob_uwd2.ge.sign).or.(prob_vwd2.ge.sign),vwd_ano2,hgt_DJF@_FillValue)      ; 显著的风场(多暴雨年-气候场)     

   uwd_ano2_pbn    = uwd_ano2
   vwd_ano2_pbn    = vwd_ano2
   uwd_ano2_pbn    = where((prob_uwd2.lt.sign).and.(prob_vwd2.lt.sign),uwd_ano2,hgt_DJF@_FillValue)     ; 不显著的风场 (多暴雨年-气候场)   
   vwd_ano2_pbn    = where((prob_uwd2.lt.sign).and.(prob_vwd2.lt.sign),vwd_ano2,hgt_DJF@_FillValue)     ; 不显著的风场 (多暴雨年-气候场)    
  print("1163")   
;----------------------------------------------------------------------------------------------------------------    
;;2.4--ttest-Quv整层水汽通量   
   prob_vint_qu    := ttest(vint_qu_ave1, var_vint_qu1, x1, vint_qu_ave2, var_vint_qu2, x2, iflag, False) ; 平均值，方差，样本数（维度同平均值）
   copy_VarCoords(vint_qu_ave1, prob_vint_qu)       
   prob_vint_qu     = where(prob_vint_qu.eq.0,10,100*(1-prob_vint_qu))
;   prob_vint_qu     = 100*(1-prob_vint_qu)  
   
   prob_vint_qv    := ttest(vint_qv_ave1, var_vint_qv1, x1, vint_qv_ave2, var_vint_qv2, x2, iflag, False) ; 平均值，方差，样本数（维度同平均值）
   copy_VarCoords(vint_qu_ave1, prob_vint_qu)   
   prob_vint_qv     = where(prob_vint_qv.eq.0,10,100*(1-prob_vint_qv))   
   
;   prob_vint_qv     = 100*(1-prob_vint_qv)     
   
;   prob_vint       = prob_vint_qu
;   prob_vint       = prob_vint_qu<prob_vint_qv                                                     ; 返回较大值

   vint_qu_ano12_pb = vint_qu_ano12  
   vint_qv_ano12_pb = vint_qv_ano12 
   vint_qu_ano12_pb = where((prob_vint_qu.ge.sign).or.(prob_vint_qv.ge.sign),vint_qu_ano12,vint_qu_ano12@_FillValue)  ; 显著的风场   
   vint_qv_ano12_pb = where((prob_vint_qu.ge.sign).or.(prob_vint_qv.ge.sign),vint_qv_ano12,vint_qu_ano12@_FillValue)

   vint_qu_ano12_pbn= vint_qu_ano12  
   vint_qv_ano12_pbn= vint_qv_ano12 
   vint_qu_ano12_pbn= where((prob_vint_qu.lt.sign).and.(prob_vint_qv.lt.sign),vint_qu_ano12,vint_qu_ano12@_FillValue)  ; 不显著的风场   
   vint_qv_ano12_pbn= where((prob_vint_qu.lt.sign).and.(prob_vint_qv.lt.sign),vint_qv_ano12,vint_qu_ano12@_FillValue)

;----------------------vint_qu-positive
  vint_qu_std1      = dim_variance_n_Wrap(vint_qu_ano({stdy1},:,:),0)                               ; 标准差
  vint_qu_std1      = sqrt(vint_qu_std1/x1)                                                         ; 标准差
  vint_qu_std1      = where(vint_qu_std1.eq.0,vwd_DJF@_FillValue,vint_qu_std1)  
  t_vint_qu1        = vint_qu_ano1/vint_qu_std1       
  prob_vint_qu1     = vint_qu_ano1
  prob_vint_qu1     = student_t(t_vint_qu1, x1-1)  
  prob_vint_qu1     = where(prob_vint_qu1.eq.0,10,100*(1-prob_vint_qu1))  
;----------------------vint_qu-negtive
  vint_qu_std2      = dim_variance_n_Wrap(vint_qu_ano({stdy2},:,:),0)                               ; 标准差
  vint_qu_std2      = sqrt(vint_qu_std2/x2)                                                         ; 标准差
  vint_qu_std2      = where(vint_qu_std2.eq.0,vwd_DJF@_FillValue,vint_qu_std2)  
  t_vint_qu2        = vint_qu_ano2/vint_qu_std2       
  prob_vint_qu2     = vint_qu_ano2
  prob_vint_qu2     = student_t(t_vint_qu2, x2-1)  
  prob_vint_qu2     = where(prob_vint_qu2.eq.0,10,100*(1-prob_vint_qu2))  
;------------------------------------------------------------------------------------------------------------------<  
;----------------------vint_qv-positive
  vint_qv_std1      = dim_variance_n_Wrap(vint_qv_ano({stdy1},:,:),0)                               ; 标准差
  vint_qv_std1      = sqrt(vint_qv_std1/x1)                                                         ; 标准差
  vint_qv_std1      = where(vint_qv_std1.eq.0,vwd_DJF@_FillValue,vint_qv_std1)  
  t_vint_qv1        = vint_qv_ano1/vint_qv_std1       
  prob_vint_qv1     = vint_qv_ano1
  prob_vint_qv1     = student_t(t_vint_qv1, x1-1)  
  prob_vint_qv1     = where(prob_vint_qv1.eq.0,10,100*(1-prob_vint_qv1))  
;----------------------vint_qv-negtive
  vint_qv_std2      = dim_variance_n_Wrap(vint_qv_ano({stdy2},:,:),0)                               ; 标准差
  vint_qv_std2      = sqrt(vint_qv_std2/x2)                                                         ; 标准差
  vint_qv_std2      = where(vint_qv_std2.eq.0,vwd_DJF@_FillValue,vint_qv_std2)  
  t_vint_qv2        = vint_qv_ano2/vint_qv_std2       
  prob_vint_qv2     = vint_qv_ano2
  prob_vint_qv2     = student_t(t_vint_qv2, x2-1)  
  prob_vint_qv2     = where(prob_vint_qv2.eq.0,10,100*(1-prob_vint_qv2))    

;----多暴雨年水汽通量场------------------------------------------------------------------------------------------------------------
  vint_qu_ano1_pb   = vint_qu_ano1                                                                  ; 显著的风场
  vint_qv_ano1_pb   = vint_qv_ano1 
  vint_qu_ano1_pb   = where((prob_vint_qu1.ge.sign).or.(prob_vint_qv1.ge.sign),vint_qu_ano1,vint_qu_ano12@_FillValue)     
  vint_qv_ano1_pb   = where((prob_vint_qu1.ge.sign).or.(prob_vint_qv1.ge.sign),vint_qv_ano1,vint_qu_ano12@_FillValue)

  vint_qu_ano1_pbn  = vint_qu_ano1                                                                  ; 不显著的风场   
  vint_qv_ano1_pbn  = vint_qv_ano1 
  vint_qu_ano1_pbn  = where((prob_vint_qu1.lt.sign).and.(prob_vint_qv1.lt.sign),vint_qu_ano2,vint_qu_ano12@_FillValue)  
  vint_qv_ano1_pbn  = where((prob_vint_qu1.lt.sign).and.(prob_vint_qv1.lt.sign),vint_qv_ano2,vint_qu_ano12@_FillValue)    
   
;----少暴雨年水汽通量场------------------------------------------------------------------------------------------------------------
  vint_qu_ano2_pb   = vint_qu_ano2                                                                  ; 显著的风场
  vint_qv_ano2_pb   = vint_qv_ano2 
  vint_qu_ano2_pb   = where((prob_vint_qu2.ge.sign).or.(prob_vint_qv2.ge.sign),vint_qu_ano2,vint_qu_ano12@_FillValue)     
  vint_qv_ano2_pb   = where((prob_vint_qu2.ge.sign).or.(prob_vint_qv2.ge.sign),vint_qv_ano2,vint_qu_ano12@_FillValue)  

  vint_qu_ano2_pbn  = vint_qu_ano2                                                                  ; 不显著的风场   
  vint_qv_ano2_pbn  = vint_qv_ano2 
  vint_qu_ano2_pbn  = where((prob_vint_qu2.lt.sign).and.(prob_vint_qv2.lt.sign),vint_qu_ano2,vint_qu_ano12@_FillValue)  
  vint_qv_ano2_pbn  = where((prob_vint_qu2.lt.sign).and.(prob_vint_qv2.lt.sign),vint_qv_ano2,vint_qu_ano12@_FillValue)    
;----------------------------------------------------------------------------------------------------------------    
;;2.5--ttest-Quv整层水汽通量散度
   prob_div_quv    := ttest(div_quv_ave1, var_div_quv1, x1, div_quv_ave2, var_div_quv2, x2, iflag, False) ; 平均值，方差，样本数（维度同平均值）
   copy_VarCoords(div_quv_ave1, prob_div_quv) 
    
   prob_div_quv     =   where(prob_div_quv.eq.0,10,100*(1-prob_div_quv))                            ; 因为数据参与rtest运算，缺失值得到结果为0 
;   prob_div_quv   =   100*(1-prob_div_quv)                                                         ; 因为数据参与rtest运算，缺失值得到结果为0 
  print("ttest-Quv")
;----------------------div_quv-positive
  div_quv_std1      = dim_variance_n_Wrap(div_quv_ano({stdy1},:,:),0)                               ; 标准差
  div_quv_std1      = sqrt(div_quv_std1/x1)                                                         ; 标准差
  div_quv_std1      = where(div_quv_std1.eq.0,div_quv_DJF@_FillValue,div_quv_std1)  
  t_div_quv1        = div_quv_ano1/div_quv_std1       
  prob_div_quv1     = div_quv_ano1
  prob_div_quv1     = student_t(t_div_quv1, x1-1)  
  prob_div_quv1     = where(prob_div_quv1.eq.0,10,100*(1-prob_div_quv1))  
;----------------------div_quv-negtive
  div_quv_std2      = dim_variance_n_Wrap(div_quv_ano({stdy2},:,:),0)                               ; 标准差
  div_quv_std2      = sqrt(div_quv_std2/x2)                                                         ; 标准差
  div_quv_std2      = sqrt(var_div_quv2/x2)                                                         ; 标准差
  div_quv_std2      = where(div_quv_std2.eq.0,div_quv_DJF@_FillValue,div_quv_std2)  
  t_div_quv2        = div_quv_ano2/div_quv_std2       
  prob_div_quv2     = div_quv_ano2
  prob_div_quv2     = student_t(t_div_quv2, x2-1)  
  prob_div_quv2     = where(prob_div_quv2.eq.0,10,100*(1-prob_div_quv2))   
  
;;2.6--ttest-div散度
   prob_div        := ttest(div_ave1, var_div1, x1, div_ave2, var_div2, x2, iflag, False)           ; 平均值，方差，样本数（维度同平均值）
   copy_VarCoords(div_ave1, prob_div)  
   prob_div         = where(prob_div.eq.0,10,100*(1-prob_div))   

;----------------------div-positive
  div_std1          = dim_variance_n_Wrap(div_ano({stdy1},:,:,:),0)                                 ; 多暴雨年的t检验
  div_std1          = sqrt(div_std1/x1)                                                            
  div_std1          = where(div_std1.eq.0,div_DJF@_FillValue,div_std1)  
  t_div1            = div_ano1/div_std1       
  prob_div1         = div_ano1
  prob_div1         = student_t(t_div1, x1-1)  
  prob_div1         = where(prob_div1.eq.0,10,100*(1-prob_div1))  
;----------------------div-negtive
  div_std2          = dim_variance_n_Wrap(div_ano({stdy2},:,:,:),0)                                 ; 少暴雨年的t检验
  div_std2          = sqrt(div_std2/x2)                                                             
  div_std2          = where(div_std2.eq.0,div_DJF@_FillValue,div_std2)  
  t_div2            = div_ano2/div_std2       
  prob_div2         = div_ano2
  prob_div2         = student_t(t_div2, x2-1)  
  prob_div2         = where(prob_div2.eq.0,10,100*(1-prob_div2))   

  print("ttest-div")
;;2.7--ttest-thq theta_se
   prob_eqt        := ttest(eqt_ave1, var_eqt1, x1, eqt_ave2, var_eqt2, x2, iflag, False)           ; 平均值，方差，样本数（维度同平均值）
   copy_VarCoords(eqt_ave1, prob_eqt)      
   prob_eqt         =   where(prob_eqt.eq.0,10,100*(1-prob_eqt))                                    ; 因为数据参与rtest运算，缺失值得到结果为0  
;   prob_eqt     =   100*(1-prob_eqt)                                                               ; 因为数据参与rtest运算，缺失值得到结果为0   

;----------------------eqt-positive
  eqt_std1          = dim_variance_n_Wrap(eqt_ano({stdy1},:,:,:),0)                                 ; 标准差
  eqt_std1          = sqrt(eqt_std1/x1)                                                             ; 标准差
  eqt_std1          = where(eqt_std1.eq.0,eqt_DJF@_FillValue,eqt_std1)  
  t_eqt1            = eqt_ano1/eqt_std1       
  prob_eqt1         = eqt_ano1
  prob_eqt1         = student_t(t_eqt1, x1-1)  
  prob_eqt1         = where(prob_eqt1.eq.0,10,100*(1-prob_eqt1))  
;----------------------eqt-negtive
  eqt_std2          = dim_variance_n_Wrap(eqt_ano({stdy2},:,:,:),0)                                 ; 标准差
  eqt_std2          = sqrt(eqt_std2/x2)                                                             ; 标准差
  eqt_std2          = where(eqt_std2.eq.0,eqt_DJF@_FillValue,eqt_std2)  
  t_eqt2            = eqt_ano2/eqt_std2       
  prob_eqt2         = eqt_ano2
  prob_eqt2         = student_t(t_eqt2, x2-1)  
  prob_eqt2         = where(prob_eqt2.eq.0,10,100*(1-prob_eqt2))    
  
 print("ttest-eqt")
;;2.8--ttest-vor涡度场
   prob_vor12       = ttest(vor_ave1, var_vor1, x1, vor_ave2, var_vor2, x2, iflag, False)           ; 平均值，方差，样本数（维度同平均值）
   copy_VarCoords(vor_ave1, prob_vor12)  
   prob_vor12       = where(prob_vor12.eq.0,10,100*(1-prob_vor12))                                  ; 暴雨多年与雨暴雨少年的差值检验

  vor_ano12_pb       = uwd_ano12
   vor_ano12_pb     = where(prob_vor12.ge.sign,vor_ano12,vor_DJF@_FillValue)      ; 显著的风场(多暴雨年-少暴雨年)   
   
;----------------------vor-positive
  vor_std1          = dim_variance_n_Wrap(vor_ano({stdy1},:,:,:),0)                                 ; 多暴雨年的t检验
  vor_std1          = sqrt(vor_std1/x1)                                                            
  vor_std1          = where(vor_std1.eq.0,vor_DJF@_FillValue,vor_std1)  
  t_vor1            = vor_ano1/vor_std1       
  prob_vor1         = vor_ano1
  prob_vor1         = student_t(t_vor1, x1-1)  
  prob_vor1         = where(prob_vor1.eq.0,10,100*(1-prob_vor1))  
  
  vor_ano1_pb       = uwd_ano1
  vor_ano1_pb       = where(prob_vor1.ge.sign,vor_ano1,vor_DJF@_FillValue)      ; 显著的风场(多暴雨年-少暴雨年)     
;----------------------vor-negtive
  vor_std2          = dim_variance_n_Wrap(vor_ano({stdy2},:,:,:),0)                                 ; 少暴雨年的t检验
  vor_std2          = sqrt(vor_std2/x2)                                                             
  vor_std2          = where(vor_std2.eq.0,vor_DJF@_FillValue,vor_std2)  
  t_vor2            = vor_ano2/vor_std2       
  prob_vor2         = vor_ano2
  prob_vor2         = student_t(t_vor2, x2-1)  
  prob_vor2         = where(prob_vor2.eq.0,10,100*(1-prob_vor2))   

  vor_ano2_pb       = uwd_ano2
  vor_ano2_pb       = where(prob_vor2.ge.sign,vor_ano2,vor_DJF@_FillValue)      ; 显著的风场(多暴雨年-少暴雨年)     
  
 print("ttest-vor")
;---------------;便于后期画图处理  
  prob_vors         = vor_anos                                                                    
  prob_vors(0,:,:,:)= prob_vor1  
  prob_vors(1,:,:,:)= prob_vor2  
  prob_vors(2,:,:,:)= prob_vor12  
     
  prob_hgts         = hgt_anos                                                                    
  prob_hgts(0,:,:,:)= prob_hgt1  
  prob_hgts(1,:,:,:)= prob_hgt2  
  prob_hgts(2,:,:,:)= prob_hgt      
     
     
     
;;2.9--ttest-Fx,Fy,Fz水平波作用通量
;;2.9.1 ttest-Fx 
   prob_Fx         := ttest(Fx_ave1, var_Fx1, x1, Fx_ave2, var_Fx2, x2, iflag, False)               ; 平均值，方差，样本数（维度同平均值）
   copy_VarCoords(Fx_ave1, prob_Fx)    
   prob_Fx          = where(prob_Fx.eq.0,10,100*(1-prob_Fx))                                        ; 因为数据参与rtest运算，缺失值得到结果为0
;   prob_Fx    = 100*(1-prob_Fx)                                                                    ; 因为数据参与rtest运算，缺失值得到结果为0     


;;2.9.2 ttest-Fy   
   prob_Fy         := ttest(Fy_ave1, var_Fy1, x1, Fy_ave2, var_Fy2, x2, iflag, False)               ; 平均值，方差，样本数（维度同平均值）
   copy_VarCoords(Fy_ave1, prob_Fy) 

   prob_Fy          = where(prob_Fy.eq.0,10,100*(1-prob_Fy))
;   prob_Fy     = 100*(1-prob_Fy)                                                                   ; 因为数据参与rtest运算，缺失值得到结果为0    
   
   Fx_ano12_pb      = Fx_ano12
   Fy_ano12_pb      = Fy_ano12   
   Fx_ano12_pb      = where((prob_Fx.ge.sign).or.(prob_Fy.ge.sign),Fx_ano12,hgt_DJF@_FillValue)         ; 显著的箭头，黑色   
   Fy_ano12_pb      = where((prob_Fx.ge.sign).or.(prob_Fy.ge.sign),Fy_ano12,hgt_DJF@_FillValue)   

   Fx_ano12_pbn     = Fx_ano12
   Fy_ano12_pbn     = Fy_ano12
   Fx_ano12_pbn     = where((prob_Fx.lt.sign).and.(prob_Fy.lt.sign),Fx_ano12,hgt_DJF@_FillValue)        ; 不显著的风场，灰色   
   Fy_ano12_pbn     = where((prob_Fx.lt.sign).and.(prob_Fy.lt.sign),Fy_ano12,hgt_DJF@_FillValue)

;;2.9.3 ttest-Fz
     
   prob_Fz         := ttest(Fz_ave1, var_Fz1, x1, Fz_ave2, var_Fz2, x2, iflag, False)               ; 平均值，方差，样本数（维度同平均值）
   copy_VarCoords(Fz_ave1, prob_Fz) 
;   prob_Fz     = where(prob_Fz.eq.0,10,100*(1-prob_Fz))
   prob_Fz          = 100*(1-prob_Fz)                                                               ; 因为数据参与rtest运算，缺失值得到结果为0    
   
   Fz_ano12_pb      = Fz_ano12
   
   Fxz_ano12_pb    := where((prob_Fx.ge.sign).or.(prob_Fz.ge.sign),Fx_ano12,hgt_DJF@_FillValue)         ; 因为数据参与rtest运算，缺失值得到结果为0   
   Fzz_ano12_pb    := where((prob_Fx.ge.sign).or.(prob_Fz.ge.sign),Fy_ano12,hgt_DJF@_FillValue)  
    
 print("ttest-fxy")   
;; 2.10--ttest 流函数
 
   prob_sf          = ttest(sf_ave1, var_sf1, x1, sf_ave2, var_sf2, x2, iflag, False)               ; 平均值，方差，样本数（维度同平均值）
   copy_VarCoords(sf_ave1, prob_sf)    
   prob_sf          =   where(prob_sf.eq.0,10,100*(1-prob_sf))                                      ; 因为数据参与rtest运算，缺失值得到结果为0
;   prob_sf         = 100*(1-prob_sf)                                                               ; 因为数据参与rtest运算，缺失值得到结果为0   

;----------------------sf-positive
  sf_std1           = dim_variance_n_Wrap(sf_ano({stdy1},:,:,:),0)                                 ; 标准差
  sf_std1           = sqrt(sf_std1/x1)                                                             ; 标准差
  sf_std1           = where(sf_std1.eq.0,sf_DJF@_FillValue,sf_std1)  
  t_sf1             = sf_ano1/sf_std1       
  prob_sf1          = sf_ano1
  prob_sf1          = student_t(t_sf1, x1-1)  
  prob_sf1          = where(prob_sf1.eq.0,10,100*(1-prob_sf1))  
;----------------------sf-negtive
  sf_std2           = dim_variance_n_Wrap(sf_ano({stdy2},:,:,:),0)                                 ; 标准差
  sf_std2           = sqrt(sf_std2/x2)                                                             ; 标准差
  sf_std2           = where(sf_std2.eq.0,sf_DJF@_FillValue,sf_std2)  
  t_sf2             = sf_ano2/sf_std2       
  prob_sf2          = sf_ano2
  prob_sf2          = student_t(t_sf2, x2-1)  
  prob_sf2          = where(prob_sf2.eq.0,10,100*(1-prob_sf2))  

 print("ttest-sf")   

;; 2.11--ttest-omega位势高度---多暴雨和少暴雨年的差异检验    
   prob_omega        := ttest(omega_ave1, var_omega1, x1, omega_ave2, var_omega2, x2, iflag, False) ; 平均值，方差，样本数（维度同平均值）
   copy_VarCoords(omega_ave1, prob_omega)    
   prob_omega         =   where(prob_omega.eq.0,10,100*(1-prob_omega))                              ; 因为数据参与rtest运算，缺失值得到结果为0
;   prob_omega         = 100*(1-prob_omega)                                                    

;----------------------ttest;;;;;;omega-positive
  omega_std1          = dim_variance_n_Wrap(omega_ano({stdy1},:,:,:),0)                             ; 标准差
  omega_std1          = sqrt(omega_std1/x1)                                                         ; 标准差
  omega_std1          = where(omega_std1.eq.0,omega_DJF@_FillValue,omega_std1)  
  t_omega1            = omega_ano1/omega_std1       
  prob_omega1         = omega_ano1
  prob_omega1         = student_t(t_omega1, x1-1)  
  prob_omega1         = where(prob_omega1.eq.0,10,100*(1-prob_omega1))  
;----------------------ttest;;;;;;omega-negtive
  omega_std2          = dim_variance_n_Wrap(omega_ano({stdy2},:,:,:),0)                             ; 标准差
  omega_std2          = sqrt(omega_std2/x2)                                                         ; 标准差
  omega_std2          = where(omega_std2.eq.0,omega_DJF@_FillValue,omega_std2)  
  t_omega2            = omega_ano2/omega_std2       
  prob_omega2         = omega_ano2
  prob_omega2         = student_t(t_omega2, x2-1)  
  prob_omega2         = where(prob_omega2.eq.0,10,100*(1-prob_omega2))       
  
;----多暴雨年风场-----------------------------------------------------------------------显著和不显著的风场-------------
   u_ano1_pb          = uwd_ano1(0:11,:,:)
   omega_ano1_pb      = omega_ano1   
   u_ano1_pb          = where((prob_uwd1(0:11,:,:).ge.sign).or.(prob_omega1.ge.sign),uwd_ano1(0:11,:,:),hgt_DJF@_FillValue) ; 显著的风场(多暴雨年-气候场)     
   omega_ano1_pb      = where((prob_uwd1(0:11,:,:).ge.sign).or.(prob_omega1.ge.sign),omega_ano1,hgt_DJF@_FillValue) ; 显著的风场(多暴雨年-气候场)     

   u_ano1_pbn         = uwd_ano1(0:11,:,:)
   omega_ano1_pbn     = omega_ano1
   u_ano1_pbn         = where((prob_uwd1(0:11,:,:).lt.sign).and.(prob_omega1.lt.sign),uwd_ano1(0:11,:,:),hgt_DJF@_FillValue)  ; 不显著的风场 (多暴雨年-气候场)   
   omega_ano1_pbn     = where((prob_uwd1(0:11,:,:).lt.sign).and.(prob_omega1.lt.sign),omega_ano1,hgt_DJF@_FillValue); 不显著的风场 (多暴雨年-气候场)   
   
;----少暴雨年风场-------------------------------------------------------------------------显著和不显著的风场-------------
   u_ano2_pb          = uwd_ano2(0:11,:,:)
   omega_ano2_pb      = omega_ano2   
   u_ano2_pb          = where((prob_uwd2(0:11,:,:).ge.sign).or.(prob_omega2.ge.sign),uwd_ano2(0:11,:,:),hgt_DJF@_FillValue) ; 显著的风场(多暴雨年-气候场)     
   omega_ano2_pb      = where((prob_uwd2(0:11,:,:).ge.sign).or.(prob_omega2.ge.sign),omega_ano2,hgt_DJF@_FillValue) ; 显著的风场(多暴雨年-气候场)     

   u_ano2_pbn         = uwd_ano2(0:11,:,:)
   omega_ano2_pbn     = omega_ano2
   u_ano2_pbn         = where((prob_uwd2(0:11,:,:).lt.sign).and.(prob_omega2.lt.sign),uwd_ano2(0:11,:,:),hgt_DJF@_FillValue)  ; 不显著的风场 (多暴雨年-气候场)   
   omega_ano2_pbn     = where((prob_uwd2(0:11,:,:).lt.sign).and.(prob_omega2.lt.sign),omega_ano2,hgt_DJF@_FillValue); 不显著的风场 (多暴雨年-气候场)     
   
   print("ttest-omega")  
   
;-------------------over;;climate winter ncepdata ----------------------------------------------------------------------------------------------------< 
 
         
;########################################### 二.大气环流#######################################################################################
;------------画图1---合成分析各种气候态
;>------气候场绘图 -------------合成正负0.6个标准差 --------------先画有阴影的，不然会被覆盖------------------------------------------------------
  minlat_l                              = -20                                                       ; 地图边界范围
  maxlat_l                              = 70
  minlon_l                              = 70
  maxlon_l                              = 170

  minlat                                = 20                                                        ; 东南地区地图边界范围（不包括海南岛）
  maxlat                                = 33
  minlon                                = 106
  maxlon                                = 123
  
;;2.1 暴雨指数合成分析
  base4                                 = new(5,"graphic")                                          ; 正标准差年合成分析
  plot4                                 = new(5,"graphic")    
  
  base5                                 = new(5,"graphic")                                          ; 负标准差年合成分析
  plot5                                 = new(5,"graphic")    
  
  base6                                 = new(5,"graphic")                                          ; 正减负标准差年合成分析
  plot6                                 = new(5,"graphic")    

  wksstd4                               = gsn_open_wks(wks_type,datadir+"pic/1-Interannual/3-circulation/1-plot-std4")    ; 标准差合成分析
  gsn_define_colormap(wksstd4,"CBR_coldhot")                                                       ; 调用rainbow+gray色板,，其它色板名称请查阅http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
  wksstd5                               = gsn_open_wks(wks_type,datadir+"pic/1-Interannual/3-circulation/1-plot-std5")    ; 标准差合成分析
  gsn_define_colormap(wksstd5,"CBR_coldhot")                                                       ; 调用rainbow+gray色板,，其它色板名称请查阅http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
  wksstd6                               = gsn_open_wks(wks_type,datadir+"pic/1-Interannual/3-circulation/1-plot-std6")    ; 标准差合成分析
  gsn_define_colormap(wksstd6,"CBR_coldhot")                                                       ; 调用rainbow+gray色板,，其它色板名称请查阅http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml

  res@mpMinLatF                         = minlat_l
  res@mpMaxLatF                         = maxlat_l
  res@mpMinLonF                         = minlon_l
  res@mpMaxLonF                         = maxlon_l
  res@tmYLMode                          = "Explicit"                                                ;对Y轴的坐标进行设置,显示自己想要的值
  res@tmYLValues                        = (/-20,0,20,40,60,80/)
  res@tmYLLabels                        = (/"20~S~o~N~S","EQ","20~S~o~N~N","40~S~o~N~N","60~S~o~N~N","80~S~o~N~N"/)
  res@tmXBMode                          = "Explicit"
  res@tmXBValues                        = (/70,90,110,130,150,170/)
  res@tmXBLabels                        = (/"70~S~o~N~E","90~S~o~N~E","110~S~o~N~E","130~S~o~N~E","150~S~o~N~E","170~S~o~N~E"/)
  print("1424")

;  ;青藏高原地形图 
;  plot_3                                :=  gsn_csm_contour(wks_850,z,resz)                          ; 5-青藏高原地形图1500m   
;  overlay(base_3,plot_3)         
  
  ;;1. 200hPa冬季标准差——急流----------------------------------------------- 
  res@gsnLeftString                     =  "Ano U200"                                               ; 子图的主标题 
  base4(0)                              =  gsn_csm_map(wksstd4,res)                                 ; ave1的三张图
  base5(0)                              =  gsn_csm_map(wksstd5,res)                                 ; ave2的三张图
  base6(0)                              =  gsn_csm_map(wksstd6,res)                                 ; ave1-ave2的三张图

;添加边框

  boxlat=(/maxlat,maxlat,minlat,minlat,maxlat/) 
  boxlon=(/minlon,maxlon,maxlon,minlon,minlon/)
  
  lnres = True  
   lnres@gsLineColor="black"                        ; gray50
   lnres@gsLineThicknessF=1.0
   lnres@gsLineDashPattern=0
  
  dum1=gsn_add_polyline(wksstd4,base4(0),boxlon,boxlat,lnres)  
  dum2=gsn_add_polyline(wksstd5,base5(0),boxlon,boxlat,lnres)  
  dum3=gsn_add_polyline(wksstd6,base6(0),boxlon,boxlat,lnres)  

;  ;;200.1纬向急流-异常场
;  resc@cnFillOn                         = True
;  resc@cnLinesOn                        = False                                                     ; 绘制等值线
;  resc@cnLevelSelectionMode             = "ExplicitLevels"
;  resc@cnLevels                        := (/-5,-4,-3,-2,-1,1,2,3,4,5/)                              ;   (/-5,-4,-3,-2,-1,1,2,3,4,5/)
;  resc@cnFillColors                    := (/2,3,4,5,6,0,8,9,10,11,12/)                  ; (/80,85,88,92,95,0,185,190,195,200,210/)
;  resc@cnLineColor                      = "Black"                                                   ; 颜色
;  resc@lbLabelBarOn                     = True                                                      ; 不画图标，统一绘制
;  resc@cnLineLabelsOn                   = False
;  resc@cnLineLabelFontHeightF           = 0.008
;  resc@pmLabelBarOrthogonalPosF         = 0.1
;  resc@lbLabelFontHeightF               = 0.01  
;  
;  print("1443")
;  
;  plot4(0)                              =  gsn_csm_c ontour(wksstd4,uwd_ano1({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)
;;  plot4(0)                             =  ColorNegDashZeroPosContour(plot4(0),"black","white","black")    
;  overlay(base4(0),plot4(0)) 
;  plot8                                 =  gsn_csm_contour(wksstd4,prob_uwd1({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)      ; 打点-显著性检验     
;  overlay(base4(0),plot8) 
;    
;  plot5(0)                              =  gsn_csm_contour(wksstd5,uwd_ano2({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc) 
;;  plot5(0)                             =  ColorNegDashZeroPosContour(plot5(0),"black","white","black")     
;  overlay(base5(0),plot5(0))  
;  plot9                                 =  gsn_csm_contour(wksstd4,prob_uwd2({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)      ; 打点-显著性检验     
;  overlay(base4(0),plot9) 
;
;  plot6(0)                              =  gsn_csm_contour(wksstd6,uwd_ano12({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
;;  plot6(0)                             =  ColorNegDashZeroPosContour(plot6(0),"black","white","black")     
;  overlay(base6(0),plot6(0))     
;  plot10                                =  gsn_csm_contour(wksstd6,prob_uwd({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)      ; 打点-显著性检验  
;  overlay(base6(0),plot10)  
;  print("1460")
  
  ;;200.2纬向急流-平均场(气候态)
  
  resc@cnFillOn                         = False
  resc@cnLinesOn                        = True                                                      ; 绘制等值线
  resc@cnLevelSelectionMode             = "ExplicitLevels"
  resc@cnLevels                        :=  (/30,50,70/)                                             ;   
  resc@cnLineColor                      = "black"                                                   ; 颜色
  resc@lbLabelBarOn                     = False                                                     ; 不画图标，统一绘制
  resc@cnLineLabelsOn                   = True
;  resc@cnLineLabelBackgroundColor      = "red"
  print("1472")
  
  plot4(0)                              =  gsn_csm_contour(wksstd4,uwd_c({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
  overlay(base4(0),plot4(0))                                                                        ; 带地图的图必须放在最下图层  
  print("1476")

  plot5(0)                              =  gsn_csm_contour(wksstd5,uwd_c({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
  overlay(base5(0),plot5(0))                                                                        ; 带地图的图必须放在最下图层 
  print("1480")

  plot6(0)                              =  gsn_csm_contour(wksstd6,uwd_c({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
  overlay(base6(0),plot6(0))                                      
   
  print("std-200")   

  ;;200.3纬向急流-平均场(极端降水多年和少年)
  
  resc@cnFillOn                         = False
  resc@cnLinesOn                        = True                                                      ; 绘制等值线
  resc@cnLevelSelectionMode             = "ExplicitLevels"
  resc@cnLevels                        :=  (/30,50,70/)                                             ;   
  resc@cnLineColor                      = "red"                                                     ; 颜色
  resc@lbLabelBarOn                     = False                                                     ; 不画图标，统一绘制
  resc@cnLineLabelsOn                   = True
;  resc@cnLineLabelBackgroundColor      = "red"
  print("1472")
  
  plot4(0)                              =  gsn_csm_contour(wksstd4,uwd_ave1({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
  overlay(base4(0),plot4(0))                                                                        ; 带地图的图必须放在最下图层  
  print("1476")

  plot5(0)                              =  gsn_csm_contour(wksstd5,uwd_ave2({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
  overlay(base5(0),plot5(0))                                                                        ; 带地图的图必须放在最下图层 
  print("1480")

  plot6(0)                              =  gsn_csm_contour(wksstd6,uwd_ano12({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
  overlay(base6(0),plot6(0))                                      
   
  print("std-200")     
  
  ;;2. 200散度场和风场异常-------------------------------------------------
  res@gsnLeftString                     =  "(a) Ano Div&UV200"                                          ; 子图的主标题 
  base4(1)                              =  gsn_csm_map(wksstd4,res)  
  base5(1)                              =  gsn_csm_map(wksstd5,res)    
  base6(1)                              =  gsn_csm_map(wksstd6,res)
  
;添加边框

  boxlat=(/maxlat,maxlat,minlat,minlat,maxlat/) 
  boxlon=(/minlon,maxlon,maxlon,minlon,minlon/)
  
  lnres = True  
   lnres@gsLineColor="black"                 ;gray50
   lnres@gsLineThicknessF=1.0
   lnres@gsLineDashPattern=0
    
  dum4=gsn_add_polyline(wksstd4,base4(1),boxlon,boxlat,lnres)  
  dum5=gsn_add_polyline(wksstd5,base5(1),boxlon,boxlat,lnres)  
  dum6=gsn_add_polyline(wksstd6,base6(1),boxlon,boxlat,lnres)    
      
  ;;2.1 200散度场异常---
  resc@cnFillOn                         = True
  resc@cnLinesOn                        = False                                                     ; 绘制等值线 AutomaticLevels
  resc@cnLevelSelectionMode             = "ExplicitLevels"
  resc@lbLabelBarOn                     = True                                                      ; 画图标，统一绘制
  resc@cnLineLabelsOn                   = False                                            
;   
  plot4(1)                              =  gsn_csm_contour(wksstd4,div_ano1({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)   
;  plot4(1)                             =  ColorNegDashZeroPosContour(plot4(1),"black","white","black")   
  overlay(base4(1),plot4(1))   
  plot8                                :=  gsn_csm_contour(wksstd4,prob_div1({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)      ; 打点-显著性检验  
  overlay(base4(1),plot8)  
;    
  plot5(1)                              =  gsn_csm_contour(wksstd5,div_ano2({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
;  plot5(1)                             =  ColorNegDashZeroPosContour(plot5(1),"black","white","black")      
  overlay(base5(1),plot5(1))   
  plot9                                :=  gsn_csm_contour(wksstd5,prob_div2({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)      ; 打点-显著性检验  
  overlay(base5(1),plot9)  
;  
  resc@cnLevels                        := (/-2.5,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5/)   
  resc@cnFillColors                    := (/2,3,4,5,6,0,8,9,10,11,12/)  ;  
  plot6(1)                              =  gsn_csm_contour(wksstd6,div_ano12({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
;  plot6(1)                             =  ColorNegDashZeroPosContour(plot6(1),"black","white","black")        
  overlay(base6(1),plot6(1))     
  plot10                               :=  gsn_csm_contour(wksstd6,prob_div({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)      ; 打点-显著性检验  
  overlay(base6(1),plot10)  

  ;;2.2 200风场异常 参考箭头10ms-1
  resv@vcRefMagnitudeF                  =  5                                                        ; 标准长度箭头所表示的大小
  resv@vcMinMagnitudeF                  =  0.5                                                      ; 要绘制箭头所表示的最小数值，即小于该数值则不绘制     
  resv@vcRefAnnoString2On               =  True                                                     ; 这里仅设定其下方的字符
  resv@vcRefAnnoString2                 =  5+" m/s"                                                 ; 参考箭头大小    

;----- 不显著的风场--灰色   
  resv@vcLineArrowColor                 =  "gray50"                                                 ; 不显著的风场--灰色 
  plot4(1)                              =  gsn_csm_vector(wksstd4,uwd_ano1_pbn({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano1_pbn({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv) 
  overlay(base4(1),plot4(1))    
  
  plot5(1)                              =  gsn_csm_vector(wksstd5,uwd_ano2_pbn({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano2_pbn({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv)  
  overlay(base5(1),plot5(1))       
  
  plot6(1)                              =  gsn_csm_vector(wksstd6,uwd_ano12_pbn({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano12_pbn({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv)   
  overlay(base6(1),plot6(1))   
;----- 显著的风场--黑色        
  resv@vcLineArrowColor                 =  "blue"                                                  ; 显著的风场--黑色
  plot4(1)                             =  gsn_csm_vector(wksstd4,uwd_ano1_pb({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano1_pb({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv) 
  overlay(base4(1),plot4(1))  
  
  plot5(1)                             =  gsn_csm_vector(wksstd5,uwd_ano2_pb({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano2_pb({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv) 
  overlay(base5(1),plot5(1)) 
  
  plot6(1)                             =  gsn_csm_vector(wksstd6,uwd_ano12_pb({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano12_pb({200},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv) 
  overlay(base6(1),plot6(1))       

 print("std-200")  

  ;;3.  整层水汽通量及其散度场 ----异常场 
  
  ;;3.1 整层水汽通量散度-异常场
  res_qu                                            = res
  res_qu@gsnLeftString             =  "(d) Ano Quv&Div_Quv"                                        ; 子图的主标题 

   res_qu@tmYLMode                  = "Explicit"                                                ; 对Y轴的坐标进行设置,显示自己想要的值
   res_qu@tmYLValues                = (/-10,0,10,20,30,40/)
   res_qu@tmYLLabels                = (/"10~S~o~S~N","EQ","10~S~o~N~N","20~S~o~N~N","30~S~o~N~N","40~S~o~N~N"/)
   res_qu@tmXBMode                  = "Explicit"
   res_qu@tmXBValues               := (/60,80,100,120,140,160/)
   res_qu@tmXBLabels               := (/"60~S~o~N~E","80~S~o~N~E","100~S~o~N~E","120~S~o~N~E","140~S~o~N~E","160~S~o~N~E"/)  
  
  base4(2)                              =  gsn_csm_map(wksstd4,res_qu) 
  base5(2)                              =  gsn_csm_map(wksstd5,res_qu)      
  base6(2)                              =  gsn_csm_map(wksstd6,res_qu)

  resc@cnFillOn                         =  True
  resc@cnLinesOn                        =  False                                                    ; 绘制等值线
  resc@cnLevelSelectionMode             =  "ExplicitLevels"
  resc@cnLevels                        :=  (/-2.5,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5/)         ;   
  resc@cnFillColors                    :=  (/2,3,4,5,6,0,8,9,10,11,12/)                 ;
;  resc@cnLineColor                      =  "Black"                                                 ; 颜色
  resc@lbLabelBarOn                     =  True                                                     ; 不画图标，统一绘制
  resc@cnLineLabelsOn                   =  False
  resc@tfPolyDrawOrder                  =  "PreDraw"
  
  plot4(2)                              =  gsn_csm_contour(wksstd4,div_quv_ano1({minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)
  overlay(base4(2),plot4(2))                                                                        ; 带地图的图必须放在最下图层
  plot8                                :=  gsn_csm_contour(wksstd4,prob_div_quv1({minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)      ; 打点-显著性检验  
  overlay(base4(2),plot8)  
;    
  plot5(2)                              =  gsn_csm_contour(wksstd5,div_quv_ano2({minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)     
  overlay(base5(2),plot5(2))  
  plot9                                :=  gsn_csm_contour(wksstd5,prob_div_quv2({minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)      ; 打点-显著性检验  
  overlay(base5(2),plot9) 
; 
  plot6(2)                              =  gsn_csm_contour(wksstd6,div_quv_ano12({minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)    
  overlay(base6(2),plot6(2))                                                                          
  plot10                               :=  gsn_csm_contour(wksstd6,prob_div_quv({minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)      ; 打点-显著性检验 prob_div_quv 
  overlay(base6(2),plot10)   
 
  ;;3.2 整层水汽通量-异常场 
  reffxa                                =  50   
  resfx_c                               =  resv
  resfx_c@vcRefMagnitudeF               =  reffxa                                                    ; 标准长度箭头所表示的大小
  resfx_c@vcMinMagnitudeF               =  5                                                         ; 要绘制箭头所表示的最小数值，即小于该数值则不绘制     
  resfx_c@vcRefAnnoString2On            =  True                                                      ; 这里仅设定其下方的字符
  resfx_c@vcRefAnnoString2              =  reffxa+" kg/(m*s)"                                        ; 参考箭头大小

;----- 不显著的风场--灰色    
  resfx_c@vcLineArrowColor              =  "gray50"                                                 ; 不显著的风场--灰色 
  plot4(2)                              =  gsn_csm_vector(wksstd4,vint_qu_ano1_pbn,vint_qv_ano1_pbn,resfx_c)   
  overlay(base4(2),plot4(2))     
  plot5(2)                              =  gsn_csm_vector(wksstd5,vint_qu_ano2_pbn,vint_qv_ano2_pbn,resfx_c)  
  overlay(base5(2),plot5(2)) 
  plot6(2)                              =  gsn_csm_vector(wksstd6,vint_qu_ano12_pbn,vint_qv_ano12_pbn,resfx_c)    
  overlay(base6(2),plot6(2))  

;----- 显著的风场--黑色 
  resfx_c@vcLineArrowColor              =  "blue"                                                  ; 显著的风场--黑色

  plot4(2)                             =  gsn_csm_vector(wksstd4,vint_qu_ano1_pb,vint_qv_ano1_pb,resfx_c)   
  overlay(base4(2),plot4(2))     
  plot5(2)                             =  gsn_csm_vector(wksstd5,vint_qu_ano2_pb,vint_qv_ano2_pb,resfx_c)  
  overlay(base5(2),plot5(2)) 
  plot6(2)                             =  gsn_csm_vector(wksstd6,vint_qu_ano12_pb,vint_qv_ano12_pb,resfx_c)  
  overlay(base6(2),plot6(2))       

;添加边框

  boxlat=(/maxlat,maxlat,minlat,minlat,maxlat/) 
  boxlon=(/minlon,maxlon,maxlon,minlon,minlon/)
  
  lnres = True  
   lnres@gsLineColor="black"
   lnres@gsLineThicknessF=1.0
   lnres@gsLineDashPattern=0  
   
  dum7=gsn_add_polyline(wksstd4,base4(2),boxlon,boxlat,lnres)  
  dum8=gsn_add_polyline(wksstd5,base5(2),boxlon,boxlat,lnres)  
  dum9=gsn_add_polyline(wksstd6,base6(2),boxlon,boxlat,lnres)    

 print("Quv")  
    
;;4. 850Hgt和风场异常-------------------------------------------------

  res@gsnLeftString                     =  "(c) Ano Div&UV850"                                          ; 子图的主标题 
  base4(3)                              =  gsn_csm_map(wksstd4,res)  
  base5(3)                              =  gsn_csm_map(wksstd5,res)       
  base6(3)                              =  gsn_csm_map(wksstd6,res) 
 
  ;;4.1 850Theta_se---
;  resc@cnFillOn                       = True
;  resc@cnLinesOn                      = False                                                      ; 绘制等值线
;  resc@cnLevelSelectionMode           = "ExplicitLevels"
;  resc@cnLevels                      := (/-5,-4,-3,-2,-1,1,2,3,4,5/)  ;   
;  resc@cnFillColors                  := (/80,85,88,92,95,0,180,185,190,195,200/)  
;  resc@lbLabelBarOn                   = True                                                       ; 画图标
;  resc@cnLineLabelsOn                 = False                                            
   
;  plot4(3)                            =  gsn_csm_contour(wksstd4,eqt_ave1_ano({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)   
;  plot4(3)                            =  ColorNegDashZeroPosContour(plot4(3),"black","white","black")   
;  overlay(base4(3),plot4(3))   
  
;  plot5(3)                            =  gsn_csm_contour(wksstd5,eqt_ave2_ano({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
;  plot5(3)                            =  ColorNegDashZeroPosContour(plot5(3),"black","white","black")      
;  overlay(base5(3),plot5(3))   
  
;  plot6(3)                            =  gsn_csm_contour(wksstd6,eqt_ano12({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
;  plot6(3)                            =  ColorNegDashZeroPosContour(plot6(3),"black","white","black")        
;  overlay(base6(3),plot6(3))    
   
;  plot10                             :=  gsn_csm_contour(wksstd6,prob_eqt({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)      ; 打点-显著性检验    
;  overlay(base6(3),plot10)   

  ;;4.1 850散度场---
  
  resc@cnFillOn                         = True
  resc@cnLinesOn                        = False                                                     ; 绘制等值线
  resc@cnLevelSelectionMode             = "ExplicitLevels"
  resc@cnLevels                        := (/-2.5,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5/)                       
  resc@cnFillColors                    := (/2,3,4,5,6,0,8,9,10,11,12/)                  
  resc@lbLabelBarOn                     = True                                                      ; 画图标，统一绘制
  resc@cnLineLabelsOn                   = False                                            
   
  plot4(3)                              =  gsn_csm_contour(wksstd4,div_ano1({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)   
;  plot4(3)                             =  ColorNegDashZeroPosContour(plot4(5),"black","white","black")   
  overlay(base4(3),plot4(3))   
  plot8                                :=  gsn_csm_contour(wksstd4,prob_div1({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)      ; 打点-显著性检验  
  overlay(base4(3),plot8)     
  
  plot5(3)                              =  gsn_csm_contour(wksstd5,div_ano2({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
;  plot5(3)                             =  ColorNegDashZeroPosContour(plot5(5),"black","white","black")      
  overlay(base5(3),plot5(3))   
  plot9                                :=  gsn_csm_contour(wksstd5,prob_div2({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)      ; 打点-显著性检验  
  overlay(base5(3),plot9)   
    
  plot6(3)                              =  gsn_csm_contour(wksstd6,div_ano12({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
;  plot6(3)                             =  ColorNegDashZeroPosContour(plot6(5),"black","white","black")        
  overlay(base6(3),plot6(3))     
  plot10                               :=  gsn_csm_contour(wksstd6,prob_div({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)       ; 打点-显著性检验  
  overlay(base6(3),plot10)      
  
  ;;4.2 850风场异常 参考箭头ms-1
  resv@vcRefMagnitudeF                  =  2                                                        ; 标准长度箭头所表示的大小
  resv@vcMinMagnitudeF                  =  0.1                                                      ; 要绘制箭头所表示的最小数值，即小于该数值则不绘制     
  resv@vcRefAnnoString2On               =  True                                                     ; 这里仅设定其下方的字符
  resv@vcRefAnnoString2                 =  2+" m/s"                                                 ; 参考箭头大小  

;---不显著的箭头为灰色  
  resv@vcLineArrowColor                 =  "gray50"                                                 ; 不显著的箭头为灰色
  plot4(3)                              =  gsn_csm_vector(wksstd4,uwd_ano1_pbn({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano1_pbn({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv) 
  overlay(base4(3),plot4(3))  
  plot5(3)                              =  gsn_csm_vector(wksstd5,uwd_ano2_pbn({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano2_pbn({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv) 
  overlay(base5(3),plot5(3))  
  plot6(3)                              =  gsn_csm_vector(wksstd6,uwd_ano12_pbn({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano12_pbn({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv)   ; create plot；不显著的风场--灰色
  overlay(base6(3),plot6(3)) 

;---显著的箭头为黑色    
  resv@vcLineArrowColor                 =  "blue"                                                  ; 显著的箭头为黑色
  plot4(3)                              =  gsn_csm_vector(wksstd4,uwd_ano1_pb({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano1_pb({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv) 
  overlay(base4(3),plot4(3))  
  plot5(3)                              =  gsn_csm_vector(wksstd5,uwd_ano2_pb({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano2_pb({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv) 
  overlay(base5(3),plot5(3))  
  plot6(3)                              =  gsn_csm_vector(wksstd6,uwd_ano12_pb({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano12_pb({850},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv) 
  overlay(base6(3),plot6(3))       
  

;  ;850.3青藏高原地形图 

  shpfn2                 = datadir+"data/Qinghai-Tibet_Plateau/Qinghai-Tibet_Plateau.shp"

  shres                  = lnres
  shres@gsLineColor      = "black"
  shres@gsLineThicknessF = 1.1
  
  shp_plot4                             = gsn_add_shapefile_polylines(wksstd4,base4(3),shpfn2,shres)      
  shp_plot5                             = gsn_add_shapefile_polylines(wksstd5,base5(3),shpfn2,shres)      
  shp_plot6                             = gsn_add_shapefile_polylines(wksstd6,base6(3),shpfn2,shres)      

;;青藏高原地形图 
;   plot4(3)                             =  gsn_csm_contour(wksstd4,z,resz)                          ; 青藏高原地形图2000m   
;   plot5(3)                             =  gsn_csm_contour(wksstd5,z,resz)                          ; 青藏高原地形图2000m   
;   plot6(3)                             =  gsn_csm_contour(wksstd6,z,resz)                          ; 青藏高原地形图2000m   
;  overlay(base4(3),plot4(3))    
;  overlay(base5(3),plot5(3))    
;  overlay(base6(3),plot6(3))   

;添加边框-----------------------------------------------------------
  boxlat=(/maxlat,maxlat,minlat,minlat,maxlat/) 
  boxlon=(/minlon,maxlon,maxlon,minlon,minlon/)
  
  lnres = True  
   lnres@gsLineColor="black"
   lnres@gsLineThicknessF=1.0
   lnres@gsLineDashPattern=0  
   
  dum10 = gsn_add_polyline(wksstd4,base4(3),boxlon,boxlat,lnres)  
  dum11 = gsn_add_polyline(wksstd5,base5(3),boxlon,boxlat,lnres)  
  dum12 = gsn_add_polyline(wksstd6,base6(3),boxlon,boxlat,lnres)   
;---------------------------------------------------------添加边框-<    
  print("std-850")  

 ;;6. 500高度场和垂直速度异常-------------------------------------------------

  res@gsnLeftString                     =  "(b) Ano Omega&UV500"                                          ; 子图的主标题 
  base4(4)                              =  gsn_csm_map(wksstd4,res)
  base5(4)                              =  gsn_csm_map(wksstd5,res)    
  base6(4)                              =  gsn_csm_map(wksstd6,res)
  
;添加边框
  boxlat=(/maxlat,maxlat,minlat,minlat,maxlat/) 
  boxlon=(/minlon,maxlon,maxlon,minlon,minlon/)
  
  lnres = True  
   lnres@gsLineColor="black"
   lnres@gsLineThicknessF=1.0
   lnres@gsLineDashPattern=0  
  dum21=gsn_add_polyline(wksstd4,base4(4),boxlon,boxlat,lnres)  
  dum22=gsn_add_polyline(wksstd5,base5(4),boxlon,boxlat,lnres)  
  dum23=gsn_add_polyline(wksstd6,base6(4),boxlon,boxlat,lnres)   
        

  ;;6.1 500垂直速度场异常---
  resc@cnFillOn                       = True
  resc@cnLinesOn                      = False                                                ; 绘制等值线
  resc@cnLevelSelectionMode           = "ExplicitLevels"                                    ; ExplicitLevels
  resc@cnLevels                      := (/-5.0,-4.0,-3.0,-2.0,-1.0,1.0,2.0,3.0,4.0,5.0/)  ;   
  resc@cnFillColors                  := (/2,3,4,5,6,0,8,9,10,11,12/)  ;
  resc@lbLabelBarOn                   = True                                                 ; 画图标，统一绘制
  resc@cnLineLabelsOn                 = False                                            
   
  plot4(4)                            =  gsn_csm_contour(wksstd4,omega_ano1({500},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)   
  overlay(base4(4),plot4(4))   
  plot9                              :=  gsn_csm_contour(wksstd4,prob_omega1({500},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)      ; 打点-显著性检验  
  overlay(base4(4),plot9)   
  
  plot5(4)                            =  gsn_csm_contour(wksstd5,omega_ano2({500},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
  overlay(base5(4),plot5(4))   
  plot10                             :=  gsn_csm_contour(wksstd5,prob_omega2({500},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)      ; 打点-显著性检验  
  overlay(base5(4),plot10)   

  
  plot6(4)                            =  gsn_csm_contour(wksstd6,omega_ano12({500},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
  overlay(base6(4),plot6(4))   
  
  plot11                             :=  gsn_csm_contour(wksstd6,prob_omega({500},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)      ; 打点-显著性检验    
  overlay(base6(4),plot11)   


;  ;;2.1 500位势高度场异常---
;  resc@cnFillOn                         = False
;  resc@cnLinesOn                        = True                                                    ; 绘制等值线
;  resc@cnLevelSelectionMode             = "AutomaticLevels"
;  resc@cnLevels                        := (/-25,-20,-15,-10,-5,5,10,15,20,25/)  ;   
;  resc@cnFillColors                    := (/2,3,4,5,6,0,8,9,10,11,12/)  ;
;  resc@lbLabelBarOn                     = False                                                      ; 画图标，统一绘制
;  resc@cnLineLabelsOn                   = True                                            
;   
;  plot9                                :=  gsn_csm_contour(wksstd4,hgt_ano1({500},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)
;  plot9                                 =  ColorNegDashZeroPosContour(plot9,"blue","white","red")      
;  overlay(base4(4),plot9)   
;  
;  plot10                               :=  gsn_csm_contour(wksstd5,hgt_ano2({500},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
;  plot10                                =  ColorNegDashZeroPosContour(plot10,"blue","white","red")      
;  overlay(base5(4),plot10)   
;;  plot9                              :=  gsn_csm_contour(wksstd5,prob_hgt2({500},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)     ; 打点-显著性检验  
;;  overlay(base4(4),plot9) 
;  
;  plot11                               :=  gsn_csm_contour(wksstd6,hgt_ano12({500},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resc)  
;  plot11                                =  ColorNegDashZeroPosContour(plot11,"blue","white","red")        
;  overlay(base6(4),plot11)     
;;  plot10                             :=  gsn_csm_contour(wksstd6,prob_hgt({500},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resd)      ; 打点-显著性检验  
;;  overlay(base6(4),plot10)  
  ;;6.2 500风场异常 参考箭头ms-1
  resv@vcRefMagnitudeF                  =  3                                                        ; 标准长度箭头所表示的大小
  resv@vcMinMagnitudeF                  =  0.1                                                      ; 要绘制箭头所表示的最小数值，即小于该数值则不绘制     
  resv@vcRefAnnoString2On               =  True                                                     ; 这里仅设定其下方的字符
  resv@vcRefAnnoString2                 =  3+" m/s"                                                 ; 参考箭头大小  

;---不显著的箭头为灰色 
  lvl                                   =  500 
  resv@vcLineArrowColor                 =  "gray50"                                                 ; 不显著的箭头为灰色
  plot4(4)                              =  gsn_csm_vector(wksstd4,uwd_ano1_pbn({lvl},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano1_pbn({lvl},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv) 
  overlay(base4(4),plot4(4))  
  plot5(4)                              =  gsn_csm_vector(wksstd5,uwd_ano2_pbn({lvl},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano2_pbn({lvl},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv) 
  overlay(base5(4),plot5(4))  
  plot6(4)                              =  gsn_csm_vector(wksstd6,uwd_ano12_pbn({lvl},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano12_pbn({lvl},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv)   ; create plot；不显著的风场--灰色
  overlay(base6(4),plot6(4)) 

;---显著的箭头为黑色    
  resv@vcLineArrowColor                 =  "blue"                                                  ; 显著的箭头为黑色
  plot4(4)                              =  gsn_csm_vector(wksstd4,uwd_ano1_pb({lvl},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano1_pb({lvl},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv) 
  overlay(base4(4),plot4(4))  
  plot5(4)                              =  gsn_csm_vector(wksstd5,uwd_ano2_pb({lvl},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano2_pb({lvl},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv) 
  overlay(base5(4),plot5(4))  
  plot6(4)                              =  gsn_csm_vector(wksstd6,uwd_ano12_pb({lvl},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),vwd_ano12_pb({lvl},{minlat_l:maxlat_l},{minlon_l:maxlon_l}),resv) 
  overlay(base6(4),plot6(4))       

 print("std-500")  
  
;;组合成一张图
  resPstd                               = True                                                      ; 绘制组图 
  resPstd@gsnPanelRowSpec               = True                                                      ; 指定每行绘制的子图的个数
  resPstd@gsnPanelCenter                = True                                                      ; 子图居中
  resPstd@gsnPanelBottom                = 0.05                                                      ; add some space at bottom 
  resPstd@gsnPanelYWhiteSpacePercent    = 5.5                                                       ; 两个相邻子图之间沿Y轴方向的空白大小     
  
  gsn_panel(wksstd4,base4(0:4),(/1,2,2/),resPstd)
  gsn_panel(wksstd5,base5(0:4),(/1,2,2/),resPstd)
  gsn_panel(wksstd6,base6(0:4),(/1,2,2/),resPstd)
  
 print("panel") 
 delete([/base4,base5,base6/])
 delete([/wksstd4,wksstd5,wksstd6/]) 
 
;######################--------------画图2---高度场剖面（Hgt/thesa/omega）-------------正负0.6个标准差-------------------------------######################-
;######################---------------画图2---高度场剖面（Hgt/thesa/omega）-------------正负0.6个标准差-------------------------------######################-

  plat                                  = 24                                                        ; 北纬24度做剖面
  lat_pmin                              = 80
  lat_pmax                              = 140
  
  wks_p                                 = gsn_open_wks(wks_type,datadir+"pic/1-Interannual/3-circulation/2-plot-prs_hgt") ; 标准差合成分析
  gsn_define_colormap(wks_p,"temp_diff_18lev")                                                      ; 调用rainbow+gray色板,，其它色板名称请查阅http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
  base_p                                = new(2,"graphic")                                          ; 正标准差年合成分析 
  
  res_p@tmXTOn                          = False                                                     ; 去除X轴上部的刻度
  res_p@tmYROn                          = True                                                      ; 不去除Y轴右侧的刻度
  res_p@tmXBMode                        = "Explicit"
  res_p@tmXBValues                      = (/80,100,120,140/)
  res_p@tmXBLabels                      = (/"80~S~o~N~E","100~S~o~N~E","120~S~o~N~E","140~S~o~N~E"/)

  res_p@tmYLMode                        = "Explicit"                                                ;对Y轴的坐标进行设置,显示自己想要的值
  res_p@tmYLValues                      = (/1000,850,700,500,300,200,100/)
  res_p@tmYLLabels                      = (/"1000","850","700","500","300","200","100"/)

;;1. Hgt（17层）----------------------------------------------- 
  res_p@cnFillOn                        = False                                                     ; 关闭等值线填色 
  res_p@cnLineThicknessF                = 0.7                                                       ; 等值线粗细  
  res_p@gsnContourZeroLineThicknessF    = 0.                                                        ; 设置0值线粗细。0则不画
  res_p@cnLineLabelsOn                  = False                                                     ; 关闭标签
  res_p@cnInfoLabelOn                   = False                                                     ; 打开图右下方的等值线信息标签
  res_p@lbLabelBarOn                    = True 
  res_p@cnSmoothingOn                   = True                                                      ; 平滑线
  res_p@cnSmoothingDistanceF            = 0.01 
  res_p@cnSmoothingTensionF             = 0.1
  res_p@cnLineLabelBackgroundColor      = -1
  res_p@cnLineLabelPlacementMode        = "Constant"
  res_p@cnLineLabelFontHeightF          = 0.010
  res_p@cnLineLabelInterval             = 1
  res_p@lbLabelFontHeightF              = 0.012                                                     ; 色标标签字体大小
  res_p@pmLabelBarOrthogonalPosF        = 0.08
  res_p@pmLabelBarHeightF               = 0.09                                                      ; labelbar高度

  res_p@cnLinesOn                       = True                                                      ; 绘制等值线
  res_p@cnLevelSelectionMode            = "ExplicitLevels"
  res_p@cnLevels                       := (/-15,-12,-9,-6,-3,3,6,9,12,15/)  
  res_p@cnFillColors                   := (/10,20,30,40,45,0,65,75,85,95,100/)  
  res_p@lbLabelBarOn                    = False                                                     ; 画图标，统一绘制
  res_p@cnLineLabelsOn                  = True                                           
  res_p@cnLineDrawOrder                 = "PostDraw"                                                ; "PostDraw" 最后画（在最上面），会覆盖 the predraw and draw "PreDraw" ;先画 就是放在最下面 包括地图下面;"Draw"  ;会覆盖the predraw 但会在the postdraw下面   
   
  base_p(0)                             =  gsn_csm_pres_hgt(wks_p,hgt_ano1(:,{plat},{lat_pmin:lat_pmax}),res_p)  ; 多暴雨年,作为底层的图决定了图的坐标，所以必须加上坐标范围
  base_p(0)                             =  ColorNegDashZeroPosContour(base_p(0),"blue","white","red")    
  
  base_p(1)                             =  gsn_csm_pres_hgt(wks_p,hgt_ano2(:,{plat},{lat_pmin:lat_pmax}),res_p)  ; 少暴雨年
  base_p(1)                             =  ColorNegDashZeroPosContour(base_p(1),"blue","white","red")    

;2---Theta_se（8层）  
  res_p@cnFillOn                        = True 
  res_p@cnLinesOn                       = False                                                     ; 绘制等值线
  res_p@cnLevelSelectionMode            = "ExplicitLevels"
  res_p@cnLevels                       := (/-1.5,-1.2,-0.9,-0.6,-0.3,0.3,0.6,0.9,1.2,1.5/)  
  res_p@cnFillColors                   := (/4,5,6,8,9,0,13,14,15,16,17/)  
  res_p@lbLabelBarOn                    = True                                                      ; 画图标，统一绘制
  res_p@cnLineLabelsOn                  = False  
  res_p@cnLineDrawOrder                 = "PreDraw"                                                 ;"PostDraw" 最后画（在最上面),会覆盖 the predraw and draw "PreDraw" :先画 就是放在最下面 包括地图下面;"Draw"  ;会覆盖the predraw 但会在the postdraw下面   
  thse_p1                               =  gsn_csm_pres_hgt(wks_p,eqt_ano1(:,{plat},:),res_p)       ; 多暴雨年
  overlay(base_p(0),thse_p1)   
  
  thse_p2                               =  gsn_csm_pres_hgt(wks_p,eqt_ano2(:,{plat},:),res_p)       ; 少暴雨年
  overlay(base_p(1),thse_p2)  
  
print("theta_se")
;3---omega（12层）  
  res_p@cnLinesOn                       = False
  res_p@cnFillOn                        = False
  res_p@cnLineLabelsOn                  = False

;箭头设置--------------------------------------------------------------------------------------------------
  res_p@vcGlyphStyle                    = "LineArrow"                                               ; 其余三种选项为“LineArrow”（默认值）、“WindBarb” 、“CurlyVector”
  res_p@vcMinDistanceF                  = 0.03                                                      ; 箭头之间的最小距离（在单位平方中）
  res_p@vcFillArrowEdgeThicknessF       = 1                                                         ; 箭头边界粗细
  res_p@vcFillArrowEdgeColor            = "white"                                                   ; 及颜色
  res_p@vcFillArrowFillColor            = "darkslateblue"                                           ; 箭头内部填充颜色
  res_p@vcFillArrowWidthF               = 0.07                                                      ; 箭头宽度
  res_p@vcFillArrowHeadXF               = 0.4                                                       ; 请参考附录中Fillarrow箭头示意图
  res_p@vcFillArrowHeadYF               = 0.10                                                      ;
  res_p@vcFillArrowHeadInteriorXF       = 0.20                                                      ;            
print("1919")    
  res_p@vcFillArrowMinFracWidthF        = 1.0 
  res_p@vcFillArrowHeadMinFracXF        = 1.0  
  res_p@vcFillArrowHeadMinFracYF        = 1.0 
print("1923")  
;0.5.2****设定参考箭头****
  res_p@vcRefAnnoOn                     =  True  
  res_p@vcRefLengthF                    =  0.040                                                    ; 标准长度箭头在单位平方中的大小
  res_p@vcRefAnnoBackgroundColor        =  -1                                                       ; 背景颜色，透明     
  res_p@vcRefAnnoPerimOn                =  False                                                    ; 关闭边框                                           
  res_p@vcRefAnnoFontHeightF            =  0.012                                                    ; 参考箭头标签字体大小         
  res_p@vcRefAnnoString1On              =  False                                                    ; 设定参考箭头上的字符        
  res_p@vcRefAnnoString2On              =  True                                                     ; 这里仅设定其下方的字符           
  res_p@vcRefAnnoSide                   =  "Top"                                                    ; 参考箭头放至图形上方
  res_p@vcRefAnnoOrthogonalPosF         =  0                                                        ; 调整箭头垂直位置(>0为向上)
  res_p@vcRefAnnoParallelPosF           =  1.00                                                     ; 调整箭头水平位置
  res_p@vcRefAnnoArrowSpaceF            =  1.3                                                      ; 调整箭头和文本间的距离 
  res_p@vcRefMagnitudeF                 =  2                                                        ; 标准长度箭头所表示的大小
  res_p@vcMinMagnitudeF                 =  0.1                                                      ; 要绘制箭头所表示的最小数值，即小于该数值则不绘制     
  res_p@vcRefAnnoString2On              =  True                                                     ; 这里仅设定其下方的字符
  res_p@vcRefAnnoString2                =  2+" m/s"                                                 ; 参考箭头大小    

;----- 不显著的风场--灰色   
  res_p@gsnLeftString                   = "(a)Ano Hgt,Theta_se,Wind"                                ; 子图的主标题 
  res_p@vcLineArrowColor                =  "gray50"                                                 ; 不显著的风场--灰色 
  
  uo1                                   =  gsn_csm_pres_hgt_vector(wks_p,u_ano1_pbn({1000:100},{plat},{80:140}),u_ano1_pbn({1000:100},{plat},{80:140}),omega_ano1_pbn({1000:100},{plat},{80:140})*-100,res_p) 
  overlay(base_p(0),uo1)    
;----- 显著的风场--黑色
  res_p@vcLineArrowColor                =  "blue"                                                  ; 显著的风场--黑色
  uo11                                  =  gsn_csm_pres_hgt_vector(wks_p,u_ano1_pb({1000:100},{plat},{80:140}),u_ano1_pb({1000:100},{plat},{80:140}),omega_ano1_pb({1000:100},{plat},{80:140})*-100,res_p) ;垂直速度量级和水平差100倍，且转换为负，箭头向上
  overlay(base_p(0),uo11)    

  res_p@gsnLeftString                   = "(b)Ano Hgt,Theta_se,Wind"                                ; 子图的主标题     
  res_p@vcLineArrowColor                =  "gray50"                                                 ; 不显著的风场--灰色 

;----- 不显著的风场--灰色   
  uo2                                   =  gsn_csm_pres_hgt_vector(wks_p,u_ano2_pbn({1000:100},{plat},{80:140}),u_ano1_pbn({1000:100},{plat},{80:140}),omega_ano2_pbn({1000:100},{plat},{80:140})*-100,res_p) 
  overlay(base_p(1),uo2)    
;----- 显著的风场--黑色        
  res_p@vcLineArrowColor                =  "blue"                                                  ; 显著的风场--黑色  
  uo22                                  =  gsn_csm_pres_hgt_vector(wks_p,u_ano2_pb({1000:100},{plat},{80:140}),u_ano2_pb({1000:100},{plat},{80:140}),omega_ano2_pb({1000:100},{plat},{80:140})*-100,res_p) 
  overlay(base_p(1),uo22)        
print("1959")  
;;组合成一张图
  resP                                  = True                                                      ; 绘制组图 
  resP@gsnPanelRowSpec                  = True                                                      ; 指定每行绘制的子图的个数
  resP@gsnPanelCenter                   = True                                                      ; 子图居中
  resP@gsnPanelBottom                   = 0.05                                                      ; add some space at bottom 
  resP@gsnPanelXWhiteSpacePercent       = 1.1                                                       ; 两个相邻子图之间沿X轴方向的空白大小  
print("1970")  
  gsn_panel(wks_p,base_p(0:1),(/1,1/),resP)           
    
print("hgt-pres")

    print("the positive years of larger than 1.0 are "+stdy1)
    print("the negetive years of smaller than -1.0 are "+stdy2) 
    
  endtime            =    get_cpu_time()
  print("   runing time: " + (endtime-startime)+"seconds")   
 
end     

;--------------------------读出到文件中------------------------------------------------------------------------------------- 
;###########################################################################################################################################

     
