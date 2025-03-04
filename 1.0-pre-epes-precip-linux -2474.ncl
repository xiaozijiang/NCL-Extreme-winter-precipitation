;####文件说明：数据的前处理，目的：节省重复运算的时间。具体做法：

; 1. 年限61-20年冬季，1961年12月-1962年2月算作是1961年冬季；

;  2. 所有冬季数据缺失超过10%的站点去除，筛选出2144个站。

;输出2个文件：

;（1）站点信息文件：data/1-Interannual/1-f_station2144.txt；

;（2）冬季降水数据：data/1-Interannual/1-f_winterain2144_61-20.txt；

;#####################################################################################################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/cnmap/cnmap.ncl"

begin

startime            =    get_cpu_time()
;########################START##############一.工作站、变量定义#######################################################################################    

   datadir        =   "/home/jupiter/jiangsh/"
;   datadir        =   "/home/tom/"

   
    syear           =   1961                                                                        ; 实际开始年份,虽然1958年很多站点不存在，但1958有极强冬季暴雨，所以从1958开始。
    eyear           =   2020                                                                        ; 实际结束年份  
    yer             =   ispan(syear,eyear,1)  
    sey             =   dimsizes(yer)     
       
    f               =   new(55285478,string)                                                        ; 读取1961-2020年站点数据的行数    
   
;#######################START############### 二.数据读取###########################################################################################################################################    

;--2.1; read  station-lat-lon  data 2479站点经纬度数据-------------------------------------------------------------------------------------
    name_station  =   asciiread(datadir+"data/2474-1951-2020-pre/NAME_STATION_2479.txt",-1,"string")   
    delim         =   str_get_tab()    
    name_num      =   stringtoint(str_get_field(name_station,1,delim))
;    name_lat      =   str_get_field(name_station,2,delim)
;    name_lon      =   stringtoint(str_get_field(name_station,3,delim))
    name_id       =   stringtoint(str_get_field(name_station,4,delim))                           ; 各省的ID
    name_wgt      =   stringtofloat(str_get_field(name_station,6,delim))                         ; 各省的面积权重
        
;    name_lat     :=   stringtoint(str_get_cols(name_lat,0,1))+ stringtoint(str_get_cols(name_lat,2,3))/60.0  ;直接转换为浮点型
;    name_lon     :=   floor(name_lon/100)+(name_lon%100)/60.0                                    ; 默认转换为浮点型  
        
;--2.2; read pricipation station  data ----1951-2019年冬季----------------------------------------------------------------------------------------------------

    diri            =   datadir+"data/2474-1951-2020-pre/PRE/"
    fils            =   systemfunc ("ls -1 "+diri+"SURF_CLI_CHN_MUL_DAY-PRE-13011-*.TXT")           ; 文件名的集合;所有路径名
    numfs           =   dimsizes(fils)
    print(numfs)                                                                                    ; 文件数 

    nn              =   0                                                                           ; 每个文件数据行数       
    do  i           =   0 , numfs - 1                                                               ; :=再赋值，改变数组大小，确定站点
        lins       :=   asciiread(fils(i),-1,"string")                                              ; 行的形式读取一个文件,lines为n*1维的数组 
        f(nn:nn+dimsizes(lins)-1)      =   lins                                                     ; 所有数据整合为一个文件。
        nn          =   nn+dimsizes(lins) 
    end do
    print(nn)
    print(f(nn-1))

    delim           =   " "                                                                         ; 分隔符：空格  
    sta_num         =   stringtoint(str_get_field(f,1,delim))  
    sta_lat         =   stringtoint(str_get_field(f,2,delim))
    sta_lon         =   stringtoint(str_get_field(f,3,delim))
    sta_high        =   stringtoint(str_get_field(f,4,delim))    
    sta_year        =   stringtoint(str_get_field(f,5,delim))  
    sta_mon         =   stringtoint(str_get_field(f,6,delim)) 
    sta_day         =   stringtoint(str_get_field(f,7,delim))    
    sta_rain        =   stringtoint(str_get_field(f,10,delim)) 
    
    delete(f)
;#######################START############### 三.数据处理###########################################################################################################################################    
  
;--3.1--------- ;经纬度处理--------------------------------------------------------------------
    sta_lat        :=   floor(sta_lat/100)+(sta_lat%100)/60.0                                       ; 直接转换为浮点型，弧度
    sta_lon        :=   floor(sta_lon/100)+(sta_lon%100)/60.0    
         
;--3.2--------- ;雨量处理--------------------------------------------------------------------
    sta_rain@_FillValue=-9999 
    sta_rain       :=   where(sta_rain.eq.32766,sta_rain@_FillValue,sta_rain)                       ; 32766为缺失值
    sta_rain       :=   where((sta_rain.ge.30000).and.(sta_rain.le.30999),sta_rain-30000,sta_rain)  ; 雪量（雨夹雪，雪暴）
    sta_rain       :=   where((sta_rain.ge.31000).and.(sta_rain.le.31999),sta_rain-31000,sta_rain)  ; 雨和雪的总量
    sta_rain       :=   where((sta_rain.ge.32000).and.(sta_rain.le.32699),sta_rain-32000,sta_rain)  ; 纯雾露霜
    sta_rain       :=   where((sta_rain.ge.32700),0,sta_rain/10.0)                                  ; 微量降水
     
;--3.3--------- ;所有的数据转化为冬季(1958-2019)------------------------------------------------------------
    sta_mon        :=   where(sta_mon.eq.12,0,sta_mon)                                              ; 将12月份转变为0月份
    sta_year       :=   where((sta_mon.eq.1).or.(sta_mon.eq.2),sta_year-1,sta_year)                 ; 将1月和2月的年份转换为前一年
    
    indys           =   ind((sta_year.ge.syear).and.(sta_year.le.eyear))                            ; 只计算1961-2020年       
    indw            =   ind(((sta_mon(indys)).eq.0).or.((sta_mon(indys)).eq.1).or.((sta_mon(indys)).eq.2)) ; 只计算冬季  

    stanum         :=   sta_num(indys(indw))
    stalat         :=   sta_lat(indys(indw))
    stalon         :=   sta_lon(indys(indw))  
;    stahigh        :=   sta_high(indys(indw)) 
    stayear        :=   sta_year(indys(indw))   
    stamon         :=   sta_mon(indys(indw))       
    staday         :=   sta_day(indys(indw))
    starain        :=   sta_rain(indys(indw))
          
    delete(sta_num)
    delete(sta_lat)
    delete(sta_lon) 
;    delete(sta_high)       
    delete(sta_year)    
    delete(sta_mon)
    delete(sta_day)
    delete(sta_rain) 
 
    print("winter data")
;#######################START############### 四.数据缺测处理###########################################################################################################################################    
    
;--4.1------START------1951-2019年每年冬季有多少非缺测站station_y------------------------------------------------    
    stn_2479        =   2479                                                                        ; 全国所有站点数，用于计算：每年非缺测站 
    stsn_2479       =   ispan(1,2479,1) 
    
    years           =   ispan(syear,eyear,1)
    station_y       =   new(sey,float)                                                              ; 每年有多少个站点，没缺测
    station_y!0     =   "years"
    station_y&years =   years  
    
    do  i           =    0   ,   eyear-syear                                                        ; 某年

        if (any(stayear.eq.(i+syear)))  then                                                        ; 判断此年是否存在
            
            ind_y  :=    ind(stayear.eq.(i+syear))                                                  ; 某一年   
            
            stas   :=    0  
            
            do  j   =    0   ,    stn_2479-1
                            
                if(any(stanum(ind_y).eq.name_num(j))) then
                
                   ind_yn :=   ind(stanum(ind_y).eq.name_num(j))                                    ; 某个站
                
                  if (any(starain(ind_y(ind_yn)).ge.0.0)) then                                      ; 判断此年是否有站点，说明有站
                
                   stas = stas+1                                                                    ; 此站存在就+1
                   
                  end if
                  
                end if
                
            end do
            
            station_y(i)= stas
            
        end if
        
        print(i)
            
    end do
         
    print(station_y)                                                                               ; 每年有多少非缺测站；；没问题
    
; 1951	1952	1953	1954	1955	1956	1957	1958	1959	1960	1961	1962	1963	1964	1965	1966	1967	1968	1969	1970	1971	1972	1973	1974	1975	1976	1977	1978	1979	1980	1981	1982	1983	1984	1985	1986	1987	1988	1989	1990	1991	1992	1993	1994	1995	1996	1997	1998	1999	2000	2001	2002	2003	2004	2005	2006	2007	2008	2009	2010	2011	2012	2013	2014	2015	2016	2017	2018	2019
; 188	  304	  387	  473	  582	  1065	1246	1730	1997	2077	2121	2106	2128	2163	2184	2213	2205	2215	2239	2279	2312	2327	2351	2362	2369	2378	2394	2402	2407	2414	2415	2415	2416	2415	2413	2413	2409	2406	2401	2394	2392	2385	2382	2371	2364	2363	2362	2388	2387	2387	2388	2396	2398	2398	2403	2415	2417	2417	2417	2421	2423	2422	2423	2422	2421	2421	2420	2420	2437  
;--4.1------OVER------1951-2019年每年冬季有多少非缺测站station_y------------------------          

;--4.2------START-----筛选站点：冬季缺测值要小于10%（最终从2479站中筛选出了2078站）------------------------------------------------------------                                                 
    stn_90          =   2144                                                                        ; 全国所有站点数（90%无缺测），用于计算：全国暴雨频次  
    stsn_90         =   ispan(1,stn_90,1) 
    st_num          =   new(stn_90,integer)                                                         ; 全国有效站点数1958-2019年冬季，总共2078个  
        
    stans           =    0
    
    do   i          =    0   ,   stn_2479-1
        
         if(any(stanum.eq.name_num(i))) then      
                
           ind_n   :=    ind(stanum.eq.name_num(i))                                                 ; 某个站
           
           stan     =    dimsizes(ind_n)                                                            ; 每个站有多少天，包括缺测值，但不包括无缺测值且无数据的站           
                
           if (any(starain(ind_n).ge.0.0)) then                                                     ; 非缺测的天数（雨量≥0.0）
                                             
               ind_n00 := ind(starain(ind_n).ge.0.0)                                                ; 非缺测的天数
                          
               n00  =    dimsizes(stanum(ind_n(ind_n00)))
               
               n10  =    n00/(sey*90+15.0)*1.0                                                      ; 1961-2020年冬季有多少天数，15为15天的闰年，15个闰年
               
               if(n10.ge.0.9)  then                                                                 ; sey为年份，90为冬季站数
                 
                 st_num(stans) = stanum(ind_n(0))                                                   ; 哪些站点符合缺测少于10%
                                               
                 stans= stans+1      

               end if                                                              
               
           end if
           
         end if
         
         print(stans)
       
    end do
  
    print(st_num)      
                                                                                 
                                                                                                                                                          
;--4.3.1--------- 处理筛选出的站点(将数据转换为2144站数据)-建立2144站站点信息文件------------------------------------------------------------ 

    idn             =   new(stn_90,integer) 
    sn              =    0
       
    do   i          =    0    ,   stn_90-1
           
        if(any(name_num.eq.st_num(i)))   then                                                       ; 筛选出2078站的基础数据（站号、id和面积权重）
        
          idn(sn)   =    ind(name_num.eq.st_num(i))           
          
          sn        =    sn + 1     
                                
        end if     
                
    end do  
   
    name_num9       =    name_num(idn)
    name_id9        =    name_id(idn)
    name_wgt9       =    name_wgt(idn)
     
;--4.3.2--------- 处理筛选出的站点(将数据转换为2144站数据)-所有冬季降水2078站数据------------------------------------------------------------
   name_lat9        =   new(stn_90,float)
   name_lon9        =   new(stn_90,float)

   nx               =   11558754

   stanum9          =   new(nx,integer)                                                       ; 2078数据的行数  
;   stahigh9         =   new(nx,float)
   stayear9         =   new(nx,integer)         
   stamon9          =   new(nx,integer)         
   staday9          =   new(nx,integer)         
   starain9         =   new(nx,float)         
   stalat9          =   new(nx,float)         
   stalon9          =   new(nx,float)    
   
   sx               =   0                                                                           ; 计数器  
            
   do    i          =   0  ,   stn_90-1
     
         if(any(stanum.eq.name_num9(i))) then
         
           indx    :=    ind(stanum.eq.name_num9(i))
           
           n        =    dimsizes(indx)
           
           stanum9(sx:sx+n-1)  = stanum(indx)
           
;           stahigh9(sx:sx+n-1) = stahigh(indx) 
          
           stayear9(sx:sx+n-1) = stayear(indx)
           
           stamon9(sx:sx+n-1)  = stamon(indx)
           
           staday9(sx:sx+n-1)  = staday(indx)
         
           starain9(sx:sx+n-1) = starain(indx)
           
           stalat9(sx:sx+n-1)  = stalat(indx)
           name_lat9(i)        = stalat(indx(0))                                                   ; 2078站站点经纬度信息    
                
           stalon9(sx:sx+n-1)  = stalon(indx)
           name_lon9(i)        = stalon(indx(0))
           
           sx       =    sx + n
           
           print(sx)
           print(n)
         
         end if
                  
   end do   
                                                                  
print("creat file")   

   alist           :=    [/stanum9,stalat9,stalon9,stayear9,stamon9,staday9,starain9/]     
   file_name       :=    datadir+"data/1-Interannual/1-f_winterain2144_61-20.txt"                ;建立冬季降水的文件
   system("/bin/rm -f "+file_name)                                                                  ;移除存在的文件                                                          
   write_table(file_name,"w",alist,"%d,%8.3f,%8.3f,%d,%d,%d,%8.1f")                
            
    delete(stanum)
    delete(stalat)
    delete(stalon)  
;    delete(stahigh)     
    delete(stayear)    
    delete(stamon)
    delete(staday)
    delete(starain)
    
;--4.3.3--------- 处理筛选出的站点(将数据转换为2144站数据)-建立2144站站点信息文件------------------------------------------------------------  

   alist            =    [/name_num9,name_lat9,name_lon9,name_id9,name_wgt9/]     
   file_name        =    datadir+"data/1-Interannual/1-f_station2144.txt"                        ;
   system("/bin/rm -f "+file_name)                                                                  ;移除存在的文件                                                          
   write_table(file_name,"w",alist,"%d,%8.3f,%8.3f,%d,%8.6f")           

   print("Game is Over!")  
   
  endtime            =    get_cpu_time()
  print("   runing time: " + (endtime-startime)/60.0+"minutes")    
    
end