function ages_isnev, specdata, justpass=justpass
; jm10dec12ucsd - return the (visually inspected) list of objects with
; significant [NeV] detections 

;  data = read_ages_gandalf(/solar)    
;  ww = where(data.nev_3426[0]/data.nev_3426[1] gt 1.0)
;  qaplot_ages_gandalf_specfit, data[ww], specfit, /solar, /nev, $
;    psfile=ages_path(/ppxf)+'qaplot_isnev_v2.0.ps'
;  niceprint, strtrim(data[ww].pass,2)+'/'+string(data[ww].aper,format='(I3.3)')
   
; original list based on v2.0 PPXF/GANDALF fitting; the full list
; includes all the objects that were examined, with the false
; detections commented out
    isnev = [$
      '101/121',$
      '101/164',$
;     '101/175',$
;     '101/295',$
;     '102/004',$
;     '102/016',$
      '102/030',$
;     '102/119',$ ; strong, but likely noise
;     '102/160',$
      '102/263',$ ; marginal
;     '102/268',$
;     '103/019',$
;     '103/032',$
;     '103/049',$
;     '103/168',$
;     '103/196',$
;     '103/207',$
      '103/252',$
;     '104/060',$
      '104/108',$ ; wow!
;     '104/230',$
;     '105/015',$
;     '105/045',$ ; strong, but right on the edge...probably noise
;     '105/098',$
;     '105/127',$
      '105/171',$
;     '107/021',$
      '107/034',$
;     '107/061',$
;     '107/122',$
;     '107/225',$
;     '107/299',$
;     '108/078',$
;     '108/171',$ ; strong, but right on the edge...probably noise
;     '109/077',$
      '111/029',$ ; wow!
;     '111/032',$ ; truncated CR
      '111/034',$ ; wow!
      '111/045',$
;     '111/182',$
;     '111/239',$
      '111/279',$ ; wow!
;     '112/114',$
;     '112/283',$
;     '113/008',$ ; maybe?
;     '113/024',$
;     '113/039',$ ; strong, but likely noise
      '114/004',$ ; wow!
;     '114/035',$ ; marginal
;     '114/056',$
;     '114/057',$
;     '114/106',$
;     '114/132',$
;     '114/273',$
;     '115/143',$
;     '115/163',$
;     '115/169',$
;     '115/185',$
;     '115/230',$
      '115/294',$
;     '201/014',$
;     '201/037',$
;     '201/046',$
;     '201/061',$
;     '201/180',$
      '201/183',$ ; strong!
;     '201/205',$
;     '201/215',$
;     '202/097',$ ; marginal
;     '202/181',$
;     '202/299',$
;     '203/069',$ ; marginal
;     '203/197',$
;     '203/272',$
;     '204/086',$
      '204/195',$ ; wow!
;     '205/001',$
;     '205/163',$
;     '205/231',$
;     '206/106',$
      '206/116',$
;     '206/175',$ ; likely noise
;     '206/244',$
;     '207/238',$
;     '207/240',$
;     '207/277',$
;     '208/184',$
;     '208/225',$
      '210/167',$ ; strong
;     '211/051',$
;     '211/185',$ ; CR
;     '212/028',$
;     '212/039',$
;     '212/187',$
;     '213/004',$ ; marginal
;     '213/024',$
;     '213/061',$
;     '213/266',$
;     '213/284',$
;     '214/278',$
;     '214/283',$
      '215/031',$ ; wow!
      '215/163',$
;     '215/237',$
;     '301/007',$
;     '302/078',$
;     '302/209',$
;     '303/011',$
;     '303/107',$
;     '303/176',$
;     '303/178',$
;     '303/209',$
      '303/241',$ 
;     '303/252',$
      '303/270',$ ; wow!
;     '304/157',$
      '304/222',$
;     '304/290',$ ; candidate, but not formally fitted
;     '306/116',$
;     '306/128',$
;     '306/237',$
;     '306/266',$
      '307/090',$ ; possibly noise?
      '307/145',$
;     '307/149',$
;     '307/227',$
;     '307/237',$
;     '308/022',$
;     '308/028',$
      '308/087',$
;     '308/098',$
      '309/057',$
;     '309/060',$
;     '309/180',$
;     '312/123',$
;     '312/132',$
;     '312/133',$
;     '312/142',$
;     '312/225',$
      '313/095',$
;     '313/159',$
      '313/172',$ ; noise?
;     '313/180',$
      '314/077',$ ; wow!
      '314/154',$
;     '314/259',$
      '314/284',$
      '314/285',$
;     '315/010',$
;     '315/034',$ ; noise?
;     '315/058',$
      '315/102',$
;     '315/118',$
      '315/123',$
;     '315/156',$
      '315/198',$
;     '315/248',$
;     '401/007',$
;     '401/021',$
;     '401/071',$ ; CR
;     '401/159',$
;     '401/185',$
;     '401/233',$
      '402/092',$ ; noise?
;     '402/109',$
;     '403/019',$
      '403/028',$
      '404/072',$
 ;    '404/255',$
;     '404/278',$ ; CR
;     '405/228',$
;     '406/012',$
      '406/049',$
;     '406/171',$
;     '406/186',$
      '406/200',$
;     '406/233',$
      '406/292',$
;     '407/165',$
;     '407/168',$
;     '408/048',$
;     '408/126',$
;     '408/141',$
;     '408/214',$
;     '410/064',$
;     '410/084',$
      '410/284',$ ; noise?
;     '411/132',$
;     '411/146',$
;     '416/075',$
;     '416/090',$
;     '416/157',$
;     '416/208',$
;     '416/274',$
;     '417/110',$
      '417/146',$
;     '417/153',$
;     '417/280',$
;     '418/030',$
;     '418/044',$
;     '418/052',$
      '418/060',$
;     '418/110',$
;     '418/145',$
;     '419/132',$
      '419/203',$ ; marginal
      '419/213',$ ; wow!
      '419/251',$
;     '420/019',$
      '420/091',$
      '420/178',$
;     '420/202',$
      '420/233',$ ; noise?
;     '420/242',$
;     '421/024',$
;     '421/048',$ ; CR
      '421/077',$
;     '421/085',$
;     '421/088',$
      '421/112',$ ; wow!
      '421/122',$
      '421/150',$
;     '421/197',$
      '421/210',$
      '421/246',$
;     '421/264',$ ; sky residuals
      '421/300',$
;     '423/292',$
;     '501/115',$
      '501/123',$
;     '501/159',$
      '501/171',$ ; noise? v. broad, but so is [OII]
      '501/267',$
      '502/031',$
;     '502/133',$
;     '502/202',$
      '502/242',$
;     '503/028',$
;     '503/238',$
      '504/039',$ ; no other lines! noise?
;     '504/181',$
;     '505/116',$
;     '505/180',$
;     '505/280',$
;     '506/009',$
;     '506/049',$
;     '506/145',$
;     '506/210',$
;     '506/277',$
;     '507/066',$
;     '507/136',$
;     '507/158',$
;     '507/298',$
;     '508/144',$
;     '508/160',$ ; CR
      '508/261',$ ; wow!
;     '510/112',$
;     '510/123',$
;     '510/138',$
;     '510/176',$
;     '510/177',$ ; CR
;     '510/184',$ ; CR
      '510/213',$ ; noise?
      '510/282',$
;     '511/041',$ ; marginal
      '511/122',$
;     '511/167',$
;     '511/245',$
;     '512/045',$
;     '512/202',$
;     '512/231',$
      '512/283',$ ; marginal
;     '524/121',$
;     '524/164',$
;     '525/135',$
      '526/028',$
;     '526/078',$
      '526/208',$
;     '526/214',$
;     '601/062',$
;     '601/176',$
;     '602/142',$
;     '602/211',$
      '602/238',$ ; wow!
;     '602/279',$ ; noise? on edge
;     '602/288',$
;     '603/003',$
;     '603/037',$ ; CR
;     '603/062',$
;     '603/213',$ ; noise?
;     '603/278',$
;     '604/130',$
;     '604/222',$ ; maybe, but telluric issues
;     '604/295',$
;     '605/034',$
;     '605/118',$
;     '605/179',$ ; sky-subtraction residuals
;     '605/189',$
;     '605/278',$
      '605/280',$
;     '606/082',$ ; CR
;     '606/086',$
      '606/142',$ ; crazy broad, but so are the other lines!
;     '606/251',$
      '606/271',$
;     '607/045',$
      '607/220',$
;     '607/222',$
;     '607/267',$
;     '608/129',$
;     '608/196',$
;     '609/018',$
;     '609/158',$
;     '609/167',$
;     '609/188',$
;     '609/209',$
;     '610/078',$ ; noise?
;     '611/058',$
;     '611/126',$
;     '611/204',$
;     '611/265',$
;     '612/019',$
;     '612/034',$
;     '612/075',$
;     '612/077',$
      '612/243',$
;     '613/261',$
;     '622/195',$ ; sky
;     '622/254',$ ; sky
;     '709/024',$
;     '709/201',$ ; marginal
;     '709/206',$
;     '710/129',$ ; marginal
;     '710/155',$
;     '710/270',$
      '712/028',$
;     '712/163',$
;     '712/166',$
;     '712/252',$
      '713/029',$
      '713/095',$
;     '713/175',$
;     '713/246',$
      '714/086',$
;     '714/260',$
      '715/037',$
;     '715/212',$
;     '715/220',$
;     '715/237',$
;     '715/242',$
;     '722/025',$
      '722/042',$
      '722/280']
    
    all = string(specdata.pass,format='(I3.3)')+'/'+$
      string(specdata.aper,format='(I3.3)')
    indx = where_array(isnev,all)
;   if (indx[0] ne -1) then out = specdata[indx] else out = -1
    
return, indx
end
    