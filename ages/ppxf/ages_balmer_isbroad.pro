function ages_balmer_isbroad, specdata, justpass=justpass
; jm09dec14ucsd - return the list of objects targetted as galaxies
;   (GSHORT>0, 0.001<z<1.0; see AGES_GET_ZABS_VDISP()) that may
;   contain broad Balmer line components; these objects must be fitted
;   with a different set of parameters in AGES_GANDALF_SPECFIT
;
; note that this list is *not* exhaustive nor rigorously derived!
;
;   ss = read_ages_gandalf(/ppxf)    
;   cc = iclassification(ss,ratios=rr,sigmacut=1000)
;   ww = where(strtrim(rr.final_class,2) eq 'AGN'
;   data = ss[ww]  & specfit = read_ages_gandalf_specfit(data)
;   qaplot_ages_gandalf_specfit, data, specfit, psfile='qaplot_balmer_isbroad_v1.0.ps'

; original list based on v1.0 PPXF/GANDALF fitting:

    isbroad = [$
      '101/044',$ ; 1200 km/s
      '101/164',$ ; 1200 km/s
      '102/229',$ ; broad Mg II (z=0.426)
      '103/012',$
      '105/228',$ ; 600 km/s
      '107/034',$
      '107/092',$
      '107/117',$ ; 2200 km/s (H-beta, z=0.69)
      '107/236',$ ; broad H-beta (Mg II on the edge) (z=0.354)
      '108/046',$ ; 2200 km/s broad!
      '109/176',$ ; 740 km/s (z=0.318)
      '111/015',$ ; 600 km/s (z=0.085)
      '111/029',$ ; 4800 km/s (z=0.677), broad Mg II
      '111/034',$
      '111/126',$
      '111/137',$
      '112/043',$ ; 465 km/s (z=0.321)
      '112/246',$ ; 7400 km/s broad!!
      '113/208',$ ; 8300 km/s (z=0.160)
      '113/209',$
      '114/004',$ ; gorgeous broad-line AGN!
      '115/083',$ ; 800 km/s (z=0.0448)
      '115/101',$ ; 900 km/s (z=0.121)
;     '115/104',$ ; not needed
      '202/064',$ ; 1200 km/s (z=0.068)
      '206/177',$ ; 740 km/s (z=0.174)
      '207/240',$ ; 1200 km/s (z=0.124)
;     '208/298',$ ; not needed
      '210/066',$ ; 1000 km/s (z=0.172)
      '210/167',$
;     '211/126',$ ; not needed
      '211/078',$ ; broad Mg II
      '211/294',$ ; 670 km/s (z=0.340)
      '302/141',$ ; broad Mg II
      '302/209',$ ; broad Mg II
      '302/283',$ ; broad Mg II
      '303/008',$ ; broad Mg II
      '303/044',$ ; broad Mg II
      '303/270',$ ; broad Mg II (z=0.662)
      '304/181',$ ; 1700 km/s broad!
      '304/290',$ ; marginal
      '308/199',$
      '308/270',$ ; broad Mg II (but low S/N)
      '309/257',$ 
      '312/132',$ ; broad Mg II
      '312/227',$ ; broad Mg II
      '313/095',$
      '313/173',$ ; broad Mg II
      '313/280',$
      '314/008',$ ; broad Mg II
      '314/030',$ ; broad Mg II
      '314/051',$ ; marginal (H-beta ~190 km/s)
      '314/146',$
      '314/154',$
      '314/284',$ ; broad Mg II, H-beta
      '315/102',$ ; broad Mg II
      '315/123',$
      '315/133',$ ; broad Mg II
      '315/149',$ ; broad Mg II
      '403/236',$ ; broad Mg II
      '404/045',$ ; broad Mg II
      '404/072',$ 
      '404/155',$ 
      '404/187',$ ; broad Mg II
      '406/280',$
      '408/093',$ ; broad Mg II
      '409/191',$ ; broad Mg II(?), but low S/N
      '410/271',$
      '416/033',$ ; not needed
;     '416/153',$ ; not needed
      '416/215',$
      '416/275',$ ; broad Mg II, H-beta
      '417/105',$ ; broad Mg II (z=0.762)
      '418/060',$ ; broad Mg II
      '418/057',$ ; maybe not?
      '418/154',$ ; broad(?) Mg II (significant?) (z=0.785)
      '419/202',$
;     '419/211',$ ; not needed (maybe?)
      '419/213',$ ; broad Mg II, but low S/N
      '421/112',$ ; broad Mg II
      '421/150',$
      '422/083',$ ; broad Mg II
      '422/242',$ ; broad Mg II
      '502/242',$ ; broad Mg II, but low S/N
      '503/104',$ ; broad Mg II, but low S/N
      '506/295',$ ; broad Mg II, H-beta
      '507/262',$
      '508/077',$ ; broad H-beta and Mg II (z=0.426)
      '508/261',$ ; broad Mg II
      '510/207',$ ; broad Mg II
      '510/281',$ ; broad Mg II
      '526/290',$ ; broad Mg II
      '602/238',$ ; broad Mg II
      '602/296',$ ; broad Mg II
      '603/213',$ ; broad Mg II
      '606/245',$ ; broad Mg II
      '606/251',$ ; broad Mg II
      '607/220',$ ; broad Mg II
      '608/036',$ ; broad Mg II
      '608/253',$ ; broad Mg II
      '611/167',$ ; broad Mg II
      '709/036',$
      '709/094',$
      '709/281',$ ; broad Mg II
      '715/136',$
      '715/153',$
      '722/144']

; just return the PASS number (see AGES_GANDALF_SPECFIT)    
    if keyword_set(justpass) then begin
       pass = strmid(isbroad,0,3)
       pass = pass[uniq(pass,sort(pass))]
       return, pass
    endif
    
    all = string(specdata.pass,format='(I3.3)')+'/'+$
      string(specdata.aper,format='(I3.3)')
    indx = where_array(isbroad,all)
;   if (indx[0] ne -1) then out = specdata[indx] else out = -1
    
return, indx
end
    
