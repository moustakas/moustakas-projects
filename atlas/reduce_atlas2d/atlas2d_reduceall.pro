pro atlas2d_reduceall, zptshift=zptshift, run=run, doplot=doplot
; jm02jan5uofa
; batch reduce all the runs, including SINGS
; RUN should be a string name
    
    if n_elements(zptshift) eq 0L then $
      zptshift = strtrim(string(0.11),2) else $
      zptshift = strtrim(string(zptshift),2)
    
    dpath = '/mount/moon1/ioannis/spectral_atlas/reduce_atlas2d/'
;   dpath = '/global/bias1/ioannis/spec2datlas/'
;   dpath = '/d1/ioannis/spec2datlas/'
;   dpath = '/d0/ioannis/kennicutt/data/'
;   dpath = '/home/ioannis/kennicutt/data/'

    subpath = [$
      '94nov/',$ ; Turner
      '95mar/',$ ; Turner
      '95oct/',$ ; Turner
      '96apr/',$ ; Turner
      '97apr/',$ ; Turner
      '98mar/',$
      '98apr/',$
      '98jun/',$
      '98oct/',$
      '99apr/',$
      '99may/',$
      '99nov/',$
      '00apr/',$
      '01nov/',$
      '01dec/',$
      '02feb/',$
      '02apr/',$
      '02may/',$
;     '02jun/',$
      '03may/',$
;     '04mar/',$
      '05apr/',$
      '06mar/',$
      '06may/']
    path = dpath+subpath
    scripts = [$
      'nov94_script',$
      'mar95_script',$
      'oct95_script',$
      'apr96_script',$
      'apr97_script',$
      'mar98_script',$
      'apr98_script',$
      'jun98_script',$
      'oct98_script',$
      'apr99_script',$
      'may99_script',$
      'nov99_script',$
      'apr00_script',$
      'nov01_script',$
      'dec01_script',$
      'feb02_script',$
      'apr02_script',$
      'may02_script',$
;     'jun02_script',$
      'may03_script',$
;     'mar04_script',$
      'apr05_script',$
      'mar06_script',$
      'may06_script']
    npath = n_elements(path)

    nrun = n_elements(run)
    if (nrun ne 0L) then begin

       indx = lonarr(nrun)
       for k = 0L, nrun-1L do indx[k] = where(strmatch(subpath,'*'+run[k]+'*') eq 1B)

       minusone = where(indx eq -1L,nminus)
       if nminus gt 0L then begin
          splog, 'RUN name not found.'
          return
       endif

    endif else indx = lindgen(npath)
    nreduce = n_elements(indx)

    if (n_elements(doplot) eq 0L) then doplot = '0' else doplot = strn(doplot)

    splog, file=dpath+'atlas2d_reduceall.log'
    t0 = systime(/seconds)
    for i = 0L, nreduce-1L do begin

       j = indx[i]
       
       splog, 'Changing directory to ', path[j]
       pushd, path[j]
       
       runit = execute(scripts[j]+', doplot='+doplot+', zptshift='+zptshift)

       popd

    endfor
    splog, 'Total time elapsed (minutes) = ', (systime(/seconds)-t0)/60.0
    splog, /close
    
return
end
