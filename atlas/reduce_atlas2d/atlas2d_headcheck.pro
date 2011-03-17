pro atlas2d_headcheck, run=run

    dpath = '/home/ioannis/kennicutt/data/'
    subpath = [$
      '98mar/','98apr/',$
      '98jun/','98oct/',$
      '99apr/','99nov/',$
      '00apr/','01nov/',$
      '01dec/','02feb/',$
      '02apr/','02may/',$
      '02jun/']
    path = dpath+subpath
    scripts = [$
      'mar98_script','apr98_script',$
      'jun98_script','oct98_script',$
      'apr99_script','nov99_script',$
      'apr00_script','nov01_script',$
      'dec01_script','feb02_script',$
      'apr02_script','may02_script',$
      'jun02_script']
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

    splog, file='temp.log'
    for i = 0L, nreduce-1L do begin

       j = indx[i]
       
       splog, 'Changing directory to ', path[j]
       pushd, path[j]

       flist = findfile('a.*.fits',count=fcount)
       for k = 0L, fcount-1L do begin
          h = headfits(flist[k])
          splog, flist[k]
;         if stregex(h,'photometric',/boolean) eq 1B then splog, flist[k]
          hgrep, h, 'photometric'
       endfor
       
       popd

    endfor
    splog, /close
    
stop    
    
return
end
