pro lampwavecalib, paramfile, arc

    ibatch, paramfile, caliblist=arc, tracename='', sensname='', extfile='', /calibrate

return
end    

function return_arcfile, paramfile
    
    paramtxt = djs_readlines(paramfile,nhead=2)
    paramtxt = paramtxt[where(strcompress(paramtxt,/remove) ne '')] ; remove blank lines

    nlines = n_elements(paramtxt)
    params = paramtxt

    comment = where(strmatch(paramtxt,'*#*') eq 1B,ncomment) ; crop trailing comments
    if ncomment ne 0L then for i = 0L, ncomment-1L do params[comment[i]] = $
      strtrim(strmid(paramtxt[comment[i]],0,strpos(paramtxt[comment[i]],'#'))) 

    arcline = where(strmatch(params,'*ARCFILE*'))
    arcfile = (strsplit(params[arcline],' ',/extract))[1]
    
return, arcfile
end

pro atlas_lampres
; jm03apr29uofa
; measure the spectral resolution of all the arc lamps

    dpath = '/home/ioannis/kennicutt/data/'
    pushd, dpath
    
    subpath = [$
      '98mar/','98apr/',$
      '98jun/','98oct/',$
      '99apr/','99nov/',$
      '00apr/','01nov/',$
      '02feb/','02apr/',$
      '02may/']
    npath = n_elements(subpath)

    for k = 0L, npath-1L do begin

       splog, 'Changing to '+subpath[k] ; change directory
       cd, subpath[k]

       batchlist = findfile('ibatch*.txt',count=nbatch)
       
       for j = 0L, nbatch-1L do begin

          arc = return_arcfile(batchlist[j])

          if file_test(arc,/regular) ne 0L then begin

;            lampwavecalib, batchlist[j], arc
             if j eq 0L then arclist = arc else arclist = [arclist,arc]

          endif

       endfor
       
       arclist = subpath[k]+'w'+arclist

       if k eq 0L then bigarclist = arclist else bigarclist = [bigarclist,arclist]
       
       cd, '../'

    endfor
    popd

    narc = n_elements(bigarclist)
    nline = 10L

;   for i = 0L, 4L do begin
    for i = 0L, narc-1L do begin
    
       results = ilampres(bigarclist[i],nline=nline,nback=nback,npoly=npoly,datapath=dpath)
       if i eq 0L then begin
          linewave = results.linewave
          linewidth = results.linewidth
          linewidth_err = results.linewidth_err
       endif else begin
          linewave = [linewave,results.linewave]
          linewidth = [linewidth,results.linewidth]
          linewidth_err = [linewidth_err,results.linewidth_err]
       endelse
       
    endfor

    cmsave, file='atlas_lampres.idlsave', /all
    
    
    stop

return
end    
