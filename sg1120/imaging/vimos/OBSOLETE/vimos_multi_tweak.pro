pro vimos_multi_tweak, findstars=findstars, sg1120=sg1120, standards=standards, _extra=extra
; jm07jan10nyu    

    datapath = vimos_path()

    if keyword_set(sg1120) then begin
       
       splog, '###########################################################################'
       splog, 'SG112'
       splog, '###########################################################################'

       if keyword_set(findstars) then begin
          imagelist = [file_search(datapath+'03dec/sg1120/sra.sg1120*_*_Q?.fits'),$
            file_search(datapath+'06feb/sg1120/sra.sg1120*_*_Q?.fits')]
          xylistlist = repstr(imagelist,'.fits','.xy.fits')
          for ii = 0L, n_elements(imagelist)-1L do begin
             if file_test(xylistlist[ii],/regular) then rmfile, xylistlist[ii]
             spawn, 'fits2xy '+imagelist[ii]
          endfor
       endif

       imagelist = [file_search(datapath+'03dec/sg1120/sra.sg1120*_*_Q?.fits[0]'),$
         file_search(datapath+'06feb/sg1120/sra.sg1120*_*_Q?.fits[0]')]
       xylistlist = repstr(imagelist,'.fits','.xy.fits')

       an_multi_tweak, imagelist, xylistlist, siporder=siporder, outdir=datapath, $
         result=sg1120_result, maxshift=1000.0, _extra=extra
       mwrfits, sg1120_result, datapath+'vimos_multi_tweak.fits', /create

    endif
       
    if keyword_set(standards) then begin

       splog, '###########################################################################'
       splog, 'Standards'
       splog, '###########################################################################'

       if keyword_set(findstars) then begin
          imagelist = [file_search(datapath+'03dec/sg1120/sra.[PG,SA]*_*_Q?.fits'),$
            file_search(datapath+'06feb/sg1120/sra.[PG,SA]*_*_Q?.fits')]
          xylistlist = repstr(imagelist,'.fits','.xy.fits')
          for ii = 0L, n_elements(imagelist)-1L do begin
             if file_test(xylistlist[ii],/regular) then rmfile, xylistlist[ii]
             spawn, 'fits2xy '+imagelist[ii]
          endfor
       endif

       imagelist = [file_search(datapath+'03dec/sg1120/sra.[PG,SA]*_*_Q?.fits'),$
         file_search(datapath+'06feb/sg1120/sra.[PG,SA]*_*_Q?.fits')]
       xylistlist = repstr(imagelist,'.fits','.xy.fits')

       an_multi_tweak, imagelist, xylistlist, siporder=siporder, outdir=datapath, $
         result=standards_result, maxshift=1000.0, _extra=extra
       mwrfits, standards_result, datapath+'standards_multi_tweak.fits', /create
       
    endif
       
return
end
