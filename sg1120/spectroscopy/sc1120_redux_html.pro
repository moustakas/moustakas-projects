pro sc1120_redux_html, _extra=extra
; jm05jan17uofa
; generate a webpage with the QAPLOT figures from each of the 4
; quadrants 
    
    html_path = sc1120_path(/web)

    quad = ['Q1','Q2','Q3','Q4']
    quadpath = [sc1120_path(/Q1),sc1120_path(/Q2),sc1120_path(/Q3),sc1120_path(/Q4)]

; for each quadrant, move the postscript files of interest into the
; public_html redux directory, then generate the web page   

    for iquad = 0L, n_elements(quadpath)-1L do begin
       
       htmlbase = 'redux_'+quad[iquad]
       spawn, ['/bin/rm '+html_path+htmlbase+'/*'], /sh
       spawn, ['/bin/cp '+quadpath[iquad]+'qaplot*.ps '+html_path+htmlbase], /sh

       pushd, html_path+htmlbase
       psfiles = file_search('qaplot*.ps',count=n_psfiles)
       popd
       
; sort the PSFILES into iSPEC2d order

       lastindx = lindgen(n_psfiles)
       psbias = where(strmatch(psfiles,'qaplot*bias*') eq 1B,npsbias)
       psresp = where(strmatch(psfiles,'qaplot*response*') eq 1B,npsresp)
       pslamp = where(strmatch(psfiles,'qaplot*arc*') eq 1B,npslamp)
       pssens = where(strmatch(psfiles,'qaplot*sens*') eq 1B,npssens)
       pstell = where(strmatch(psfiles,'qaplot*telluric*') eq 1B,npstell)

       newpsfiles = ''
       if (npsbias ne 0L) then newpsfiles = [newpsfiles,psfiles[psbias]]
       if (npsresp ne 0L) then newpsfiles = [newpsfiles,psfiles[psresp]]
       if (npslamp ne 0L) then newpsfiles = [newpsfiles,psfiles[pslamp]]
       if (npssens ne 0L) then newpsfiles = [newpsfiles,psfiles[pssens]]
       if (npstell ne 0L) then newpsfiles = [newpsfiles,psfiles[pstell]]

; subscript any additional files

       if (n_elements(newpsfiles)-1L lt n_elements(psfiles)) then begin

          remove, [psbias,psresp,pslamp,pssens,pstell], lastindx
          newpsfiles = [newpsfiles,psfiles[lastindx]]
          
       endif

       psfiles = newpsfiles[1L:n_elements(newpsfiles)-1L]
       im_ps2html, htmlbase, html_path=html_path, pslist=psfiles, cleanpng=0, $
         npscols=3, _extra=extra

    endfor
       
    stop
    
return
end
    
