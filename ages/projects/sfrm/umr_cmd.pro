pro umr_cmd
; jm10mar30ucsd - u-r CMD for Richard

    sample = read_sfrm_sample()
    zbins = sfrm_zbins(nzbins)
    
    psfile = 'umr_cmd.ps'
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.6, $
      height=3.0*[1,1,1]

    xrange = [-16,-25]
    yrange = [0.7,3.6]

    magaxis = im_array(-24.5,-16.2,0.1)
    
    for ii = 0, nzbins-1 do begin
       if odd(ii) then begin
          ytitle = ''
          ytickname = replicate(' ',10)
       endif else begin
          ytitle = textoidl('^{0.1}(u-r)')
          delvarx, ytickname 
       endelse

       if (ii lt 4) then begin
          xtitle = ''
          xtickname = replicate(' ',10)
       endif else begin
          xtitle = textoidl('M_{0.1r}')
          delvarx, xtickname 
       endelse
       
       these = where((sample.z gt zbins[ii].zlo) and $
         (sample.z lt zbins[ii].zup),nthese)
       xx = sample[these].ugriz_absmag[2]
       yy = sample[these].ugriz_absmag[0]-sample[these].ugriz_absmag[2]
       
       hogg_scatterplot, xx, yy, noerase=(ii gt 0), $
         /internal, /outliers, position=pos[*,ii], xsty=1, ysty=1, $
         yrange=yrange, xrange=xrange, xtickname=xtickname, $
         ytickname=ytickname, xtitle=xtitle, ytitle=ytitle
       djs_oplot, magaxis, poly(magaxis,[1.0,-0.08]), line=0
       legend, 'z='+strtrim(string(zbins[ii].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[ii].zup,format='(F12.2)'),2), $
         /right, /top, box=0, margin=0, charsize=1.5
    endfor

    im_plotconfig, /psclose

return
end
