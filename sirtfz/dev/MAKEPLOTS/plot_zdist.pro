pro plot_zdist, postscript=postscript
; jm01may13uofa
; generate redshift distribution plots based on SIMULATIONS runs

; overplot redshift distribution with strong luminosity evolution 
   
    common sirtf_simulations, sirtf

    colortable2 & plotsym, 0, 1, /fill
    path = (sirtf_datapath())[0]
    ppath = path+'plots/predictions/'

    sim0 = 'catalogs/SIMULATIONS/ZDIST/zdist_60_00_3.6_4.5_5.8_8_24_70_160.idlsave'
;   sim1 = 'catalogs/SIMULATIONS/ZDIST/zdist_60_10_3.6_4.5_5.8_8_24_70_160.idlsave'
    sim3 = 'catalogs/SIMULATIONS/ZDIST/zdist_60_30_3.6_4.5_5.8_8_24_70_160.idlsave'

    cmrestore, path+sim0, dndz0, zarray, filters, mininames, flimits, $
      names=['dndz','zarray','filters','mininames','flimits']
;   cmrestore, path+sim1, dndz1, names=['dndz']
    cmrestore, path+sim3, dndz3, names=['dndz']
    
    for k = 0L, n_elements(filters)-1L do begin

       if keyword_set(postscript) then begin
          ps_open, ppath+'zdist_'+mininames[k], /ps_fonts
          device, /inches, /times
       endif else window, 0, xs=450, ys=450

       plot, zarray, dndz0[*,k]/max(dndz0[*,k]), color=4, xsty=1, ysty=1, $
         xtit='Redshift', ytit='Fraction', charsize=1.8, xthick=2.0, ythick=2.0, $
         charthick=2.0, tit=filters[k]+textoidl('\mu')+'m Redshift Distribution', $
         thick=2.5, ps=10, yrange=[0,1.05]
;      oplot, zarray, dndz1[*,k]/max(dndz1[*,k]), color=5, line=1, ps=10, thick=2.5
       oplot, zarray, dndz3[*,k]/max(dndz3[*,k]), color=3, line=2, ps=10, thick=2.5
       legend, [textoidl('S_{lim}')+' = '+strn(2.5*1E3*flimits[k],format='(G0.0)')+$
                ' '+textoidl('\mu')+'Jy (2.5'+textoidl('\sigma')+')'], box=0, $
         /right, /top, color=16, charsize=1.2, charthick=2.0
       
       if keyword_set(postscript) then ps_close

    endfor

return
end
