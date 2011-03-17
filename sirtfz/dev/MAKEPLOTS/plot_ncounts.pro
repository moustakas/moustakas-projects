pro plot_ncounts, postscript=postscript
; jm01may13uofa
; generate number counts plots based on SIMULATIONS runs
    
    common sirtf_simulations, sirtf

    colortable2 & plotsym, 0, 1, /fill
    ppath = (sirtf_datapath())[0]+'plots/predictions/'

    fname = (sirtf_datapath())[0]+'catalogs/SIMULATIONS/NCOUNTS/counts_60_00_3.6_4.5_5.8_8_24_70_160.idlsave'
    cmrestore, fname, dnds, flim, filters, mininames

    for k = 0L, n_elements(filters)-1L do begin

       if keyword_set(postscript) then begin
          ps_open, ppath+'counts_'+mininames[k], /ps_fonts
          device, /inches, /times
       endif else window, 0, xs=450, ys=450

       nozero = where(cumulate(dnds[*,k]) gt 0.0,nr)
       if nr eq 0L then nozero = lindgen(n_elements(flim))

       plot, alog10(flim[nozero]), alog10(cumulate(dnds[nozero,k])), color=4, xsty=3, ysty=3, $
         xtit='log S (mJy)', ytit='log N(>S)', charsize=2.0, xthick=2.0, ythick=2.0, $
         charthick=2.0, tit=filters[k]+' Integral Counts', thick=2.5

       if keyword_set(postscript) then ps_close

    endfor

return
end
