pro hdole_seds, ps=ps
; jm01sep10uofa

    readfast, 'plot_galaxyspectrum_template.dat', sedata, skip=12, /double, nlines=nlines

    wave = sedata[0,*]*1E6      ; [micron]
    srtwave = sort(wave)
    
    sednames = 'hdole_'+['9.00','10.0','11.0','12.0','13.0']+'.dat'
    nseds = n_elements(sednames)
    
    colortable2, color
    plotfaves

    if keyword_set(ps) then begin
       ps_open, 'hdole_figure', /ps_fonts, /portrait;, /color
       device, /inches, /times, xsize=7, ysize=7
    endif else window, 0, xs=450, ys=450

    plot, [0], [0], xr=[min(wave),max(wave)], yr=[4,13], /xlog, $
      /nodata, charsize=1.5, charthick=2.0, xsty=3, ysty=3, $
      ytitle=textoidl('log \nu L_{\nu} (L_{sun})'), $
      xtitle=textoidl('\lambda (\mu')+'m)'

    for i = 0L, nseds-1L do begin

       openw, lun, sednames[i], /get_lun
       for j = 0L, nlines-1L do printf, lun, wave[srtwave[j]], sedata[i+1L,srtwave[j]], $
         format='(2x,E9.3,2x,E10.4)'

       free_lun, lun

       nu = 2.99793D14/wave[srtwave]         ; [Hz]
       lum = nu*sedata[i+1,srtwave]/3.826D26 ; L_sun

       oplot, wave[srtwave], alog10(lum), line=i, thick=2.0

    endfor

    if keyword_set(ps) then ps_close

return
end
