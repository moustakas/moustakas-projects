pro highzea_spec1d_postscript, nsmooth=nsmooth, postscript=postscript
; jm06jul07uofa - make one big postscript file

    datapath = highzea_path(/spec1d)
    analysis_path = highzea_path(/analysis)

    specfile = file_basename(file_search(datapath+'*.ms.fits'))
    ngalaxy = n_elements(specfile)

    if (n_elements(nsmooth) eq 0L) then nsmooth = 3L
    
    if keyword_set(postscript) then begin
       dfpsplot, analysis_path+'highzea_spec1d.ps', /color, xsize=8.5, ysize=5.5;, /square
       postthick = 4.0
       postthick2 = 2.0
    endif else begin
       im_window, 0, xratio=0.5, /square
       postthick = 2.0
       postthick2 = 1.0
    endelse
    
    pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, xmargin=[1.1,0.4], ymargin=[0.4,1.1], $
      xpage=8.5, ypage=5.5, width=7.0, height=4.0, position=pos, /normal
    
    scale = 1E17
    xtitle = 'Wavelength [\AA]'
    ytitle='Flux [10^{-17} '+flam_units()+']'

    for i = 0L, ngalaxy-1L do begin

       scube = rd1dspec(specfile[i],datapath=datapath)
       wave = scube.wave
       flux = smooth(scube.spec,nsmooth)
       header = scube.header
       galaxy = sxpar(header,'GALAXY')

; S/N statistics

       goodsig = where(scube.sigspec gt 0.0,ngoodsig,ncomp=nbadsig)
       if (ngoodsig ne 0L) then begin
          snr = scube.spec[goodsig]/scube.sigspec[goodsig]
          snrmean = mean(snr) & snrmedian = median(snr) & snrsig = stddev(snr)
          snrstr = 'S/N = '+strtrim(string(snrmedian,format='(F12.1)'),2)
       endif else snrstr = ''

; make the plot       

       yrange = minmax(scale*flux)*[0.9,1.1]
       
       djs_plot, wave, scale*flux, xsty=3, ysty=3, ps=10, xrange=xrange, yrange=yrange, $
         xthick=postthick, ythick=postthick, charsize=1.5, charthick=postthick, $
         thick=postthick2, xtitle=xtitle, ytitle=ytitle, position=pos[*,0];, color='grey'
       legend, [galaxy,snrstr], /right, /top, box=0, charsize=1.4, charthick=postthick, $
         clear=keyword_set(postscript)

       if (not keyword_set(postscript)) then cc = get_kbrd(1)
       icleanup, scube

    endfor
    
    if keyword_set(postscript) then begin
       dfpsclose
       spawn, ['gzip -f '+analysis_path+'highzea_spec1d.ps'], /sh
    endif

return
end
    
