pro plotbcgsfhs_sbprofiles, pdf=pdf
; jm13oct19siena - plot the SB profiles

    if keyword_set(pdf) then begin
       paperpath = bcgsfhs_path()+'qaplots/'
       suffix = '.ps'
    endif else begin
       paperpath = bcgsfhs_path(/paper)
       suffix = '.eps'
    endelse

    sample = read_bcgsfhs_sample(/noa2261)
    ncl = n_elements(sample)

    pixscale = 0.065
    
; plot these bands    
    thesefilt = ['f160w','f814w','f475w']
    color = ['orange','firebrick','navy']
    color2 = ['orange red','orange','dodger blue']
    line = [0,3,5]
    nthese = n_elements(thesefilt)

; make the plot
;   xtitle = 'Semi-Major Axis (kpc)'    
;   ytitle = '\mu (mag arcsec^{-2})'
    xrange = [0.5,200]
;   xrange = [3*pixscale*arcsec2kpc,150]

    psfile = paperpath+'bcg_sbprofiles'+suffix
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.3
    pos = im_getposition(nx=5,ny=3,/landscape,yspace=0.1,xspace=0.1,$
      xmargin=[1.0,0.4])
    
;   for ic = 0, 0 do begin
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       print & splog, cluster, sample[ic].z
       datapath = bcgsfhs_path(/bcg)+cluster+'/'

       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
       
       galphot = mrdfits(datapath+cluster+'-ellipse-image.fits.gz',1,/silent)
       modphot = mrdfits(datapath+cluster+'-ellipse-model.fits.gz',1,/silent)
       nfilt = n_elements(modphot)

       reffilt = where(strtrim(galphot.band,2) eq 'f160w') ; reference filter
       
;      xrange = [3*pixscale*arcsec2kpc,max(galphot[reffilt].radius_kpc)]
;      xrange = [pixscale*arcsec2kpc,max(galphot.radius_kpc)*1.5]
;      yrange = [max(modphot.sb0fit),min(modphot.sb0fit)]
       yrange = [30,16]

       if ic gt 9 then begin
          delvarx, xtickname
       endif else begin
          xtickname = replicate(' ',10)
       endelse

       if (ic mod 5) eq 0 then begin
          delvarx, ytickname
       endif else begin
          ytickname = replicate(' ',10)
       endelse
       
       djs_plot, [0], [0], /nodata, position=pos[*,ic], noerase=ic gt 0, $
         xsty=1, ysty=1, yrange=yrange, xrange=xrange, /xlog, $
         xtitle=xtitle, ytitle=ytitle, ytickname=ytickname, xtickname=xtickname, $
         ytickinterval=4.0
       im_legend, strupcase(cluster), /right, /top, box=0, margin=0, charsize=1.0

       if ic eq 0 then begin
          im_legend, ['F160W','F814W','F475W'], /left, /bottom, box=0, $
            line=line, color=color, pspacing=1.6, margin=-0.2, charsize=0.9
       endif
       
;      for ii = 0, nfilt-1 do begin
;         modgood = where(modphot[ii].sb0fit gt 10^(-0.4*modphot[ii].sblimit))
;         djs_oplot, modphot[ii].radius_kpc[modgood], $
;           -2.5*alog10(modphot[ii].sb0fit[modgood])
;      endfor

       for ii = 0, nthese-1 do begin
          this = where(thesefilt[ii] eq strtrim(galphot.band,2))

;         notzero = lindgen(n_elements(galphot[this].sb0fit_ivar))
          galgood = where(galphot[this].sb0fit gt 0)
          djs_oplot, galphot[this].radius_kpc[galgood], $
            -2.5*alog10(galphot[this].sb0fit[galgood]), $
            line=0, color=cgcolor('medium gray') ; color=cgcolor(color2[ii])
;         oploterror, galphot[this].radius_kpc[notzero], galphot[this].sb0fit[notzero], $
;           1.0/sqrt(galphot[this].sb0fit_ivar[notzero]), color=cgcolor('dodger blue'), $
;           psym=symcat(16), symsize=0.4

          modgood = where(modphot[this].sb0fit gt 10^(-0.4*modphot[this].sblimit))
          splog, thesefilt[ii], max(modphot[this].radius_kpc[modgood])
          djs_oplot, modphot[this].radius_kpc[modgood], $
            -2.5*alog10(modphot[this].sb0fit[modgood]), $
            color=cgcolor(color[ii]), line=line[ii]
;         djs_oplot, 10^!x.crange, galphot[this].sblimit*[1,1], line=5, $
;           color=cgcolor(color[ii])
       endfor
;      cc = get_kbrd(1)
    endfor

    xyouts, min(pos[0,*])-0.05, (max(pos[3,*])-min(pos[1,*]))/2.0+min(pos[1,*]), $
      textoidl('\mu (mag arcsec^{-2})'), orientation=90, align=0.5, charsize=1.4, /norm
    xyouts, (max(pos[2,*])-min(pos[0,*]))/2.0+min(pos[0,*]), min(pos[1,*])-0.09, $
      textoidl('Semi-Major Axis (kpc)'), align=0.5, charsize=1.4, /norm

    im_plotconfig, psfile=psfile, /psclose, pdf=pdf

stop    
    
return
end
    
