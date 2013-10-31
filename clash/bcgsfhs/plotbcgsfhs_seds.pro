pro plotbcgsfhs_seds, pdf=pdf
; jm13oct25siena - plot the radial surface brightness profiles

    if keyword_set(pdf) then begin
       paperpath = bcgsfhs_path()
       suffix = '.ps'
    endif else begin
       paperpath = bcgsfhs_path(/paper)
       suffix = '.eps'
    endelse

    sample = read_bcgsfhs_sample()
    ncl = n_elements(sample)

    ellpath = bcgsfhs_path()+'ellipse/'
    sersicpath = bcgsfhs_path()+'sersic/'
    
; plot these bands    
    thesefilt = ['f160w','f814w','f475w']
    color1 = ['tomato','powder blue','medium grey']
    color2 = ['firebrick','navy','black']

;   color = ['orange','firebrick','navy']
;   color2 = ['orange red','orange','dodger blue']
    line = [0,3,5]
    nthese = n_elements(thesefilt)

; make the plot
;   xtitle = 'Semi-Major Axis (kpc)'    
;   ytitle = '\mu (mag arcsec^{-2})'
    xrange = [0.5,120]
;   xrange = [3*pixscale*arcsec2kpc,150]
    yrange = [28,16]

    psfile = paperpath+'bcg_seds'+suffix
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.3
    pos = im_getposition(nx=5,ny=3,/landscape,yspace=0.1,xspace=0.1,$
      xmargin=[1.0,0.4])

;   for ic = 0, 0 do begin
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       print & splog, cluster, sample[ic].z
;      datapath = bcgsfhs_path(/bcg)+cluster+'/'

       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
       
       sersic = mrdfits(sersicpath+cluster+'-sersic.fits.gz',1,/silent)
       galphot = mrdfits(ellpath+cluster+'-ellipse-image.fits.gz',1,/silent)
       modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
       nfilt = n_elements(modphot)

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
            line=line, color=color2, pspacing=1.8, margin=-0.2, charsize=0.8
       endif
       
       for ii = 0, nthese-1 do begin
          this = where(thesefilt[ii] eq strtrim(galphot.band,2))

          modgood = where(modphot[this].sb0fit gt 10^(-0.4*modphot[this].sblimit) and $
            modphot[this].sb0fit_ivar gt 0 and modphot[this].radius_kpc lt 100.0) ; r<100 kpc
          splog, thesefilt[ii], max(modphot[this].radius_kpc[modgood])

          rr = modphot[this].radius_kpc[modgood]
          sb = -2.5*alog10(modphot[this].sb0fit[modgood])
          sberr = 2.5/(modphot[this].sb0fit[modgood]*sqrt(modphot[this].sb0fit_ivar[modgood])*alog(10))

; just for the plot make the uncertainty have a minimum floor so that
; the SB profile shows up!          
          sberr = sberr>0.05

          polyfill, [rr,reverse(rr)], [sb-sberr,reverse(sb+sberr)], /fill, $
            color=cgcolor(color1[ii])
       endfor
;      cc = get_kbrd(1)
    endfor

    xyouts, min(pos[0,*])-0.05, (max(pos[3,*])-min(pos[1,*]))/2.0+min(pos[1,*]), $
      textoidl('\mu (mag arcsec^{-2})'), orientation=90, align=0.5, charsize=1.4, /norm
    xyouts, (max(pos[2,*])-min(pos[0,*]))/2.0+min(pos[0,*]), min(pos[1,*])-0.09, $
      textoidl('Equivalent Radius (kpc)'), $
;     textoidl('Equivalent Radius r=a\sqrt{1-\epsilon} (kpc)'), $
      align=0.5, charsize=1.4, /norm

    im_plotconfig, psfile=psfile, /psclose, pdf=pdf

stop    
    
return
end
    
