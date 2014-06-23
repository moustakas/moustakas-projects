pro plotbcgmstar_ellipse
; jm14may18siena - generate a paper plot of just one cluster showing
; how we do the ellipse-fitting (see BCGMSTAR_ELLIPSE)

    ellpath = bcgmstar_path(/ellipse)
    sersicpath = bcgmstar_path(/sersic)
    paperpath = bcgmstar_path(/paper)

; ---------------------------------------------------------------------------    
; show all the ellipticity and position angle profiles

    

; ---------------------------------------------------------------------------    
; ellipse example plot
    cluster = 'a611'
    filt = 'f160w'
    sample = read_bcgmstar_sample()
    sample = sample[where(strtrim(sample.shortname,2) eq cluster)]

    pixscale = 0.065

    plotgalphot = 1
    
;   xtitle = 'Semi-Major Axis (kpc)'    
;   ytitle = '\mu (mag arcsec^{-2})'

;   yrange = [30,16]
    xrange = [0.7,100]
;   xrange = [0.6,300]
;   xrange = [3*pixscale*arcsec2kpc,150]

    galcolor = 'dodger blue' ; 'dark gray'
    galline = 0
    modline = 0
    modcolor = 'black'
    refmodcolor = 'firebrick'
    refmodline = 5
    galsym = 16
    galsymsize = 0.7

    psfile = paperpath+'ellipse_'+cluster+'.eps'
    im_plotconfig, 4, pos, psfile=psfile, charsize=1.4, $
      width=4.0, height=2.0*[1,1,1]

;   pos = im_getposition(nx=1,ny=3,yspace=0.0,$
;     xmargin=[1.6,0.9],width=6.0)
    arcsec2kpc = dangular(sample.z,/kpc)/206265D ; [kpc/arcsec]
           
    sersic = mrdfits(sersicpath+cluster+'-sersic.fits.gz',1,/silent)
    modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
    thisfilt = where(strtrim(modphot.band,2) eq filt)
    modphot = modphot[thisfilt]
    sersic = sersic[thisfilt]
    modgood = where(modphot.sb0fit gt 0 and modphot.majora lt sersic.amax_pixels)
    refmodgood = where(modphot[0].sb0fit gt 0)

    if plotgalphot then begin
       galphot = mrdfits(ellpath+cluster+'-ellipse-refimage.fits.gz',1,/silent)
       galgood = where(galphot.sb0fit gt 0 and galphot.majora lt sersic.amax_pixels)
    endif   

;   amax = max(pp.sma)
;   splog, '   ', band, amax
    sbmin = modphot.sblimit
    rcut = 150

; SB profile (data + model + Sersic fit)
;   if plotgalphot then begin
;      yrange = [min(modphot.sb0fit[modgood])<min(modphot[0].sb0fit[refmodgood])<$
;        min(galphot.sb0fit[galgood]),$
;        max(modphot.sb0fit[modgood])>max(modphot[0].sb0fit[refmodgood])>$
;        max(galphot.sb0fit[galgood])]
;   endif
;   yrange = -2.5*alog10(yrange)
    yrange = [26,17]
    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=5, ysty=5, $
      yrange=yrange, xrange=xrange, /xlog
;   polyfill, [rcut,10^!x.crange[1],10^!x.crange[1],rcut], $
;     [!y.crange[1],!y.crange[1],!y.crange[0],!y.crange[0]], $
;     /fill, color=cgcolor('light grey')
;   polyfill, 10^[!x.crange,!x.crange[1],!x.crange[0]], $
;     [sbmin,sbmin,!y.crange[0],!y.crange[0]], $
;     /fill, color=cgcolor('light grey')

     djs_plot, [0], [0], /nodata, /noerase, position=pos[*,0], xsty=1, ysty=1, $
       yrange=yrange, xrange=xrange, /xlog, xtickname=replicate(' ',10), $
       ytitle='\mu_{F160W} (mag arcsec^{-2})'

     xyouts, 25, 19, 'Abell 611', charsize=1.4, align=0.0
     
;   djs_oplot, [3.0,10^!x.crange[1]], sbmin*[1,1], $
;     line=0, color=cgcolor('grey')
    
;   im_legend, ['F160W/Model','F160W/Data'], /left, /bottom, box=0, $
;     margin=0, line=[modline,galline], pspacing=1.7, $
;     color=cgcolor([modcolor,galcolor]), thick=8, charsize=1.0
    
     if plotgalphot then begin
        aa = galphot.majora[galgood]*pixscale*arcsec2kpc
        up = -2.5*alog10(galphot.sb0fit[galgood]+1.0/sqrt(galphot.sb0_ivar[galgood]))
        lo = -2.5*alog10(galphot.sb0fit[galgood]-1.0/sqrt(galphot.sb0_ivar[galgood]))
;       polyfill, [aa,reverse(aa)], -2.5*alog10([lo,reverse(up)]), $
;         /fill, color=cgcolor(galcolor)
       djs_oplot, aa, -2.5*alog10(galphot.sb0fit[galgood]), $
         color=cgcolor(galcolor), line=galline, psym=-symcat(galsym,thick=2), symsize=galsymsize        
;       oploterror, aa, -2.5*alog10(galphot.sb0fit[galgood]), up+2.5*alog10(galphot.sb0fit[galgood]), $
;         psym=-symcat(galsym,thick=2), symsize=galsymsize, /up, $
;         color=cgcolor(galcolor), errcolor=cgcolor(galcolor)
;       oploterror, aa, -2.5*alog10(galphot.sb0fit[galgood]), -2.5*alog10(galphot.sb0fit[galgood])-lo, $
;         psym=-symcat(galsym,thick=2), symsize=galsymsize, /lo, $
;         color=cgcolor(galcolor), errcolor=cgcolor(galcolor)
     endif
     djs_oplot, modphot.majora[modgood]*pixscale*arcsec2kpc, $
       -2.5*alog10(modphot.sb0fit[modgood]), thick=6, $
       color=cgcolor(modcolor), line=modline ;, psym=-symcat(16,thick=2), symsize=galsymsize        

; ellipticity (model)
     djs_plot, [0], [0], /nodata, position=pos[*,1], /noerase, xsty=5, ysty=5, $
       yrange=[0.0,0.7], xrange=xrange, /xlog
;   polyfill, [rcut,10^!x.crange[1],10^!x.crange[1],rcut], $
;     [!y.crange[1],!y.crange[1],!y.crange[0],!y.crange[0]], $
;     /fill, color=cgcolor('light grey')
    
     djs_plot, [0], [0], /nodata, position=pos[*,1], /noerase, xsty=1, ysty=1, $
       yrange=[0.0,0.7], xrange=xrange, /xlog, xtickname=replicate(' ',10), $
       ytickinterval=0.2, ytitle='Ellipticity', yminor=2
     
     if plotgalphot then begin
;       aa = galphot.majora[galgood]*pixscale*arcsec2kpc
;       up = galphot.ellipticityfit[galgood]+1.0/sqrt(galphot.ellipticity_ivar[galgood])
;       lo = galphot.ellipticityfit[galgood]-1.0/sqrt(galphot.ellipticity_ivar[galgood])
;       polyfill, [aa,reverse(aa)], [lo,reverse(up)], $
;         /fill, color=cgcolor(galcolor)
        djs_oplot, galphot.majora[galgood]*pixscale*arcsec2kpc, galphot.ellipticityfit[galgood], $
          color=cgcolor(galcolor), line=galline, psym=-symcat(galsym,thick=2), symsize=galsymsize
     endif
     djs_oplot, modphot.majora[modgood]*pixscale*arcsec2kpc, modphot.ellipticityfit[modgood], $
       color=cgcolor(modcolor), line=modline, thick=6
     
;   im_legend, 'Abell 611', /left, /top, box=0, charsize=1.8
     
; position angle (model)
    xtitle = 'Semi-Major Axis (kpc)'
    delvarx, xtickname
    
    djs_plot, [0], [0], /nodata, position=pos[*,2], /noerase, xsty=5, ysty=5, $
      yrange=[0.0,180], xrange=xrange, /xlog, xtitle=xtitle
;   polyfill, [rcut,10^!x.crange[1],10^!x.crange[1],rcut], $
;     [!y.crange[1],!y.crange[1],!y.crange[0],!y.crange[0]], $
;     /fill, color=cgcolor('light grey')
    
    djs_plot, [0], [0], /nodata, position=pos[*,2], /noerase, xsty=1, ysty=1, $
      yrange=[0.0,180], xrange=xrange, /xlog, xtickname=xtickname, xtitle=xtitle, $
      ytitle='PA (degree)';, yminor=2

    if plotgalphot then begin
;      aa = galphot.majora[galgood]*pixscale*arcsec2kpc
;      up = (galphot.pafit[galgood]+1.0/sqrt(galphot.pa_ivar[galgood]))*!radeg+90
;      lo = (galphot.pafit[galgood]-1.0/sqrt(galphot.pa_ivar[galgood]))*!radeg+90
;      polyfill, [aa,reverse(aa)], [lo,reverse(up)], $
;        /fill, color=cgcolor(galcolor)
       djs_oplot, galphot.majora[galgood]*pixscale*arcsec2kpc, galphot.pafit[galgood]*!radeg+90, $
         color=cgcolor(galcolor), line=galline, psym=-symcat(galsym,thick=2), symsize=galsymsize
    endif
    djs_oplot, modphot.majora[modgood]*pixscale*arcsec2kpc, modphot.pafit[modgood]*!radeg+90, $
      color=cgcolor(modcolor), line=modline, thick=6

    im_plotconfig, psfile=psfile, /psclose, /pdf;, /pskeep

stop    
    
return
end
    
