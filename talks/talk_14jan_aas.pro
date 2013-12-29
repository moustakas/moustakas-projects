pro talk_14jan_aas, pdf=pdf
; jm13dec28siena - generate some plots for my poster for the 223 AAS
; meeting in DC in Jan/2014

    ellpath = bcgsfhs_path()+'ellipse/'
    sersicpath = bcgsfhs_path()+'sersic/'
    talkpath = getenv('IM_RESEARCH_DIR')+'/talks/2014/14jan_aas/'
    
; --------------------------------------------------
; plot: radial SB profile for A611
    cluster = 'a611'
    sample = read_bcgsfhs_sample()
    sample = sample[where(strupcase(strtrim(sample.shortname,2)) eq cluster)]
    arcsec2kpc = dangular(sample.z,/kpc)/206265D ; [kpc/arcsec]

    pixscale = 0.065
    
; plot these bands    
    thesefilt = ['f475w','f814w','f160w']
    color1 = ['medium grey','powder blue','tomato']
    color2 = ['black','navy','firebrick']
    line = [2,3,0]
    nthese = n_elements(thesefilt)

    xrange = [0.5,100]
    yrange = [27,17]

    psfile = talkpath+'a611_sbprofile.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, height=4.7, width=5.0
       
    sersic = mrdfits(sersicpath+cluster+'-sersic.fits.gz',1,/silent)
    galphot = mrdfits(ellpath+cluster+'-ellipse-image.fits.gz',1,/silent)
    modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
    nfilt = n_elements(modphot)
       
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      yrange=yrange, xrange=xrange, /xlog, $
      xtitle='Equivalent Radius (kpc)', ytitle='\mu (mag arcsec^{-2})', $
      xthick=8, ythick=8, charthick=4
    im_legend, strupcase(cluster), /right, /top, box=0, margin=0, $
      charsize=2.0, charthick=4.0

;   im_legend, ['F475W','F814W','F160W'], /left, /bottom, box=0, $
;     line=line, color=color2, pspacing=1.8, margin=-0.2, charsize=2.0

    for ii = 0, nthese-1 do begin
       this = where(thesefilt[ii] eq strtrim(galphot.band,2))

       modgood = where(modphot[this].majora*pixscale*arcsec2kpc le sersic[this].amax_kpc and $
         modphot[this].sb0fit gt 0 and modphot[this].sb0fit_ivar gt 0)

       rr = modphot[this].radius_kpc[modgood]
       sb = -2.5*alog10(modphot[this].sb0fit[modgood])
       sberr = 2.5/(modphot[this].sb0fit[modgood]*sqrt(modphot[this].sb0fit_ivar[modgood])*alog(10))

; just for the plot make the uncertainty have a minimum floor so that
; the SB profile shows up!          
       sberr = sberr>0.1
       polyfill, [rr,reverse(rr)], [sb-sberr,reverse(sb+sberr)], /fill, $
         color=cgcolor(color1[ii]), noclip=0
       
; overplot the single Sersic function
       rrser = [0,range(1E-3,200.0,500,/log)]
       djs_oplot, rrser, bcgsfhs_sersic_func(rrser,params=sersic[this]), $
         line=line[ii], thick=6; color=cgcolor(color2[ii]), 
    endfor

    xyouts, 7.0, 19.5, 'F160W', align=0.0, /data, charsize=1.6, charthick=4.0
    xyouts, 3.0, 21.5, 'F814W', align=0.0, /data, charsize=1.6, charthick=4.0
    xyouts, 1.5, 23.2, 'F475W', align=0.0, /data, charsize=1.6, charthick=4.0
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
return
end
    
