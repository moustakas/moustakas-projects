pro talk_14jan_aas, pdf=pdf
; jm13dec28siena - generate some plots for my poster for the 223 AAS
; meeting in DC in Jan/2014

    ellpath = bcgsfhs_path()+'ellipse/'
    sersicpath = bcgsfhs_path()+'sersic/'
    talkpath = getenv('IM_RESEARCH_DIR')+'/talks/2014/14jan_aas/'

    prefix = 'bcgsfhs'
    isedfit_dir = bcgsfhs_path(/isedfit)
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

    filt = bcgsfhs_filterlist(pivotwave=weff,width=hwhm)
    weff = weff/1D4
    hwhm = hwhm/1D4*0
    
; --------------------------------------------------
; BCG SED-fitting example
    cluster = 'rxj2248'
    outprefix = prefix+'_'+cluster

    phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
    rad = phot[0].photradius_kpc
    ised = read_isedfit(isedfit_paramfile,outprefix=outprefix,$
      thissfhgrid=3,/silent,index=8,/getmodels,isedfit_post=post) ; R=5 kpc
;     thissfhgrid=4,/silent,index=0,/getmodels) ; integrated

    xtitle = textoidl('Wavelength (\mu'+'m)')
    ytitle = textoidl('AB Magnitude')

;   yrange = [23,13] ; integrated
    yrange = [24,17.5]
    xrange = [0.3,1.7]
    col = 'forest green'
    
    psfile = talkpath+'rxj2248_sed.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=4.0

    ticks = loglevels(xrange,/fine)
    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
      xsty=1, ysty=1, xlog=0, position=pos, xtitle=xtitle, ytitle=ytitle, $
      ytickname=ytickname, xtickname=xtickname, xtickinterval=0.5, $
      xthick=8, ythick=8, charthick=5, ytickinterval=2           ;, xtickv=ticks, xticks=n_elements(ticks)-1
    djs_oplot, ised.wave/1D4, ised.flux, color=im_color('grey60')

; model photometry
    notzero = where(ised.bestmaggies gt 0.0 and ised.maggies gt 0)
    bestmab = -2.5*alog10(ised.bestmaggies[notzero])
    djs_oplot, weff[notzero], bestmab, psym=symcat(6,thick=6), $
      symsize=3.2;, color=cgcolor('grey')
; overplot the data; distinguish between three different cases, based
; on the input photometry
    mab = maggies2mag(ised.maggies,nsigma=0.0,$
      ivar=ised.ivarmaggies,magerr=maberr,$
      lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
    used = where(mab gt -90.0,nused)
    upper = where(mab lt -90.0 and mabupper gt -90,nupper)
    if (nused ne 0L) then begin
       oploterror, weff[used], mab[used], hwhm[used], $
         mabhierr[used], psym=symcat(16), $
         symsize=2.0, color=cgcolor(col), /hibar, $
         errcolor=cgcolor(col), errthick=5
       oploterror, weff[used], mab[used], hwhm[used], $
         mabloerr[used], psym=3, color=cgcolor(col), /lobar, $
         errcolor=cgcolor(col), errthick=5
    endif
    if (nupper ne 0) then begin
       djs_oplot, weff[upper], mabupper[upper], $
         psym=symcat(11,thick=6), symsize=3.0, color=cgcolor('forest green')
    endif
;   djs_oplot, weff, -2.5*alog10(ised.bestmaggies), psym=7
;   djs_oplot, weff, -2.5*alog10(ised.maggies), psym=6, color='red'
    
    im_legend, 'RXJ2248 at R=10 kpc', /left, /top, box=0, $
      charsize=1.7, charthick=5, margin=0

; inset with P(age)
    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=[0,1.05], $
      xrange=[0.0,9.0], position=[0.53,0.34,0.88,0.64]
    im_plothist, post.age, bin=0.4, /peak, /overplot, /fill, $
      fcolor=im_color('grey')
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.05], $
      xrange=[0.0,9.0], position=[0.53,0.34,0.88,0.64], ytitle='Probability', $
      xtitle='Age (Gyr)', xtickinterval=2.0, ytickname=replicate(' ',10), $
      charsize=1.5, xthick=8, ythick=8, charthick=5
    
    im_plotconfig, psfile=psfile, /psclose, /pdf
    
; --------------------------------------------------
; plot age and stellar mass vs radius for two example clusters 
    agerange = [1.0,9.0]
    massrange = [10.0,12.2]

    xrange = [1,50]
    ticks = loglevels(xrange)
    
    psfile = talkpath+'agemass.ps'
    im_plotconfig, 6, pos, psfile=psfile, charsize=2.0, width=6.3, $
      xmargin=[1.1,1.1], height=[4.7,4.7], yspace=0.1
       
; MACS1149
    cluster = 'macs0744'
    outprefix = prefix+'_'+cluster

    phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
    rad = phot[0].photradius_kpc
    ised = read_isedfit(isedfit_paramfile,outprefix=outprefix,$
      thissfhgrid=4,/silent,index=lindgen(n_elements(rad))+1) ; no integrated
    good = where(ised.chi2 lt 1E6 and total(ised.ivarmaggies gt 0,1) ge 6) ; minimum of 5 bands 

    djs_plot, [0], [0], /nodata, /xlog, xsty=1, ysty=9, position=pos[*,0], $
      xrange=xrange, yrange=agerange, ytitle='', $
      xthick=8, ythick=8, charthick=5, xtitle='', xtickname=replicate(' ',10), $
      xtickv=ticks, xticks=n_elements(ticks)-1
    polyfill, [rad[good],reverse(rad[good])], [ised[good].age_50+ised[good].age_err/2,$
      reverse(ised[good].age_50-ised[good].age_err/2)], /fill, color=cgcolor('purple')

    xyouts, 18.0, 6.8, 'Mass', align=0.5, color=cgcolor('dodger blue'), $
      /data, charthick=5.0
    xyouts, 18.0, 2.5, 'Age', align=0.5, color=cgcolor('purple'), /data, charthick=5.0
    
    axis, /yaxis, yrange=massrange, ysty=1, /save, ytitle='', $
      ytickinterval=1.0, ythick=8, charthick=5
    cumumass = alog10(total(10.0^ised[good].mstar_50,/cumu))
    polyfill, [rad[good],reverse(rad[good])], [cumumass+0.06,reverse(cumumass-0.06)], $
      /fill, color=cgcolor('dodger blue')
;   djs_oplot, rad[good], alog10(total(10.0^ised[good].mstar_50,/cumu)), psym=8
    im_legend, 'MACS0744', /left, /top, box=0, charthick=5

; RXJ2248
    cluster = 'rxj2248'
    outprefix = prefix+'_'+cluster

    phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
    rad = phot[0].photradius_kpc
    ised = read_isedfit(isedfit_paramfile,outprefix=outprefix,$
      thissfhgrid=4,/silent,index=lindgen(n_elements(rad))+1) ; no integrated
    good = where(ised.chi2 lt 1E6 and total(ised.ivarmaggies gt 0,1) ge 6) ; minimum of 5 bands 

    djs_plot, [0], [0], /nodata, /xlog, xsty=1, ysty=9, position=pos[*,1], /noerase, $
      xrange=xrange, yrange=agerange, ytitle='', $
      xthick=8, ythick=8, charthick=5, xtitle='Radius (kpc)', $
      xtickv=ticks, xticks=n_elements(ticks)-1
    polyfill, [rad[good],reverse(rad[good])], [ised[good].age_50+ised[good].age_err/2,$
      reverse(ised[good].age_50-ised[good].age_err/2)], /fill, color=cgcolor('purple')
;   djs_oplot, rad[good], ised[good].age, psym=8

    xyouts, 18.0, 6.8, 'Mass', align=0.5, color=cgcolor('dodger blue'), $
      /data, charthick=5.0
    xyouts, 18.0, 4.1, 'Age', align=0.5, color=cgcolor('purple'), /data, charthick=5.0
    
    axis, /yaxis, yrange=massrange, ysty=1, /save, ytitle='', $
      ytickinterval=1.0, ythick=8, charthick=5
    cumumass = alog10(total(10.0^ised[good].mstar_50,/cumu))
    polyfill, [rad[good],reverse(rad[good])], [cumumass+0.06,reverse(cumumass-0.06)], $
      /fill, color=cgcolor('dodger blue')
;   djs_oplot, rad[good], alog10(total(10.0^ised[good].mstar_50,/cumu)), psym=8
    im_legend, 'RXJ2248', /left, /top, box=0, charthick=5

    xyouts, pos[0,0]-0.08, pos[1,0], 'Onset of Star Formation (Gyr Ago)', $
      align=0.5, orientation=90.0, charthick=5, /norm, charsize=2.2
    xyouts, pos[2,0]+0.1, pos[1,0], textoidl('Log(Enclosed Stellar Mass, M_{\odot})'), $
      align=0.5, orientation=90.0, charthick=5, /norm, charsize=2.2
    
    im_plotconfig, psfile=psfile, /psclose, /pdf
    
stop    
    
; --------------------------------------------------
; plot: radial SB profile for A611
    cluster = 'a611'
    nicename = 'Abell 611'
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
      xtitle='Radius (kpc)', ytitle='\mu (mag arcsec^{-2})', $
      xthick=8, ythick=8, charthick=5
    im_legend, nicename, /right, /top, box=0, margin=0, $
      charsize=2.0, charthick=5

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
    
