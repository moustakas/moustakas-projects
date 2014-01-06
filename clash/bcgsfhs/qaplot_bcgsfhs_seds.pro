pro qaplot_bcgsfhs_seds
; jm14jan01siena - generate a QAplot of the SEDs

    showmodel = 1
    
    prefix = 'bcgsfhs'

    qapath = bcgsfhs_path()
    ancpath = bcgsfhs_path(/ancillary)
    sersicpath = bcgsfhs_path(/sersic)
    isedfit_dir = bcgsfhs_path(/isedfit)
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

; read the sample
    sample = read_bcgsfhs_sample()
    ncl = n_elements(sample)

; --------------------------------------------------
; read the ancillary data
; SDSS
    dr9 = mrdfits(ancpath+'ancillary_sdss_dr9.fits.gz',1)
    sdss_to_maggies, sdssmm, sdssii, calib=dr9
    sdssweff = k_lambda_eff()

; WISE
    wise = mrdfits(ancpath+'bcgsfhs_wise.fits.gz',1)
    wise_to_maggies, wise, wisemm, wiseii, /mpro
    wiseweff = k_lambda_eff(filterlist=wise_filterlist())
    
; make the plot    
    psfile = qapath+'qa_seds.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8, $
      height=5.0
    
    xrange = [0.3,15.0]
    yrange = [24,13]
    ticks = loglevels(xrange)
    
;   for ic = 0, 3 do begin
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
       nrad = n_elements(phot[0].photradius_kpc)
       nfilt = n_elements(phot)

       outprefix = prefix+'_'+cluster

       ised = read_isedfit(isedfit_paramfile,outprefix=outprefix,$
         thissfhgrid=1,index=0,getmodels=showmodel)
       
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xrange=xrange, yrange=yrange, xtitle='Wavelength \lambda (\mu'+'m)', $
         ytitle='Magnitude (AB)', /xlog, title=strupcase(cluster), $
         xtickv=ticks, xticks=n_elements(ticks)-1
       if showmodel then djs_oplot, ised.wave/1D4, ised.flux, color=cgcolor('grey')
       
; integrated light       
       mab = maggies2mag(phot.maggies_int,ivarmaggies=phot.ivarmaggies_int,$
         magerr=maberr,lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper,$
         nsigma=1.0)
       used = where(mab gt -90.0,nused)
       upper = where(mab lt -90.0 and mabupper gt -90,nupper)

       oploterror, phot[used].weff/1D4, mab[used], mabhierr[used], $
         psym=symcat(16), symsize=1.5, color=cgcolor('dodger blue'), $
         /hibar, errcolor=cgcolor('dodger blue')
       oploterror, phot[used].weff/1D4, mab[used], mabloerr[used], psym=3, $
         color=cgcolor('dodger blue'), /lobar, errcolor=cgcolor('dodger blue')
       if (nupper ne 0) then djs_oplot, [phot[upper].weff/1D4], [mabupper[upper]], $
         psym=symcat(11,thick=6), symsize=3.0, color=im_color('dodger blue')

; SDSS/DR9       
       this = where(cluster eq strtrim(dr9.shortname,2))
       if this[0] ne -1 then begin
          mab = maggies2mag(sdssmm[*,this],ivarmaggies=sdssii[*,this],magerr=maberr,$
            lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
;         djs_oplot, sdssweff/1D4, mab, psym=symcat(6), color=cgcolor('forest green')
          oploterror, sdssweff/1D4, mab, mabhierr, $
            psym=symcat(16), symsize=1.5, color=cgcolor('forest green'), $
            /hibar, errcolor=cgcolor('forest green')
          oploterror, sdssweff/1D4, mab, mabloerr, psym=3, $
            color=cgcolor('forest green'), /lobar, errcolor=cgcolor('forest green')
       endif
       
; WISE
       this = where(cluster eq strtrim(wise.shortname,2))
       if this[0] ne -1 then begin
          mab = maggies2mag(wisemm[*,this],ivarmaggies=wiseii[*,this],magerr=maberr,$
            lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
;         djs_oplot, wiseweff/1D4, mab, psym=symcat(6), color=cgcolor('orange')
          oploterror, wiseweff/1D4, mab, mabhierr, $
            psym=symcat(16), symsize=1.5, color=cgcolor('orange'), $
            /hibar, errcolor=cgcolor('orange')
          oploterror, wiseweff/1D4, mab, mabloerr, psym=3, $
            color=cgcolor('orange'), /lobar, errcolor=cgcolor('orange')
       endif
       
;; radial bins       
;       for ir = 0, nrad-1 do begin
;          mab = maggies2mag(phot.maggies[ir],ivarmaggies=phot.ivarmaggies[ir],$
;            magerr=maberr,lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
;          used = where(mab gt -90.0,nused)
;          upper = where(mab lt -90.0 and mabupper gt -90,nupper)
;          if (nused ne 0) then begin
;             oploterror, phot[used].weff/1D4, mab[used], mabhierr[used], $
;               psym=-symcat(16), symsize=1.5, color=cgcolor('dodger blue'), $
;               /hibar, errcolor=cgcolor('dodger blue')
;             oploterror, phot[used].weff/1D4, mab[used], mabloerr[used], psym=3, $
;               color=cgcolor('dodger blue'), /lobar, errcolor=cgcolor('dodger blue')
;          endif
;          if (nupper ne 0) then begin
;             djs_oplot, [phot[upper].weff/1D4], [mabupper[upper]], $
;               psym=symcat(11,thick=6), symsize=2.0, color=cgcolor('forest green')
;          endif
;       endfor 
    endfor 
    im_plotconfig, psfile=psfile, /psclose, /pdf
    
return
end
    
