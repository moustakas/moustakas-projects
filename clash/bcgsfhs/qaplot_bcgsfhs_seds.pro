pro qaplot_bcgsfhs_seds
; jm13oct29siena - generate a QAplots of the SEDs

    sample = read_bcgsfhs_sample()
    ncl = n_elements(sample)

    sersicpath = bcgsfhs_path()+'sersic/'
    psfile = sersicpath+'qa_seds.ps'

    im_plotconfig, 0, pos, psfile=psfile, charsize=1.8, $
      height=5.0
    
; make the plot
    xrange = [0.3,2.0]
    yrange = [26,10]
    
    for ic = 0, 3 do begin
;   for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       phot = mrdfits(sersicpath+cluster+'-sersic.fits.gz',1,/silent)
       nrad = n_elements(phot[0].photradius_kpc)
       nfilt = n_elements(phot)

       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xrange=xrange, yrange=yrange, xtitle='Wavelength \lambda (\mu'+'m)', $
         ytitle='Magnitude (AB)', /xlog, title=strupcase(cluster)

; integrated light       
       mab = maggies2mag(phot.maggies_int,ivarmaggies=phot.ivarmaggies_int,$
         magerr=maberr,lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
       used = where(mab gt -90.0,nused)
       upper = where(mab lt -90.0 and mabupper gt -90,nupper)

       oploterror, phot[used].weff/1D4, mab[used], mabhierr[used], $
         psym=-symcat(16), symsize=1.5, color=cgcolor('firebrick'), $
         /hibar, errcolor=cgcolor('firebrick')
       oploterror, phot[used].weff/1D4, mab[used], mabloerr[used], psym=3, $
         color=cgcolor('firebrick'), /lobar, errcolor=cgcolor('firebrick')

; radial bins       
       for ir = 0, nrad-1 do begin
          mab = maggies2mag(phot.maggies[ir],ivarmaggies=phot.ivarmaggies[ir],$
            magerr=maberr,lomagerr=mabloerr,himagerr=mabhierr,magnsigma=mabupper)
          used = where(mab gt -90.0,nused)
          upper = where(mab lt -90.0 and mabupper gt -90,nupper)

          if (nused ne 0) then begin
             oploterror, phot[used].weff/1D4, mab[used], mabhierr[used], $
               psym=-symcat(16), symsize=1.5, color=cgcolor('dodger blue'), $
               /hibar, errcolor=cgcolor('dodger blue')
             oploterror, phot[used].weff/1D4, mab[used], mabloerr[used], psym=3, $
               color=cgcolor('dodger blue'), /lobar, errcolor=cgcolor('dodger blue')
          endif

          if (nupper ne 0) then begin
             djs_oplot, [phot[upper].weff/1D4], [mabupper[upper]], $
               psym=symcat(11,thick=6), symsize=2.0, color=cgcolor('forest green')
          endif
       endfor 
    endfor 

    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
return
end
    
