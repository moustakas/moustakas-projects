pro write_sdss_alpha_sample, sdssdust, sdssancillary, write=write, $
  dr6=dr6, nopsfile=nopsfile
; jm07oct15nyu - written based on WRITE_DEEP2_ALPHA_SAMPLE
; jm08mar25nyu - updated for apr08 observing run

    if keyword_set(dr6) then dr = '_dr6' else dr = ''
    
    run = '08apr'+dr
    nfinal = 100L ; 500L

    red, h100=0.7, omega0=0.3, omega_lambda=0.7
    h100 = redh100()
    lsun = 3.826D33

    alphapath = deep2_path(/projects)+'alpha/'

    stime0 = systime(1)
    if keyword_set(write) then begin
       splogfile = alphapath+'write_sdss_alpha_sample_'+run+'.log'
       splog, filename=splogfile
       splog, 'Log file '+splogfile+' opened '+systime()
       splog, 'IDL version: ' + string(!version,format='(99(A," "))')
       spawn, 'uname -a', uname
       splog, 'UNAME: '+uname[0]
    endif
       
    if keyword_set(write) then begin
       postscript = 1L
       postthick1 = 5.0
       postthick2 = 3.0
    endif else begin
       postthick1 = 1.0
       postthick2 = 1.0
;      im_window, 0, xratio=0.45, /square
    endelse

; ---------------------------------------------------------------------------
; read the emission-line data
; ---------------------------------------------------------------------------

    if (n_elements(sdssdust) eq 0L) then begin
       splog, 'Reading SDSS data...'
       if keyword_set(dr6) then begin
          sdsscat = read_sdss_main(/mpacat,/dr6)
          sdssancillary1 = read_sdss_main(/kcorr,/dr6)
          sdssdust1 = read_sdss_main(/ispec,/dr6)
          these = where(strtrim(sdsscat.release,2) eq 'dr5' or $
            strtrim(sdsscat.release,2) eq 'dr6')
;         these = where(strmatch(sdsscat.release,'*dr6*',/fold),nthese)
          sdssdust = (temporary(sdssdust1))[these]
          sdssancillary = (temporary(sdssancillary1))[these]
       endif else begin
          sdssancillary = read_sdss_main(/kcorr)
          sdssdust = read_sdss_main(/ispec)
       endelse
    endif
    ngalaxy = n_elements(sdssdust)

; ---------------------------------------------------------------------------
; select the parent sample
; ---------------------------------------------------------------------------

    snrcut = 10.0
    minmagerr = 0.2
    sample_zmin = 0.01
    sample_zmax = 0.20

    alpha = where((sdssancillary.z gt sample_zmin) and (sdssancillary.z lt sample_zmax) and $
      (sdssancillary.ra/15.0 gt 8.0) and (sdssancillary.ra/15.0 lt 20.0) and (sdssancillary.dec lt 15.0) and $
;     (sdssancillary.kcorr_mobs_ab[1] gt -900.0) and (sdssancillary.kcorr_mobs_ab[2] gt -900.0) and $
;     (sdssancillary.kcorr_mobs_ab[3] gt -900.0) and (sdssancillary.kcorr_chi2 lt 5.0) and $
;     (abs(sdssancillary.kcorr_mobs_ab[2]-sdssancillary.kcorr_mobs_ab_synth[2]) lt 0.2) and $
;     (sdssancillary.kcorr_mobs_ab_ivar[1] gt 1.0/minmagerr^2) and $
;     (sdssancillary.kcorr_mobs_ab_ivar[2] gt 1.0/minmagerr^2) and $
;     (sdssancillary.kcorr_mobs_ab_ivar[3] gt 1.0/minmagerr^2) and $
      (sdssdust.oiii_4959_ew[0] gt 0.0) and (sdssdust.oiii_4959_ew[1] gt 0.0) and $ ; EW>0; EW_err>0
      (sdssdust.oiii_4959_ew[0]/sdssdust.oiii_4959_ew[1] gt snrcut) and $           ; hard S/N cut
      (sdssdust.oiii_5007_ew[0] gt 0.0) and (sdssdust.oiii_5007_ew[1] gt 0.0) and $ ; EW>0; EW_err>0
      (sdssdust.oiii_5007_ew[0]/sdssdust.oiii_5007_ew[1] gt snrcut),nalpha)         ; hard S/N cut

; throw away objects that are within LAMBDA_TOL of a night sky line

    lambda_tol = 1.0 ; wavelength tolerance for proximity to a sky line
    
    lambda_oiii_4959 = 4958.911
    lambda_oiii_5007 = 5006.843

    readcol, getenv('ISPEC_DIR')+'/etc/skylist_o2_oh.dat', lambda_sky, format='D', /silent
    lambda_sky = lambda_sky[where((lambda_sky gt 4950.0) and (lambda_sky lt 6010.0))]
    nsky = n_elements(lambda_sky)

    allgood = bytarr(nsky,nalpha)
    for ii = 0L, nalpha-1L do begin
       for isky = 0L, nsky-1L do begin         
          allgood[isky,ii] = (abs(lambda_oiii_4959*(1.0+sdssancillary[alpha[ii]].z)-lambda_sky[isky]) gt lambda_tol) and $
            (abs(lambda_oiii_5007*(1.0+sdssancillary[alpha[ii]].z)-lambda_sky[isky]) gt lambda_tol)
       endfor
    endfor

    good = where((nsky-total(allgood,1) eq 0.0),nalpha)
    alpha = alpha[good]

    splog, 'alpha sample: '+string(nalpha,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+' ('+$
      strtrim(string(100.0*nalpha/(ngalaxy),format='(F12.1)'),2)+'%).'

    zstats = im_stats(sdssancillary[alpha].z)
    splog, '   Redshift: ['+strtrim(string(zstats.min,format='(F12.2)'),2)+'-'+strtrim(string(zstats.max,format='(F12.2)'),2)+'] '+$
      strtrim(string(zstats.median,format='(F12.2)'),2)+' ('+strtrim(string(zstats.mean,format='(F12.2)'),2)+'+/-'+$
      strtrim(string(zstats.sigma,format='(F12.2)'),2)+')'

; ---------------------------------------------------------------------------    
; write out
; ---------------------------------------------------------------------------    
    
    if keyword_set(write) then begin

       splog, 'Computing luminosity distances'
       dlum = dluminosity(sdssancillary[alpha].z,/cm)
       
       splog, 'Appending additional tags'
       moretags = replicate({oiii_4959_flux: [-999.0,-999.0], oiii_4959_lum: [-999.0,-999.0], $
         oiii_5007_flux: [-999.0,-999.0], oiii_5007_lum: [-999.0,-999.0], $
         magu: -999.0, magg: -999.0, magr: -999.0, magi: -999.0, magz: -999.0, galaxy: ''},nalpha)
       moretags.oiii_4959_flux[0] = sdssancillary[alpha].cflux_4959*sdssdust[alpha].oiii_4959_ew[0]
       moretags.oiii_4959_flux[1] = sdssancillary[alpha].cflux_4959*sdssdust[alpha].oiii_4959_ew[1]
       moretags.oiii_4959_lum[0] = moretags.oiii_4959_flux[0]*4.0*!dpi*dlum*dlum/lsun
       moretags.oiii_4959_lum[1] = moretags.oiii_4959_flux[1]*4.0*!dpi*dlum*dlum/lsun

       moretags.oiii_5007_flux[0] = sdssancillary[alpha].cflux_5007*sdssdust[alpha].oiii_5007_ew[0]
       moretags.oiii_5007_flux[1] = sdssancillary[alpha].cflux_5007*sdssdust[alpha].oiii_5007_ew[1]
       moretags.oiii_5007_lum[0] = moretags.oiii_5007_flux[0]*4.0*!dpi*dlum*dlum/lsun
       moretags.oiii_5007_lum[1] = moretags.oiii_5007_flux[1]*4.0*!dpi*dlum*dlum/lsun

       moretags.magu = -2.5*alog10(sdssancillary[alpha].abmaggies[0])
       moretags.magg = -2.5*alog10(sdssancillary[alpha].abmaggies[1])
       moretags.magr = -2.5*alog10(sdssancillary[alpha].abmaggies[2])
       moretags.magi = -2.5*alog10(sdssancillary[alpha].abmaggies[3])
       moretags.magz = -2.5*alog10(sdssancillary[alpha].abmaggies[4])
       
       sdssalpha1 = struct_addtags(struct_addtags(sdssancillary[alpha],$
         struct_trimtags(sdssdust[alpha],except=['RA','DEC','Z'])),$
         moretags)
       sdssalpha1.galaxy = 'sdss.'+string(sdssalpha1.mjd,format='(I0)')+'.'+string(sdssalpha1.plateid,format='(I4.4)')+'.'+$
         string(sdssalpha1.fiberid,format='(I3.3)')

; select the NFINAL galaxies with the highest [O III] luminosities

       final = (reverse(sort(sdssalpha1.oiii_5007_lum[0])))[0L:nfinal-1L]
       sdssalpha = sdssalpha1[final]
       sdssalpha_dennis = struct_trimtags(sdssalpha,select=['GALAXY','Z',$
         'RA','DEC','MAGU','MAGG','MAGR','MAGI','MAGZ','OIII_5007_LUM'])
       
; k-correct preliminaries

       vname = 'default.nolines'
       light = 2.99792458D18    ; speed of light [A/s]

       k_load_vmatrix, vmatrix, lambda, vfile=vfile, vpath=vpath, lfile=lfile, vname=vname
       plotsym, 0, 1.5, /fill
       
       filters = 'sdss_'+['u0','g0','r0','i0','z0']+'.par'
       filtinfo = im_filterspecs(filterlist=filters,/verbose)
       nfilter = n_elements(filtinfo)

; plot the sdss spectra

       if keyword_set(nopsfile) then begin
          sdssspectra = struct_addtags(sdssalpha_dennis,replicate({wave: fltarr(3891), $
             flux: fltarr(3891), ivar: fltarr(3891)},nfinal))
       endif else begin

          readspec, sdssalpha.plateid, sdssalpha.fiberid, $
            mjd=sdssalpha.mjd, flux=allflux, flerr=allferr, $
            invvar=allivar, wave=allwave, /align
          npix = n_elements(allwave)
          sdssspectra = struct_addtags(sdssalpha_dennis,replicate({wave: float(allwave), $
            flux: fltarr(npix), ivar: fltarr(npix)},nfinal))
          sdssspectra.flux = allflux*1D-17
          sdssspectra.ivar = allivar

          psname = alphapath+'sdss_alpha_'+run+'.ps'
          dfpsplot, psname, /color, /landscape

          pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, ymargin=[0.5,1.1], $
            xmargin=[1.1,0.2], position=pos, /normal

          for igal = 0L, n_elements(sdssalpha)-1L do begin

             galaxy = strtrim(sdssalpha[igal].galaxy,2)
             z = sdssalpha[igal].z
             
             wave = sdssspectra[igal].wave
             flux = sdssspectra[igal].flux
             ivar = sdssspectra[igal].ivar

             good = where(ivar gt 0,ngood)
             wave = wave[good]
             flux = 1D17*flux[good]
             ivar = ivar[good]
             ferr = 1.0/sqrt(ivar)
             npix = n_elements(wave)

             restwave = wave/(1.0+z)
             restflux = flux*(1.0+z)
             restivar = ivar/(1.0+z)^2.0
             restferr = ferr*(1.0+z)

; make the plot

             title = galaxy+', z = '+strtrim(string(z,format='(F12.4)'),2)+$
               ', L[O III] = '+strtrim(string(sdssalpha[igal].oiii_5007_lum[0],format='(E10.2)'),2)+' L_{\odot}'

             xrange = [4800.0,5050.0]
             get_element, restwave, xrange, xx
             yrange = minmax(restflux[xx[0]:xx[1]])
             yrange[0] = yrange[0]*0.9
             
             djs_plot, [0], [0], /nodata, xsty=1, ysty=1, position=pos, $
               xthick=postthick1, ythick=postthick1, charthick=postthick2, $
               xtitle='Rest Wavelength (\AA)', charsize=1.4, xrange=xrange, $
               yrange=yrange, ytitle='Flux (10^{-17} '+flam_units()+')', title=title
             djs_oplot, restwave, restflux, thick=postthick1, color='grey', ps=10

; make an inset with the SED fit

             restwave = lambda
             restflux = vmatrix#sdssalpha[igal].coeffs ; f_lambda
             restflux_fnu = restflux*restwave^2.0/light*10.0^(0.4*48.6) ; f_nu
             restflux_ab = -2.5*alog10(restflux_fnu>1D-32)

             wave = restwave*(1.0+z)
             flux = restflux/(1.0+z)
             flux_fnu = flux*wave^2.0/light*10.0^(0.4*48.6) ; f_nu
             flux_ab = -2.5*alog10(flux_fnu>1D-32)

             abgood = where((sdssalpha[igal].abmaggies_ivar gt 0.0),nabgood)
             mab = -2.5*alog10(sdssalpha[igal].abmaggies[abgood])
             mab_err = 1.0/sqrt(sdssalpha[igal].abmaggies_ivar[abgood])/$
               sdssalpha[igal].abmaggies[abgood]/alog(10.0)

             xrange = [min(filtinfo[abgood].weff-1.5*filtinfo[abgood].fwhm),$
               max(filtinfo[abgood].weff+1.5*filtinfo[abgood].fwhm)] ;/(1.0+zobj)
             xrange[0] = xrange[0]>90.0
             get_element, wave, xrange, xx
             
             yrange = [min(mab-1.5*mab_err),max(mab+1.5*mab_err)]
             yrange[0] = yrange[0]<min(flux_ab[xx[0]:xx[1]])
             yrange[1] = yrange[1]>max(flux_ab[xx[0]:xx[1]])
             yrange = reverse(yrange)
             
             xhoff = 0.1 & yhoff = 0.37 & xsize = 0.5 & ysize = 0.4
             insetpos = [pos[0]+xhoff,pos[1]+yhoff,pos[0]+xhoff+xsize,pos[1]+yhoff+ysize]

             djs_plot, [0], [0], /nodata, /noerase, xthick=postthick1, ythick=postthick1, $
               xtitle='Observed Wavelength (\AA)', ytitle='m_{AB}', charsize=1.3, $
               charthick=postthick2, xsty=3, ysty=1, xrange=xrange, position=insetpos[*,0], $
               yrange=yrange
             oplot, wave, flux_ab, line=0.1
             oploterror, filtinfo[abgood].weff, mab, filtinfo[abgood].fwhm, mab_err, ps=8, thick=2.0, $
               color=djs_icolor('red'), errcolor=djs_icolor('red'), errthick=2.0
             djs_oplot, 5007*(1+z)*[1,1], !y.crange, line=2, thick=postthick2, color='blue'

          endfor
          
          dfpsclose
          spawn, 'gzip -f '+psname, /sh

       endelse
          
; write out the spectra

;      splog, 'Writing '+alphapath+'sdss_alpha_specdata.fits.gz'
;      mwrfits, sdssalpha, alphapath+'sdss_alpha_specdata.fits', /create
;      spawn, 'gzip -f '+alphapath+'sdss_alpha_specdata.fits', /sh

       alphafile = alphapath+'sdss_alpha_'+run+'.fits'
       splog, 'Writing '+alphafile
       mwrfits, sdssspectra, alphafile, /create
       spawn, 'gzip -f '+alphafile, /sh
       
;      splog, 'Writing '+alphapath+'sdss_alpha.fits'
;      mwrfits, sdssalpha_dennis, alphapath+'sdss_alpha.fits', /create
;      splog, 'Writing '+alphapath+'sdss_alpha.dat'
;      struct_print, sdssalpha_dennis, filename=alphapath+'sdss_alpha.dat'
       
       splog, /close

    endif

return
end
