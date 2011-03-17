pro write_deep2_alpha_sample, ispec, ancillary, write=write, nopsfile=nopsfile
; jm07oct15nyu - written based on WRITE_DEEP2_OIILF_SAMPLE
; jm08mar31nyu - updated for apr08 run
; jm08sep04nyu - updated for sep08 run (major overhaul and new catalogs)

    light = 2.99792458D5        ; speed of light [km/s]

    run = '08sep' ; '08apr' ; '07nov'
    nfinal = 100L
    
    red, h100=0.7, omega0=0.3, omega_lambda=0.7
    h100 = redh100()
    lsun = 3.826D33

    specfitpath = deep2_path(/specfit)
    analysis_path = deep2_path(/analysis)
    alphapath = deep2_path(/projects)+'alpha/'

    stime0 = systime(1)
    if keyword_set(write) then begin
       splogfile = alphapath+'write_deep2_alpha_sample_'+run+'.log'
       splog, filename=splogfile
       splog, 'Log file '+splogfile+' opened '+systime()
       splog, 'IDL version: ' + string(!version,format='(99(A," "))')
       spawn, 'uname -a', uname
       splog, 'UNAME: '+uname[0]
    endif
       
; ---------------------------------------------------------------------------
; read the emission-line data
; ---------------------------------------------------------------------------

; code corresponding to 08sep run
    if (n_elements(ancillary) eq 0L) then ancillary = read_deep2(/kcorr)
    if (n_elements(ispec) eq 0L) then begin
       ispec1 = read_deep2(/ispec)
       ispecmore = {$
;        oiii_4959_flux: [-999.0,-999.0], $
;        oiii_4959_lum:  [-999.0,-999.0], $
;        oiii_4959_peak:          -999.0, $
         oiii_5007_flux: [-999.0,-999.0], $
         oiii_5007_lum:  [-999.0,-999.0], $
         oiii_5007_peak:          -999.0}
       ispecmore = replicate(ispecmore,n_elements(ispec1))
       ispec = struct_addtags(ispec1,ispecmore)
    endif

;; code corresponding to 08apr run
;   if (n_elements(ispec) eq 0L) then begin
;      specdatafile = specfitpath+'deep2_specdata_'+deep2_version(/specdata)+'.fits.gz'
;      ancillaryfile = analysis_path+'deep2_kcorr_'+deep2_version(/kcorr)+'.fits.gz'
;      splog, 'Reading '+specdatafile
;      ispec = mrdfits(specdatafile,1,/silent)
;      splog, 'Reading '+ancillaryfile
;      ancillary = mrdfits(ancillaryfile,1,/silent)
;   endif
;; code corresponding to 07nov run
;   if (n_elements(ispec) eq 0L) then begin
;      splog, 'Reading '+specfitpath+'deep2_specdata_'+deep2_version(/specdata)+'.fits.gz'
;      ispec = mrdfits(specfitpath+'deep2_specdata.fits.gz',1,/silent)
;      splog, 'Reading '+analysis_path+'deep2.dr3.kcorr.fits.gz'
;      ancillary = mrdfits(analysis_path+'deep2.dr3.kcorr.fits.gz',1,/silent)
;   endif
    ngalaxy = n_elements(ispec)

; ---------------------------------------------------------------------------
; select the parent sample
; ---------------------------------------------------------------------------

    snrcut = 5.0
    minmagerr = 0.2
    sample_zmin = 0.4
    sample_zmax = 0.8

    alpha = where((ancillary.z gt sample_zmin) and (ancillary.z lt sample_zmax) and $
      (ispec.oiii_4959_ew[0] gt 0.0) and (ispec.oiii_4959_ew[1] gt 0.0) and $ ; EW>0; EW_err>0
      (ispec.oiii_4959_ew[0]/ispec.oiii_4959_ew[1] gt snrcut) and $           ; hard S/N cut
      (ispec.oiii_5007_ew[0] gt 0.0) and (ispec.oiii_5007_ew[1] gt 0.0) and $ ; EW>0; EW_err>0
      (ispec.oiii_5007_ew[0]/ispec.oiii_5007_ew[1] gt snrcut),nalpha)         ; hard S/N cut

    alpha_ispec = ispec[alpha]
    alpha_ancillary = ancillary[alpha]
    
; throw away objects that are within LAMBDA_TOL of a night sky line

    lambda_tol = 1.0 ; wavelength tolerance for proximity to a sky line
    
    lambda_oiii_4959 = 4958.911
    lambda_oiii_5007 = 5006.843

    readcol, getenv('ISPEC_DIR')+'/etc/skylist_o2_oh.dat', lambda_sky, format='D', /silent
    nsky = n_elements(lambda_sky)

    splog, 'Rejecting sky-contaminated objects'
    allgood = bytarr(nsky,nalpha)
    for ii = 0L, nalpha-1L do begin
       for isky = 0L, nsky-1L do begin         
          allgood[isky,ii] = (abs(lambda_oiii_4959*(1.0+alpha_ancillary[ii].z)-lambda_sky[isky]) gt lambda_tol) and $
            (abs(lambda_oiii_5007*(1.0+alpha_ancillary[ii].z)-lambda_sky[isky]) gt lambda_tol)
       endfor
    endfor

; put everything together
    
    good = where((nsky-total(allgood,1) eq 0.0),nalpha)
    alpha = alpha[good]

    alpha_ispec = ispec[alpha]
    alpha_ancillary = ancillary[alpha]
;   alpha_ispec = struct_addtags(ancillary[alpha],ispec[alpha])

; some statistics    
    
    splog, 'alpha sample: '+string(nalpha,format='(I0)')+'/'+$
      string(ngalaxy,format='(I0)')+' ('+$
      strtrim(string(100.0*nalpha/(ngalaxy),format='(F12.1)'),2)+'%).'

    zstats = im_stats(alpha_ispec.z)
    splog, '   Redshift: ['+strtrim(string(zstats.min,format='(F12.2)'),2)+'-'+$
      strtrim(string(zstats.max,format='(F12.2)'),2)+'] '+$
      strtrim(string(zstats.median,format='(F12.2)'),2)+' ('+$
      strtrim(string(zstats.mean,format='(F12.2)'),2)+'+/-'+$
      strtrim(string(zstats.sigma,format='(F12.2)'),2)+')'

; compute the integrated flux and luminosity

    dlum = dluminosity(alpha_ispec.z,/cm)
    
;   alpha_ispec.oiii_4959_flux[0] = alpha_ancillary.cflux_4959*alpha_ispec.oiii_4959_ew[0]
;   alpha_ispec.oiii_4959_flux[1] = alpha_ancillary.cflux_4959*alpha_ispec.oiii_4959_ew[1]
;   alpha_ispec.oiii_4959_lum[0]  = alpha_ispec.oiii_4959_flux[0]*4.0*!dpi*dlum*dlum/lsun
;   alpha_ispec.oiii_4959_lum[1]  = alpha_ispec.oiii_4959_flux[1]*4.0*!dpi*dlum*dlum/lsun

    alpha_ispec.oiii_5007_flux[0] = alpha_ancillary.cflux_5007*alpha_ispec.oiii_5007_ew[0]
    alpha_ispec.oiii_5007_flux[1] = alpha_ancillary.cflux_5007*alpha_ispec.oiii_5007_ew[1]
    alpha_ispec.oiii_5007_lum[0]  = alpha_ispec.oiii_5007_flux[0]*4.0*!dpi*dlum*dlum/lsun
    alpha_ispec.oiii_5007_lum[1]  = alpha_ispec.oiii_5007_flux[1]*4.0*!dpi*dlum*dlum/lsun

; compute the 4959 and 5007 peak-luminosity

    restwave = findgen(200)+4900.0
;   restwave = findgen(500)+4600.0
    linemodel = irestore_speclinefit(alpha_ispec,outwave=restwave,$
      linename='OIII_5007',/nocontinuum)
    for jj = 0L, nalpha-1L do alpha_ispec[jj].oiii_5007_peak = $
      max(linemodel[jj].restflux)
;   plot, linemodel[0].restwave, linemodel[0].restflux, ps=10

    final = (reverse(sort(alpha_ispec.oiii_5007_peak)))[0L:nfinal-1L]
    alpha_final = struct_trimtags(alpha_ispec[final],select=['GALAXY','Z',$
      'RA','DEC','MAGB','MAGR','MAGI','OIII_5007_LUM','OIII_5007_PEAK'])
    alpha_final.galaxy = repstr(alpha_final.galaxy,'spec1d.','deep2.')
    
; ---------------------------------------------------------------------------    
; write out
; ---------------------------------------------------------------------------    

    if keyword_set(write) then begin

; k-correct preliminaries

;      vname = 'default'
       vname = 'default.nolines'
       light = 2.99792458D18    ; speed of light [A/s]

       k_load_vmatrix, vmatrix, lambda, vname=vname
       plotsym, 0, 1.5, /fill
       
       filters = 'deep_'+['B','R','I']+'.par'
       filtinfo = im_filterspecs(filterlist=filters,/verbose)
       nfilter = n_elements(filtinfo)

; plot the deep2 spectra

       alpha_spectra = struct_addtags(alpha_final,replicate({wave: fltarr(8192), $
         flux: fltarr(8192), ivar: fltarr(8192)},nfinal))

       if (not keyword_set(nopsfile)) then begin
          
          datapath = deep2_path(/dr3)
          psname = alphapath+'deep2_alpha_'+run+'.ps'
          dfpsplot, psname, /color, /landscape
          im_plotfaves, /postscript

          pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, ymargin=[0.5,1.1], $
            xmargin=[1.1,0.2], position=pos, /normal

          spec = read_deep2_specfit(alpha_ispec[final].file,/obswave)
          alpha_spectra.wave = reform(spec[*,0,*])
          alpha_spectra.flux = reform(spec[*,1,*])
          alpha_spectra.ivar = reform(spec[*,5,*])
          
          for igal = 0L, n_elements(alpha_final)-1L do begin
             
             galaxy = strtrim(alpha_final[igal].galaxy,2)
             z = alpha_final[igal].z
             
             restwave = alpha_spectra[igal].wave/(1.0+z)
             restflux = alpha_spectra[igal].flux*(1.0+z)
             restivar = alpha_spectra[igal].ivar/(1.0+z)^2.0

; make the plot

             title = galaxy+', z = '+strtrim(string(z,format='(F12.4)'),2)+$
               ', L[O III] = '+strtrim(string(alpha_final[igal].oiii_5007_lum[0],$
               format='(E10.2)'),2)+' L_{\odot}'

             xrange = [4800.0,5050.0]
             get_element, restwave, xrange, xx
             yrange = minmax(restflux[xx[0]:xx[1]])
             yrange[0] = yrange[0]*0.9
             
             djs_plot, [0], [0], /nodata, xsty=1, ysty=1, position=pos, $
               xtitle='Rest Wavelength (\AA)', charsize=1.4, xrange=xrange, $
               yrange=yrange, ytitle='Flux (ADU)', title=title
             djs_oplot, restwave, restflux, color='grey', ps=10

; make an inset with the SED fit

             restwave = k_lambda_to_centers(lambda)
             restflux = vmatrix#alpha_ancillary[final[igal]].coeffs        ; f_lambda
             restflux_fnu = restflux*restwave^2.0/light*10.0^(0.4*48.6) ; f_nu
             restflux_ab = -2.5*alog10(restflux_fnu>1D-32)

             wave = restwave*(1.0+z)
             flux = restflux/(1.0+z)
             flux_fnu = flux*wave^2.0/light*10.0^(0.4*48.6) ; f_nu
             flux_ab = -2.5*alog10(flux_fnu>1D-32)

             abgood = where((alpha_ancillary[final[igal]].abmaggies_ivar gt 0.0),nabgood)
             mab = -2.5*alog10(alpha_ancillary[final[igal]].abmaggies[abgood])
             mab_err = 1.0/sqrt(alpha_ancillary[final[igal]].abmaggies_ivar[abgood])/$
               alpha_ancillary[final[igal]].abmaggies[abgood]/alog(10.0)
             
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

             djs_plot, [0], [0], /nodata, /noerase, xtitle='Observed Wavelength (\AA)', $
               ytitle='m_{AB}', charsize=1.3, xsty=3, ysty=1, xrange=xrange, $
               position=insetpos[*,0], yrange=yrange
             oplot, wave, flux_ab, line=0.1
             oploterror, filtinfo[abgood].weff, mab, filtinfo[abgood].fwhm, mab_err, ps=8, thick=2.0, $
               color=djs_icolor('red'), errcolor=djs_icolor('red'), errthick=2.0
             djs_oplot, 5007*(1+z)*[1,1], !y.crange, line=2, thick=postthick2, color='blue'

          endfor
          
          dfpsclose
          spawn, 'gzip -f '+psname, /sh
          im_plotfaves

       endif 
          
; write out the spectra

       alphafile = alphapath+'deep2_alpha_'+run+'.fits'
       splog, 'Writing '+alphafile
       mwrfits, alpha_spectra, alphafile, /create
       spawn, 'gzip -f '+alphafile, /sh

       splog, /close

    endif

return
end
