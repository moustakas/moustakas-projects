pro write_ediscs_alpha_sample, ediscsdust, ediscsancillary, write=write
; jm07oct15nyu - written based on WRITE_EDISCS_OIILF_SAMPLE

    run = '08apr'
    nfinal = 100L
    
    red, h100=0.7, omega0=0.3, omega_lambda=0.7
    h100 = redh100()
    lsun = 3.826D33

    specfitpath = ediscs_path(/specfit)
    analysis_path = ediscs_path(/analysis)
    datapath = ediscs_path(/spec1d)
    alphapath = deep2_path(/projects)+'alpha/'

    stime0 = systime(1)
    if keyword_set(write) then begin
       splogfile = alphapath+'write_ediscs_alpha_sample_'+run+'.log'
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

    if (n_elements(ediscsdust) eq 0L) then begin
       specdatafile = specfitpath+'ediscs_specdata_'+ediscs_version(/specfit)+'.fits.gz'
       ancillaryfile = analysis_path+'ediscs_ancillary_data_'+ediscs_version(/ancillary)+'.fits.gz'
       splog, 'Reading '+specdatafile
       ediscsdust = mrdfits(specdatafile,1,/silent)
       splog, 'Reading '+ancillaryfile
       ediscsancillary = mrdfits(ancillaryfile,1,/silent)
    endif
    ngalaxy = n_elements(ediscsdust)

; ---------------------------------------------------------------------------
; select the parent sample
; ---------------------------------------------------------------------------

    snrcut = 5.0
    minmagerr = 0.2
    sample_zmin = 0.25
    sample_zmax = 0.90

    alpha = where((ediscsancillary.z gt sample_zmin) and (ediscsancillary.z lt sample_zmax) and $
;     (ediscsancillary.kcorr_mobs_ab[0] gt -900.0) and (ediscsancillary.kcorr_mobs_ab[1] gt -900.0) and $
;     (ediscsancillary.kcorr_mobs_ab[2] gt -900.0) and (ediscsancillary.kcorr_chi2 lt 5.0) and $
;     (abs(ediscsancillary.kcorr_mobs_ab[1]-ediscsancillary.kcorr_mobs_ab_synth[1]) lt 0.2) and $
;     (ediscsancillary.kcorr_mobs_ab_ivar[0] gt 1.0/minmagerr^2) and $
;     (ediscsancillary.kcorr_mobs_ab_ivar[1] gt 1.0/minmagerr^2) and $
;     (ediscsancillary.kcorr_mobs_ab_ivar[2] gt 1.0/minmagerr^2) and $
      (ediscsdust.oiii_4959_ew[0] gt 0.0) and (ediscsdust.oiii_4959_ew[1] gt 0.0) and $ ; EW>0; EW_err>0
      (ediscsdust.oiii_4959_ew[0]/ediscsdust.oiii_4959_ew[1] gt snrcut) and $           ; hard S/N cut
      (ediscsdust.oiii_5007_ew[0] gt 0.0) and (ediscsdust.oiii_5007_ew[1] gt 0.0) and $ ; EW>0; EW_err>0
      (ediscsdust.oiii_5007_ew[0]/ediscsdust.oiii_5007_ew[1] gt snrcut),nalpha)         ; hard S/N cut

; throw away objects that are within LAMBDA_TOL of a night sky line

    lambda_tol = 1.0 ; wavelength tolerance for proximity to a sky line
    
    lambda_oiii_4959 = 4958.911
    lambda_oiii_5007 = 5006.843

    readcol, getenv('ISPEC_DIR')+'/etc/skylist_o2_oh.dat', lambda_sky, format='D', /silent
    nsky = n_elements(lambda_sky)
    
    allgood = bytarr(nsky,nalpha)
    for ii = 0L, nalpha-1L do begin
       for isky = 0L, nsky-1L do begin         
          allgood[isky,ii] = (abs(lambda_oiii_4959*(1.0+ediscsancillary[alpha[ii]].z)-lambda_sky[isky]) gt lambda_tol) and $
            (abs(lambda_oiii_5007*(1.0+ediscsancillary[alpha[ii]].z)-lambda_sky[isky]) gt lambda_tol)
       endfor
    endfor

    good = where((nsky-total(allgood,1) eq 0.0),nalpha)
    alpha = alpha[good]

    splog, 'alpha sample: '+string(nalpha,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+' ('+$
      strtrim(string(100.0*nalpha/(ngalaxy),format='(F12.1)'),2)+'%).'

    zstats = im_stats(ediscsancillary[alpha].z)
    splog, '   Redshift: ['+strtrim(string(zstats.min,format='(F12.2)'),2)+'-'+$
      strtrim(string(zstats.max,format='(F12.2)'),2)+'] '+$
      strtrim(string(zstats.median,format='(F12.2)'),2)+' ('+$
      strtrim(string(zstats.mean,format='(F12.2)'),2)+'+/-'+$
      strtrim(string(zstats.sigma,format='(F12.2)'),2)+')'

; ---------------------------------------------------------------------------    
; write out
; ---------------------------------------------------------------------------    
    
    if keyword_set(write) then begin

       splog, 'Computing luminosity distances'
       dlum = dluminosity(ediscsancillary[alpha].z,/cm)
       
       splog, 'Appending additional tags'
       moretags = replicate({oiii_4959_flux: [-999.0,-999.0], oiii_4959_lum: [-999.0,-999.0], $
         oiii_5007_flux: [-999.0,-999.0], oiii_5007_lum: [-999.0,-999.0]},nalpha)
       moretags.oiii_4959_flux[0] = ediscsancillary[alpha].cflux_4959*ediscsdust[alpha].oiii_4959_ew[0]
       moretags.oiii_4959_flux[1] = ediscsancillary[alpha].cflux_4959*ediscsdust[alpha].oiii_4959_ew[1]
       moretags.oiii_4959_lum[0] = moretags.oiii_4959_flux[0]*4.0*!dpi*dlum*dlum/lsun
       moretags.oiii_4959_lum[1] = moretags.oiii_4959_flux[1]*4.0*!dpi*dlum*dlum/lsun
;      moretags.oiii_4959_lum[0] = alog10(moretags.oiii_4959_flux[0]*4.0*!dpi*dlum*dlum/lsun)
;      moretags.oiii_4959_lum[1] = moretags.oiii_4959_flux[1]/moretags.oiii_4959_flux[0]/alog(10.0)

       moretags.oiii_5007_flux[0] = ediscsancillary[alpha].cflux_5007*ediscsdust[alpha].oiii_5007_ew[0]
       moretags.oiii_5007_flux[1] = ediscsancillary[alpha].cflux_5007*ediscsdust[alpha].oiii_5007_ew[1]
       moretags.oiii_5007_lum[0] = moretags.oiii_5007_flux[0]*4.0*!dpi*dlum*dlum/lsun
       moretags.oiii_5007_lum[1] = moretags.oiii_5007_flux[1]*4.0*!dpi*dlum*dlum/lsun
;      moretags.oiii_5007_lum[0] = alog10(moretags.oiii_5007_flux[0]*4.0*!dpi*dlum*dlum/lsun)
;      moretags.oiii_5007_lum[1] = moretags.oiii_5007_flux[1]/moretags.oiii_5007_flux[0]/alog(10.0)

       ediscsalpha1 = struct_addtags(struct_addtags(ediscsdust[alpha],$
         struct_trimtags(ediscsancillary[alpha],except=['GALAXY','SPECFILE'])),$
         moretags)

; select the NFINAL galaxies with the highest [O III] luminosities
       
       final = (reverse(sort(ediscsalpha1.oiii_5007_lum[0])))[0L:nfinal-1L]
       ediscsalpha = ediscsalpha1[final]
       ediscsalpha_dennis = struct_trimtags(ediscsalpha,select=['GALAXY','Z',$
         'RA','DEC','PHOT_B','PHOT_V','PHOT_R','PHOT_I','PHOT_J','PHOT_K',$
         'OIII_5007_LUM'])
       
; k-correct preliminaries

       vname = 'default.nolines'
       light = 2.99792458D18    ; speed of light [A/s]

       k_load_vmatrix, vmatrix, lambda, vfile=vfile, vpath=vpath, lfile=lfile, vname=vname
       plotsym, 0, 1.5, /fill
       
       filters = [['FORS_B','FORS_V','FORS2_R','FORS_I']+'_ccd',['SOFI_J','SOFI_Ks']]+'_atm.par'
       filtinfo = im_filterspecs(filterlist=filters,/verbose)
       nfilter = n_elements(filtinfo)

; plot the ediscs spectra

       ediscsspectra = struct_addtags(ediscsalpha_dennis,replicate({wave: fltarr(2201), $
         flux: fltarr(2201), ivar: fltarr(2201)},nfinal))
       
       psname = alphapath+'ediscs_alpha_'+run+'.ps'
       dfpsplot, psname, /color, /landscape

       pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, ymargin=[0.5,1.1], $
         xmargin=[1.1,0.2], position=pos, /normal

       for igal = 0L, n_elements(ediscsalpha)-1L do begin
       
          galaxy = strtrim(ediscsalpha[igal].galaxy,2)
          z = ediscsalpha[igal].z

          spec = rd1dspec(ediscsalpha[igal].specfile,datapath=datapath)

          wave = spec.wave
          flux = 1D17*spec.spec
          ivar = 1.0/spec.sigspec^2.0
          
          good = where(ivar gt 0.0,ngood)
          wave = wave[good]
          flux = flux[good]
          ivar = ivar[good]
          ferr = 1.0/sqrt(ivar)
          npix = n_elements(wave)
       
          ediscsspectra[igal].wave = wave
          ediscsspectra[igal].flux = flux
          ediscsspectra[igal].ivar = ivar

          restwave = wave/(1.0+z)
          restflux = flux*(1.0+z)
          restivar = ivar/(1.0+z)^2.0
          restferr = ferr*(1.0+z)

; make the plot

          title = galaxy+', z = '+strtrim(string(z,format='(F12.4)'),2)+$
            ', L[O III] = '+strtrim(string(ediscsalpha[igal].oiii_5007_lum[0],format='(E10.2)'),2)+' L_{\odot}'
;           ', L[O III] = '+strtrim(string(10.0^ediscsalpha[igal].oiii_5007_lum[0],format='(E10.2)'),2)+' L_{\odot}'
;         title = galaxy+', z = '+strtrim(string(z,format='(F12.4)'),2)+', !c '+$
;           '(\alpha,\delta) = '+strtrim(string(ediscsalpha[igal].ra,format='(F15.8)'),2)+', '+$
;           strtrim(string(ediscsalpha[igal].dec,format='(F15.8)'),2)

          xrange = [4800.0,5050.0]
          get_element, restwave, xrange, xx
          yrange = [0.8*min(smooth(restflux[xx[0]:xx[1]],20)),1.05*max(restflux[xx[0]:xx[1]])]
;         yrange = minmax(restflux[xx[0]:xx[1]])
          
          djs_plot, [0], [0], /nodata, xsty=1, ysty=1, position=pos, $
            xthick=postthick1, ythick=postthick1, charthick=postthick2, $
            xtitle='Rest Wavelength (\AA)', charsize=1.4, xrange=xrange, $
            yrange=yrange, ytitle='Flux (10^{-17} '+flam_units()+')', title=title
          djs_oplot, restwave, restflux, thick=postthick1, color='grey', ps=10

; make an inset with the SED fit

          restwave = lambda
          restflux = vmatrix#ediscsalpha[igal].coeffs ; f_lambda
          restflux_fnu = restflux*restwave^2.0/light*10.0^(0.4*48.6) ; f_nu
          restflux_ab = -2.5*alog10(restflux_fnu>1D-32)

          wave = restwave*(1.0+z)
          flux = restflux/(1.0+z)
          flux_fnu = flux*wave^2.0/light*10.0^(0.4*48.6) ; f_nu
          flux_ab = -2.5*alog10(flux_fnu>1D-32)

          abgood = where((ediscsalpha[igal].abmaggies_ivar gt 0.0),nabgood)
          mab = -2.5*alog10(ediscsalpha[igal].abmaggies[abgood])
          mab_err = 1.0/sqrt(ediscsalpha[igal].abmaggies_ivar[abgood])/$
            ediscsalpha[igal].abmaggies[abgood]/alog(10.0)
          
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

;         label = ['log (M/M'+sunsymbol()+') = '+strtrim(string(ediscsalpha[igal].kcorr_mass,format='(F12.2)'),2),$
;           '\chi^{2}_{\nu} = '+strtrim(string(ediscsalpha[igal].kcorr_chi2,format='(F12.2)'),2)]
;         legend, textoidl(label), /left, /top, box=0, charsize=1.3, charthick=2.0

;         cc = get_kbrd(1)
          
       endfor

       dfpsclose
       spawn, 'gzip -f '+psname, /sh

; write out the spectra

;      splog, 'Writing '+alphapath+'ediscs_alpha_specdata.fits.gz'
;      mwrfits, ediscsalpha, alphapath+'ediscs_alpha_specdata.fits', /create
;      spawn, 'gzip -f '+alphapath+'ediscs_alpha_specdata.fits', /sh

       alphafile = alphapath+'ediscs_alpha_'+run+'.fits'
       splog, 'Writing '+alphafile
       mwrfits, ediscsspectra, alphafile, /create
       spawn, 'gzip -f '+alphafile, /sh
       
;      splog, 'Writing '+alphapath+'ediscs_alpha.fits'
;      mwrfits, ediscsalpha_dennis, alphapath+'ediscs_alpha.fits', /create
;      splog, 'Writing '+alphapath+'ediscs_alpha.dat'
;      struct_print, ediscsalpha_dennis, filename=alphapath+'ediscs_alpha.dat'

       splog, /close

    endif

stop    
    
return
end
