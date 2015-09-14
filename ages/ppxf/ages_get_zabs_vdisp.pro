pro ages_get_zabs_vdisp, pass1, firstpass=firstpass, lastpass=lastpass, $
  doplot=doplot, index=index1, test=test
; jm09nov18ucsd - the point of this routine is to get reliable
;   estimates of the absorption-line redshift and velocity dispersion;
;   therefore, we do two iterations of PPXF; in the first iteration we
;   fit over a large wavelength range with lots of additive and
;   multiplicative degrees of freedom, which we subsequently use to
;   correct the spectrum shape for any calibration problems; we then
;   refit centered on the Balmer/4000-A break to lock down ZABS and
;   VDISP 

; jm10dec06ucsd - just use the solar metallicity templates
    
    common ages_catalogs, allcodes

    light = 2.99792458D5        ; speed of light [km/s]

; path names and emission-line file name, for masking
    version = ages_version(/ppxf_specfit)
    zabsvdisppath = ages_path(/spec1d)+'fluxed/zabs_vdisp/'+version+'/'
    linefile = ages_path(/ppxf)+'gandalf_elinelist_all.dat'

; figure out which passes we are going to fit
    if (n_elements(pass1) eq 0) then $
      allpass = ages_allpasses(/fluxed) else $ 
      allpass = string(pass1,format='(I3.3)')
    npass = n_elements(allpass)

    if (n_elements(firstpass) eq 0) then firstpass = 0
    if (n_elements(lastpass) eq 0) then lastpass = npass-1

; read the *solar* templates (see BUILD_AGES_PPXF_TEMPLATES); resample
; and convolve to AGES pixel size and instrumental resolution
    velscale = ages_ppxf_velscale()
    inst_vdisp = ages_ppxf_instvdisp()

    tempflux = read_ages_ppxf_templates(tempwave,ntemp=ntemp,$
      velscale=velscale,inst_vdisp=inst_vdisp,/solar)

; initial fitting parameters    
    vmaxshift = 300.0D ; maximum velocity tweak [km/s]
    sigmamax = 350.0D  ; maximum velocity dispersion [km/s]
    
    vdisp_guess = 150.0D
    degree = 5 ; additive Legendre polynomials
;   mdegree = 3 ; multiplicative Legendre polynomials

    min_fitwave2 = 3660.0 ; rest-frame
    max_fitwave2 = 4600.0 ; 4200.0
;   max_fitwave2 = 6600.0 ; 4200.0
    
; read the AGES codes
    if (n_elements(allcodes) eq 0) then allcodes = mrdfits($
      ages_path(/cat)+'catalog.codes.fits.gz',1)

; initialize the output data structure
    nfinalpix = 5000
    result1 = {$
      ages_id:      0L, $
      pass:          0, $
      aper:          0, $
      z:           0.0, $
      zabs:        0.0, $
      zabs_err:    0.0, $
      vdisp:       0.0, $
      vdisp_err:   0.0, $
      chi2:        0.0, $ ; reduced chi^2
      snr:         0.0}   ; mean S/N per pixel

; fit each plate separately
    t0 = systime(1)
    for ipass = firstpass, lastpass do begin

       suffix = string(allpass[ipass],format='(I0)')       
       if keyword_set(test) then suffix = suffix+'_test'
       outfile = zabsvdisppath+'zabs_vdisp_'+suffix+'.fits'

; read all the spectra on this plate but just fit the non-quasars at
; 0<z<1 
       spec1d = read_ages_spec1d(allpass[ipass],/fix)
       codes = allcodes[spec1d.ages_id]
       if (n_elements(index1) ne 0) then index = index1 else begin
;         index = lindgen(n_elements(spec1d.z))
          index = where((codes.gshort gt 0) and $
            (spec1d.z gt 0.001) and (spec1d.z lt 1.0))
       endelse 
       nobj = n_elements(index)
       if (nobj eq 0) then message, 'Problem here!'

       result = replicate(result1,nobj)
       result.ages_id = spec1d.ages_id[index]
       result.pass = spec1d.pass
       result.aper = spec1d.aper[index]
       result.z = spec1d.z[index]
       
; fit each object using PPXF
       t0 = systime(1)
       if (keyword_set(doplot) eq 0) then $
         psfile = zabsvdisppath+'qaplot_zabs_vdisp_'+suffix+'.ps'
       im_plotconfig, 8, pos, psfile=psfile, ymargin=[0.8,1.1], $
         charsize=1.6, thick=1
;      for iobj = 2, 2 do begin
       for iobj = 0, nobj-1 do begin
;         print, format='("Fitting object ",I0,"/",I0,A10,$)', $
;           iobj, nobj-1, string(13b)
          splog, allpass[ipass], iobj

          zobj = spec1d.z[index[iobj]]
          vsys = alog(zobj+1.0D)*light ; systemic velocity

; interpolate the galaxy spectrum logarithmically in wavelength
          good = where((spec1d.wave ge spec1d.minwave[index[iobj]]) and $
            (spec1d.wave le spec1d.maxwave[index[iobj]]),ngood)
          wave = spec1d.wave[good]
          flux = spec1d.flux[good,index[iobj]]
          ferr = spec1d.ferr[good,index[iobj]]

          lnflux = im_log_rebin(wave,flux,var=ferr^2,$
            outwave=lnwave,outvar=lnvar,vsc=velscale)
;         log_rebin, minmax(wave), flux, lnflux, lnwave, velscale=velscale
;         log_rebin, minmax(wave), ferr^2, lnvar, velscale=velscale
          lnferr = sqrt(abs(lnvar))

          lnrestwave = lnwave - alog(zobj+1.0D)
          lnrestflux = lnflux*(zobj+1.0D)
          lnrestferr = lnferr*(zobj+1.0D)
          npix = n_elements(lnrestwave)

; initial guesses: compute the *rest-frame* velocity offset between
; the templates and the data and restrict the wavelength range of the
; templates
          keep = where((exp(tempwave) gt exp(min(lnrestwave))-10.0) and $
            (exp(tempwave) lt exp(max(lnrestwave))+10.0),nkeep)
          if (nkeep eq 0) then message, 'Problem here!'
          fit_tempwave = tempwave[keep]
          fit_tempflux = tempflux[keep,*]
          voffset = -(min(lnrestwave)-min(fit_tempwave))*light
          start = [vsys*0.0+voffset,vdisp_guess]

; iteration 1: call PPFX with lots of polynomial freedom to take out
; any gross shape errors
          min_fitwave1 = min(wave)+100.0
          max_fitwave1 = 8300.0<(max(wave)-200.0)
          linepars = read_gandalf_elinelist(linefile)
          goodpixels = gandalf_mask_emission_lines(npix,vsys,$
            linepars,velscale,alog(wave[0]),velscale/light,$
            l_rf_range=[min_fitwave1,max_fitwave1]/(zobj+1.0D),$
            sigma=300.0)
;         djs_plot, exp(lnrestwave), lnrestflux, xsty=3, ysty=3
;         djs_oplot, exp(lnrestwave[goodpixels]), lnrestflux[goodpixels], psym=4, color='red'

          ages_ppxf, fit_tempflux, lnrestflux, lnrestferr, velscale, start, $
            sol, goodpixels=goodpixels, plot=doplot, moments=2, degree=degree, $
            weights=weights, polyweights=polyweights, mdegree=mdegree, $
            bestfit=bestfit, clean=0, /quiet, vmaxshift=vmaxshift, $
            sigmamax=sigmamax

; reconstruct the additive and multiplicative polynomials and correct
; the input spectrum accordingly
          xx = range(-1D,1D,npix)
;         mpoly = 1D & for jj = 1, mdegree do mpoly += legendre(xx,jj)*sol[6+jj]
          apoly = 0D & for jj = 0, degree do apoly += legendre(xx,jj)*polyweights[jj]

          lnrestflux_cor = lnrestflux-apoly
          lnrestferr_cor = lnrestferr

;         lnrestflux_cor = (lnrestflux-apoly)/(mpoly+(mpoly eq 0.0))*(mpoly ne 0.0)
;         lnrestferr_cor = lnrestferr/(mpoly+(mpoly eq 0.0))*(mpoly ne 0.0)
;         zero = where(mpoly le 0.0,nzero)
;         if (nzero ne 0) then lnrestferr_cor[zero] = 1E6
;         djs_plot, exp(lnrestwave), lnrestflux, xsty=3, ysty=3
;         djs_oplot, exp(lnrestwave), lnrestflux_cor, color='red'

; iteration 2: fit over a restricted wavelength range          

; just fit the Ca doublet region; also mask emission lines 
          goodpixels = gandalf_mask_emission_lines(npix,vsys,$
            linepars,velscale,alog(wave[0]),velscale/light,$
            l_rf_range=[min_fitwave2,max_fitwave2],sigma=300.0)
;         dfpsclose
;         djs_plot, exp(lnrestwave), lnrestflux_cor, xsty=3, ysty=3
;         djs_oplot, exp(lnrestwave[goodpixels]), lnrestflux_cor[goodpixels], psym=4, color='red'

          ages_ppxf, fit_tempflux, lnrestflux_cor, lnrestferr_cor, velscale, $
            start, sol, goodpixels=goodpixels, plot=doplot, $
            moments=2, degree=degree, error=err, bestfit=bestfit, $
            /clean, /quiet, vmaxshift=vmaxshift, sigmamax=sigmamax
          err = err*sqrt(sol[6]) ; scale by chi^2/dof

          zabs = exp((vsys+sol[0]-voffset)/light)-1.0D ; new absorption-line redshift
          zabs_err = zabs*(err[0]/light)
;         print, zobj, zabs, zabs_err
;
;         dfpsclose
;         djs_plot, exp(lnrestwave), lnrestflux_cor, ps=10, xrange=[3700,4050]
;         djs_oplot, exp(lnrestwave), bestfit, ps=10, color='green'
;         cc = get_kbrd(1)

          result[iobj].zabs = zabs
          result[iobj].zabs_err = zabs_err
          result[iobj].vdisp = sol[1]
          result[iobj].vdisp_err = err[1]
          result[iobj].chi2 = sol[6]

; get the mean S/N          
          inrange = where((flux gt 0.0) and (wave gt 6400.0) and $
            (wave lt 6600.0) and (ferr gt 0.0))
          result[iobj].snr = djs_median(flux[inrange]/ferr[inrange])

; make a QAplot
          xrange1 = [min_fitwave2,max_fitwave2]
          inrange = where((exp(lnrestwave) gt xrange1[0]) and $
            (exp(lnrestwave) lt xrange1[1]),ninrange)
          if (ninrange ne 0) then begin
             stats = im_stats(lnrestflux_cor[inrange],sigrej=5.0)
             yrange = 1E17*[-0.15,1.2]*stats.maxrej
          endif
          djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=3, ysty=3, $
            xrange=xrange1, yrange=yrange, xtitle='Rest Wavelength (\AA)', $
            ytitle='Flux (10^{-17} '+flam_units()+')', title='Pass '+$
            repstr(suffix,'_',' ')
          djs_oplot, exp(lnrestwave), 1E17*(lnrestflux_cor-bestfit), psym=10, color='light grey'
          djs_oplot, exp(lnrestwave), 1E17*lnrestflux_cor, psym=10
          djs_oplot, exp(lnrestwave), 1E17*bestfit, psym=10, color='red', thick=2
          djs_oplot, exp(lnrestwave[goodpixels]), 1E17*bestfit[goodpixels], psym=6, $
            symsize=0.3, color='blue'
          djs_oplot, !x.crange, [0,0], line=0, thick=1

          im_legend, [$
            'Fiber '+string(result[iobj].aper,format='(I3.3)'),$
            'z_{AGES}='+strtrim(string(result[iobj].z,format='(F12.5)'),2),$
            'z_{PPXF}='+strtrim(string(result[iobj].zabs,format='(F12.5)'),2),$
            '\sigma='+strtrim(string(result[iobj].vdisp,format='(F12.1)'),2)+'\pm'+$
            strtrim(string(result[iobj].vdisp_err,format='(F12.1)'),2)+' km s^{-1}', $
            '\chi^{2}_{\nu}='+strtrim(string(result[iobj].chi2,format='(F12.3)'),2)], $
            /left, /top, box=0, charsize=1.3, margin=0
       endfor
       if (keyword_set(doplot) eq 0) then $
         im_plotconfig, psfile=psfile, /gzip, /psclose
       splog, 'Total time = ', (systime(1)-t0)/60.0, ' minutes'
       im_mwrfits, result, outfile, /clobber
    endfor 

return
end
