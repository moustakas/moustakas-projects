;+
; NAME:
;       NFGS_STELLAR_MASS_KCORRECT
;
; PURPOSE:
;       Compute stellar masses for the NFGS using Blanton's
;       k-correct.  
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2006 Apr 30, U of A - written, based on
;         AGES_PHOTO_KCORR 
;-

pro nfgs_stellar_mass_kcorrect, nfgs, nfgsinfo, bigresult, debug=debug, write=write

    outpath = nfgs_path(/analysis)
    outfile = 'nfgs_stellar_mass_kcorrect.fits'
    
    if (n_elements(nfgs) eq 0L) then nfgs = read_nfgs()
    if (n_elements(nfgsinfo) eq 0L) then nfgsinfo = nfgs_read_info()
    ngalaxy = n_elements(nfgs)
    nnfgsinfo = n_elements(nfgsinfo)

; cosmology    
    
    red, omega0=0.3, omegalambda=0.7, h100=0.70
    omega0 = redomega0() & omegal = redomegal() & h100 = redh100()

    imf_chabrier_to_salpeter = +0.25
    light = 2.99792458D5        ; speed of light [km/s]
    light2 = 2.99792458D18      ; speed of light [A/s]

; define the filters for which we have observed-frame photometry  

    filterlist = [$
      'bessell_U',$
      'bessell_B',$
      'bessell_V',$
      'bessell_R',$
      'twomass_J',$
      'twomass_H',$
      'twomass_Ks'$
    ]+'.par'
       
    filtinfo = im_filterspecs(filterlist=filterlist,/verbose)
    nfilt = n_elements(filtinfo)

; Vega --> AB conversions

    U_vega2ab  = (k_vega2ab(filterlist='bessell_U.par',/kurucz))[0]
    B_vega2ab  = (k_vega2ab(filterlist='bessell_B.par',/kurucz))[0]
    V_vega2ab  = (k_vega2ab(filterlist='bessell_V.par',/kurucz))[0]
    R_vega2ab  = (k_vega2ab(filterlist='bessell_R.par',/kurucz))[0]
    J_vega2ab  = (k_vega2ab(filterlist='twomass_J.par',/kurucz))[0]
    H_vega2ab  = (k_vega2ab(filterlist='twomass_H.par',/kurucz))[0]
    Ks_vega2ab = (k_vega2ab(filterlist='twomass_Ks.par',/kurucz))[0]

; initialize the output data structure

    result = {$
      galaxy:                      ' ', $
      kcorr_chi2:               -999.0, $
      kcorr_mass:               -999.0, $
      kcorr_mets:               -999.0, $
      kcorr_intsfh:             -999.0, $
      kcorr_maggies:     fltarr(nfilt), $ ; input photometry
      kcorr_maggies_err: fltarr(nfilt), $ ; error in above
      kcorr_coeffs:          fltarr(5)  $ ; eigentemplate coefficients
      }
    bigresult = replicate(result,nnfgsinfo)
    result = replicate(result,ngalaxy)

    result.galaxy = nfgs.galaxy
    bigresult.galaxy = nfgsinfo.galaxy

    redshift = nfgs.distance*h100*100.0/light

; setup the k-correction photometry vectors

    mags_vega = transpose([ $
      [nfgs.u], $
      [nfgs.b], $
      [nfgs.v], $
      [nfgs.r], $
      [nfgs.twomass_j], $
      [nfgs.twomass_h], $
      [nfgs.twomass_ks] $
      ])
    mags_vega_err = transpose([ $
      [nfgs.u_err], $
      [nfgs.b_err], $
      [nfgs.v_err], $
      [nfgs.r_err], $
      [nfgs.twomass_j_err], $
      [nfgs.twomass_h_err], $
      [nfgs.twomass_ks_err] $
      ])

    maggies = mags_vega*0.0D ; needs to be double
    maggies_err = mags_vega*0.0D
    maggies_ivar = mags_vega*0.0D

    vega2ab = [u_vega2ab,b_vega2ab,v_vega2ab,r_vega2ab,j_vega2ab,h_vega2ab,ks_vega2ab]
    
    nphot = (size(mags_vega,/dimension))[0]
    for iphot = 0L, nphot-1L do begin

       good = where(mags_vega[iphot,*] gt -900.0,comp=flag,ncomp=nflag)
       maggies[iphot,good] = 10.0^(-0.4*(mags_vega[iphot,good] + vega2ab[iphot]))
       maggies_err[iphot,good] = 0.4*alog(10.0)*maggies[iphot,good]*mags_vega_err[iphot,good]

       if (nflag ne 0L) then begin
          maggies[iphot,flag] = 0.0D
          maggies_err[iphot,flag] = 1D16
       endif
       
    endfor

    maggies_ivar = 1.0/maggies_err^2.0

; compute the k-corrections      

    splog, 'Fitting k-correct templates.'
    t0 = systime(1)
    kcorrect, maggies, maggies_ivar, redshift, kcorrect, $
      band_shift=band_shift, filterlist=filterlist, $
      filterpath=filterpath, rmatrix=rmatrix, zvals=zvals, $
      lambda=lambda, vmatrix=vmatrix, coeffs=coeffs, chi2=chi2, $
      maxiter=maxiter, omega0=omega0, omegal0=omegal, mass=mass, $
      intsfh=intsfh, mets=mets, b300=b300, b1000=b1000, mtol=mtol
    splog, 'Total time = '+string((systime(1)-t0)/60.0,format='(G0.0)')+' minutes.'

    result.kcorr_chi2 = chi2

; store the basic results; convert the stellar mass to my IMF and
; Hubble constant; store the reduced chi2 below

    good = where((mass gt 0.0) and finite(mass),ngood)
    if (ngood ne 0L) then begin

       result[good].kcorr_mass = alog10(mass[good]) - alog10(h100) + imf_chabrier_to_salpeter
       result[good].kcorr_mets = mets[good]
       result[good].kcorr_intsfh = alog10(intsfh[good])

       result[good].kcorr_maggies = maggies[*,good]
       result[good].kcorr_maggies_err = maggies_err[*,good]
       result[good].kcorr_coeffs = coeffs[*,good]
       
    endif
    
; compute k-corrections and rest-frame apparent magnitudes; m_R = m_Q + K_QR

    plotsym, 0, 2.0, /fill, color=djs_icolor('red')

    restwave = lambda
       
    splog, 'Computing k-corrections.'
    t0 = systime(1)
    for i = 0L, ngalaxy-1L do begin

       if (result[i].kcorr_mass gt 0.0) then begin
   
          good = where((result[i].kcorr_maggies gt 0.0) and (result[i].kcorr_maggies_err lt 1E15),ngood)
          if (ngood ne 0L) then result[i].kcorr_chi2 = chi2[i] / float(ngood)

          wave = lambda*(1.0+redshift[i])
          restflux = vmatrix#coeffs[*,i] ; f_lambda
          flux = restflux/(1.0+redshift[i])
          flux_fnu = flux*wave^2.0/light2*10.0^(0.4*48.6) ; f_nu

          if keyword_set(debug) then begin

             if (ngood ne 0L) then begin

;               maggies_flam = maggies*10^(-0.4*48.6)*light2/rebin(filtinfo.weff,nfilter,ngalaxy)^2
                
                xrange = [min(filtinfo[good].weff-1.5*filtinfo[good].fwhm),$
                  max(filtinfo[good].weff+1.5*filtinfo[good].fwhm)] ;/(1.0+zobj)
                xrange[0] = xrange[0]>90.0
                get_element, wave, xrange, xx

                yrange = [min(maggies[good,i]-maggies_err[good,i]),max(maggies[good,i]+maggies_err[good,i])]
                yrange[0] = yrange[0]<min(flux_fnu[xx[0]:xx[1]])
                yrange[1] = yrange[1]>max(flux_fnu[xx[0]:xx[1]])
                
                plot, [0], [0], /nodata, xsty=3, ysty=3, charsize=1.5, charthick=2.0, xthick=2.0, $
                  ythick=2.0, xtitle='Observed Wavelength [\AA]', ytitle=textoidl('f_{\nu}'), yrange=yrange, $
                  xrange=xrange, title=strtrim(nfgs[i].galaxy,2)
                oploterror, filtinfo[good].weff, maggies[good,i], filtinfo[good].fwhm, maggies_err[good,i], ps=8
                oplot, wave, flux_fnu, line=0.1
                cc = get_kbrd(1)

             endif
         
          endif

       endif

    endfor
       
; copy RESULT into BIGRESULT and write out    
    
    doit = match_string(result.galaxy,bigresult.galaxy,/exact,findex=index)
;   niceprint, result.galaxy, bigresult[index].galaxy

    bigresult[index] = result

    if keyword_set(write) then begin
       splog, 'Writing '+outpath+outfile+'.'
       mwrfits, bigresult, outpath+outfile, /create
       spawn, ['gzip -f '+outpath+outfile], /sh
    endif
    
return
end
    
