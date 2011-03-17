pro ediscs_sdss_vmax
; jm07aug16nyu - 

    if (n_elements(vname) eq 0L) then vname = 'default'
    if (n_elements(band_shift) eq 0L) then band_shift = 0.0

    red, omega0=0.3, omegalambda=0.7, h100=0.70
    omega0 = redomega0() & omegal = redomegal() & h100 = redh100()
    imf_chabrier_to_salpeter = +0.25

    survey_filter = 'sdss_r0.par'
    sample_zmin = 0.02
    sample_zmax = 0.09
    dr4_area = 4783.0*!radeg^2 ; [sr]
    sdss_mlimit = [14.5,17.77]

    path = ediscs_path(/sdss)
    sdss = mrdfits(path+'EDisCS_comparison_photometry_SED_table.fits',1,/silent)
    sdss = sdss[1000:2000]
    ngalaxy = n_elements(sdss)

; ---------------------------------------------------------------------------    
;;  photo = mrdfits(sdss_path(/dr4)+'gal_photo_trim_dr4_v5_1c.fit',1,/silent)
;   photo = mrdfits(sdss_path()+'sdss_main_dr4.fits.gz',1,/silent)
;;  casid_extract, sdss.objid, run, rerun, camcol, field, sdssid, sky_version=sky, first_field=first
;   photoid = casid(photo.run,photo.rerun,photo.camcol,photo.field,photo.id,$
;     sky_version=1L,first_field=0L)
;   help, cmset_op(photoid,'AND',sdss.objid,/index), sdss
    
; ---------------------------------------------------------------------------    
    
; define the observed and rest-frame filters

    splog, 'Defining filter functions with observed photometry:'
    obsfilters = 'sdss_'+['u0','g0','r0','i0','z0']+'.par'
    obsbands = ['u','g','r','i','z']
    nobsfilter = n_elements(obsfilters)

    splog, 'Loading rest-frame filters:'
    restfilters = ['bessell_'+['U','B','V','R','I'],'sdss_'+['u0','g0','r0','i0','z0']]+'.par'
    restbandpasses = ['U','B','V','R','I','sdss_u','sdss_g','sdss_r','sdss_i','sdss_z']
    restvega = [1,1,1,1,1,0,0,0,0,0]
    nrestfilter = n_elements(restfilters)

    restvega2ab = k_vega2ab(filterlist=restfilters,/kurucz)
    notvega = where((restvega eq 0L),nnotvega)
    if (nnotvega ne 0L) then restvega2ab[notvega] = 0.0

; convert the observed photometry to maggies and store

    splog, 'Converting observed photometry to maggies.'
    obsmaggies = convert_to_maggies(sdss,obsbands,fltarr(nobsfilter),$
      maggies_invvar=obsmaggies_ivar,mags=mobs_ab,err_mags=mobs_ab_err,$
      ivar_mags=mobs_ab_ivar,err_suffix='err',badvalue=-999.0,/nanomaggies)
    setzero = where(obsmaggies_ivar eq 1D-32,nsetzero)
    if (nsetzero ne 0L) then obsmaggies_ivar[setzero] = 0.0

; compute k-corrections

;   splog, 'Fitting k-correct templates.'
;   t0 = systime(1)
;   kcorrect, obsmaggies, obsmaggies_ivar, sdss.redshift, kcorrect, $
;     band_shift=band_shift, filterlist=obsfilters, filterpath=filters_path, $
;     rmatrix=rmatrix, zvals=zvals, lambda=lambda, vmatrix=vmatrix, $
;     coeffs=coeffs, rmaggies=obsmaggies_synth, chi2=chi2, maxiter=maxiter, omega0=omega0, $
;     omegal0=omegal, mass=mass, intsfh=intsfh, mets=mets, b300=b300, $
;     b1000=b1000, mtol=mtol, vname=vname, /silent
;   splog, 'Total time = '+string((systime(1)-t0)/60.0,format='(G0.0)')+' minutes.'

; compute the k-corrections

    t0 = systime(1)
    splog, 'Computing UBVRI k-corrections.'
    ubvri_kcorrect = sdss2bessell(sdss.redshift,nmgy=obsmaggies,ivar=obsmaggies_ivar,$
      band_shift=band_shift,chi2=chi2,coeffs=ubvri_coeffs,rmaggies=rmaggies,$
      omaggies=omaggies,oivar=oivar,vname=vname,mass=mass,mtol=mtol,absmag=ubvri_absmag,$
      amivar=ubvri_amivar,omega0=omega0,omegal0=omegal,mets=mets,b300=b300,b1000=b1000,$
      intsfh=intsfh,/vega)      ; note VEGA

    splog, 'Computing ugriz k-corrections.'
    ugriz_kcorrect = sdss_kcorrect(sdss.redshift,nmgy=obsmaggies,ivar=obsmaggies_ivar,$
      band_shift=band_shift,chi2=chi2,coeffs=coeffs,rmaggies=rmaggies,omaggies=omaggies,$
      oivar=oivar,vname=vname,mass=mass,mtol=mtol,absmag=ugriz_absmag,$
      amivar=ugriz_amivar,omega0=omega0,omegal0=omegal,$
      mets=mets,b300=b300,b1000=b1000,intsfh=intsfh)
    splog, 'Total time = '+string((systime(1)-t0)/60.0,format='(G0.0)')+' minutes.'

    ugriz_appmag_obs       = fltarr(5,ngalaxy)
    ugriz_appmag_obs_ivar  = fltarr(5,ngalaxy)
    ugriz_appmag_obs_synth = fltarr(5,ngalaxy)

    for iobs = 0L, 4L do begin
       pos = where((oivar[iobs,*] gt 0.0) and (omaggies[iobs,*] gt 0.0),npos)
       if (npos ne 0L) then begin
          ugriz_appmag_obs[iobs,pos]      = reform(-2.5*alog10(omaggies[iobs,pos]))                    ; maggies used in the fit
          ugriz_appmag_obs_ivar[iobs,pos] = (0.4*alog(10.0))^2.0*reform(oivar[iobs,pos]*omaggies[iobs,pos]^2.0)  
       endif
       rpos = where((rmaggies[iobs,*] gt 0.0),nrpos)
       if (nrpos ne 0L) then ugriz_appmag_obs_synth[iobs,rpos] = reform(-2.5*alog10(rmaggies[iobs,rpos])) ; maggies reconstructed from the fit
    endfor
    
    ubvri_appmag_rest = ubvri_absmag + rebin(reform(lf_distmod(sdss.z,omega0=omega0,omegal0=omegal0),1,ngalaxy),5,ngalaxy) ; apparent magnitudes
    ugriz_appmag_rest = ugriz_absmag + rebin(reform(lf_distmod(sdss.z,omega0=omega0,omegal0=omegal0),1,ngalaxy),5,ngalaxy)

    ubvri_absmag = ubvri_absmag + 5.0*alog10(h100) ; h=1 --> h=0.7    
    ugriz_absmag = ugriz_absmag + 5.0*alog10(h100)

; store everything

    result = {$
      zmin:                                   -999.0, $
      zmax:                                   -999.0, $
      vmax:                                   -999.0, $
      vol:                                    -999.0, $
      kcorr_mass:                             -999.0, $
      kcorr_mets:                             -999.0, $
      kcorr_intsfh:                           -999.0, $
      kcorr_coeffs:                        fltarr(5), $ ; eigentemplate coefficients
      kcorr_chi2:                             -999.0, $ ; reduced chi2
      kcorr_kcorr:         fltarr(nrestfilter)-999.0, $ ; k-correction
      kcorr_mobs_ab:        fltarr(nobsfilter)-999.0, $ ; "actual observed photometry"
      kcorr_mobs_ab_ivar:         fltarr(nobsfilter), $ ; inverse variance in MOBS_AB (=0 if synthesized)
      kcorr_mobs_ab_synth:        fltarr(nobsfilter), $ ; "synthesized observed photometry"
      kcorr_mrest:         fltarr(nrestfilter)-999.0, $ ; k-corrected, rest frame
      kcorr_mrest_ivar:          fltarr(nrestfilter)}   ; inverse variance in KCORR_MREST (=0 if synthesized)

    for i = 0L, nrestfilter-1L do result = create_struct(result,$
      'M_'+restbandpasses[i],-999.0,'M_'+restbandpasses[i]+'_ivar',0.0)
    result = replicate(result,ngalaxy)
    
    good = where((mass gt 0.0) and finite(mass),ngood)
    if (ngood ne 0L) then begin
       result[good].kcorr_mass   = alog10(mass[good]) - 2.0*alog10(h100) + imf_chabrier_to_salpeter
       result[good].kcorr_mets   = mets[good]
       result[good].kcorr_intsfh = alog10(intsfh[good])
       result[good].kcorr_coeffs = coeffs[*,good]
       result[good].kcorr_chi2   = chi2[good]
    endif
    
    result.kcorr_mobs_ab       = ugriz_appmag_obs
    result.kcorr_mobs_ab_ivar  = ugriz_appmag_obs_ivar
    result.kcorr_mobs_ab_synth = ugriz_appmag_obs_synth
    
    result.kcorr_kcorr      = [ubvri_kcorrect,ugriz_kcorrect]
    result.kcorr_mrest      = [ubvri_appmag_rest,ugriz_appmag_rest]
    result.kcorr_mrest_ivar = [ubvri_amivar,ugriz_amivar]

    result.M_U = reform(ubvri_absmag[0,*]) ; Vega, h=0.7
    result.M_B = reform(ubvri_absmag[1,*]) ; Vega, h=0.7
    result.M_V = reform(ubvri_absmag[2,*]) ; Vega, h=0.7
    result.M_R = reform(ubvri_absmag[3,*]) ; Vega, h=0.7
    result.M_I = reform(ubvri_absmag[4,*]) ; Vega, h=0.7
    
    result.M_U_ivar = reform(ubvri_amivar[0,*]) ; Vega, h=0.7
    result.M_B_ivar = reform(ubvri_amivar[1,*]) ; Vega, h=0.7
    result.M_V_ivar = reform(ubvri_amivar[2,*]) ; Vega, h=0.7
    result.M_R_ivar = reform(ubvri_amivar[3,*]) ; Vega, h=0.7
    result.M_I_ivar = reform(ubvri_amivar[4,*]) ; Vega, h=0.7
    
    result.M_SDSS_U = reform(ugriz_absmag[0,*]) ; AB, h=0.7
    result.M_SDSS_G = reform(ugriz_absmag[1,*]) ; AB, h=0.7
    result.M_SDSS_R = reform(ugriz_absmag[2,*]) ; AB, h=0.7
    result.M_SDSS_I = reform(ugriz_absmag[3,*]) ; AB, h=0.7
    result.M_SDSS_Z = reform(ugriz_absmag[4,*]) ; AB, h=0.7

    result.M_SDSS_U_ivar = reform(ugriz_amivar[0,*]) ; AB, h=0.7
    result.M_SDSS_G_ivar = reform(ugriz_amivar[1,*]) ; AB, h=0.7
    result.M_SDSS_R_ivar = reform(ugriz_amivar[2,*]) ; AB, h=0.7
    result.M_SDSS_I_ivar = reform(ugriz_amivar[3,*]) ; AB, h=0.7
    result.M_SDSS_Z_ivar = reform(ugriz_amivar[4,*]) ; AB, h=0.7

; compute Vmax    

    splog, 'Computing Vmax'
    lf_calc_vmax, result.kcorr_mobs_ab[2], result.m_sdss_r+5.0*alog10(1.0/h100), $ ; h=0.7-->h=1
      result.kcorr_coeffs, survey_filter, dr4_area, sdss_mlimit[0], sdss_mlimit[1], sample_zmin, $
      sample_zmax, actual_z=sdss.redshift, vmax=vmax, zmin=zmin, zmax=zmax;, vname='default.nolines'
    
    result.zmin = zmin
    result.zmax = zmax
    result.vmax = vmax/h100^3.0 ; h=1-->h=0.7 [Mpc^3]
    result.vol  = (dr4_area/3.0)*(lf_comvol(sdss.redshift)-(lf_comvol(sample_zmin))[0])/h100^3.0 ; h=1-->h=0.7 [Mpc^3]

stop    
    
return
end
