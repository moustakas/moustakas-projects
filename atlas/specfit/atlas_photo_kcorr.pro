;+
; NAME:
;       ATLAS_PHOTO_KCORR
;
; PURPOSE:
;       Compute k-corrections and stellar masses for the ATLAS using
;       Blanton's k-correct. 
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;      From Blanton in response to "How can I overplot the observed
;      fluxes on the best-fitting template?" 
; 
;      "I claim in the documentation that the units are erg/S/cm2/A  
; 
;         fluxfactor= 3.631D-29*(3D18)/(central)^2 ; AB nMgy to f_lambda
;         flambda= nmgys*fluxfactor
; 
;      Takes you from nanomaggies to f_lambda. 
; 
;      Then remember that the vmatrix is in the restframe, so you have
;      to multiply the wavelengths by (1+z) and divide the spectrum by
;      (1+z) to compare to the observed fluxes.
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2006 Apr 30, U of A - written, based on
;         AGES_PHOTO_KCORR 
;-

pro atlas_photo_kcorr, atlas, atlasinfo, mass, write=write

; define the output data path and the cosmology
    
    analysis_path = atlas_path(/analysis)

    red, omega0=0.3, omegalambda=0.7, h100=0.70
    omega0 = redomega0() & omegal = redomegal() & h100 = redh100()

    imf_chabrier_to_salpeter = +0.25

; Vega --> AB conversions

    U_vega2ab  = (k_vega2ab(filterlist='bessell_U.par',/kurucz))[0]
    B_vega2ab  = (k_vega2ab(filterlist='bessell_B.par',/kurucz))[0]
    V_vega2ab  = (k_vega2ab(filterlist='bessell_V.par',/kurucz))[0]
    R_vega2ab  = (k_vega2ab(filterlist='bessell_R.par',/kurucz))[0]
    I_vega2ab  = (k_vega2ab(filterlist='bessell_I.par',/kurucz))[0]
    H_vega2ab  = (k_vega2ab(filterlist='twomass_H.par',/kurucz))[0]
    K_vega2ab  = (k_vega2ab(filterlist='dwfs_onis_K.par',/kurucz))[0]

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
       
    filt = im_filterspecs(filterlist=filterlist)

; read the spectral atlas

    if (n_elements(atlas) eq 0L) then atlas = read_integrated()
    if (n_elements(atlasinfo) eq 0L) then atlasinfo = atlas_read_info()
    ngalaxy = n_elements(atlas)
    natlasinfo = n_elements(atlasinfo)

; initialize the output data structure

    result = {$
      ages_id:                 0L, $
      galaxy:                 ' ', $
      z:                   -999.0, $

      kcorr_Bw_U:          -999.0, $ ; derived k-corrections
;     kcorr_Bw_B:          -999.0, $
      kcorr_Bw_u_sdss:     -999.0, $

      kcorr_R_B:           -999.0, $
      kcorr_R_V:           -999.0, $
      kcorr_R_g_sdss:      -999.0, $

      kcorr_I_R:           -999.0, $
      kcorr_I_r_sdss:      -999.0, $
;     kcorr_I_I:           -999.0, $
;     kcorr_I_i_sdss:      -999.0, $

;     kcorr_g_sdss_U:      -999.0, $ 
;     kcorr_g_sdss_u_sdss: -999.0, $

;     kcorr_z_sdss_I:      -999.0, $
;     kcorr_z_sdss_i_sdss: -999.0, $

;     kcorr_K_K:           -999.0, $
;     kcorr_K_H:           -999.0  $

      kcorr_chi2:          -999.0, $
      kcorr_mass:          -999.0, $
      kcorr_mets:          -999.0, $
      kcorr_intsfh:        -999.0  $
      }
    result = replicate(result,ngalaxy)

    redshift = agescats.z
    
    result.ages_id = agescats.ages_id
    result.galaxy = agescats.galaxy
    result.z = redshift

; setup the k-correction photometry vectors: BwRIK

    bwrik_mags_vega = transpose([ $
      [agescats.phot_bw],$
      [agescats.phot_r], $
      [agescats.phot_i], $
      [agescats.phot_k]  $
      ])
    bwrik_mags_vega_err = transpose([ $
      [agescats.phot_bw_err],$
      [agescats.phot_r_err], $
      [agescats.phot_i_err], $
      [agescats.phot_k_err]  $
      ])

    bwrik_mags_ab = bwrik_mags_vega
    bwrik_mags_ab_err = bwrik_mags_vega_err
    bwrik_vega2ab = [bw_vega2ab,r_vega2ab,i_vega2ab,k_vega2ab]
    
    nphot = (size(bwrik_mags_vega,/dimension))[0]
    for iphot = 0L, nphot-1L do begin

       good = where(bwrik_mags_ab[iphot,*] gt -900.0,comp=flag,ncomp=nflag)
       bwrik_mags_ab[iphot,good] = bwrik_mags_vega[iphot,good] + bwrik_vega2ab[iphot]
       bwrik_mags_ab_err[iphot,good] = bwrik_mags_vega_err[iphot,good]

       if (nflag ne 0L) then begin
          bwrik_mags_ab[iphot,flag] = 0.0
          bwrik_mags_ab_err[iphot,flag] = 1E32
       endif
       
    endfor

; subset of the available photometry    
    
    bwri_mags_vega = bwrik_mags_vega[0:2,*]
    bwri_mags_vega_err = bwrik_mags_vega_err[0:2,*]
    bwri_mags_ab = bwrik_mags_ab[0:2,*]
    bwri_mags_ab_err = bwrik_mags_ab_err[0:2,*]
    
;;; SDSS ugriz
;;
;;    sdss_mags_ab = transpose([ $
;;      [agescats.phot_sdss_u],$
;;      [agescats.phot_sdss_g],$
;;      [agescats.phot_sdss_r],$
;;      [agescats.phot_sdss_i],$
;;      [agescats.phot_sdss_z] $
;;      ])
;;    sdss_mags_ab_err = transpose([ $
;;      [agescats.phot_sdss_u_err],$
;;      [agescats.phot_sdss_g_err],$
;;      [agescats.phot_sdss_r_err],$
;;      [agescats.phot_sdss_i_err],$
;;      [agescats.phot_sdss_z_err] $
;;      ])
;;
;;    nphot = (size(sdss_mags_ab,/dimension))[0]
;;    for iphot = 0L, nphot-1L do begin
;;
;;       good = where(sdss_mags_ab[iphot,*] gt -900.0,comp=flag,ncomp=nflag)
;;       sdss_mags_ab_err[iphot,good] = sdss_mags_ab_err[iphot,good]
;;
;;       if (nflag ne 0L) then begin
;;          sdss_mags_ab[iphot,flag] = 0.0
;;          sdss_mags_ab_err[iphot,flag] = 1E32
;;       endif
;;    endfor
    
; concatenate all the photometry and convert to maggies; for now,
; daniel recommends against mixing the SDSS and NDWFS photometry; so
; listen to him!  and he also advises against using the K-band, since
; it only exists for about half the sample and is of considerably
; lower quality than the optical photometry (jm05aug04uofa)

    if keyword_set(usekband) then begin
       mags_ab = bwrik_mags_ab
       mags_err = bwrik_mags_ab_err
    endif else begin
       mags_ab = bwri_mags_ab
       mags_err = bwri_mags_ab_err
    endelse
;   mags_ab = [bwri_mags_ab,sdss_mags_ab]
;   mags_err = [bwri_mags_ab_err,sdss_mags_ab_err]

    maggies = mags_ab*0.0
    maggies_err = mags_ab*0.0
    maggies_ivar = mags_ab*0.0
    
    good = where(mags_err ne 1E32,ngood,comp=flag,ncomp=nflag)

    maggies[good] = 10.0^(-0.4*mags_ab[good])                       ; flux density
    maggies_err[good] = 0.4*alog(10.0)*maggies[good]*mags_err[good] ; error
    maggies_ivar[good] = 1.0/maggies_err[good]^2                    ; inverse variance

    if (nflag ne 0L) then maggies_ivar[flag] = 0.0 ; zero weight

; compute the k-corrections      

    splog, 'Computing k-corrections.'
    t0 = systime(1)
    kcorrect, maggies, maggies_ivar, redshift, kcorrect, $
      band_shift=band_shift, filterlist=filterlist, $
      filterpath=filterpath, rmatrix=rmatrix, zvals=zvals, $
      lambda=lambda, vmatrix=vmatrix, coeffs=coeffs, chi2=chi2, $
      maxiter=maxiter, omega0=omega0, omegal0=omegal, mass=mass, $
      intsfh=intsfh, mets=mets, b300=b300, b1000=b1000, mtol=mtol
    splog, 'Total time = '+string((systime(1)-t0)/60.0,format='(G0.0)')+' minutes.'

    result.kcorr_chi2 = chi2

; store the stellar mass; convert the IMF and Hubble constant    

    good = where(mass gt 0.0,ngood)
    if (ngood ne 0L) then begin
       result[good].kcorr_mass = alog10(mass[good]) - alog10(h100) + imf_chabrier_to_salpeter
       result[good].kcorr_mets = mets[good]
       result[good].kcorr_intsfh = alog10(intsfh[good])
    endif
    
; adopted k-corrections: g-->Uu; R-->BVg; I-->Rr; z-->Ii; K-->H 
    
; UBVRIH at z = 0

    k_reconstruct_maggies, coeffs, redshift*0.0, maggies_U_rest, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='bessell_U.par', /silent
    mag_U_rest = -2.5*alog10(reform(maggies_U_rest)) - U_vega2ab ; Vega!
    
    k_reconstruct_maggies, coeffs, redshift*0.0, maggies_B_rest, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='bessell_B.par', /silent
    mag_B_rest = -2.5*alog10(reform(maggies_B_rest)) - B_vega2ab ; Vega

    k_reconstruct_maggies, coeffs, redshift*0.0, maggies_V_rest, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='bessell_V.par', /silent
    mag_V_rest = -2.5*alog10(reform(maggies_V_rest)) - V_vega2ab ; Vega

    k_reconstruct_maggies, coeffs, redshift*0.0, maggies_R_rest, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='bessell_R.par', /silent
    mag_R_rest = -2.5*alog10(reform(maggies_R_rest)) - R_vega2ab ; Vega

    k_reconstruct_maggies, coeffs, redshift*0.0, maggies_I_rest, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='bessell_I.par', /silent
    mag_I_rest = -2.5*alog10(reform(maggies_I_rest)) - I_vega2ab ; Vega

    k_reconstruct_maggies, coeffs, redshift*0.0, maggies_H_rest, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='twomass_H.par', /silent
    mag_H_rest = -2.5*alog10(reform(maggies_H_rest)) - H_vega2ab ; Vega

; SDSS ugriz at z = 0
    
    k_reconstruct_maggies, coeffs, redshift*0.0, maggies_sdss_u_rest, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='sdss_u0.par', /silent
    mag_sdss_u_rest = -2.5*alog10(reform(maggies_sdss_u_rest)) ; AB

    k_reconstruct_maggies, coeffs, redshift*0.0, maggies_sdss_g_rest, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='sdss_g0.par', /silent
    mag_sdss_g_rest = -2.5*alog10(reform(maggies_sdss_g_rest)) ; AB

    k_reconstruct_maggies, coeffs, redshift*0.0, maggies_sdss_r_rest, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='sdss_r0.par', /silent
    mag_sdss_r_rest = -2.5*alog10(reform(maggies_sdss_r_rest)) ; AB

    k_reconstruct_maggies, coeffs, redshift*0.0, maggies_sdss_i_rest, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='sdss_i0.par', /silent
    mag_sdss_i_rest = -2.5*alog10(reform(maggies_sdss_i_rest)) ; AB

    k_reconstruct_maggies, coeffs, redshift*0.0, maggies_sdss_z_rest, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='sdss_z0.par', /silent
    mag_sdss_z_rest = -2.5*alog10(reform(maggies_sdss_z_rest)) ; AB

; BwRIK at z = zobj (just do BwRI for now - jm05aug04uofa)
    
    k_reconstruct_maggies, coeffs, redshift, maggies_Bw_obs, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='dwfs_Bw.par', /silent
    mag_Bw_obs = -2.5*alog10(reform(maggies_Bw_obs)) - Bw_vega2ab ; Vega

    k_reconstruct_maggies, coeffs, redshift, maggies_R_obs, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='bessell_R.par', /silent
    mag_R_obs = -2.5*alog10(reform(maggies_R_obs)) - R_vega2ab ; Vega

    k_reconstruct_maggies, coeffs, redshift, maggies_I_obs, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='bessell_I.par', /silent
    mag_I_obs = -2.5*alog10(reform(maggies_I_obs)) - I_vega2ab ; Vega
    
    k_reconstruct_maggies, coeffs, redshift, maggies_K_obs, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='dwfs_onis_K.par', /silent
    mag_K_obs = -2.5*alog10(reform(maggies_K_obs)) - K_vega2ab ; Vega

; SDSS ugriz at z = zobj
    
    k_reconstruct_maggies, coeffs, redshift, maggies_sdss_u_obs, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='sdss_u0.par', /silent
    mag_sdss_u_obs = -2.5*alog10(reform(maggies_sdss_u_obs)) ; AB

    k_reconstruct_maggies, coeffs, redshift, maggies_sdss_g_obs, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='sdss_g0.par', /silent
    mag_sdss_g_obs = -2.5*alog10(reform(maggies_sdss_g_obs)) ; AB

    k_reconstruct_maggies, coeffs, redshift, maggies_sdss_r_obs, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='sdss_r0.par', /silent
    mag_sdss_r_obs = -2.5*alog10(reform(maggies_sdss_r_obs)) ; AB

    k_reconstruct_maggies, coeffs, redshift, maggies_sdss_i_obs, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='sdss_i0.par', /silent
    mag_sdss_i_obs = -2.5*alog10(reform(maggies_sdss_i_obs)) ; AB

    k_reconstruct_maggies, coeffs, redshift, maggies_sdss_z_obs, $
      band_shift=band_shift, zvals=zvals, rmatrix=rmatrix, $
      vmatrix=vmatrix, lambda=lambda, filterpath=filterpath, $
      filterlist='sdss_z0.par', /silent
    mag_sdss_z_obs = -2.5*alog10(reform(maggies_sdss_z_obs)) ; AB

; the following bit of code compares my output from kcorrect to
; daniel's 
    
;   kcorr = mrdfits(analysis_path+'catalog.kcorr.v3.fits.gz',1,/silent)
;   kid = lindgen(n_elements(kcorr))
;   match, result.ages_id, kid, indx1, indx2
;   plot, mag_bw_obs[indx1], mag_bw_obs[indx1]-kcorr[indx2].kcorr_bw_model, ps=4
;   plot, mag_r_obs[indx1], mag_r_obs[indx1]-kcorr[indx2].kcorr_r_model, ps=4
;   plot, mag_i_obs[indx1], mag_i_obs[indx1]-kcorr[indx2].kcorr_i_model, ps=4
;
;   junk1 = im_stats(mag_bw_obs[indx1]-kcorr[indx2].kcorr_bw_model,/verbose)
;   junk2 = im_stats(mag_r_obs[indx1]-kcorr[indx2].kcorr_r_model,/verbose)
;   junk3 = im_stats(mag_i_obs[indx1]-kcorr[indx2].kcorr_i_model,/verbose)
    
; store the k-corrections of interest

    result.kcorr_Bw_U = (mag_U_rest-mag_Bw_obs) ; Bw-->U
;   result.kcorr_Bw_B = (mag_B_rest-mag_Bw_obs) ; Bw-->B
    result.kcorr_Bw_u_sdss = (mag_sdss_u_rest-mag_Bw_obs) ; Bw-->u_sdss

    result.kcorr_R_B = (mag_B_rest-mag_R_obs) ; R-->B
    result.kcorr_R_V = (mag_V_rest-mag_R_obs) ; R-->V
    result.kcorr_R_g_sdss = (mag_sdss_g_rest-mag_R_obs) ; R-->g_sdss

    result.kcorr_I_R = (mag_R_rest-mag_I_obs)           ; I-->R
    result.kcorr_I_r_sdss = (mag_sdss_r_rest-mag_I_obs) ; I-->r_sdss
;   result.kcorr_I_I = (mag_I_rest-mag_I_obs)           ; I-->I
;   result.kcorr_I_i_sdss = (mag_sdss_i_rest-mag_I_obs) ; I-->i

;   result.kcorr_g_sdss_U = (mag_U_rest-mag_sdss_g_obs)           ; g-->U
;   result.kcorr_g_sdss_u_sdss = (mag_sdss_u_rest-mag_sdss_g_obs) ; g-->u

;   result.kcorr_z_sdss_I = (mag_I_rest-mag_sdss_z_obs)           ; z-->I
;   result.kcorr_z_sdss_i_sdss = (mag_sdss_i_rest-mag_sdss_z_obs) ; z-->i

;   result.kcorr_K_H = (mag_H_rest-mag_K_obs) ; K-->H
;   result.kcorr_K_K = (mag_K_rest-mag_K_obs) ; K-->K

; write out    
    
    if keyword_set(write) then begin

       if keyword_set(usekband) then suffix = '_withk' else suffix = ''

       splog, 'Writing '+analysis_path+'atlas_kcorr'+suffix+'.fits.gz'
       mwrfits, result, analysis_path+'atlas_kcorr'+suffix+'.fits', /create
       spawn, ['gzip -f '+analysis_path+'atlas_kcorr'+suffix+'.fits'], /sh
       
    endif

return
end
    
