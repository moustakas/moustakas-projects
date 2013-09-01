;
; NAME:
;   BUILD_AGES_KCORRECT
;
; PURPOSE:
;   Compute K-corrections, rest-frame luminosities, and stellar 
;   masses for AGES using K-correct.
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   result - stellar masses and K-corrections
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2007 Jun 20, NYU - excised from WRITE_SDSS_ANCILLARY
;   jm08mar24nyu - major re-write (see WRITE_SDSS_MAIN)  
;   jm08jun18nyu - bug fix in computing masses; remove Bell & de
;     Jong mass estimate
;   jm08aug08nyu - another bit of a rewrite to now use
;     IM_KCORRECT() 
;   jm09mar03nyu - use AGES_TO_MAGGIES
;   jm10may01ucsd - renamed BUILD_AGES_KCORRECT and cleaned up 
;-

function init_ages_kcorrect, filterlist=filterlist
; output data structure
    nfilt = n_elements(filterlist)
    result = {$
      ages_id:                           0L,$
      z:                             -999.0,$
      filterlist:                filterlist,$
      maggies:                fltarr(nfilt),$
      ivarmaggies:            fltarr(nfilt),$
      bestmaggies:            fltarr(nfilt),$
      mass:                          -999.0,$
      coeffs:                     fltarr(5),$
      chi2:                          -999.0,$
      cflux:                fltarr(5)-999.0,$
      uvflux:               fltarr(2)-999.0,$

      galex_absmag_00:         fltarr(2)-999.0,$
      galex_absmag_ivar_00:    fltarr(2)-999.0,$
      galex_kcorrect_00:       fltarr(2)-999.0,$
      galex_absmag_01:         fltarr(2)-999.0,$
      galex_absmag_ivar_01:    fltarr(2)-999.0,$
      galex_kcorrect_01:       fltarr(2)-999.0,$

      ugriz_absmag_00:         fltarr(5)-999.0,$
      ugriz_absmag_ivar_00:    fltarr(5)-999.0,$
      ugriz_kcorrect_00:       fltarr(5)-999.0,$
      ugriz_absmag_01:         fltarr(5)-999.0,$
      ugriz_absmag_ivar_01:    fltarr(5)-999.0,$
      ugriz_kcorrect_01:       fltarr(5)-999.0,$

      ubvrijhk_absmag_00:      fltarr(8)-999.0,$
      ubvrijhk_absmag_ivar_00: fltarr(8)-999.0,$
      ubvrijhk_kcorrect_00:    fltarr(8)-999.0,$
      ubvrijhk_absmag_01:      fltarr(8)-999.0,$
      ubvrijhk_absmag_ivar_01: fltarr(8)-999.0,$
      ubvrijhk_kcorrect_01:    fltarr(8)-999.0}
return, result
end

pro build_ages_kcorrect, result, nowrite=nowrite, clobber=clobber

    mycatpath = ages_path(/mycatalogs)
    kcorr_version = ages_version(/kcorrect)
    kcorrfile = mycatpath+'ages_kcorrect_'+kcorr_version+'.fits'
    if file_test(kcorrfile+'.gz',/regular) and $
      (keyword_set(clobber) eq 0) and $
      (keyword_set(nowrite) eq 0) then begin
       splog, 'FITS file '+kcorrfile+'.gz exists; use /CLOBBER'
       return
    endif

; read the AGES photometry and then fit the subset of the sample with
; GSHORT>0, 0.001<z<1.0, and good BwRI photometry
    phot = read_ages(/phot)
    ages_to_maggies, phot, maggies, ivarmaggies, /itot, $
      filterlist=filterlist, use_aper='04', /totalmag
    
    bwindx = where(strmatch(filterlist,'*_Bw*'))
    rindx = where(strmatch(filterlist,'*_R*'))
    iindx = where(strmatch(filterlist,'*_I*'))
    index = where((phot.gshort gt 0) and (phot.z gt 0.001) and $ ; (phot.z lt 1.0) and $
      (ivarmaggies[rindx,*] gt 0.0) and (ivarmaggies[iindx,*] gt 0.0),ngal)
;   index = where((phot.gshort gt 0) and (phot.z gt 0.001) and $
;     (phot.z lt 1.0) and (ivarmaggies[bwindx,*] gt 0.0) and $
;     (ivarmaggies[rindx,*] gt 0.0) and (ivarmaggies[iindx,*] gt 0.0),ngal)
    splog, 'Computing K-corrections for '+strtrim(ngal,2)+' galaxies'

; final photometry; do not use ch2-4
    in_redshift = phot[index].z
    in_maggies = maggies[*,index]
    in_ivarmaggies = ivarmaggies[*,index]

    toss = where(strmatch(filterlist,'*_ch[2-4].par*',/fold),ntoss)
    if (ntoss ne 0) then in_ivarmaggies[toss,*] = 0.0
    
; initialize the output data structure
    result = init_ages_kcorrect(filterlist=filterlist)
    result = replicate(result,n_elements(phot))
    result.ages_id            = phot.ages_id
    result.z                  = phot.z
    result[index].maggies     = in_maggies
    result[index].ivarmaggies = in_ivarmaggies
    
    vname = 'default.nolines' ; note!

; compute k-corrections    
    splog, 'Computing ugriz K-corrections'
    ugriz_kcorrect_00 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,sdss_filterlist(),band_shift=0.0,chi2=chi2,mass=mass,$
      coeffs=coeffs,bestmaggies=bestmaggies,absmag=ugriz_absmag_00,$
      ivarabsmag=ugriz_absmag_ivar_00,clineflux=cflux,uvflux=uvflux,$
      /silent,vname=vname) ; AB, band_shift=0.0
    ugriz_kcorrect_01 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,sdss_filterlist(),band_shift=0.1,absmag=ugriz_absmag_01,$
      ivarabsmag=ugriz_absmag_ivar_01,/silent,vname=vname) ; AB, band_shift=0.1

    splog, 'Computing NUV,FUV K-corrections'
    galex_kcorrect_00 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,galex_filterlist(),band_shift=0.0,absmag=galex_absmag_00,$
      ivarabsmag=galex_absmag_ivar_00,/silent,vname=vname) ; AB, band_shift=0.0
    galex_kcorrect_01 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,galex_filterlist(),band_shift=0.0,absmag=galex_absmag_01,$
      ivarabsmag=galex_absmag_ivar_01,/silent,vname=vname) ; AB, band_shift=0.1

    splog, 'Computing UBVRIJHKs K-corrections'
    kcorrect_00 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,[bessell_filterlist(),twomass_filterlist()],band_shift=0.0,$
      absmag=absmag_00,ivarabsmag=absmag_ivar_00,/silent,vname=vname) ; AB, band_shift=0.0
    kcorrect_01 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
      filterlist,[bessell_filterlist(),twomass_filterlist()],band_shift=0.0,$
      absmag=absmag_01,ivarabsmag=absmag_ivar_01,/silent,vname=vname) ; AB, band_shift=0.1

; pack the structure and write out   
    result[index].bestmaggies = bestmaggies
    result[index].mass        = alog10(mass)   ; Chabrier IMF, h=0.7
    result[index].coeffs      = coeffs
    result[index].chi2        = chi2
    result[index].uvflux      = uvflux ; flux density at 1500, 2800 A
    result[index].cflux       = cflux  ; line-continuum (see IM_KCORRECT)

    result[index].galex_absmag_00         = galex_absmag_00
    result[index].galex_absmag_ivar_00    = galex_absmag_ivar_00
    result[index].galex_kcorrect_00       = galex_kcorrect_00
    result[index].galex_absmag_01         = galex_absmag_01
    result[index].galex_absmag_ivar_01    = galex_absmag_ivar_01
    result[index].galex_kcorrect_01       = galex_kcorrect_01

    result[index].ugriz_absmag_00         = ugriz_absmag_00
    result[index].ugriz_absmag_ivar_00    = ugriz_absmag_ivar_00
    result[index].ugriz_kcorrect_00       = ugriz_kcorrect_00
    result[index].ugriz_absmag_01         = ugriz_absmag_01
    result[index].ugriz_absmag_ivar_01    = ugriz_absmag_ivar_01
    result[index].ugriz_kcorrect_01       = ugriz_kcorrect_01

    result[index].ubvrijhk_absmag_00      = absmag_00
    result[index].ubvrijhk_absmag_ivar_00 = absmag_ivar_00
    result[index].ubvrijhk_kcorrect_00    = kcorrect_00
    result[index].ubvrijhk_absmag_01      = absmag_01
    result[index].ubvrijhk_absmag_ivar_01 = absmag_ivar_01
    result[index].ubvrijhk_kcorrect_01    = kcorrect_01

    if (keyword_set(nowrite) eq 0) then begin
       im_mwrfits, result, kcorrfile, clobber=clobber
;      psfile = repstr(kcorrfile,'.fits','.ps')
;      kcorrect_qaplot, result[0:100], filterlist, psfile=psfile, /clobber
    endif
       
return
end
