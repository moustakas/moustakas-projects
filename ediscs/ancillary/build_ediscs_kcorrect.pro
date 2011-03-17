;
; NAME:
;   BUILD_EDISCS_KCORRECT
;
; PURPOSE:
;   Compute K-corrections and rest-frame luminosities (and stellar
;   masses) for EDisCS using K-correct.
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;   nowrite - do not write out the data structure
;
; OUTPUTS:
;   allediscsphot - photometry data structure
;   allresult   - stellar masses and K-corrections
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2008 Mar 25, NYU - stolen from EDISCS_KCORRECT_ALL
;   jm09feb22nyu - use IM_KCORRECT()
;   jm10may03ucsd - renamed BUILD_EDISCS_KCORRECT and rewritten to
;     match BUILD_AGES_KCORRECT
;-

function init_ediscs_kcorrect, filterlist=filterlist
; output data structure
    nfilt = n_elements(filterlist)
    result = {$
      ediscs_id:                          0,$
      galaxy:                            '',$
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

;     galex_absmag_00:         fltarr(2)-999.0,$
;     galex_absmag_ivar_00:    fltarr(2)-999.0,$
;     galex_kcorrect_00:       fltarr(2)-999.0,$
;     galex_absmag_01:         fltarr(2)-999.0,$
;     galex_absmag_ivar_01:    fltarr(2)-999.0,$
;     galex_kcorrect_01:       fltarr(2)-999.0,$

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

pro build_ediscs_kcorrect, result, nowrite=nowrite, clobber=clobber

    mycatpath = ediscs_path(/mycatalogs)
    kcorr_version = ediscs_version(/kcorrect)
    kcorrfile = mycatpath+'ediscs_kcorrect_'+kcorr_version+'.fits'
    if file_test(kcorrfile+'.gz',/regular) and $
      (keyword_set(clobber) eq 0) and $
      (keyword_set(nowrite) eq 0) then begin
       splog, 'FITS file '+kcorrfile+'.gz exists; use /CLOBBER'
       return
    endif

; read the EDISCS photometry and spectroscopic redshifts; only compute
; k-corrections for objects with photometry in at least three
; bandpasses (24/1739 objects have *no* photometry)
    spec1d = read_ediscs(/spec1d)
    phot = read_ediscs(/phot)
    ediscs_to_maggies, phot, maggies, ivarmaggies, $
      filterlist=filterlist
    index = where((strmatch(spec1d.galaxy,'*NONAME*',/fold) eq 0),ngal)
    splog, 'Computing K-corrections for '+strtrim(ngal,2)+' galaxies'

    in_redshift = spec1d[index].z
    in_maggies = maggies[*,index]
    in_ivarmaggies = ivarmaggies[*,index]

; initialize the output data structure
    result = replicate(init_ediscs_kcorrect(filterlist=filterlist),n_elements(phot))
    result.ediscs_id          = spec1d.ediscs_id
    result.galaxy             = spec1d.galaxy
    result.z                  = spec1d.z
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

;   splog, 'Computing NUV,FUV K-corrections'
;   galex_kcorrect_00 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
;     filterlist,galex_filterlist(),band_shift=0.0,absmag=galex_absmag_00,$
;     ivarabsmag=galex_absmag_ivar_00,/silent,vname=vname) ; AB, band_shift=0.0
;   galex_kcorrect_01 = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
;     filterlist,galex_filterlist(),band_shift=0.0,absmag=galex_absmag_01,$
;     ivarabsmag=galex_absmag_ivar_01,/silent,vname=vname) ; AB, band_shift=0.1

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

;   result[index].galex_absmag_00         = galex_absmag_00
;   result[index].galex_absmag_ivar_00    = galex_absmag_ivar_00
;   result[index].galex_kcorrect_00       = galex_kcorrect_00
;   result[index].galex_absmag_01         = galex_absmag_01
;   result[index].galex_absmag_ivar_01    = galex_absmag_ivar_01
;   result[index].galex_kcorrect_01       = galex_kcorrect_01

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
