;+
; NAME:
;   build_vagc_kcorrect
; PURPOSE:
;   builds kcorrect dir for vagc
; CALLING SEQUENCE:
;   build_kcorrect
; INPUTS:
; OPTIONAL INPUTS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; DATA DEPENDENCIES:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   03-Nov-2002  Written by Mike Blanton, NYU
;   10-Mar-2010 J. Moustakas, UCSD - modified version of the
;     VAGC's BUILD_KCORRECT to have more control over the
;     rest-frame bandpasses; many of the bells and whistles have also
;     been removed
;-

pro build_vagc_kcorrect, collision_type=collision_type, $
  flux_type=flux_type, vname=vname, outfile=outfile

; echo "build_vagc_kcorrect" | idl > & log &     
    
    common com_vagc_kcorrect, im, twomass
    
    if (n_elements(flux_type) eq 0) then flux_type = 'model' ; = 'petro'
    if (n_elements(collision_type) eq 0) then collision_type = 'nearest'
    if (n_elements(vname) eq 0) then vname = 'default.nolines'

; read in object_sdss_imaging, twomass and collisions
    object_sdss_imaging_file = vagc_name('object_sdss_imaging')+'.gz'
    if (n_elements(im) eq 0L) then begin
       im = hogg_mrdfits(object_sdss_imaging_file,1,$
         nrowchunk=75000,columns=['object_position',$
         'ra','dec',flux_type+'flux',flux_type+'flux_ivar'])
    endif
       
    iobject_position = tag_indx(im[0],'object_position')
    if (iobject_position eq -1) then $
      object_position = lindgen(n_elements(im)) $
    else object_position = im.(iobject_position)

    if (n_elements(twomass) eq 0L) then begin
       if (file_test(vagc_name('object_twomass')+'.gz')) then begin
          twomass = hogg_mrdfits(vagc_name('object_twomass')+'.gz',1,$
            nrowchunk=75000, columns=['ra', $
            'decl','twomass_tag','k_m_ext','k_msig_ext', $
            'h_m_ext','h_msig_ext','j_m_ext','j_msig_ext'])
          twomass = twomass[object_position]
       endif
    endif

; figure out which collisions corrections to use, if any
    words = strsplit(collision_type,'_',/extr)
    if (words[0] eq 'lss') then begin
       collisions = mrdfits(vagc_name('lss_index',sample=words[1]),1)
       collisions = collisions[object_position]
    endif else begin 
       collisions = mrdfits(vagc_name('collisions',collision_type=collision_type)+'.gz',1)
       collisions = collisions[object_position]
    endelse

; compute K-corrections
;   test_im = im[0:1000]
;   test_twomass = twomass[0:1000]
;   test_collisions = collisions[0:1000]
    
    ngal = n_elements(im)
    filterlist = [sdss_filterlist(),twomass_filterlist()]
    nfilt = n_elements(filterlist)
    
    kcorr = {$
      ra:                              0.0D,$
      dec:                             0.0D,$
      z:                                0.0,$
      maggies:                fltarr(nfilt),$
      ivarmaggies:            fltarr(nfilt),$
      bestmaggies:            fltarr(nfilt),$
      mass:                          -999.0,$
      coeffs:                     fltarr(5),$
      chi2:                          -999.0,$
      cflux:                fltarr(5)-999.0,$
;     uvflux:               fltarr(2)-999.0,$
;     galex_absmag:         fltarr(2)-999.0,$
;     galex_absmag_ivar:    fltarr(2)-999.0,$
;     galex_kcorrect:       fltarr(2)-999.0,$
      ugriz_absmag:         fltarr(5)-999.0,$
      ugriz_absmag_ivar:    fltarr(5)-999.0,$
      ugriz_kcorrect:       fltarr(5)-999.0,$
      ubvrijhk_absmag:      fltarr(8)-999.0,$
      ubvrijhk_absmag_ivar: fltarr(8)-999.0,$
      ubvrijhk_kcorrect:    fltarr(8)-999.0}
    kcorr = replicate(kcorr,ngal)

    sdss_to_maggies, smaggies, sivarmaggies, calib=im, flux=flux_type
    twomass_to_maggies, twomass, tmaggies, tivarmaggies
    
    kcorr.ra = collisions.ra
    kcorr.dec = collisions.dec
    kcorr.z = collisions.z
    kcorr.maggies = [smaggies,tmaggies]
    kcorr.ivarmaggies = [sivarmaggies,tivarmaggies]

    good = where((total(kcorr.ivarmaggies gt 0,1) ge 3) and $
      (kcorr.z gt 0.0),ngood)
    splog, 'Number of galaxies with good redshifts', ngood
    
    splog, 'Computing ugriz K-corrections'
    t0 = systime(1)
    ugriz_kcorrect = im_kcorrect(kcorr[good].z,kcorr[good].maggies,kcorr[good].ivarmaggies,$
      filterlist,sdss_filterlist(),band_shift=0.1,chi2=chi2,mass=mass,$
      coeffs=coeffs,bestmaggies=bestmaggies,absmag=ugriz_absmag,$ ; AB
      ivarabsmag=ugriz_absmag_ivar,clineflux=cflux,/silent,vname=vname);,uvflux=uvflux)
    splog, 'Total time = ', (systime(1)-t0)/60.0, ' minutes'
    
;   galex_kcorrect = im_kcorrect(in_redshift,in_maggies,in_ivarmaggies,$
;     filterlist,galex_filterlist(),band_shift=0.0,absmag=galex_absmag,$
;     ivarabsmag=galex_absmag_ivar,/silent,vname=vname) ; AB

    splog, 'Computing UBVRIJHKs K-corrections'
    t0 = systime(1)
    kcorrect = im_kcorrect(kcorr[good].z,kcorr[good].maggies,kcorr[good].ivarmaggies,$ 
      filterlist,[bessell_filterlist(),twomass_filterlist()],$
      band_shift=0.0,absmag=absmag,ivarabsmag=absmag_ivar,$
      /vega,/silent,vname=vname) ; Vega           
    splog, 'Total time = ', (systime(1)-t0)/60.0, ' minutes'
    
    kcorr[good].bestmaggies = bestmaggies
    kcorr[good].mass = alog10(mass) ; h=0.7, Chabrier
    kcorr[good].coeffs = coeffs
    kcorr[good].chi2 = chi2
    kcorr[good].cflux = cflux
;   kcorr[good].uvflux = uvflux
    kcorr[good].ugriz_absmag = ugriz_absmag
    kcorr[good].ugriz_absmag_ivar = ugriz_absmag_ivar
    kcorr[good].ugriz_kcorrect = ugriz_kcorrect
;   kcorr[good].galex_absmag = galex_absmag
;   kcorr[good].galex_absmag_ivar = galex_absmag_ivar
;   kcorr[good].galex_kcorrect = galex_kcorrect
    kcorr[good].ubvrijhk_absmag = absmag
    kcorr[good].ubvrijhk_absmag_ivar = absmag_ivar
    kcorr[good].ubvrijhk_kcorrect = kcorrect

; write out
    spawn, 'mkdir -p '+getenv('VAGC_REDUX')+'/kcorrect'
    if (n_elements(outfile) eq 0) then $
      outfile = getenv('VAGC_REDUX')+'/kcorrect/kcorrect.'+collision_type+'.'+flux_type+'.fits'
;     outfile = vagc_name('kcorrect',collision_type=collision_type,$
;     flux_type=flux_type,band_shift=band_shift)
    
;    hdr=['']
;    sxaddpar, hdr, 'VAGCTIME',systime(),'Time of creation'
;    sxaddpar, hdr, 'VAGCVERS',vagc_version(),'Version of vagc used'
;    sxaddpar, hdr, 'KCORRVER',k_version(),'version of kcorrect'
;    sxaddpar, hdr, 'COLLTYPE',collision_type,'type of collision correction'
;    sxaddpar, hdr, 'OMEGA0',omega0,'matter density for absmag'
;    sxaddpar, hdr, 'OMEGAL0',omegal0,'cosmo const for absmag'
;    sxaddpar, hdr, 'FLUXTYPE',flux_type,'type of flux used'
;    sxaddpar, hdr, 'BANDSHFT',string(band_shift),'band blueshift used'
;    sxaddpar, hdr, 'PHOTO_RESOLVE',getenv('PHOTO_RESOLVE'),'type of resolve used'
;    sxaddpar, hdr, 'PHOTO_CALIB',getenv('PHOTO_CALIB'),'type of calib used'
;    if (n_elements(twomass) gt 0) then $
;      sxaddpar,hdr,'TWOMASS','yes','2mass data included'

    im_mwrfits, kcorr, outfile, /clobber

stop    
    
return    
end
