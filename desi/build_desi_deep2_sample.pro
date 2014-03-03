pro build_desi_deep2_sample
; jm13dec18siena - build a sample of DEEP2 galaxies for DESI; this is
; a super-sample of objects we are actually probably interested in 

    outpath = getenv('IM_PROJECTS_DIR')+'/desi/deep2/'

    zcat = read_deep2_zcat(photo=phot)
    kcorr = read_deep2(/kcorr)
    ngal = n_elements(zcat)
    
    deep2_to_maggies, phot, maggies, ivarmaggies, filterlist=filt
    nfilt = n_elements(filt)

    zcat = struct_addtags(zcat,replicate({maggies: fltarr(nfilt), $
      ivarmaggies: fltarr(nfilt), source: ''},n_elements(zcat)))
    zcat.maggies = maggies
    zcat.ivarmaggies = ivarmaggies
    zcat.source = phot.source

    ispec = read_deep2(/ppxf)

;   oiiflux = fltarr(
;   both = where(ispec.oii_3727_1[1] gt 0 and ispec.oii_3727_2[1] gt 0,nboth)
;   if nboth ne 0L then begin
;   endif

    these = where($
      zcat.zbest gt 0.7 and zcat.zbest lt 1.6 and $
      total(maggies gt 0,1) ge 3 and $
      -2.5*alog10(maggies[1,*]) lt 24.0,ngal)
;     ispec.oii_3727_1_ew[0] gt 1.0 and $
;     ispec.oii_3727_2_ew[0] gt 1.0 and $
;     ispec.oii_3727_1[0]/ispec.oii_3727_1[1] gt 5.0 and $
;     ispec.oii_3727_2[0]/ispec.oii_3727_2[1] gt 5.0,ngal)
    splog, 'Sample ', ngal

; build a poor-mans desi sample by applying on [OII] flux simple cut
; >5e-17 cgs 
    im_mwrfits, zcat[these], outpath+'deep2_zcat.fits', /clobber
    im_mwrfits, kcorr[these], outpath+'deep2_kcorr.fits', /clobber
    im_mwrfits, ispec[these], outpath+'deep2_ispec.fits', /clobber
    
return
end
