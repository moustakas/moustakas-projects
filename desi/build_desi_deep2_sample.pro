pro build_desi_deep2_sample
; jm13dec18siena - build a sample of DEEP2 galaxies for DESI

    outpath = getenv('IM_PROJECTS_DIR')+'/desi/deep2/'

    zcat = read_deep2_zcat(photo=phot)
    deep2_to_maggies, phot, maggies, ivarmaggies, filterlist=filt
    nfilt = n_elements(filt)
    zcat = struct_addtags(zcat,replicate({maggies: fltarr(nfilt), $
      ivarmaggies: fltarr(nfilt), source: ''},n_elements(zcat)))
    zcat.maggies = maggies
    zcat.ivarmaggies = ivarmaggies
    zcat.source = phot.source

    ispec = read_deep2(/ispec)

    these = where(zcat.zbest gt 0.7 and zcat.zbest lt 1.6 and $
      total(maggies gt 0,1) ge 3 and $
      ispec.oii_3727_1[0]/ispec.oii_3727_1[1] gt 1.0 and $
      ispec.oii_3727_2[0]/ispec.oii_3727_2[1] gt 1.0,ngal)

    im_mwrfits, zcat[these], outpath+'deep2_zcat.fits', /clobber
    im_mwrfits, ispec[these], outpath+'deep2_ispec.fits', /clobber
    
return
end
