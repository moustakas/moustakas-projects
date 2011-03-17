pro ages_isedfit_ugriz, models=models, isedfit=isedfit, qaplot=qaplot, $
  doplots=doplots, clobber=clobber
; jm09may20nyu - fit the AGES galaxies with just the SDSS ugriz
;   photometry  

    analysis_path = ages_path(/analysis)
    iopath = ages_path(/isedfit)
    paramfile = iopath+'ages_isedfit_ugriz.par'

    if keyword_set(models) then isedfit_models, $
      paramfile, iopath=iopath, /noigm

    ages = mrdfits(analysis_path+'ages_merged_catalogs.fits.gz',1,/silent)
    these = where((ages.qshort eq 0) and (ages.phot_bw gt 0.0) and $
      (ages.phot_r gt 0.0) and (ages.phot_i gt 0.0) and $
      (ages.phot_z gt 0.0) and (ages.phot_k gt 0.0) and $
      (ages.phot_flamj gt 0.0) and (ages.phot_flamk gt 0.0) and $
      (ages.z ge 0.01) and (ages.z le 1.0) and (ages.sdss_match eq 1),ngal)

    if keyword_set(isedfit) then begin
       ages_to_maggies, ages, maggies, ivarmaggies, /sdss
       isedfit, paramfile, maggies, ivarmaggies, ages.z, result, $
         iopath=iopath, outprefix=outprefix, nminphot=nminphot, $
         clobber=clobber, debug=0, index=these
    endif

; build a QAplot
    if keyword_set(qaplot) then begin
       isedfit_qaplot, paramfile, iopath=iopath, outprefix=outprefix, $
         galaxy=ages[these].ndwfs_galaxy, clobber=clobber, index=these
    endif

    
return
end
