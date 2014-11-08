pro decamlegacy_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber
; jm14sep07siena

    prefix = 'decals'
    isedfit_dir = getenv('DECALS_DIR')+'/isedfit/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    
; gather the photometry
;   cat = mrdfits(isedfit_dir+'phot-spec.fits.gz',1)

    filterlist = decamlegacy_filterlist()
    nfilt = n_elements(filterlist)

; --------------------------------------------------
; build the parameter files
    if keyword_set(write_paramfile) then begin
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, zminmax=[0.001,2.0], nzz=5, $
         spsmodels='fsps_v2.4_miles', imf='chab', redcurve='calz', igm=0, $
         sfhgrid=1, nmodel=5000L, age=[0.01,0.6], AV=[0.0,0.0], tau=[0.0,0.6], $
         Zmetal=[0.0008,0.019], oiiihb=[0.0,1.0], /nebular, /delayed, $
         clobber=clobber
    endif

; --------------------------------------------------
; build the Monte Carlo grids    
    if keyword_set(build_grids) then begin
       isedfit_montegrids, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, clobber=clobber
    endif

; --------------------------------------------------
; calculate the model photometry 
    if keyword_set(model_photometry) then begin
       isedfit_models, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber
    endif

; --------------------------------------------------
; do the fitting!
    if keyword_set(isedfit) then begin
       hff_to_maggies, cat, maggies, ivarmaggies, /nJy, filterlist=filt
       if keyword_set(photoz) eq 0 then z = cat.z ; cat.bpz
       if keyword_set(noirac) then begin
          isirac = where(strmatch(filt,'*irac*'))
          ivarmaggies[isirac,*] = 0
       endif
       isedfit, isedfit_paramfile, maggies, ivarmaggies, z, thissfhgrid=thissfhgrid, $
         isedfit_dir=isedfit_dir, outprefix=outprefix, isedfit_results=ised, $
         isedfit_post=isedpost, clobber=clobber, photoz=photoz, index=index
    endif 

; --------------------------------------------------
; compute K-corrections
    if keyword_set(kcorrect) then begin
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         absmag_filterlist=galex_filterlist(), band_shift=0.0, $
         clobber=clobber
    endif 

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, galaxy=cat.galaxy, yrange=[32,23], $
         xrange=[0.3,7.0]*1D4, nsigma=2.0, index=index, outprefix=outprefix
    endif
    
return
end
