pro a2744_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber, photoz=photoz
; jm12sep17siena

    if keyword_set(photoz) then prefix = 'a2744_photoz' else prefix = 'a2744'
    isedfit_dir = getenv('CLASH_PROJECTS')+'/hff/a2744/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    
; gather the photometry
    cat = read_a2744(photoz=photoz)

    filterlist = hff_filterlist(/useirac)
    nfilt = n_elements(filterlist)

; --------------------------------------------------
; build the parameter files
    if keyword_set(write_paramfile) then begin
       spsmodels = 'fsps_v2.4_miles'
       imf = 'chab'
       if keyword_set(photoz) then begin
          nmodel = 10000
          redcurve = 'charlot'
          write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
            prefix=prefix, filterlist=filterlist, zminmax=[0.1,10.0], nzz=75, $
            spsmodels=spsmodels, imf=imf, redcurve=redcurve, /igm, pburst=0.2, $
            sfhgrid=1, nmodel=nmodel, age=[0.01,12.0], AV=[0.0,3.0], tau=[0.01,5.0], $
            Zmetal=[0.0008,0.03], oiiihb=[0.0,1.0], /nebular, /delayed, /flatav, $
            modelchunksize=1000L, clobber=clobber
       endif else begin
          nmodel = 20000L
          write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
            prefix=prefix, filterlist=filterlist, zminmax=[6.95,10.0], nzz=20, $
            spsmodels=spsmodels, imf=imf, redcurve=redcurve, /igm, $
            sfhgrid=1, nmodel=nmodel, age=[0.01,0.75], AV=[0.0,0.0], tau=[0.01,1.0], $
            Zmetal=[0.0008,0.019], oiiihb=[0.0,1.0], /nebular, /delayed, $
            clobber=clobber
       endelse
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
       if keyword_set(photoz) eq 0 then z = cat.z_b_1
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
         clobber=clobber, /xlog, galaxy=cat.galaxy, yrange=[30.5,23], $
         xrange=[0.3,7.0]*1D4, nsigma=2.0, index=index
    endif
    
return
end
