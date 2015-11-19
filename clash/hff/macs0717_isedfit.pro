pro macs0717_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber, photoz=photoz
; jm15sep11siena

    if keyword_set(photoz) then prefix = 'macs0717_photoz' else prefix = 'macs0717'
    isedfit_dir = getenv('CLASH_PROJECTS')+'/hff/macs0717/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    
; gather the photometry
    date = '15oct10'
    cat = rsex(isedfit_dir+'macs0717.'+date)

    filterlist = hff_filterlist(/useirac)
    nfilt = n_elements(filterlist)

; --------------------------------------------------
; build the parameter files
    if keyword_set(write_paramfile) then begin
       spsmodels = 'fsps_v2.4_miles'
       imf = 'chab'
       if keyword_set(photoz) then begin
          nmodel = 20000L
          redcurve = 'charlot'
          write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
            prefix=prefix, filterlist=filterlist, zminmax=[0.1,12.0], nzz=150, $
            spsmodels=spsmodels, imf=imf, redcurve=redcurve, /igm, pburst=0.2, $
            sfhgrid=1, nmodel=nmodel, age=[0.01,12.0], AV=[0.0,3.0], tau=[0.0,5.0], $
            Zmetal=[0.0008,0.03], oiiihb=[-0.3,1.0], /nebular, /delayed, /flatav, $
            modelchunksize=1000L, clobber=clobber
       endif else begin
          nmodel = 20000L
          write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
            prefix=prefix, filterlist=filterlist, zminmax=[6.9,12.0], zbin=0.05, $
            spsmodels=spsmodels, imf=imf, redcurve='none', /igm, $
            sfhgrid=1, nmodel=nmodel, age=[0.01,0.75], AV=[0.0,0.0], tau=[0,1.0], $
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
       zobj = cat.bpz
       isedfit, isedfit_paramfile, maggies, ivarmaggies, zobj, thissfhgrid=thissfhgrid, $
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
       galaxy = strtrim(cat.id,2)
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, galaxy=galaxy, yrange=[30.5,23], $
         xrange=[0.3,7.0]*1D4, nsigma=2.0, index=index, outprefix=outprefix
    endif
    
return
end
