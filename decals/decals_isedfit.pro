pro decals_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber, nodecals=nodecals
; jm14sep07siena

    prefix = 'decals'
    isedfit_dir = getenv('DECALS_DIR')+'/isedfit/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    
; gather the photometry
    cat = mrdfits(isedfit_dir+'decals_edr.fits.gz',1)
    these = where(cat.z ge 0.05 and cat.z le 0.7 and cat.zwarning eq 0 and $
      cat.brick_primary eq 'T' and strtrim(cat.class,2) eq 'GALAXY',ngal)
    cat = cat[these]
    splog, 'Number of galaxies ', ngal
    
    decals_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist
    galaxy = cat.brickname+'-'+string(cat.objid,format='(I5.5)')
    zobj = cat.z

; do not use DECaLS photometry in the fitting    
    if keyword_set(nodecals) then begin
       outprefix = 'nodecals'
       ivarmaggies[0:2,*] = 0   ; no weight!
    endif

; --------------------------------------------------
; build the parameter files
    if keyword_set(write_paramfile) then begin
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, zminmax=[0.05,0.7], nzz=50, $
         spsmodels='bc03_basel', imf='chab', redcurve='calzetti', igm=0, $
         sfhgrid=1, nmodel=5000L, age=[0.01,13.0], AV=[0.35,2.0], tau=[0.0,5.0], $
         Zmetal=[0.0004,0.05], oiiihb=[-1.0,1.0], /nebular, /delayed, $
         pburst=0.1, clobber=clobber
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
       index = lindgen(30)
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, nrandom=40, galaxy=galaxy, $
         index=index, outprefix=outprefix
    endif
    
return
end
