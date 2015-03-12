pro decals_edr_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber, nodecals=nodecals
; jm14sep07siena

    prefix = 'decals'
    isedfit_dir = getenv('DECALS_DIR')+'/isedfit/'
    montegrids_dir = isedfit_dir+'montegrids/'
    edr_dir = getenv('HOME')+'/edr/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    
; define the catalog
    tractor = mrdfits(edr_dir+'edr-specz-dr10.fits',1)
    specz = mrdfits(edr_dir+'edr-specz-dr10.fits',2)
    phot = mrdfits(edr_dir+'edr-specz-dr10.fits',3)
    these = where(specz.z ge 0.05 and specz.z le 0.7 and specz.zwarning eq 0 and $
      strtrim(specz.class,2) eq 'GALAXY' and tractor.brick_primary eq 'T',ngal)
    tractor = tractor[these]
    specz = specz[these]
    phot = phot[these]
    splog, 'Number of galaxies ', ngal

    galaxy = tractor.brickname+'-'+string(tractor.objid,format='(I5.5)')

; gather the photometry    
    decals_to_maggies, tractor, dmaggies, divarmaggies, filterlist=filterlist

;   sdss_to_maggies, smaggies, sivarmaggies, calib=phot, flux='model'
    sdss_to_maggies, modelmaggies, modelivarmaggies, calib=phot, flux='model'
    sdss_to_maggies, cmodelmaggies, cmodelivarmaggies, calib=phot, flux='cmodel'
    ratio = cmodelmaggies[2,*]/modelmaggies[2,*]
    neg = where(modelmaggies[2,*] le 0)
    if (neg[0] ne -1L) then ratio[neg] = 1.0

    factor = rebin(ratio,5,ngal)
    smaggies = modelmaggies*factor
    sivarmaggies = modelivarmaggies/factor^2

    maggies = [dmaggies,smaggies]
    ivarmaggies = [divarmaggies,sivarmaggies]
    zobj = specz.z

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
