pro nsa_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber
; jm15jun04siena

    prefix = 'nsa_v1_2'
    isedfit_dir = getenv('IM_DATA_DIR')+'/nsa/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    
; define the catalog
    nsa = read_nsa()
    ngal = n_elements(nsa)
    splog, 'Number of galaxies ', ngal

; gather the photometry    
    maggies = nsa.nmgy*1D-9
    ivarmaggies = nsa.nmgy_ivar*1D18
    zobj = nsa.zdist
    filterlist = [galex_filterlist(),sdss_filterlist()]

; --------------------------------------------------
; build the parameter files
    if keyword_set(write_paramfile) then begin
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, zminmax=[0.0001,0.056], nzz=30, /zlog, $
         spsmodels='fsps_v2.4_miles', imf='chab', redcurve='charlot', igm=0, $
         sfhgrid=1, nmodel=20000L, age=[0.01,13.0], AV=[0.35,2.0], tau=[0.0,5.0], $
         Zmetal=[0.004,0.03], oiiihb=[-1.0,1.0], /nebular, /delayed, $
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
         isedfit_post=isedpost, clobber=clobber, photoz=photoz, index=index, $
         ra=nsa.ra, dec=nsa.dec
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
       nrandom = 100
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, nrandom=nrandom, galaxy=galaxy, $
         outprefix=outprefix; index=index, 
    endif
    
return
end
