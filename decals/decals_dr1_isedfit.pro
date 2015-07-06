pro decals_dr1_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber, nodecals=nodecals, nosdss=nosdss
; jm14sep07siena

; echo "decals_dr1_isedfit, /ised, /qaplot_sed, /cl" | /usr/bin/nohup idl > & logall.log & 
; echo "decals_dr1_isedfit, /ised, /qaplot_sed, /nosdss, /cl" | /usr/bin/nohup idl > & lognosdss.log & 
; echo "decals_dr1_isedfit, /ised, /qaplot_sed, /nodecals, /cl" | /usr/bin/nohup idl > & lognodecals.log & 

    prefix = 'decals_dr1'
    dr1_dir = getenv('DECALS_DIR_DR1')+'/'
    isedfit_dir = getenv('IM_ARCHIVE_DIR')+'/decam/isedfit/dr1/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

; define the catalog
    tractor = mrdfits(dr1_dir+'decals-specz.fits',1)
    specz = mrdfits(dr1_dir+'decals-specz.fits',2)
    phot = mrdfits(dr1_dir+'decals-specz.fits',3)
    these = where(specz.z ge 0.05 and specz.z le 0.7 and specz.zwarning eq 0 and $
      strtrim(specz.class,2) eq 'GALAXY' and tractor.decam_nobs[1] ge 1 and $
      tractor.decam_nobs[2] ge 1 and tractor.decam_nobs[4] ge 1,ngal)
;   these = these[0:50] & ngal = n_elements(these)
    tractor = tractor[these]
    specz = specz[these]
    phot = phot[these]
    splog, 'Number of galaxies ', ngal
    
    galaxy = tractor.brickname+'-'+string(tractor.objid,format='(I5.5)')
       
; gather the photometry    
    decals_to_maggies, tractor, dmaggies, divarmaggies, $
      filterlist=dfilterlist, /shortwise, /decam_grz
    
;   sdss_to_maggies, smaggies, sivarmaggies, calib=phot, flux='model'
    sdss_to_maggies, modelmaggies, modelivarmaggies, calib=phot, flux='model'
    sdss_to_maggies, cmodelmaggies, cmodelivarmaggies, calib=phot, flux='cmodel'
    ratio = cmodelmaggies[2,*]/(modelmaggies[2,*]+(modelmaggies[2,*] eq 0))
    neg = where(modelmaggies[2,*] le 0)
    if (neg[0] ne -1L) then ratio[neg] = 1.0
    
    factor = rebin(ratio,5,ngal)
    smaggies = modelmaggies*factor
    sivarmaggies = modelivarmaggies/factor^2
    
    maggies = [dmaggies,smaggies]
    ivarmaggies = [divarmaggies,sivarmaggies]
    filterlist = [dfilterlist,sdss_filterlist()]
    zobj = specz.z
       
; do not use DECaLS photometry in the fitting    
    if keyword_set(nodecals) then begin
       outprefix = 'nodecals'
       ivarmaggies[0:2,*] = 0   ; no weight!
    endif
       
; do not use SDSS photometry in the fitting    
    if keyword_set(nosdss) then begin
       outprefix = 'nosdss'
       ivarmaggies[5:9,*] = 0   ; no weight!
    endif
       
; --------------------------------------------------
; build the parameter files
    if keyword_set(write_paramfile) then begin
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, zminmax=[0.05,0.7], nzz=30, $
         spsmodels='fsps_v2.4_miles', imf='chab', redcurve='charlot', igm=0, $
         sfhgrid=1, nmodel=10000L, age=[0.01,13.0], AV=[0.35,2.0], tau=[0.0,5.0], $
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
         ra=tractor.ra, dec=tractor.dec
    endif 

; --------------------------------------------------
; compute K-corrections
    if keyword_set(kcorrect) then begin
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         absmag_filterlist=galex_filterlist(), band_shift=0.0, $
         clobber=clobber, outprefix=outprefix
    endif 

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
       index = lindgen(50)
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, nrandom=40, galaxy=galaxy, $
         index=index, outprefix=outprefix
    endif
    
return
end
