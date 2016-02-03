pro desi_ages_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, isedfit=isedfit, kcorrect=kcorrect, $
  build_oiiflux=build_oiiflux, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber
; jm15apr23siena - fit the parent sample of AGES galaxies for the DESI
; project

;   echo "desi_ages_isedfit, /write_param, /build_grids, /model_phot, /isedfit, /cl" | /usr/bin/nohup idl > & ~/desi-ages-isedfit.log & 
    
    version = desi_bgs_templates_version(/isedfit)

    prefix = 'desi_ages'
    isedfit_dir = getenv('IM_ARCHIVE_DIR')+'/projects/desi/templates/'+$
      'bgs_templates/isedfit/'+version+'/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

; select the sample    
    zminmax = [0.01,0.85]
    nzz = 100
    
    phot = read_ages(/photo)
    index = where((phot.imain eq 1) and (phot.z ge zminmax[0]) and $
      (phot.z le zminmax[1]),ngal)
;   index = index[0:30] & ngal = n_elements(index)

    ages_to_maggies, phot, maggies, ivarmaggies, /totalmag, $
      filterlist=filterlist, /itot, use_aper='04'

; don't use the GALEX, Y-band, or Spitzer/ch3-4 photometry
    keep = where(strmatch(filterlist,'*galex*') eq 0 and $
      strmatch(filterlist,'*yfilter*') eq 0 and $
      strmatch(filterlist,'*ch3*') eq 0 and $
      strmatch(filterlist,'*ch4*') eq 0)
    maggies = maggies[keep,*]
    ivarmaggies = ivarmaggies[keep,*]
    filterlist = filterlist[keep]

; --------------------------------------------------
; write the parameter file
    if keyword_set(write_paramfile) then begin
       spsmodels = 'fsps_v2.4_miles'
       imf = 'chab'
       redcurve = 'charlot'
       Zmetal = [0.004,0.03]
       age = [0.1,12.5]
       tau = [0.0,10]
       nmodel = 30000L
       pburst = 0.2
       interval_pburst = 2.0

; without emission lines
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels=spsmodels, $
         imf=imf, redcurve=redcurve, igm=0, zminmax=zminmax, nzz=nzz, $
         nmodel=nmodel, age=age, tau=tau, Zmetal=Zmetal, $
         pburst=pburst, interval_pburst=interval_pburst, clobber=clobber
; with emission lines
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels=spsmodels, $
         imf=imf, redcurve=redcurve, igm=0, zminmax=zminmax, nzz=nzz, $
         nmodel=nmodel, age=age, tau=tau, Zmetal=Zmetal, $
         pburst=pburst, interval_pburst=interval_pburst, /nebular, /append
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
; fit!
    if keyword_set(isedfit) then begin
       isedfit, isedfit_paramfile, maggies, ivarmaggies, phot.z, ra=phot.ra, $
         dec=phot.dec, isedfit_dir=isedfit_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, index=index, outprefix=outprefix
    endif 

; --------------------------------------------------
; compute K-corrections
    if keyword_set(kcorrect) then begin
;      index = [9,23,29]
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         absmag_filterlist=sdss_filterlist(), band_shift=0.0, $
         clobber=clobber, index=index, outprefix=outprefix
    endif 

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
       galaxy = phot.galaxy
       these = shuffle_indx(ngal,num=50)
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, galaxy=galaxy, index=index[these]
    endif

return
end
