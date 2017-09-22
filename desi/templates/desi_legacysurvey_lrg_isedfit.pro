pro desi_legacysurvey_lrg_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, isedfit=isedfit, kcorrect=kcorrect, $
  build_oiiflux=build_oiiflux, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber
; jm15apr23siena - fit the parent sample of AGES galaxies for the DESI
; project

;   echo "desi_ages_isedfit, /write_param, /build_grids, /model_phot, /isedfit, /cl" | /usr/bin/nohup idl > & ~/desi-ages-isedfit.log & 
;   echo "desi_ages_isedfit, /kcorrect, thissfhgrid=2, /cl" | /usr/bin/nohup idl > ~/desi-ages-kcorrect.log 2>&1 & 
    
    version = 'v1.0'
;   version = desi_lrg_templates_version(/isedfit)
    splog, 'Fixing the iSEDfit version number to '+version+'!'

    prefix = 'legacysurvey_lrg'
    isedfit_dir = getenv('IM_ARCHIVE_DIR')+'/data/desi-archive/templates/'+$
      'lrg_templates/isedfit/'+version+'/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

; select the sample    
    zminmax = [0.2,1.2]
    nzz = 51 ; 101

    catfile = isedfit_dir+'targets-dr3.1-0.14.0-lrg-rf-photoz-0.2.fits'
    print, 'Reading '+catfile
    cat = mrdfits(catfile, 1)

    legacysurvey_to_maggies, cat, allmaggies, allivarmaggies, filterlist=filterlist
    nfilt = n_elements(filterlist)
    nobsmin = 3

    these = where( (cat.z_spec gt zminmax[0]) and (cat.z_spec lt zminmax[1]) and $
      (cat.nobs_g ge nobsmin) and (cat.nobs_r ge nobsmin) and (cat.nobs_z ge nobsmin) and $
      (total(allivarmaggies gt 0, 1) eq nfilt),ngal)

    maggies = allmaggies[*, these]
    ivarmaggies = allivarmaggies[*, these]
    zobj = cat[these].z_spec
    ra = cat[these].ra
    dec = cat[these].dec

;   gr = -2.5*alog10(maggies[0, *] / maggies[1, *])
;   rz = -2.5*alog10(maggies[1, *] / maggies[2, *])
;   plot, rz, gr, psym=3, xsty=3, ysty=3
    
; --------------------------------------------------
; write the parameter file
    if keyword_set(write_paramfile) then begin
       spsmodels = 'ckc14z'
       imf = 'kroupa01'
;      spsmodels = 'fsps_v2.4_miles' ; v1.0
;      imf = 'chab'
       redcurve = 'charlot'
       Zmetal = [0.004,0.03]
       age = [1.0,11.05]
       tau = [0.0,6]
       nmodel = 5000L

; no emission lines, no bursts
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels=spsmodels, $
         imf=imf, redcurve=redcurve, igm=0, zminmax=zminmax, nzz=nzz, $
         nmodel=nmodel, age=age, tau=tau, Zmetal=Zmetal, $
         /delayed, nebular=0, clobber=clobber
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
;      splog, 'Temporary hack!'
;      toss = where(filterlist eq 'spitzer_irac_ch2.par')
;      ivarmaggies[toss,*] = 0
       isedfit, isedfit_paramfile, maggies, ivarmaggies, zobj, ra=ra, $
         dec=dec, isedfit_dir=isedfit_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, index=index, outprefix=outprefix
    endif 

; --------------------------------------------------
; compute K-corrections
    if keyword_set(kcorrect) then begin
;      index = [9,23,29]
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         absmag_filterlist=legacysurvey_filterlist(), band_shift=0.0, $
         clobber=clobber, index=index, outprefix=outprefix
    endif 

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
       these = lindgen(50)
;      these = lindgen(ngal)
;      these = shuffle_indx(ngal,num=50)
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, index=these, /xlog
    endif

return
end
