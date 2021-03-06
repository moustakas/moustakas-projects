pro legacysurvey_lrgs_dr8_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, isedfit=isedfit, kcorrect=kcorrect, $
  build_oiiflux=build_oiiflux, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber
; jm19nov05siena - fit the DR8 LRGs with existing spectroscopic redshifts

;   echo "legacysurvey_lrgs_dr8_isedfit, /write_param, /build_grids, /model_phot, /cl" | /usr/bin/nohup idl > & ~/desi-lrg-isedfit.log & 
;   echo "legacysurvey_lrgs_dr8_isedfit, /isedfit, /kcorrect, /qaplot_sed, /cl" | /usr/bin/nohup idl > ~/lrg-isedfit-v2.0.log 2>&1 & 
    
    prefix = 'lrgs_dr8'
    isedfit_dir = getenv('IM_ARCHIVE_DIR')+'/projects/legacysurvey/lrgs/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

; select the sample    
    zminmax = [0.2,1.2]
    nzz = 51 ; 101

    catfile = isedfit_dir+'lrgs-dr8.fits'
    print, 'Reading '+catfile
    cat = mrdfits(catfile, 1);, range=[0, 99])

    legacysurvey_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist
    nfilt = n_elements(filterlist)
    nobsmin = 3

    zobj = cat.z
    ra = cat.ra
    dec = cat.dec

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
