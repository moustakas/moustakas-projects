pro ediscs_lss_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, isedfit=isedfit, kcorrect=kcorrect, $
  build_oiiflux=build_oiiflux, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber, cfht=cfht
; jm14jul22siena - fit the EDisCS sample for Pascale
    
    isedfit_dir = getenv('IM_PROJECTS_DIR')+'/ediscs/lss/'
    montegrids_dir = isedfit_dir+'montegrids/'

    zminmax = [0.05,1.37]
    zbin = 0.05

    if keyword_set(cfht) then begin
       prefix = 'cfht_lss'
       phot = mrdfits(isedfit_dir+'CFHT_photometry_all_spec.fits',1)

       filterlist = cfhtls_filterlist()
       maggies = transpose([[phot.u],[phot.g],[phot.r],[phot.i],[phot.z]])
       ivarmaggies = transpose([[phot.u_ivar],[phot.g_ivar],[phot.r_ivar],[phot.i_ivar],[phot.z_ivar]])
       fix = where(finite(ivarmaggies) eq 0)
       ivarmaggies[fix] = 0
       
       zobj = phot.z_spec
       ra = phot.ra
       dec = phot.dec
       index = where(zobj ge zminmax[0] and zobj le zminmax[1])
    endif else begin
       prefix = 'ediscs_lss'
       phot = read_ediscs(/phot)
       spec1d = read_ediscs(/spec1d)
       zobj = spec1d.z
       ra = spec1d.ra
       dec = spec1d.dec
       ediscs_to_maggies, phot, maggies, ivarmaggies, filterlist=filterlist
    endelse
    ngal = n_elements(phot)

    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    
; --------------------------------------------------
; write the parameter file
    if keyword_set(write_paramfile) then begin
       spsmodels = 'fsps_v2.4_miles'
       imf = 'chab'
       redcurve = 'charlot'
       Zmetal = [0.004,0.03]
       age = [0.1,12.5]
       tau = [0.0,10]
       nmodel = 20000L
       pburst = 0.2
       interval_pburst = 2.0

       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels=spsmodels, $
         imf=imf, redcurve=redcurve, igm=0, zminmax=zminmax, nzz=nzz, $
         zbin=zbin, nmodel=nmodel, age=age, tau=tau, Zmetal=Zmetal, $
         pburst=pburst, interval_pburst=interval_pburst, /nebular, /cl
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
       isedfit, isedfit_paramfile, maggies, ivarmaggies, zobj, ra=ra, $
         dec=dec, isedfit_dir=isedfit_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, index=index, outprefix=outprefix
    endif 

; --------------------------------------------------
; compute K-corrections
    if keyword_set(kcorrect) then begin
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         absmag_filterlist=sdss_filterlist(), band_shift=0.0, $
         clobber=clobber, index=index, outprefix=outprefix
    endif 

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
       if keyword_set(cfht) then begin
          these = shuffle_indx(ngal,num=20)
          isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
            montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
            clobber=clobber, /xlog, galaxy='Obj'+strtrim(index[these],2), $
            index=index[these], outprefix=outprefix
       endif else begin
          these = shuffle_indx(ngal,num=20)
          isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
            montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
            clobber=clobber, /xlog, galaxy=spec1d.galaxy, index=these, $
            outprefix=outprefix
       endelse
    endif

return
end
