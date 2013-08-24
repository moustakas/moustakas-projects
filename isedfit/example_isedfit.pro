;+
; NAME:
;   EXAMPLE_ISEDFIT
;
; PURPOSE:
;   Simple script demonstrating one possible way to call iSEDfit.
;
; INPUTS:
;   None required.
;
; OPTIONAL INPUTS:
;   thissfhgrid - 
;
; KEYWORD PARAMETERS:
;   write_paramfile - 
;   build_grids - 
;   model_photometry - 
;   qaplot_models - 
;   isedfit - 
;   kcorrect - 
;   qaplot_sed - 
;   clobber - overwrite existing files of the same name
;
; OUTPUTS:
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2013 Aug 15, Siena
;-

pro example_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, clobber=clobber, $
  thissfhgrid=thissfhgrid

    isedfit_dir = './'  ; /path/to/example/data
    montegrids_dir = isedfit_dir+'montegrids/'

    prefix = 'example'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

; read the data (not needed to do everything, but it's convenient to
; have this piece of code here)
    datafile = isedfit_dir+'exampledata.fits.gz'
    if file_test(datafile) eq 0 then begin
       splog, 'Example data file '+datafile+' not found!'
       return
    endif
    example = mrdfits(datafile,1)

; specify the filters    
    filterlist = [$
      'galex_FUV.par',$
      'galex_NUV.par',$
      'sdss_u0.par',$
      'sdss_g0.par',$
      'sdss_r0.par',$
      'sdss_i0.par',$
      'sdss_z0.par',$
      'wise_w1.par',$
      'wise_w2.par']
    
; --------------------------------------------------
; (mandatory) choose your priors: write the iSEDfit parameter file 
    if keyword_set(write_paramfile) then begin
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels='fsps_v2.4_miles', $
         imf='chab', redcurve='charlot', /igm, zminmax=[0.05,0.2], zbin=0.005, $
         nmodel=25000L, age=[0.1,13.0], tau=[0.1,5.0], Zmetal=[0.004,0.03], $
         AV=[0.35,2.0], mu=[0.1,4.0], pburst=0.2, interval_pburst=2.0, $
         tburst=[0.1,13.0], /delayed, galchunksize=2500L, clobber=clobber
    endif

; --------------------------------------------------
; (mandatory) build the Monte Carlo grids    
    if keyword_set(build_grids) then begin
       isedfit_montegrids, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, clobber=clobber
    endif

; --------------------------------------------------
; (mandatory) calculate the model photometry 
    if keyword_set(model_photometry) then begin
       isedfit_models, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber
    endif

; --------------------------------------------------
; (optional) generate the model photometry QAplots
    if keyword_set(qaplot_models) then begin
       thesefilters = ['galex_NUV','sdss_g0','sdss_r0','sdss_i0','wise_w1']
       isedfit_qaplot_models, isedfit_paramfile, example.maggies, $
         example.ivarmaggies, example.z, isedfit_dir=isedfit_dir, $
         thissfhgrid=thissfhgrid, thesefilters=thesefilters, clobber=clobber
    endif
    
; --------------------------------------------------
; (mandatory) fit!
    if keyword_set(isedfit) then begin
       isedfit, isedfit_paramfile, example.maggies, example.ivarmaggies, $
         example.z, ra=example.ra, dec=example.dec, isedfit_dir=isedfit_dir, $
         thissfhgrid=thissfhgrid, clobber=clobber
    endif 

; --------------------------------------------------
; (optional) compute K-corrections
    if keyword_set(kcorrect) then begin
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber
    endif 

; --------------------------------------------------
; (optional) generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
       isedfit_qaplot_sed, isedfit_paramfile, nrandom=50, $
         isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, $
         thissfhgrid=thissfhgrid, clobber=clobber, /xlog
    endif
    
return
end
