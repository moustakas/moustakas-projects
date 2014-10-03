function deep2_get_oiiflux, ppxf, cflux_3727_rest=cflux_3727_rest
; jm14mar11siena - given a pPXF-style emission-line catalog for DEEP2
; and the rest-frame continuum flux at 3727 A, return a data structure
; which contains a DESI-style catalog of [OII] fluxes
;
; recall that all the measured quantities are in the rest-frame,
; including OII_3727_EW [A], OII_3727_AMP [erg/s/cm2/A],
; OII_3727_CONTINUUM [erg/s/cm2/A]
;
; the integrated flux is redshift-independent (other than the distance
; modulus dependence) OII_3727 [erg/s/cm2]

    ngal = n_elements(ppxf)

    oii = struct_trimtags(ppxf,select=['oii_3727*'],except=$
      ['*wave*','*continuum*','*linez*','*sigma*'])
    oii = struct_addtags(im_struct_trimtags(ppxf,select='sigma_forbidden',$
      newtags='sigma_kms'),temporary(oii))

    zobj = ppxf.z
    cflux_3727_obs = cflux_3727_rest/(1.0+zobj) ; [erg/s/cm2/A]
    
; [OII] 3726    
    good = where(oii.oii_3727_1[1] gt 0,ngood,comp=crap,ncomp=ncrap)
    if ngood ne 0L then begin
       Fc_3727_rest = cflux_3727_rest[good] ; *rest-frame* [erg/s/cm2/A]
       Fc_3727_obs = cflux_3727_obs[good]   ; *observed-frame* [erg/s/cm2/A]

       oii[good].oii_3727_1 = oii[good].oii_3727_1_ew*rebin(reform(Fc_3727_rest,1,ngood),2,ngood) ; [erg/s/cm2]
       oii[good].oii_3727_1_limit = oii[good].oii_3727_1_ew_limit*Fc_3727_rest
       
       oii[good].oii_3727_1_amp = oii[good].oii_3727_1_amp*$ ; *rest-frame* [erg/s/cm2/A]
         rebin(reform(Fc_3727_rest,1,ngood),2,ngood)/$ 
         rebin(reform(ppxf[good].oii_3727_1_continuum[0],1,ngood),2,ngood)
    endif

; [OII] 3729
    good = where(oii.oii_3727_2[1] gt 0,ngood,comp=crap,ncomp=ncrap)
    if ngood ne 0L then begin
       Fc_3727_rest = cflux_3727_rest[good] ; *rest-frame* [erg/s/cm2/A]
       Fc_3727_obs = cflux_3727_obs[good]   ; *observed-frame* [erg/s/cm2/A]

       oii[good].oii_3727_2 = oii[good].oii_3727_2_ew*rebin(reform(Fc_3727_rest,1,ngood),2,ngood) ; [erg/s/cm2]
       oii[good].oii_3727_2_limit = oii[good].oii_3727_2_ew_limit*Fc_3727_rest
       
       oii[good].oii_3727_2_amp = oii[good].oii_3727_2_amp*$ ; *rest-frame* [erg/s/cm2/A]
         rebin(reform(Fc_3727_rest,1,ngood),2,ngood)/$ 
         rebin(reform(ppxf[good].oii_3727_2_continuum[0],1,ngood),2,ngood)
    endif

; [OII] 3726,3729
    good = where(oii.oii_3727[1] gt 0,ngood,comp=crap,ncomp=ncrap)
    if ngood ne 0L then begin
       Fc_3727_rest = cflux_3727_rest[good] ; *rest-frame* [erg/s/cm2/A]
       Fc_3727_obs = cflux_3727_obs[good]   ; *observed-frame* [erg/s/cm2/A]

       oii[good].oii_3727 = oii[good].oii_3727_ew*rebin(reform(Fc_3727_rest,1,ngood),2,ngood) ; [erg/s/cm2]
       oii[good].oii_3727_limit = oii[good].oii_3727_ew_limit*Fc_3727_rest
       
       oii[good].oii_3727_amp = oii[good].oii_3727_amp*$ ; *rest-frame* [erg/s/cm2/A]
         rebin(reform(Fc_3727_rest,1,ngood),2,ngood)/$ 
         rebin(reform(ppxf[good].oii_3727_continuum[0],1,ngood),2,ngood)
    endif

; unpack the arrays in the catalog
    oiicat = struct_trimtags(ppxf,select=['objno','ra','dec','z'])
    oiicat = struct_addtags(oiicat,struct_trimtags(oii,select=['sigma_kms']))
    oiicat = struct_addtags(oiicat,replicate($
      {$
      oii_3726:         0.0,$
      oii_3726_err:    -2.0,$
      oii_3729:         0.0,$
      oii_3729_err:    -2.0,$
      oii_3727:         0.0,$
      oii_3727_err:    -2.0,$
      oii_3727_ew:      0.0,$
      oii_3727_ew_err: -2.0},ngal))

    good = where(cflux_3727_rest gt 0.0,comp=bad)
    oiicat[good].oii_3726 = oii[good].oii_3727_1[0]
    oiicat[good].oii_3726_err = oii[good].oii_3727_1[1]
    oiicat[good].oii_3729 = oii[good].oii_3727_2[0]
    oiicat[good].oii_3729_err = oii[good].oii_3727_2[1]
    oiicat[good].oii_3727 = oii[good].oii_3727[0]
    oiicat[good].oii_3727_err = oii[good].oii_3727[1]
    oiicat[good].oii_3727_ew = oii[good].oii_3727_ew[0]
    oiicat[good].oii_3727_ew_err = oii[good].oii_3727_ew[1]

; assign upper limits as if they were fluxes
    limit = where(cflux_3727_rest gt 0.0 and oiicat.oii_3727 eq 0.0 and $
      oii.oii_3727_limit gt 0,nlimit)
    oiicat[limit].oii_3726 = oii[limit].oii_3727_1_limit
    oiicat[limit].oii_3726_err = oii[limit].oii_3727_1_limit ; set uncertainty=flux
    oiicat[limit].oii_3729 = oii[limit].oii_3727_2_limit
    oiicat[limit].oii_3729_err = oii[limit].oii_3727_2_limit
    oiicat[limit].oii_3727 = oii[limit].oii_3727_limit
    oiicat[limit].oii_3727_err = oii[limit].oii_3727_limit
    oiicat[limit].oii_3727_ew = oii[limit].oii_3727_ew_limit
    oiicat[limit].oii_3727_ew_err = oii[limit].oii_3727_ew_limit
    
return, oiicat
end

pro desi_deep2_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, isedfit=isedfit, kcorrect=kcorrect, $
  build_oiiflux=build_oiiflux, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber
; jm13dec18siena - fit the parent sample of DEEP2 galaxies for the
; DESI project; also build the final catalog of [OII] fluxes for the
; DEEP2 sample 
    
    version = desi_elg_templates_version()

    prefix = 'desi_deep2'
    isedfit_dir = getenv('IM_ARCHIVE_DIR')+'/projects/desi/templates/'+$
      'elg_templates/isedfit/'+version+'/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

;   filterlist = deep2_filterlist()
;   cat = mrdfits(isedfit_dir+'deep2_zcat.fits.gz',1)
    
; fit everything in DR4 so that I can use the iSEDfit results for both
; targeting tests and template simulations
    cat = read_deep2_zcat(photo=phot)
    deep2_to_maggies, phot, maggies, ivarmaggies, /unwise, $
      filterlist=filterlist

    zminmax = [0.1,2.0]
    nzz = 140

    index = where(cat.zbest ge zminmax[0] and cat.zbest le zminmax[1])
    ngal = n_elements(cat)

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
         imf=imf, redcurve=redcurve, /igm, zminmax=zminmax, nzz=nzz, $
         nmodel=nmodel, age=age, tau=tau, Zmetal=Zmetal, $
         pburst=pburst, interval_pburst=interval_pburst, clobber=clobber
; with emission lines
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels=spsmodels, $
         imf=imf, redcurve=redcurve, /igm, zminmax=zminmax, nzz=nzz, $
         nmodel=nmodel, age=age, tau=tau, Zmetal=Zmetal, $
         pburst=pburst, interval_pburst=interval_pburst, /nebular, /append
    endif

; --------------------------------------------------
; build the Monte Carlo grids    
    if keyword_set(build_grids) then begin
       isedfit_montegrids, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, clobber=clobber;, $
;        montefile=montegrids_dir+'sfhgrid02/fsps_v2.4_miles/charlot/desi_deep2_fsps_v2.4_miles_chab_montegrid.fits.gz'
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
;      outprefix = 'unwise'
;      index = where(phot.w1_nanomaggies_ivar ne 0 and cat.zbest ge zminmax[0] and $
;        cat.zbest le zminmax[1])
       isedfit, isedfit_paramfile, maggies, ivarmaggies, cat.zbest, ra=cat.ra, $
         dec=cat.dec, isedfit_dir=isedfit_dir, thissfhgrid=thissfhgrid, $
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
; build the [OII] flux catalog using the EWs from pPXF and the
; continuum flux from iSEDfit based on the models which include
; nebular emission (because they have slightly lower chi2)
    if keyword_set(build_oiiflux) then begin
;      kised = mrdfits(isedfit_dir+'desi_deep2_fsps_v2.4_miles_'+$
;        'chab_charlot_sfhgrid01_kcorr.z0.0.fits.gz',1)
       kised = mrdfits(isedfit_dir+'desi_deep2_fsps_v2.4_miles_'+$
         'chab_charlot_sfhgrid02_kcorr.z0.0.fits.gz',1)

; assume a fixed doublet ratio       
       ppxf = read_deep2(/ppxf,/fixoii)
       oiicat = deep2_get_oiiflux(ppxf,cflux_3727_rest=kised.cflux_3727)

; supplement with [OII] measurements where one line was dropped by
; MPFIT...
       ppxf = read_deep2(/ppxf)
       oiicat1 = deep2_get_oiiflux(ppxf,cflux_3727_rest=kised.cflux_3727)

       these = where(oiicat1.oii_3727 gt 0 and oiicat.oii_3727 eq 0,nthese)
       oiicat[these] = oiicat1[these]

; ...but distribute the flux according to the assumed ratio
       doubletratio = 0.73      ; =3726/3729
       factor = [doubletratio/(1.0+doubletratio),1.0/(1.0+doubletratio)]
       oiicat[these].oii_3726 = factor[0]*oiicat[these].oii_3727
       oiicat[these].oii_3729 = factor[1]*oiicat[these].oii_3727
       oiicat[these].oii_3726_err = sqrt(factor[0])*oiicat[these].oii_3727_err
       oiicat[these].oii_3729_err = sqrt(factor[1])*oiicat[these].oii_3727_err

       ww = where(oiicat.oii_3727_err eq -2)
       im_plothist, oiicat.z, bin=0.05
       im_plothist, oiicat[ww].z, bin=0.05, /over, /fill

       im_mwrfits, oiicat, isedfit_dir+'deep2_oiicat_'+version+'.fits', clobber=clobber
    endif

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
;      outprefix = 'unwise'
;      index = (where(phot.w1_nanomaggies_ivar ne 0 and cat.zbest ge zminmax[0] and $
;        cat.zbest le zminmax[1]))[0:30]
       galaxy = 'DEEP2/'+strtrim(cat.objno,2);+'/'+strtrim(cat.source,2)
       these = shuffle_indx(ngal,num=25)
;      these = where(cat[index].objno eq 12024524) & yrange = [24,20]
;      these = where(cat[index].objno eq 12024078)
;      these = where(cat[index].objno eq 12101118)
;      these = where(cat[index].objno eq 12015944)
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, galaxy=galaxy, index=index[these];, yrange=yrange
;        yrange=[26,15] ;, outprefix=outprefix
    endif

return
end
