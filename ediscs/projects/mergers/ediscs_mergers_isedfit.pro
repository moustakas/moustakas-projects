pro ediscs_mergers_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, clobber=clobber, thissfhgrid=thissfhgrid, $
  use_clusterz=use_clusterz, use_specz=use_specz
  
; jm17nov21siena - 

; echo "ediscs_mergers_isedfit, /write_param, /build_grids, /model_phot, /cl" | /usr/bin/nohup idl > ~/mergers.log 2>&1 &
; echo "ediscs_mergers_isedfit, /isedfit, /qaplot_sed, /cl" | /usr/bin/nohup idl > ~/mergers-isedfit.log 2>&1 &
    
    prefix = 'ediscs_mergers'
    isedfit_dir = ediscs_path(/projects)+'mergers/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    ediscs_mergers_dir = ediscs_path(/projects)+'mergers/'

    if keyword_set(use_specz) then begin
       phot = read_ediscs(/phot)
       ediscs_to_maggies, phot, maggies, ivarmaggies, filterlist=filterlist

       spec1d = read_ediscs(/spec1d)
       zobj = spec1d.z

       outprefix = prefix+'_specz'
    endif else begin
       phot = read_ediscs(/allphot)
       ediscs_to_maggies, phot, maggies, ivarmaggies, filterlist=filterlist
       
       if keyword_set(use_clusterz) then begin
          info = rsex(ediscs_path(/cat)+'ediscs_clusters.sex')
          shortcl = strarr(n_elements(info))
          for ii = 0, n_elements(info)-1 do begin
             cl = strtrim(info[ii].cluster_fullname,2)
             shortcl[ii] = strmid(cl,0,6)+'-'+strmid(cl,9,4)
          endfor
          niceprint, shortcl, info.cluster_fullname, info.z
          match2, strtrim(phot.cluster,2), shortcl, m1, m2

          zobj = fltarr(n_elements(phot))-1.0
          good = where(m1 ne -1) ; no photometry for cl1122
          zobj[good] = info[m1[good]].z

          index = where(zobj gt 0 and phot.starflag eq 0 and (total(ivarmaggies gt 0,1) ge 3))
          
          outprefix = prefix+'_zcl'
       endif else begin
          ; individual photoz
          index = where(phot.bestz gt 0.05 and phot.bestz lt 1.5 and $
            phot.starflag eq 0 and (total(ivarmaggies gt 0,1) ge 3))
          zobj = phot.bestz
          outprefix = prefix+'_photoz'
       endelse
    endelse

; --------------------------------------------------
; write the iSEDfit parameter file 
    if keyword_set(write_paramfile) then begin
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels='fsps_v2.4_miles', $
         imf='chab', redcurve='charlot', igm=0, zminmax=[0.05,1.5], nzz=75, $
         nmodel=25000L, age=[0.1,12.5], tau=[0.1,10.0], Zmetal=[0.004,0.03], $
         pburst=0.2, /nebular, /delayed, clobber=clobber
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

;; --------------------------------------------------
;; generate the model photometry QAplots
;    if keyword_set(qaplot_models) then begin
;       thesefilters = ['galex_NUV','sdss_g0','sdss_r0','sdss_i0','wise_w1']
;       phot = mrdfits(ediscs_mergers_dir+'ediscs_mergers_'+ver+'_phot.fits.gz',1)
;       if keyword_set(bcgs) then begin
;          outprefix = 'bcgs'
;          phot = phot[where(phot.isbcg)]
;       endif
;       isedfit_qaplot_models, isedfit_paramfile, phot.maggies, phot.ivarmaggies, $
;         phot.z, isedfit_dir=isedfit_dir, outprefix=outprefix, thissfhgrid=thissfhgrid, $
;         thesefilters=thesefilters, clobber=clobber
;    endif
    
; --------------------------------------------------
; do the fitting!
    if keyword_set(isedfit) then begin
       isedfit, isedfit_paramfile, maggies, ivarmaggies, zobj, $
         isedfit_dir=isedfit_dir, outprefix=outprefix, isedfit_results=ised, $
         isedfit_post=isedpost, clobber=clobber, thissfhgrid=thissfhgrid, $
         photoz=photoz, index=index, ra=phot.ra, dec=phot.dec
    endif 

; --------------------------------------------------
; compute K-corrections
    if keyword_set(kcorrect) then begin
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, outprefix=outprefix
    endif 

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
       isedfit_qaplot_sed, isedfit_paramfile, nrandom=50, outprefix=outprefix, $
         isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, $
         thissfhgrid=thissfhgrid, clobber=clobber, galaxy=galaxy
    endif
    
return
end


