pro clash_spitzer_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  kcorrect=kcorrect, qaplot_sed=qaplot_sed, clobber=clobber, $
  thissfhgrid=thissfhgrid
; jm13aug22siena - 

    clash = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')

;; parse Leonardo's Spitzer catalogs
;    path = '/moustakas-archive/clash-archive/spitzer/merged_catalogs_temp/'
;;   for ii = 10, 10 do begin
;    for ii = 0, n_elements(clash)-1 do begin
;       catfile = path+strtrim(clash[ii].dirname,2)+'_final.cat'
;       if file_test(catfile) then begin
;          splog, 'Reading '+catfile
;          cc = rsex(catfile)
;          keep1 = where(cc.redshift gt 0.0 and cc.redshift lt 100.0,nkeep1)
;          splog, file_basename(catfile), nkeep1
;          if nkeep1 gt 0 then begin
;             clash_to_maggies, cc[keep1], mm, ivar, /useirac
;             keep2 = where(total(ivar[17:18,*] gt 0,1) ge 1.0 and $ ; at least one IRAC
;               total(ivar[0:16,*] gt 0,1) ge 3,ngal)
;             splog, file_basename(catfile), ngal
;             if ngal gt 0 then begin
;                for jj = 0, n_tags(cc)-1 do begin
;                   fix = where(finite(cc[keep1[keep2]].(jj)) eq 0,nfix)
;                   if nfix ne 0L then cc[keep1[keep2[fix]]].(jj) = -99.0
;                endfor
;                outfile = getenv('IM_PROJECTS_DIR')+'/clash/spitzer/'+$
;                  strtrim(clash[ii].dirname,2)+'_final_specz.cat'
;                wsex, cc[keep1[keep2]], outfile=outfile
;             endif
;          endif
;       endif else splog, 'No catalog '+catfile
;    endfor

    isedfit_dir = getenv('IM_RESEARCH_DIR')+'/projects/clash/spitzer/'
    montegrids_dir = isedfit_dir+'montegrids/'

    prefix = 'spitzer'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

; read and stack all the clusters together
    for ii = 0, n_elements(clash)-1 do begin
       catfile = isedfit_dir+strtrim(clash[ii].dirname,2)+'_final_specz.cat'
       if file_test(catfile) then begin
          cat1 = rsex(catfile)
          cat1 = struct_addtags(replicate({cluster: strtrim(clash[ii].shortname,2)},$
            n_elements(cat1)),cat1)
          if ii eq 0 then cat = cat1 else cat = [cat,cat1]
       endif
    endfor
    
    clash_to_maggies, cat, maggies, ivarmaggies, /useirac, $
      filterlist=filterlist
    zz = cat.redshift

    zminmax = minmax(zz)
    zbin = 0.03

; --------------------------------------------------
; (mandatory) choose your priors: write the iSEDfit parameter file 
    if keyword_set(write_paramfile) then begin
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels='fsps_v2.4_miles', $
         imf='salp', redcurve='charlot', /igm, zminmax=zminmax, zbin=zbin, $
         nmodel=1000L, age=[0.1,11.0], tau=[0.1,10.0], Zmetal=[0.004,0.03], $
         pburst=0.1, interval_pburst=1.5, tburst=[0.1,11.0], /delayed, $
         /nebular, clobber=clobber
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
       thesefilters = ['clash_wfc3_f390w.par','clash_acs_f814w',$
         'clash_wfc3_f110w','clash_wfc3_f160w']
       isedfit_qaplot_models, isedfit_paramfile, maggies, ivarmaggies, $
         zz, isedfit_dir=isedfit_dir, thissfhgrid=thissfhgrid, $
         thesefilters=thesefilters, clobber=clobber
    endif
    
; --------------------------------------------------
; (mandatory) fit!
    if keyword_set(isedfit) then begin
       isedfit, isedfit_paramfile, maggies, ivarmaggies, zz, $
         ra=cat.ra, dec=cat.dec, isedfit_dir=isedfit_dir, $
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
