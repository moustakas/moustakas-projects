pro ediscs_mergers_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  merge_isedfit=merge_isedfit, kcorrect=kcorrect, qaplot_sed=qaplot_sed, clobber=clobber, $
  thissfhgrid=thissfhgrid, bcgs=bcgs, firstchunk=firstchunk, lastchunk=lastchunk
; jm17nov21siena - 

; echo "ediscs_mergers_isedfit, /write_param, /build_grids, /model_phot, /cl" | /usr/bin/nohup idl > ~/mergers.log 2>&1 &
; echo "ediscs_mergers_isedfit, /isedfit, /qaplot_sed, /cl" | /usr/bin/nohup idl > ~/mergers-isedfit.log 2>&1 &
    
    prefix = 'ediscs_mergers'
    isedfit_dir = ediscs_path(/projects)+'mergers/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    ediscs_mergers_dir = ediscs_path(/projects)+'mergers/'

    filterlist = ediscs_filterlist()

;   phot = mrdfits(ediscs_mergers_dir+'ediscs_mergers_v1.fits',1)
;   ediscs_mergers_to_maggies, phot, maggies, ivarmaggies
;   zobj = phot.zredmagic
    
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
; fit! but split the members catalog into chunks since it's so big 
    if keyword_set(junk_isedfit) then begin
       if keyword_set(bcgs) then begin
          phot = mrdfits(ediscs_mergers_dir+'ediscs_mergers_'+ver+'_phot.fits.gz',1)

          outprefix = 'bcgs'
          phot = phot[where(phot.isbcg)]
;         phot = phot[0:100]

          isedfit, isedfit_paramfile, phot.maggies, phot.ivarmaggies, $
            phot.z, ra=phot.ra, dec=phot.dec, isedfit_dir=isedfit_dir, $
            outprefix=outprefix, thissfhgrid=thissfhgrid, clobber=clobber
       endif else begin
          splog, 'Hard-coding NGAL to speed things up!'
          ngal = 11235932L
          chunksize = ceil(ngal/float(nchunk))
          if n_elements(firstchunk) eq 0 then firstchunk = 0
          if n_elements(lastchunk) eq 0 then lastchunk = nchunk-1

          for ii = firstchunk, lastchunk do begin
             splog, 'Working on CHUNK '+strtrim(ii,2)+', '+strtrim(lastchunk,2)
             splog, im_today()
             t0 = systime(1)
             outprefix = 'ediscs_mergers_chunk'+string(ii,format='(I3.3)')
             these = lindgen(chunksize)+ii*chunksize
             these = these[where(these lt ngal)]
;            these = these[0:99] ; test!

             phot = mrdfits(ediscs_mergers_dir+'ediscs_mergers_'+ver+'_phot.fits.gz',1,rows=these)
             isedfit, isedfit_paramfile, phot.maggies, phot.ivarmaggies, $
               phot.z, ra=phot.ra, dec=phot.dec, isedfit_dir=isedfit_dir, $
               outprefix=outprefix, thissfhgrid=thissfhgrid, clobber=clobber
             splog, 'Total time (min) = '+strtrim((systime(1)-t0)/60.0,2)
          endfor                    
       endelse
    endif 

; --------------------------------------------------
; merge the iSEDfit results
    if keyword_set(merge_isedfit) then begin
       params = read_isedfit_paramfile(isedfit_paramfile,thissfhgrid=thissfhgrid)
       fp = isedfit_filepaths(params,isedfit_dir=isedfit_dir)
;      for ii = 0, 10 do begin
       for ii = 0, nchunk-1 do begin
          outprefix = 'ediscs_mergers_chunk'+string(ii,format='(I3.3)')
          ised1 = read_isedfit(isedfit_paramfile,isedfit_dir=isedfit_dir,$
            outprefix=outprefix)
          if ii eq 0 then ised = temporary(ised1) else ised = [temporary(ised),temporary(ised1)]
       endfor
;      im_mwrfits, ised, ediscs_mergers_dir+fp.isedfit_outfile, clobber=clobber
       im_mwrfits, ised, fp.isedfit_dir+fp.isedfit_outfile, clobber=clobber
    endif

; --------------------------------------------------
; compute K-corrections
    if keyword_set(kcorrect) then begin
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber
    endif 

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
;      phot = mrdfits(ediscs_mergers_dir+'ediscs_mergers_'+ver+'_phot.fits.gz',1)
;      galaxy = 'BCG '+string(phot.mem_match_id,format='(I6.6)')
;      allcl = phot.mem_match_id
;      cl = allcl[uniq(allcl,sort(allcl))]

       isedfit_qaplot_sed, isedfit_paramfile, nrandom=50, outprefix=outprefix, $
         isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, $
         thissfhgrid=thissfhgrid, clobber=clobber, galaxy=galaxy
    endif
    
return
end


