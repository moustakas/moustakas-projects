pro ages_mgfe_isedfit, sdss=sdss, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, qaplot_models=qaplot_models, isedfit=isedfit, $
  merge_isedfit=merge_isedfit, kcorrect=kcorrect, qaplot_sed=qaplot_sed, clobber=clobber, $
  thissfhgrid=thissfhgrid, firstchunk=firstchunk, lastchunk=lastchunk

    isedfit_dir = getenv('IM_ARCHIVE_DIR')+'/projects/ages/mgfe/'
    mgfepath = getenv('IM_PROJECTS_DIR')+'/ages/mgfe/'
    montegrids_dir = isedfit_dir+'montegrids/'

    if keyword_set(sdss) then begin
       prefix = 'sdss_mgfe'
       filterlist = [galex_filterlist(),sdss_filterlist(),(wise_filterlist())[0:1]]
    endif else begin
       prefix = 'ages_mgfe'
       filterlist = ages_filterlist()
       filterlist = filterlist[where(strmatch(filterlist,'*ch3*') eq 0 and $
         strmatch(filterlist,'*ch4*') eq 0)]
    endelse

    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

    nchunk = 8
    band_shift = 0.1
    absmag_filterlist = sdss_filterlist()
    
; --------------------------------------------------
; write the iSEDfit parameter file 
    if keyword_set(write_paramfile) then begin
       if keyword_set(sdss) then begin
          write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
            prefix=prefix, filterlist=filterlist, spsmodels='fsps_v2.4_miles', $
            imf='chab', redcurve='charlot', /igm, zminmax=[0.05,0.2], zbin=0.005, $
            nmodel=25000L, age=[0.1,13.0], tau=[0.1,5.0], Zmetal=[0.004,0.03], $
            AV=[0.35,2.0], mu=[0.1,4.0], pburst=0.2, interval_pburst=2.0, $
            tburst=[0.1,13.0], /delayed, /nebular, galchunksize=2500L, clobber=clobber
       endif else begin
          write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
            prefix=prefix, filterlist=filterlist, spsmodels='fsps_v2.4_miles', $
            imf='chab', redcurve='charlot', /igm, zminmax=[0.05,0.75], zbin=0.02, $
            nmodel=25000L, age=[0.1,13.0], tau=[0.1,5.0], Zmetal=[0.004,0.03], $
            AV=[0.35,2.0], mu=[0.1,4.0], pburst=0.2, interval_pburst=2.0, $
            tburst=[0.1,13.0], /delayed, /nebular, galchunksize=2500L, clobber=clobber
       endelse
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
;  generate the model photometry QAplots
    if keyword_set(qaplot_models) then begin
       if keyword_set(sdss) then begin
          cat = mrdfits(isedfit_dir+'sdss_mgfe_parent.fits.gz',1)
          thesefilters = ['galex_NUV','sdss_g0','sdss_r0','sdss_i0','wise_w1']
       endif else begin
          cat = mrdfits(isedfit_dir+'ages_mgfe_parent.fits.gz',1)
;         thesefilters = 
       endelse
       isedfit_qaplot_models, isedfit_paramfile, cat.maggies, $
         cat.ivarmaggies, cat.z, isedfit_dir=isedfit_dir, $
         thissfhgrid=thissfhgrid, thesefilters=thesefilters, clobber=clobber
    endif
    
; --------------------------------------------------
; fit!
    if keyword_set(isedfit) then begin
       if keyword_set(sdss) then begin
          cat = mrdfits(isedfit_dir+'sdss_mgfe_parent.fits.gz',1,$
            columns=['ra','dec','z','maggies','ivarmaggies'])
          ngal = n_elements(cat)

          chunksize = ceil(ngal/float(nchunk))
          if n_elements(firstchunk) eq 0 then firstchunk = 0
          if n_elements(lastchunk) eq 0 then lastchunk = nchunk-1

          for ii = firstchunk, lastchunk do begin
             splog, 'Working on CHUNK '+strtrim(ii,2)+'/'+strtrim(lastchunk-firstchunk+1,2)
             splog, im_today()
             t0 = systime(1)
             outprefix = prefix+'_chunk'+string(ii,format='(I0)')
             these = lindgen(chunksize)+ii*chunksize
             these = these[where(these lt ngal)]
;            these = these[0:99] ; test!
             isedfit, isedfit_paramfile, cat[these].maggies, cat[these].ivarmaggies, $
               cat[these].z, ra=cat[these].ra, dec=cat[these].dec, isedfit_dir=isedfit_dir, $
               outprefix=outprefix, thissfhgrid=thissfhgrid, clobber=clobber
             splog, 'Total time (min) = '+strtrim((systime(1)-t0)/60.0,2)
          endfor                    
          
       endif else begin
          cat = mrdfits(isedfit_dir+'ages_mgfe_parent.fits.gz',1)
          isedfit, isedfit_paramfile, cat.maggies, cat.ivarmaggies, $
            cat.z, ra=cat.ra, dec=cat.dec, isedfit_dir=isedfit_dir, $
            thissfhgrid=thissfhgrid, clobber=clobber
; *copy* the output to MGFEPATH so that we can still make QAplots here 
          params = read_isedfit_paramfile(isedfit_paramfile,thissfhgrid=thissfhgrid)
          fp = isedfit_filepaths(params,isedfit_dir=isedfit_dir)
          if im_file_test(isedfit_dir+fp.isedfit_outfile+'*',clobber=clobber) eq 0 then begin
             file_copy, isedfit_dir+fp.isedfit_outfile+'.gz', $
               mgfepath+fp.isedfit_outfile, /overwrite
             file_copy, isedfit_dir+fp.post_outfile+'.gz', $
               mgfepath+fp.isedfit_outfile, /overwrite
          endif
       endelse
    endif 

; --------------------------------------------------
; merge the SDSS iSEDfit results but not the posteriors (too large)
; and then make a copy for MGFEPATH
    if keyword_set(sdss) and keyword_set(merge_isedfit) then begin
       params = read_isedfit_paramfile(isedfit_paramfile,thissfhgrid=thissfhgrid)
       fp = isedfit_filepaths(params,isedfit_dir=isedfit_dir)
       if im_file_test(isedfit_dir+fp.isedfit_outfile+'*',clobber=clobber) eq 0 then begin
          for ii = 0, nchunk-1 do begin
             outprefix = prefix+'_chunk'+string(ii,format='(I0)')
             ised1 = read_isedfit(isedfit_paramfile,isedfit_dir=isedfit_dir,$
               outprefix=outprefix)
             if ii eq 0 then ised = temporary(ised1) else $
               ised = [temporary(ised),temporary(ised1)]
          endfor
          im_mwrfits, ised, isedfit_dir+fp.isedfit_outfile, /clobber
          file_copy, isedfit_dir+fp.isedfit_outfile+'.gz', mgfepath, /overwrite
       endif
    endif

; --------------------------------------------------
;  compute K-corrections
    if keyword_set(kcorrect) then begin
       params = read_isedfit_paramfile(isedfit_paramfile,thissfhgrid=thissfhgrid)
       fp = isedfit_filepaths(params,isedfit_dir=isedfit_dir,band_shift=band_shift)
       if im_file_test(mgfepath+fp.kcorr_outfile+'*',clobber=clobber) eq 0 then begin
          isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
            montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
            absmag_filterlist=absmag_filterlist, band_shift=band_shift, $
            clobber=clobber
; copy the output to MGFEPATH
          file_copy, isedfit_dir+fp.kcorr_outfile+'.gz', mgfepath, /overwrite
       endif
    endif 

; --------------------------------------------------
;  generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
       isedfit_qaplot_sed, isedfit_paramfile, nrandom=20, $
         isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, $
         thissfhgrid=thissfhgrid, clobber=clobber, /xlog
    endif
    
return
end
