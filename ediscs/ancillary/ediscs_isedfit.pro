pro ediscs_isedfit, models=models, isedfit=isedfit, qaplot=qaplot, $
  measure=measure, clobber=clobber, index=index, merge=merge, $
  debug=debug
; jm07mar25nyu
    
; ssh -X evolve    
; echo "ediscs_isedfit, /isedfit, nmonte=0L" | idl > & ltest.log &
; echo "ediscs_isedfit, /isedfit, nmonte=0L, /maxold" | idl > & ltest_maxold.log &

    phot = read_ediscs(/phot)
    spec1d = read_ediscs(/spec1d)

    cluster = strtrim(spec1d.cluster,2)
    allcluster = cluster[uniq(cluster,sort(cluster))]
    ncluster = n_elements(allcluster)
    ngal = n_elements(phot)

    iopath = ediscs_path(/isedfit)
    paramfile = iopath+'ediscs_isedfit.par'

; --------------------------------------------------
; build the model grids    
    if keyword_set(models) then isedfit_models, $
      paramfile, iopath=iopath, clobber=clobber

; --------------------------------------------------
; do the fitting    
    if keyword_set(isedfit) then begin
       for ic = 0, ncluster-1 do begin
          splog, '##################################################'
          splog, 'Fitting cluster '+allcluster[ic]
          outprefix = allcluster[ic];+'_BVRIJK'
          these = where(allcluster[ic] eq cluster,ngal)
;         index = where((spec1d[these].galaxy eq 'EDCSNJ1037539-1247248') or $
;           (spec1d[these].galaxy eq 'EDCSNJ1037532-1247270'))
          ediscs_to_maggies, phot[these], maggies, ivarmaggies
          isedfit, paramfile, maggies, ivarmaggies, spec1d[these].z, $
            result, iopath=iopath, outprefix=outprefix, nminphot=nminphot, $
            clobber=clobber, debug=debug, index=index
       endfor          
    endif

; --------------------------------------------------
; build a QAplot for each cluster    
    if keyword_set(qaplot) then begin
;      for ic = 1, 1 do begin
       for ic = 0, ncluster-1 do begin
          outprefix = allcluster[ic];+'_BVRIJK'
          these = where(allcluster[ic] eq cluster,ngal)
          isedfit_qaplot, paramfile, iopath=iopath, outprefix=outprefix, $
            galaxy=spec1d[these].galaxy, clobber=clobber, index=index
       endfor
    endif

; --------------------------------------------------
; measure rest-frame quantities 
    if keyword_set(measure) then begin
;      for ic = 1, 1 do begin
       for ic = 0, ncluster-1 do begin
          outprefix = allcluster[ic];+'_BVRIJK'
          these = where(allcluster[ic] eq cluster,ngal)
          isedfit_measure, paramfile, measure, isedfit, iopath=iopath, $
            index=index, abmag=abmag, clobber=clobber, outprefix=outprefix
       endfor
    endif

; --------------------------------------------------
; merge all the results together    
    if keyword_set(merge) then begin
       params1 = read_isedfit_paramfile(paramfile)
       nsfhgrid = n_elements(params1.sfhgrid)
       if (nsfhgrid ge 1) then begin
          for ii = 0, nsfhgrid-1 do begin
             params = struct_trimtags(params1,except='sfhgrid')
             params = struct_addtags(params,{sfhgrid: params1.sfhgrid[ii]})
             fp = isedfit_filepaths(params,iopath=iopath)
;            isedfiles = file_search(iopath+'*_'+$
;              fp.isedfit_outfile+'.gz',count=nised)
             isedfiles = file_search(iopath+'cl*sfhgrid'+$
               string(params.sfhgrid,format='(I2.2)')+'.fits.gz',count=nised)

;            measurefiles = file_search(iopath+'*_'+$
;              fp.measure_outfile+'.gz')
             for jj = 0, nised-1 do begin
                print, format='("Merging cluster ",I0,"/",I0,".",A1,$)', $
                  jj+1, nised, string(13b)
                ised1 = mrdfits(isedfiles[jj],1,/silent)
;               measure1 = mrdfits(measurefiles[jj],1,/silent)
                if (jj eq 0) then ised = ised1 else ised = [ised,ised1]
;               if (jj eq 0) then measure = measure1 else $
;                 measure = [measure,measure1]
             endfor
             ised.isedfit_id = lindgen(ngal)
             outfile = iopath+fp.isedfit_outfile
             im_mwrfits, ised, outfile, clobber=clobber

;            measure.isedfit_id = ised.isedfit_id
;            outfile = iopath+fp.measure_outfile
;            im_mwrfits, measure, outfile, clobber=clobber
          endfor
       endif
    endif
    
return
end
