pro ediscs_merge_specfit, specdata, test=test, write=write
; jm06sep27nyu - merge all the individual EDISCS_SPECFIT results into
;                one data structure, that can then be parsed
; jm07apr24nyu - also merge all the SPECFIT results together
; jm09may22nyu - major rewrite

    version = ediscs_version(/specfit)

    basepath = ediscs_path(/specfit)
    specfitpath = basepath+version+'/'

    if keyword_set(test) then prefix = 'test_' else prefix = ''
    specdatafiles = file_basename(file_search(specfitpath+$
      '?????_ediscs_*'+prefix+'*specdata.fits.gz',count=fcount))
    specfitfiles = file_basename(file_search(specfitpath+$
      '?????_ediscs_*'+prefix+'*specfit.fits.gz',count=fcount))
    niceprint, specdatafiles, specfitfiles
    if (fcount gt 19L) then begin
       splog, 'Remove old SPECDATAFILES!!'
       return
    endif

    mjdall = strmid(specdatafiles,0,5)
    mjd = string(max(mjdall),format='(I0)')

    outfile = 'ediscs_specdata_'+version+'.fits'
    outfile_specfit = 'ediscs_specfit_'+version+'.fits'
;   outfile = mjd+'_ediscs_specdata.fits'

    specdata = mrdfits(specfitpath+specdatafiles[0],1,/silent)
    for i = 1L, fcount-1L do begin
       print, format='("Merging cluster ",I0,"/",I0,".",A1,$)', i+1, fcount, string(13b)
       specdata = struct_append(specdata,mrdfits(specfitpath+specdatafiles[i],1,/silent))
    endfor

;   specfit = mrdfits(specfitpath+specfitfiles[0],0,/silent)
    for i = 0, fcount-1 do begin
       print, format='("Merging cluster ",I0,"/",I0,".",A1,$)', i+1, fcount, string(13b)
       thisspecfit = mrdfits(specfitpath+specfitfiles[i],0,/silent)
       sz = size(thisspecfit,/dim)
;      print, sz[0]
       specfit1 = fltarr(2000,6,sz[2])
       specfit1[0:sz[0]-1,*,*] = thisspecfit
       if (i eq 0) then specfit = specfit1 else specfit = [[[specfit]],[[specfit1]]]
    endfor

    if keyword_set(write) then begin
       splog, 'Writing '+outfile
       im_mwrfits, specdata, basepath+outfile
       splog, 'Writing '+outfile_specfit
       im_mwrfits, specfit, basepath+outfile_specfit
    endif

return    
end
