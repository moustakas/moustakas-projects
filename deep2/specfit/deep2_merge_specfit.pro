pro deep2_merge_specfit, specdata, test=test, write=write
; jm07sep28nyu - based on AGES_MERGE_SPECFIT

    catpath = deep2_path(/catalogs)
    specfitpath = deep2_path(/specfit,/dr4)

    if keyword_set(test) then prefix = 'test_' else prefix = ''
    specdatafiles = file_basename(file_search(specfitpath+'?????_????_'+prefix+$
      'specdata.fits.gz',count=fcount))

    allmaskname = strmid(specdatafiles,6,4)
    srt = sort(allmaskname)
    allmaskname = allmaskname[srt]
    specdatafiles = specdatafiles[srt]
    niceprint, allmaskname, specdatafiles

    mjdall = strmid(specdatafiles,0,5)
    mjd = string(max(mjdall),format='(I0)')

    specdata = mrdfits(specfitpath+specdatafiles[0],1,/silent)
    for i = 1L, fcount-1L do begin
       print, format='("Merging plate ",I0,"/",I0,".",A1,$)', i+1, fcount, string(13b)
       specdata1 = mrdfits(specfitpath+specdatafiles[i],1,/silent)
       specdata = [temporary(specdata),specdata1]
    endfor

    outfile = catpath+'speclinefit.dr4.goodspec1d.Q34.fits'
    if keyword_set(write) then begin
       splog, 'Writing '+outfile+'.'
       mwrfits, specdata, outfile, /create
       spawn, 'gzip -f '+outfile, /sh
    endif

return    
end
