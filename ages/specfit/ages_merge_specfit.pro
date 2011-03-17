pro ages_merge_specfit, specdata, unfluxed=unfluxed, notweak=notweak, $
  test=test, write=write
; jm06feb17uofa - merge all the individual AGES_SPECFIT results into
;                 one data structure, that can then be parsed
; jm08aug18nyu - updated to the latest data model

    base_specfitpath = ages_path(/specfit)
    base_spec1dpath = ages_path(/spec1d)

    if keyword_set(unfluxed) then begin
       version = ages_version(/unfluxed_specfit)
       specfitpath = base_specfitpath+'unfluxed/'+version+'/'
       spec1dpath = base_spec1dpath+'unfluxed/after_skysubpca/'
       suffix = '_unfluxed'
    endif else begin
       version = ages_version(/ispec_specfit)
       if keyword_set(notweak) then begin
          specfitpath = base_specfitpath+'fluxed/notweak/'+version+'/'
          spec1dpath = base_spec1dpath+'fluxed/after_skysubpca/'
          suffix = '_ispec_notweak'
       endif else begin
          specfitpath = base_specfitpath+'fluxed/tweak/'+version+'/'
          spec1dpath = base_spec1dpath+'fluxed/tweak/'
          suffix = '_ispec_tweak'
       endelse
    endelse 

    if keyword_set(test) then prefix = 'test_' else prefix = ''
    specdatafiles = file_basename(file_search(specfitpath+$
      '?????_???_'+prefix+'specdata.fits.gz',count=fcount))
    
    passfield = strmid(specdatafiles,6,3) ; strmid(specdatafiles,11,3)
    srt = sort(passfield)
    passfield = passfield[srt]
    specdatafiles = specdatafiles[srt]
;   refitbadfiles = repstr(specdatafiles,'specdata.fits.gz',$ ; NOTE!
;     'refitbad_specdata.fits.gz') 

    specpassfiles = file_basename(file_search(spec1dpath+$
      'ages_'+passfield+'.fits.gz'))
    niceprint, passfield, specdatafiles, specpassfiles

    mjdall = strmid(specdatafiles,0,5)
    mjd = string(max(mjdall),format='(I0)')

    for ii = 0L, fcount-1L do begin
       print, format='("Merging plate ",I0,"/",I0,".",A1,$)', ii+1, $
         fcount, string(13b)
       specdata1 = mrdfits(specfitpath+specdatafiles[ii],1,/silent)
       refitbadfile = '?????_'+passfield[ii]+'_refitbad_specdata.fits.gz'
       if file_test(refitbadfile,/regular) then begin
          refitbad = mrdfits(specfitpath+refitbadfile,1,/silent)
          these = lonarr(n_elements(refitbad))
          for jj = 0L, n_elements(refitbad)-1L do these[jj] = $
            where(refitbad[jj].ages_id eq specdata1.ages_id)
; copy the structure over; generally, REFITBAD is fit just with
; the solar-metallicity templates, so we have to do this trick to get
; the TEMPLATEFILE array right
          junk = specdata1[these]
          struct_assign, refitbad, junk, /nozero
          specdata1[these] = junk
       endif
;      print, specdata1[0].pass
;      if (n_elements(specdata1) ne n_elements(where(strmid(ii.galaxy,5,3) eq $
;        passfield[ii]))) then stop
       if (ii eq 0L) then specdata = specdata1 else $
         specdata = [temporary(specdata),specdata1]
    endfor

    outfile = base_specfitpath+'ages_'+prefix+'specdata'+$
      suffix+'_'+version+'.fits'

    if keyword_set(write) then begin
       splog, 'Writing '+outfile
       mwrfits, specdata, outfile, /create
       spawn, 'gzip -f '+outfile, /sh
    endif

return    
end
