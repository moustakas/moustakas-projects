function read_ages_specfit, galaxy, unfluxed=unfluxed, notweak=notweak, $
  test=test, obswave=obswave, silent=silent
; jm04jun05uofa
; jm04dec17uofa    
; jm05aug03uofa - updated for release 2.0 
; jm06feb20uofa - major re-write to handle the new way I'm treating
;                 the plates in AGES_SPECFIT
; jm08aug20nyu - another major update

    ngalaxy = n_elements(galaxy)
    if (ngalaxy eq 0L) then begin
       print, 'specfit = read_ages_specfit(galaxy,_extra=extra)'
       return, -1L
    endif

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
    
    pass = strmid(galaxy,5,3)
    upass = pass[uniq(pass,sort(pass))]
    npass = n_elements(upass)

; sort through the various SPECFILES and find the most recent one

    out_specfit = fltarr(4650,6,ngalaxy)
    
    for ipass = 0L, npass-1L do begin

       print, format='("Reading pass ",I0,"/",I0,".",A1,$)', $
         ipass+1, npass, string(13b)

       specdatafiles = file_basename(file_search(specfitpath+'?????_'+$
         upass[ipass]+prefix+'_specdata.fits.gz',count=fcount))

       if (fcount eq 0L) then begin
          splog, 'No matching files found in path '+specfitpath
          return, -1.0
       endif
       specdatafile = specdatafiles[(reverse(sort(specdatafiles)))[0]]
       specfitfile = repstr(specdatafile,'_specdata','_specfit')

       specdata = mrdfits(specfitpath+specdatafile,1,/silent)
       specfit = mrdfits(specfitpath+specfitfile,0,/silent)
       npix = (size(specfit,/dim))[0]

       these = where(upass[ipass] eq pass,nthese)
       indx = lonarr(nthese)
       for ii = 0L, nthese-1L do indx[ii] = where(galaxy[these[ii]] eq specdata.galaxy)
       flag = where(indx eq -1L,nflag)
       if (nflag ne 0L) then message, 'Problem here!'

       out_specfit[0:npix-1,*,these] = specfit[*,*,indx]
       if keyword_set(obswave) then begin
          for ii = 0L, nthese-1L do begin
             zplus1 = (1.0 + specdata[indx[ii]].z_abs)
             out_specfit[*,0,these[ii]] = out_specfit[*,0,these[ii]]*zplus1
             out_specfit[*,1,these[ii]] = out_specfit[*,1,these[ii]]/zplus1
             out_specfit[*,2,these[ii]] = out_specfit[*,2,these[ii]]/zplus1
             out_specfit[*,3,these[ii]] = out_specfit[*,3,these[ii]]/zplus1
             out_specfit[*,4,these[ii]] = out_specfit[*,4,these[ii]]/zplus1
             out_specfit[*,5,these[ii]] = out_specfit[*,5,these[ii]]*zplus1^2.0
          endfor
       endif
       
    endfor
    if (not keyword_set(silent)) then print
    
return, out_specfit
end
