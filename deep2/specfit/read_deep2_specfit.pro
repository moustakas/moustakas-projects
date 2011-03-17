function read_deep2_specfit, file, test=test, obswave=obswave, silent=silent
; jm07apr16nyu - read the best-fitting deep2 spectra; sort by mask 
; jm08sep04nyu - updated to the new data model

    nfile = n_elements(file)
    if (nfile eq 0L) then begin
       print, 'specfit = read_deep2_specfit(file,_extra=extra)'
       return, -1L
    endif

    version = deep2_version(/ispec)
    base_specfitpath = deep2_path(/specfit)
    specfitpath = base_specfitpath+version+'/'

    if keyword_set(test) then prefix = 'test_' else prefix = ''

    file = strtrim(file,2)
    mask = strmid(file,0,4)
    umask = mask[uniq(mask,sort(mask))]
    nmask = n_elements(umask)

    npix = 4096*2
    out_specfit = fltarr(npix,6,nfile)
    
    for imask = 0L, nmask-1L do begin

       if (not keyword_set(silent)) then $
         print, format='("Reading mask ",I0,"/",I0,".",A1,$)', $
         imask+1, nmask, string(13b)

       specdatafiles = file_basename(file_search(specfitpath+'?????_'+$
         umask[imask]+prefix+'_specdata.fits.gz',count=fcount))

       if (fcount eq 0L) then begin
          splog, 'No matching files found in path '+specfitpath
          return, -1.0
       endif
       specdatafile = specdatafiles[(reverse(sort(specdatafiles)))[0]]
       specfitfile = repstr(specdatafile,'_specdata','_specfit')

       specdata = mrdfits(specfitpath+specdatafile,1,/silent)
       specfit = mrdfits(specfitpath+specfitfile,0,/silent)
       npix = (size(specfit,/dim))[0]

       these = where(umask[imask] eq mask,nthese)
       indx = lonarr(nthese)
       for ii = 0L, nthese-1L do indx[ii] = where(file[these[ii]] eq strtrim(specdata.file,2))
       flag = where(indx eq -1L,nflag)
       if (nflag ne 0L) then message, 'Problem here!'

       out_specfit[0:npix-1,*,these] = specfit[*,*,indx]
       
       if keyword_set(obswave) then begin
          for ii = 0L, nthese-1L do begin
             zplus1 = (1.0 + specdata[indx[ii]].z_obj)
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

;;function read_deep2_specfit, file, mask, silent=silent
;;; jm07apr16nyu - read the best-fitting Deep2 spectra; sort by mask 
;;
;;    nfile = n_elements(file)
;;    if (nfile eq 0L) then begin
;;       print, 'specfit = read_deep2_specfit(file,_extra=extra)'
;;       return, -1L
;;    endif
;;
;;    datapath = deep2_path(/specfit)
;;
;;    file = strcompress(file,/remove)
;;    mask = strcompress(mask,/remove)
;;    umask = mask[uniq(mask,sort(mask))]
;;    nmask = n_elements(umask)
;;
;;; sort through the various FILES and find the most recent one
;;
;;    specfit = fltarr(2201,5,nfile)
;;    
;;    for imask = 0L, nmask-1L do begin
;;
;;       specdatafiles = file_basename(file_search(datapath+'?????_'+umask[imask]+'_specdata.fits.gz',count=fcount))
;;
;;       if (fcount eq 0L) then message, 'Problem here!'
;;       specdatafile = specdatafiles[(reverse(sort(specdatafiles)))[0]]
;;       specfitfile = repstr(specdatafile,'specdata','specfit')
;;
;;       print, datapath+specdatafile
;;       specdata = mrdfits(datapath+specdatafile,1,/silent)
;;
;;       these = where((umask[imask] eq mask),nthese)
;;       for ii = 0L, nthese-1L do begin
;;
;;          if (not keyword_set(silent)) then print, format='("Reading mask ",I0,"/",I0," and file ",I0,"/",I0,".    ",A1,$)', $
;;            imask+1, nmask, ii+1, nthese, string(13b)
;;
;;          indx = where(strtrim(file[these[ii]],2) eq strtrim(specdata.file,2),nindx)
;;          if (nindx eq 0L) then message, 'Problem here!'
;;
;;          ext_no = indx[0]+2
;;          specfit1 = mrdfits(datapath+specfitfile,ext_no,/silent)
;;          specfit[0L:n_elements(reform(specfit1[*,0]))-1L,*,these[ii]] = specfit1
;;
;;       endfor
;;
;;    endfor
;;    if (not keyword_set(silent)) then print
;;    
;;return, specfit
;;end
