;+
; NAME:
;       SINGS_COMBINE_REPEATERS
;
; PURPOSE:
;       Combine repeat observations.
;
; CALLING SEQUENCE:
;       sings_combine_repeaters, /wfits
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jul 22, U of A
;-

pro sings_combine_repeaters, wfits=wfits
    
    analysis_path = sings_path(/analysis)
    outpath = sings_path(/spec1d);+'repeaters/temp/'
    repeatpath = sings_path(/spec1d)+'repeaters/'

; read the text file    
    
    repeatersfile = 'sings_combine_repeaters.txt'
    if (file_test(analysis_path+repeatersfile,/regular) eq 0L) then begin
       splog, 'File '+analysis_path+repeatersfile+' not found.'
       return
    endif

    repeaters = djs_readlines(analysis_path+repeatersfile)
    keep = where((strmatch(repeaters,'*#*') eq 0B) and (strcompress(repeaters,/remove) ne ''))
    repeaters = repeaters[keep]
    ngalaxy = n_elements(repeaters)

; loop on each unique object    
    
    for igalaxy = 0L, ngalaxy-1L do begin

; parse the text file

       line = strsplit(repeaters[igalaxy],' ',/extract)
       nrepeat = long(line[0]) ; number of repeat observations to be combined
       repeatlist = strtrim(line[1L:1L+nrepeat-1L],2)
       outfile = strtrim(line[1L+nrepeat],2)
       refnumber = long(line[2L+nrepeat])

; combine the repeaters; this code is lifted from
; SINGS_PLOT_REPEATERS

       s1 = rd1dspec(repeatlist[0],/silent,datapath=repeatpath)

       if (nrepeat eq 1L) then begin
       
          if keyword_set(wfits) then begin
             splog, 'Writing '+outpath+outfile+'.'
             wrt1dspec, outfile, s1.spec, s1.sigspec, s1.sky, $
               s1.mask, s1.header, datapath=outpath
          endif 

          continue
          
       endif

       s2 = rd1dspec(repeatlist[1],/silent,datapath=repeatpath)
       if (nrepeat eq 3L) then s3 = rd1dspec(repeatlist[2],/silent,datapath=repeatpath)

; interpolate both spectra onto a common wavelength vector to prevent
; extrapolation

       dwave = s1.wave[1]-s1.wave[0] ; assume a common linear dispersion

       minwave = min(s1.wave)>min(s2.wave) 
       maxwave = max(s1.wave)<max(s2.wave)
       if (nrepeat eq 3L) then begin
          minwave = minwave>min(s3.wave)
          maxwave = maxwave<max(s3.wave)
       endif

       finalwave = findgen((maxwave-minwave)/dwave+1L)*dwave+minwave
       newloglam = alog10(finalwave)
       nfinalpix = n_elements(finalwave)

       dlogwave = newloglam[1]-newloglam[0]
;      binsz = min((newloglam-shift(newloglam,1))[1L:nfinalpix-1L])
       binsz = dlogwave

; interpolate the spectra       

       newflux1 = interpol(s1.spec,s1.wave,finalwave)
       newivar1 = interpol(1.0/s1.sigspec^2.0,s1.wave,finalwave)
       newsky1 = interpol(s1.sky,s1.wave,finalwave)
       
       newflux2 = interpol(s2.spec,s2.wave,finalwave)
       newivar2 = interpol(1.0/s2.sigspec^2.0,s2.wave,finalwave)
       newsky2 = interpol(s2.sky,s2.wave,finalwave)

       if (nrepeat eq 3L) then begin
          newflux3 = interpol(s3.spec,s3.wave,finalwave)
          newivar3 = interpol(1.0/s3.sigspec^2.0,s3.wave,finalwave)
          newsky3 = interpol(s3.sky,s3.wave,finalwave)       
       endif

; compute the normalization constants       
       
       junk = im_normalize(newflux1,finalwave,normwave=5500.0,binsize=50.0,const=s1norm)
       junk = im_normalize(newflux2,finalwave,normwave=5500.0,binsize=50.0,const=s2norm)

       if (nrepeat eq 3L) then begin
          junk = im_normalize(newflux3,finalwave,normwave=5500.0,binsize=50.0,const=s3norm)
          maxnorm = max([s1norm,s2norm,s3norm],normindx)
;         print, repeatlist[0], s1norm, s2norm, s3norm
       endif else begin
          maxnorm = max([s1norm,s2norm],normindx)
;         print, repeatlist[0], s1norm, s2norm
       endelse

       objflux = [ [newflux1*(maxnorm/s1norm)], [newflux2*(maxnorm/s2norm)] ]
       objivar = [ [newivar1*(s1norm/maxnorm)^2], [newivar2*(s2norm/maxnorm)^2] ]
       skyflux = [ [newsky1*(maxnorm/s1norm)], [newsky2*(maxnorm/s1norm)] ]

       if (nrepeat eq 3L) then begin
          objflux = [ [objflux], [newflux3*(maxnorm/s3norm)] ]
          objivar = [ [objivar], [newivar3*(s3norm/maxnorm)^2] ]
          skyflux = [ [skyflux], [newsky3*(maxnorm/s3norm)] ]
       endif
          
; now combine the two spectra with appropriate variance weighting

       finalflux = total(objivar*objflux,2)/total(objivar,2)
       finalivar = total(objivar,2)
       finalsky = total(objivar*skyflux,2)/total(objivar,2)

       finalsigspec = 1.0/sqrt(finalivar)
       finalmedsnr = median(finalflux/finalsigspec)
       nfinalpix = n_elements(finalflux)
       finalmask = bytarr(nfinalpix) ; <-- NOT CORRECT!

; finalize the output FITS header

       case refnumber of
          0L: outheader = s1.header
          1L: outheader = s2.header
          2L: outheader = s3.header
          else: begin
             splog, 'REFNUMBER out of range.'
             return
          endelse
       endcase
       
       sxaddpar, outheader, 'NAXIS1', nfinalpix
       sxaddpar, outheader, 'CRVAL1', float(minwave), ' wavelength at CRPIX1', before='HISTORY'
       sxaddpar, outheader, 'CRPIX1', float(1.0), ' reference pixel number', before='HISTORY'
       sxaddpar, outheader, 'CD1_1', float(dwave), ' dispersion [Angstrom/pixel]', before='HISTORY'
       sxaddpar, outheader, 'CDELT1', float(dwave), ' dispersion [Angstrom/pixel]', before='HISTORY'
       sxaddpar, outheader, 'CTYPE1', 'LINEAR', ' projection type', before='HISTORY'
       sxaddhist, "'Weighted average of "+string(nrepeat,format='(I0)')+$
         " repeat observations "+im_today()+"'", outheader
          
; write out

       if keyword_set(wfits) then begin
          splog, 'Writing '+outpath+outfile+'.'
          wrt1dspec, outfile, finalflux, finalsigspec, finalsky, $
            finalmask, outheader, datapath=outpath
       endif 
       
    endfor 

return
end
    
