pro gto_getdss, gto, object, dsswidth=dsswidth
; jm07oct23nyu 

; use the QUERYDSS Goddard routine to grab DSS images from STScI

    dsspath = gto_path(/dss)

    if (n_elements(gto) eq 0L) then gto = gto_read_ancillary()
    
    if (n_elements(object) ne 0L) then begin
       doit = match_string(object,gto.galaxy,findex=match,/exact) ; match by galaxy name
       if (match[0] eq -1L) then message, 'Object '+object+' not found!'
       gto1 = gto[match]
    endif else begin
       gto1 = gto
       if file_test(dsspath,/directory) ne 0L then begin
          splog, 'Would you like to remove all DSS files from '+dsspath+' [Y/N]?'
          cc = get_kbrd(1)
          if strupcase(cc) eq 'Y' then spawn, ['/bin/rm -f '+dsspath+'*.fits.gz'], /sh
       endif
    endelse
    ngto = n_elements(gto1)

    sizes = [2.5,5.0,10.0,15.0,20.0,25.0,30.0] ; possible image sizes [arcmin]
    
    for j = 0L, ngto-1L do begin ; loop on each object

       galaxy = strcompress(strlowcase(gto1[j].galaxy),/remove)
       nedgalaxy = strcompress(strlowcase(gto1[j].ned_galaxy),/remove)
       radec = [15.0*im_hms2dec(gto1[j].ra),im_hms2dec(gto1[j].dec)]

; compute the required FITS image size based on the size of the galaxy 

       if (gto[j].d25_maj lt -900.0) then imsize = sizes[1] else $
         imsize = ceil(2.0*sqrt(gto[j].d25_maj*gto[j].d25_min)) ; [arcmin]

       if (n_elements(dsswidth) ne 0L) then imsize = dsswidth

       splog, 'Querying '+strupcase(galaxy)+', ', gto[j].d25_maj, imsize
;      querydss, nedgalaxy, im, h, imsize=imsize, /ned; survey='2r'
       querydss, radec, im, h, imsize=imsize, /ned, /eso;, survey='2r'

       splog, 'Writing '+dsspath+galaxy+'.fits'
       writefits, dsspath+galaxy+'.fits', im, h
       spawn, 'gzip -f '+dsspath+galaxy+'.fits', /sh
       
    endfor

return
end
