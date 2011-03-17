pro kenn92_html_spectra, object
; jm03sep30uofa
; write ASCII, FITS spectra and tarballs for the HTML atlas and ASCII
; data tables  
    
; initialize path names
    
    webpath = atlas_path(/web)+'kenn92/'     ; web path
    kenn92path = atlas_path(/kenn92)+'data/' ; kenn92 path

    asciipath = webpath+'ascii/'
    fitspath = webpath+'fits/'

; restore all the fitting results

    kenn92 = read_kenn92(/silent)
    
; to examine just one object crop the list (this function should
; become obsolete when the web page is up and running)
    
    if n_elements(object) ne 0L then begin

       doit = match_string(object,kenn92.galaxy,index=match) ; match by galaxy name
       if match[0] eq -1L then begin
          splog, 'Object '+object+' not found!'
          return
       endif
       kenn92 = kenn92[match] 
       
    endif

    srt = sort(kenn92.galaxy) ; alphabetize
    kenn92 = kenn92[srt]
    
    galaxy = strupcase(strcompress(kenn92.galaxy,/remove))
    ngalaxy = n_elements(galaxy)

; default output file names
    
    asciiname = strlowcase(galaxy)+'.spec'
    fitsname = strlowcase(galaxy)+'.fits'
    fitstablename = strlowcase(galaxy)+'_table.fits'
    asciitablename = strlowcase(galaxy)+'.dat'

    for k = 0L, ngalaxy-1L do begin

       print, format='("Galaxy ",I3,"/",I3,".",A1,$)', k+1, ngalaxy, string(13b)

; ASCII spectrum

       specfile = strcompress(kenn92[k].specfile,/remove)
       scube = rd1dspec(specfile,datapath=kenn92path,/silent)
       openw, lun, asciipath+asciiname[k], /get_lun
       printf, lun, '# '+galaxy[k]
       printf, lun, '# All quantities are in the observed frame, corrected for '
       printf, lun, '# foreground Galactic extinction but not internal reddening.'
       printf, lun, '# '
       printf, lun, '# (1) Wavelength [A]; (2) Spectrum [erg/s/cm2/A]; '
       printf, lun, '# (3) Error spectrum [erg/s/cm2/A]; (4) Sky spectrum [erg/s/cm2/A]'
       printf, lun, '# '+im_today()
       printf, lun, '# '
       for j = 0L, scube.npix-1L do printf, lun, scube.wave[j], scube.spec[j], $
         scube.sigspec[j], scube.sky[j], format='(F7.2,3x,E0.0,3x,E0.0,3x,E0.0)'
       free_lun, lun

; FITS spectrum

       spawn, ['ln -fs '+kenn92path+specfile+' '+fitspath+fitsname[k]], /sh
;      spawn, ['gzip -fc '+kenn92path+kenn92[k].specfile+' > '+fitspath+fitsname[k]], /csh

; ASCII data table

;      openw, lun, asciipath+asciitablename[k], /get_lun
;      printf, lun, '# Not written yet!'
;      free_lun, lun

; FITS data table

       mwrfits, kenn92[k], fitspath+fitstablename[k], /create, /silent
       spawn, ['gzip -f '+fitspath+fitstablename[k]], /sh

    endfor    
    print

; generate and write tar-balls for the whole kenn92 atlas

    splog, 'Writing the FITS tarball.'
    pushd, kenn92path
    spawn, ['tar czvf kenn92_spectra_fits.tar.gz *.fits*'], /sh
    spawn, ['mv -f kenn92_spectra_fits.tar.gz '+webpath], /sh
    popd    
 
    splog, 'Writing the ASCII tarball.'
    pushd, webpath+'ascii'
    spawn, ['tar czf kenn92_spectra_ascii.tar.gz *.spec'], /sh
    spawn, ['mv -f kenn92_spectra_ascii.tar.gz '+webpath], /sh
    popd    
    
    splog, 'All done!'

return
end    
