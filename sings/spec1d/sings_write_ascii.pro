pro sings_doit, scube, specinfo, asciifile
; the headers are IRSA compliant
    
    asciipath = sings_path(/ascii)

    rastr = strtrim(string(15.0*hms2dec(specinfo.ra),format='(F12.5)'),2)
    decstr = strtrim(string(hms2dec(specinfo.dec),format='(F12.5)'),2)

;\ ngc_1566_nuclear_002.5.ascii
;\ RA = 65.01000 [degree]
;\ DEC = -54.93894 [degree]
;\ Aperture = 2.5 x 2.5 [sq. arcsec]
;\ Pa = 60 [degree, N-->E]
;| wave  |     flux    |     error     |      sky     |
;|  real |     real    |      real     |     real     |
; 3548.85    2.66390E-15    1.36783E-16    1.39966E-16

    openw, lun, asciipath+asciifile, /get_lun
    printf, lun, '\ '+asciifile
    printf, lun, '\ RA = '+rastr+' [degree, J2000]'
    printf, lun, '\ DEC = '+decstr+' [degree, J2000]'
    printf, lun, '\ APERTURE = '+strtrim(specinfo.strap,2)+' [sq. arcsec]'
    printf, lun, '\ POSANGLE = '+strtrim(specinfo.posangle,2)+' [degree, N-->E]'
    printf, lun, '\ WAVE = observed wavelength [Angstrom]'
    printf, lun, '\ FLUX = observed spectrum [erg/s/cm^2/Angstrom]'
    printf, lun, '\ ERROR = one-sigma error spectrum [erg/s/cm^2/Angstrom]'
    printf, lun, '\ SKY = sky spectrum [erg/s/cm^2/Angstrom]'
    printf, lun, '\'
    printf, lun, '|  wave   |     flux     |     error    |     sky      |'
    printf, lun, '|  real   |     real     |     real     |     real     |'

;   openw, lun, asciipath+asciifile, /get_lun
;   printf, lun, '## '+asciifile
;   printf, lun, '## ASCII spectrum generated '+im_today()
;   printf, lun, '## '
;   printf, lun, '## Slit configuration:'
;   printf, lun, '## '
;   printf, lun, '## RA = '+rastr+' [degree]'
;   printf, lun, '## DEC = '+decstr+' [degree]'
;   printf, lun, '## Aperture = '+strtrim(specinfo.strap,2)+' [sq. arcsec]'
;   printf, lun, '## Pa = '+strtrim(specinfo.pa,2)+' [degree, N-->E]'
;   printf, lun, '## '
;   printf, lun, '## All quantities are in the observed frame, corrected for '
;   printf, lun, '## foreground Galactic extinction but not internal reddening.'
;   printf, lun, '## '
;   printf, lun, '# 1 Wavelength [Angstrom]'
;   printf, lun, '# 2 Spectrum   [erg/s/cm2/A]'
;   printf, lun, '# 3 Error      [erg/s/cm2/A]'
;   printf, lun, '# 4 Sky        [erg/s/cm2/A]'
;   printf, lun, '## '

    for j = 0L, scube.npix-1L do printf, lun, scube.wave[j], scube.spec[j], $
      scube.sigspec[j], scube.sky[j], format='(F9.2,3x,E12.5,3x,E12.5,3x,E12.5)'
    free_lun, lun

return
end    

pro sings_write_ascii
; jm04sep27uofa
; jm05apr08uofa - updated
; write ASCII spectra for the SINGS one-dimensional spectra

    datapath = sings_path(/analysis)
    spec1dpath = sings_path(/spec1d)

    sings = sings_read_info()
    ngalaxy = n_elements(sings)
    
    for k = 0L, ngalaxy-1L do begin

       print, format='("Galaxy ",I0,"/",I0,".",A1,$)', k+1, ngalaxy, string(13b)
;      match, strlowcase(strtrim(sings[k].galaxy,2)), strlowcase(strtrim(sings.galaxy,2)), indx1, indx2

; nuclear       

       if sings[k].nuclear then begin

          asciifile = repstr(repstr(repstr(strtrim(sings[k].nuclear_file,2),'.fits','.txt'),'.ms',''),'.gz','')
          scube = rd1dspec(sings[k].nuclear_file,datapath=spec1dpath,/silent)
          specinfo = im_struct_trimtags(sings[k],select=['NUCLEAR_RA','NUCLEAR_DEC',$
            'NUCLEAR_STRAP','NUCLEAR_POSANGLE'],newtags=['RA','DEC','STRAP','POSANGLE'])

          specinfo = struct_trimtags(specinfo,select=tag_names(specinfo),$
            format=['A0','A0','A0','I0'])

          sings_doit, scube, specinfo, asciifile

       endif
       
; drift20

       if sings[k].drift20 then begin

          asciifile = repstr(repstr(repstr(strtrim(sings[k].drift20_file,2),'.fits','.txt'),'.ms',''),'.gz','')
          scube = rd1dspec(sings[k].drift20_file,datapath=spec1dpath,/silent)
          specinfo = im_struct_trimtags(sings[k],select=['DRIFT20_RA','DRIFT20_DEC',$
            'DRIFT20_STRAP','DRIFT20_POSANGLE'],newtags=['RA','DEC','STRAP','POSANGLE'])

          specinfo = struct_trimtags(specinfo,select=tag_names(specinfo),$
            format=['A0','A0','A0','I0'])

          sings_doit, scube, specinfo, asciifile

       endif
       
; drift56

       if sings[k].drift56 then begin

          asciifile = repstr(repstr(repstr(strtrim(sings[k].drift56_file,2),'.fits','.txt'),'.ms',''),'.gz','')
          scube = rd1dspec(sings[k].drift56_file,datapath=spec1dpath,/silent)
          specinfo = im_struct_trimtags(sings[k],select=['DRIFT56_RA','DRIFT56_DEC',$
            'DRIFT56_STRAP','DRIFT56_POSANGLE'],newtags=['RA','DEC','STRAP','POSANGLE'])

          specinfo = struct_trimtags(specinfo,select=tag_names(specinfo),$
            format=['A0','A0','A0','I0'])

          sings_doit, scube, specinfo, asciifile

       endif
       
    endfor    
    print

return
end    
