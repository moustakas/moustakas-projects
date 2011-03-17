pro atlas_doit, scube, specinfo, asciifile
    
    asciipath = atlas_path(/ascii)

    rastr = strtrim(string(15.0*hms2dec(specinfo.ra),format='(F12.5)'),2)
    decstr = strtrim(string(hms2dec(specinfo.dec),format='(F12.5)'),2)

;   openw, lun, asciipath+asciifile, /get_lun
;   printf, lun, '\ '+asciifile
;   printf, lun, '\ RA = '+rastr+' [degree, J2000]'
;   printf, lun, '\ DEC = '+decstr+' [degree, J2000]'
;   printf, lun, '\ APERTURE = '+strtrim(specinfo.strap,2)+' [sq. arcsec]'
;   printf, lun, '\ PA = '+strtrim(specinfo.pa,2)+' [degree, N-->E]'
;   printf, lun, '| wave  |   flux    |   error   |   sky     |'
;   printf, lun, '| real  |   real    |   real    |   real    |'

    openw, lun, asciipath+asciifile, /get_lun
    printf, lun, '# '+asciifile
;   printf, lun, '# ASCII spectrum generated '+im_today()
    printf, lun, '# RA = '+rastr+' [degree, J2000]'
    printf, lun, '# DEC = '+decstr+' [degree, J2000]'
    printf, lun, '# APERTURE = '+strtrim(specinfo.strap,2)+' [sq. arcsec]'
    printf, lun, '# POSANGLE = '+strtrim(specinfo.posangle,2)+' [degree, N-->E]'
    printf, lun, '# '
    printf, lun, '# 1 Wavelength [Angstrom]'
    printf, lun, '# 2 Spectrum   [erg/s/cm2/A]'
    printf, lun, '# 3 Error      [erg/s/cm2/A]'
    printf, lun, '# 4 Sky        [erg/s/cm2/A]'
    printf, lun, '# '

    for j = 0L, scube.npix-1L do printf, lun, scube.wave[j], scube.spec[j], $
      scube.sigspec[j], scube.sky[j], format='(F9.2,3x,E12.5,3x,E12.5,3x,E12.5)'
    free_lun, lun

return
end    

pro atlas1d_write_ascii
; jm05may03uofa
; write ASCII spectra

    datapath = atlas_path(/analysis)
    atlaspath = atlas_path(/atlas1d)

    atlas = mrdfits(datapath+'atlas1d_info.fits.gz',1,/silent)
    ngalaxy = n_elements(atlas)
    
    for k = 0L, ngalaxy-1L do begin

       print, format='("Galaxy ",I0,"/",I0,".",A1,$)', k+1, ngalaxy, string(13b)

; nuclear       

       if atlas[k].nuclear then begin

          asciifile = repstr(repstr(repstr(strtrim(atlas[k].nuclear_file,2),'.fits','.txt'),'.ms',''),'.gz','')
          scube = rd1dspec(atlas[k].nuclear_file,datapath=atlaspath,/silent)
          specinfo = im_struct_trimtags(atlas[k],select=['NUCLEAR_RA','NUCLEAR_DEC',$
            'NUCLEAR_STRAP','NUCLEAR_POSANGLE'],newtags=['RA','DEC','STRAP','POSANGLE'])

          specinfo = struct_trimtags(specinfo,select=tag_names(specinfo),$
            format=['A0','A0','A0','I0'])

          atlas_doit, scube, specinfo, asciifile

       endif
       
; integrated

       if atlas[k].drift then begin

          asciifile = repstr(repstr(repstr(strtrim(atlas[k].drift_file,2),'.fits','.txt'),'.ms',''),'.gz','')
          scube = rd1dspec(atlas[k].drift_file,datapath=atlaspath,/silent)
          specinfo = im_struct_trimtags(atlas[k],select=['DRIFT_RA','DRIFT_DEC',$
            'DRIFT_STRAP','DRIFT_POSANGLE'],newtags=['RA','DEC','STRAP','POSANGLE'])

          specinfo = struct_trimtags(specinfo,select=tag_names(specinfo),$
            format=['A0','A0','A0','I0'])

          atlas_doit, scube, specinfo, asciifile

       endif
       
    endfor    
    print

return
end    
