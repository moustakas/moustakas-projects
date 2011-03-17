pro atlas_getdss, object, atlas=atlas, dsswidth=dsswidth
; jm03apr29uofa
; jm05may16uofa - updated
; jm05jul25uofa - updated

; use the QUERYDSS Goddard routine to grab DSS images from STScI

    dsspath = atlas_path(/dss)
    atlaspath = atlas_path(/analysis)

    if (n_elements(atlas) eq 0L) then atlas = atlas_read_info()
    
    if (n_elements(object) ne 0L) then begin
       doit = match_string(object,atlas.galaxy,findex=match,/exact) ; match by galaxy name
       if (match[0] eq -1L) then message, 'Object '+object+' not found!'
       atlas1 = atlas[match]
    endif else begin
       atlas1 = atlas
       if file_test(dsspath,/directory) ne 0L then begin
          splog, 'Would you like to remove all DSS files from '+dsspath+' [Y/N]?'
          cc = get_kbrd(1)
          if strupcase(cc) eq 'Y' then spawn, ['/bin/rm -f '+dsspath+'*.fits.gz'], /sh
       endif
    endelse

    natlas = n_elements(atlas1)

    sizes = [2.5,5.0,10.0,15.0,20.0,25.0,30.0] ; possible image sizes [arcmin]
    
    pushd, dsspath
    for j = 0L, natlas-1L do begin ; loop on each object

       galaxy = strcompress(strlowcase(atlas1[j].galaxy),/remove)
       nedgalaxy = strcompress(strlowcase(atlas1[j].ned_galaxy),/remove)
       radec = [15.0*im_hms2dec(atlas1[j].ra),im_hms2dec(atlas1[j].dec)]

;      case galaxy of
;         'mrk0960': nedgalaxy = [15.0,1.0]*im_hms2dec(['00:48:35.440','-12:42:59.9'])
;         'ic1623b': nedgalaxy = [15.0,1.0]*im_hms2dec(['01:07:48.179','-17:30:23.89'])
;         'ugc05028': nedgalaxy = [15.0,1.0]*im_hms2dec(['09:27:50.13','+68:24:42.7'])
;         'ugca225': nedgalaxy = [15.0,1.0]*im_hms2dec(['11:04:58.34','+29:08:17.5'])
;         'ugc06541': nedgalaxy = [15.0,1.0]*im_hms2dec(['11:33:29.12','+49:14:17.4'])
;         'arp244': nedgalaxy = [15.0,1.0]*im_hms2dec(['12:01:52.48','-18:52:02.9'])
;         'ngc5954': nedgalaxy = [15.0,1.0]*im_hms2dec(['15:34:35.07','+15:12:00.4'])
;         'ngc7580': nedgalaxy = [15.0,1.0]*im_hms2dec(['23 17 36.61','+14 00 03.0'])
;         'ngc7592a': nedgalaxy = [15.0,1.0]*im_hms2dec(['23:18:21.73','-04:24:57.4'])
;      else: 
;      endcase

; compute the required FITS image size based on the size of the galaxy 

       if (atlas[j].d25_maj lt -900.0) then imsize = sizes[1] else $
         imsize = ceil(1.4*atlas[j].d25_maj) ; [arcmin]

;      dmaj = ceil(atlas1[j].d25_maj) ; [arcmin]
;      if (dmaj lt -900.0) then imsize = sizes[1]
;      if (dmaj gt   0.0) and (dmaj lt  2.0) then imsize = sizes[0]
;      if (dmaj ge   2.0) and (dmaj lt  5.0) then imsize = sizes[1]
;      if (dmaj ge   5.0) and (dmaj lt 10.0) then imsize = sizes[2]
;      if (dmaj ge  10.0) and (dmaj lt 15.0) then imsize = sizes[3]
;      if (dmaj ge  15.0) and (dmaj lt 20.0) then imsize = sizes[4]
;      if (dmaj ge  20.0) and (dmaj lt 25.0) then imsize = sizes[5]
;      if (dmaj ge  25.0) and (dmaj lt 30.0) then imsize = sizes[6]

       if (n_elements(dsswidth) ne 0L) then imsize = dsswidth

       splog, 'Querying '+strupcase(galaxy)+', ', atlas[j].d25_maj, imsize
;      stop
;      querydss, nedgalaxy, im, h, imsize=imsize, /ned; survey='2r'
       querydss, radec, im, h, imsize=imsize, /ned, /eso;, survey='2r'

       writefits, galaxy+'.fits', im, h
       spawn, ['gzip -f '+galaxy+'.fits'], /sh
       
    endfor
    popd

return
end
