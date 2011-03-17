pro unpack_ldss3, unpack=unpack
; jm06nov06nyu - written
; jm07jan15nyu - make this routine do the cross-talk correction

    originalpath = ldss3_path(/original)
    rawpath = ldss3_path(/feb06)+'raw/'

    info = ldss3_forage(file_search(originalpath+'*.fits'))
    nfits = n_elements(info)
    struct_print, info

    ignore = ['0128','200?','20[4-5]?','2060','2061',$ ; junk and/or spectroscopy
      '2062',$ ; mis-labeled header (bias, not a twiflat)
      '2063']  ; bright (saturated?) star

; /bin/rm ccd30[1-3][1-9]c?.fits

    for i = 0L, nfits-1L do begin

       file = strcompress(info[i].file)
       
       skip = 0
       for j = 0L, n_elements(ignore)-1L do if strmatch(file,'*'+ignore[j]+'*') then skip = 1

       if (skip eq 0L) then begin

          if (strmatch(file,'*c1*') eq 0L) and (strmatch(file,'*c2*') eq 0L) then message, 'Problem here!'

          if strmatch(file,'*c1*') then suffix = 'c1'
          if strmatch(file,'*c2*') then suffix = 'c2'

          object = strlowcase(info[i].object)
          filter = strcompress(strlowcase(info[i].filter),/remove)
          speed = strcompress(strlowcase(info[i].speed),/remove)

          prefix = strmid(file,3,4)
          if strmatch(object,'*bias*',/fold) or strmatch(object,'*dark*',/fold) then begin
             if strmatch(object,'*bias*',/fold) then outfile = 'a.'+prefix+'_'+object+'_'+speed+'_'+suffix+'.fits'
             if strmatch(object,'*dark*',/fold) then outfile = 'a.'+prefix+'_'+strjoin(strsplit(object,' ',/extract),'_')+'_'+$
               suffix+'.fits'
          endif else begin
             if strmatch(object,'*sdss*') then object = repstr(object,' 1','_1')
             if (strpos(object,' '))[0] ne -1L then object = strcompress(strmid(object,0,strpos(object,' ')),/remove)
             outfile = 'a.'+prefix+'_'+object+'_'+filter+'_'+suffix+'.fits'
          endelse
          
          if keyword_set(unpack) and (strmatch(outfile,'*slow*',/fold) eq 0L) then begin
;            if file_test(rawpath+outfile,/regular) then rmfile, rawpath+outfile
;            spawn, ['ln -s '+'../original/'+file+' '+rawpath+outfile], /sh
; do the cross-talk correction; blow away the BSCALE and BZERO header keywords
             if (suffix eq 'c1') then file2 = repstr(file,'c1','c2') else file2 = repstr(file,'c2','c1')
             im1 = readfits(originalpath+file,h1,/silent)
             im2 = readfits(originalpath+file2,h2,/silent)
             im1 = temporary(im1) - 0.000162*im2
             sxaddhist, "'Cross-talk correction applied "+hogg_iso_date()+"'", h1
             sxdelpar, h1, 'BSCALE'
             sxdelpar, h1, 'BZERO'
             sxdelpar, h1, 'COMMENT'
             sxaddpar, h1, 'FILTER', filter ; rename the FILTER parameter
; ---------------------------------------------------------------------------
; update the header with preliminary astrometry
             orientation = 270.0 ; NOTE!
             imsize = size(im1,/dim)
             ra = sxpar(h1,'RA',count=racount)
             dec = sxpar(h1,'DEC',count=deccount)
             pixscale = sxpar(h1,'SCALE',count=scalecount)
             if ((racount+deccount+scalecount) ne 3L) then begin
                splog, 'Unexpected LDSS3 header!'
                stop
             endif
             racen = 15.0*hms2dec(ra)
             deccen = hms2dec(dec)
             dra1 = imsize[0]*pixscale/3600.0
             ddec1 = imsize[1]*pixscale/3600.0
             astr = hogg_make_astr(racen,deccen,dra1,ddec1,$
               pixscale=pixscale/3600.0,orientation=orientation)
             putast, h1, astr
             h1 = h1[where((strmatch(strcompress(h1,/remove),'*HISTORY*PUTAST*',/fold) eq 0B))]
             sxaddhist, "'Initial WCS astrometry added "+hogg_iso_date()+"'", h1
; ---------------------------------------------------------------------------
             splog, 'Writing '+rawpath+outfile
;            print, outfile, outfile2
; reflect chip 1 to account for the way it was read out
;            if strmatch(outfile,'*c1*') then im1 = reverse(temporary(im1))
             writefits, rawpath+outfile, float(im1), h1
          endif else print, outfile
          
       endif
       
    endfor

return
end
    
