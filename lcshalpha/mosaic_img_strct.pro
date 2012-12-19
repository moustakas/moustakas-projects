;+ 
; NAME:
;   MOSAIC_IMG_STRCT
;
; PURPOSE:
;    Creates and outputs a structure for a series of astronomy camp
;    images.  This structure organizes the data for the night and is
;    used  to run most of the programs in the MOSAIC imaging reduction
;    package.
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;   struct     -  IDL structure 
;
; OPTIONAL KEYWORDS:
;   LIST       - Image list:  e.g.  'gd_files.lst'
;              Default is 'Raw/*.fits'
;   MKDIR      - Make directories
;   OUTFIL     - Name of fits output file
;   
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mosaic_img_strct, nght1_strct, /MKDIR
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   11-Jul-2002 Written by JXP
;   29-Jan-2003 Polished by JXP
;-

pro mosaic_img_strct, struct, ccd=ccd, tel=tel, LIST=list, $
  NOMKDIR=nomkdir, NOFILE=nofile, OUTFIL=outfil, IMG=img, $
  clobber=clobber

;  Optional Keywords
    if n_elements(OUTFIL) eq 0 then outfil = 'imgstrct.fits'

    if file_test(outfil) and (keyword_set(clobber) eq 0) then begin
       splog, 'Summary structure '+outfil+' exists'
       if arg_present(struct) then begin
          splog, 'Reading '+outfil
          struct = mrdfits(outfil,1)
       endif
       return
    endif

; List
    if not keyword_set( IMG ) then begin
       if not keyword_set( LIST ) then begin
          img = findfile('Raw/*.fits*',count=nimg) 
          if nimg EQ 0 then begin
             print, 'mosaic_img_strct: No images in Raw!'
             stop
             return
          endif
       endif else begin
          if x_chkfil(list[0]) NE 1 then begin
             print, 'mosaic_img_strct: Trouble with filename ', list
             return
          endif
          readcol, list, img, FORMAT='A'
          nimg = n_elements(img)
       endelse
    endif else nimg = n_elements(img)


; Make directories
    if not keyword_set( NOMKDIR ) then begin
       a = findfile('OV/..', count=count)
       if count EQ 0 then file_mkdir, 'OV'
       a = findfile('Final/..', count=count)
       if count EQ 0 then file_mkdir, 'Final'
       a = findfile('Flats/..', count=count)
       if count EQ 0 then file_mkdir, 'Flats'
       a = findfile('Bias/..', count=count)
       if count EQ 0 then file_mkdir, 'Bias'
       a = findfile('Std/..', count=count)
       if count EQ 0 then file_mkdir, 'Std'
    endif

;  Header Keywords
    expcrd = 'EXPTIME'
    utcrd = 'DATE-OBS'
    racrd = 'RA'
    deccrd = 'DEC'
    eqxcrd = 'EQUINOX'
    objcrd = 'OBJECT'
    amcrd = 'AIRMASS'
    datecrd = 'DATE-OBS'
    gaincrd = 'GAIN'
    readnocrd = 'RDNOISE'
    bincrd = 'CCDSUM'
    filtcrd = 'FILTER'
    
;  Create the Structure
    tmp = {mosaic_imgstrct}
    struct = replicate(tmp,nimg) 
; Default all files to be analysed
    struct.flg_anly = 1
    struct.tel = tel
    struct.ccd = ccd

;  Loop on Indv images
    for q = 0, nimg-1 do begin
;      splog, 'Reading ', img[q]
       head = xheadfits(img[q], /silent,ext=1) ; NOTE!
       
;  Parse the Header
       struct[q].naxis1 = sxpar(head, 'NAXIS1')
       struct[q].naxis2 = sxpar(head, 'NAXIS2')
       if keyword_set(expcrd) then struct[q].exp = sxpar(head,expcrd)
       if keyword_set(gaincrd) then struct[q].gain = sxpar(head,gaincrd)
       if keyword_set(readnocrd) then struct[q].readno = sxpar(head,readnocrd)
       if keyword_set(eqxcrd) then struct[q].equinox = sxpar(head,eqxcrd)
       if keyword_set(utcrd) then struct[q].UT = strtrim(sxpar(head,utcrd),2)
       if keyword_set(racrd) then struct[q].RA = 15D*hms2dec(strtrim(sxpar(head,racrd),2))
       if keyword_set(deccrd) then struct[q].DEC = hms2dec(strtrim(sxpar(head,deccrd),2))
       if keyword_set(filtcrd) then struct[q].filter = strtrim(sxpar(head,filtcrd),2)
       
; name
       if keyword_set(objcrd) then begin
          struct[q].Obj = sxpar(head,objcrd)
          struct[q].Obj = repchr(struct[q].Obj, '"',' ')
          struct[q].Obj = strcompress(struct[q].Obj, /remove_all)
       endif
       
; set image names and root path
       struct[q].img_root = repstr(file_basename(img[q]),'.gz','')
       struct[q].rootpth = file_dirname(img[q])+'/'

       if keyword_set(amcrd) then struct[q].am = sxpar(head,amcrd)
       if keyword_set(datecrd) then begin
          date = strtrim(sxpar(head,datecrd),2)
          if strlen(strtrim(date,2)) LT 8 then begin
             splog, 'No date in header.  Telescope not talking to you..'
             splog, 'Hopefully this is a calibration file'
             splog, 'Continue if it is or fix the header otherwise'
             stop
             date = '1999-01-01'
             struct[q].UT = '00:00:00.1'
          endif
;         splog, 'FIX UT PROBLEM!!!'
;         struct[q].date = x_setjdate(date,struct[q].ut)
       endif

; binning
       if keyword_set(bincrd) then begin
          binc = strcompress(sxpar(head, bincrd),/rem)
          struct[q].cbin = long(strmid(binc, 0, 1))
          struct[q].rbin = long(strmid(binc, 1, 1))
       endif

; get the pixel scale and frame number
       struct[q].pixscale = 0.6D*(struct[q].cbin>1) ; [arcsec/pixel]
       struct[q].frame = strmid(repstr(struct[q].img_root,'.fits',''),3,/reverse)
stop
       
             struct[q].pixscale = 0.45D*(struct[q].cbin>1) ; [arcsec/pixel]

; image type
       case ccd of
          'S2KB' : begin
             imtype = strtrim(sxpar(head, 'IMAGETYP'),2)
             case imtype of 
                'dark' : struct[q].type = 'DRK'   ; Dark
                'flat' : struct[q].type = 'DFT'  ; Dome Flat
                'zero' : struct[q].type = 'ZRO'   ; Bias
                else : begin
                   struct[q].type = 'OBJ' ; Object
                end 
             endcase  
          end 
          'SBIG' : begin
             imtype = strtrim(sxpar(head, 'IMAGETYP'),2)
             case imtype of 
                'dark' : struct[q].type = 'DRK'   ; Dark
                'flat' : struct[q].type = 'DFT'  ; Dome Flat
                'zero' : struct[q].type = 'ZRO'   ; Bias
                else : begin
                   struct[q].type = 'OBJ' ; Object
                end 
             endcase 
          end 
       endcase 

; catch zero frame
       if (struct[q].exp LT 0.1) then struct[q].type = 'ZRO' ; Bias
    endfor 

;; try to find standard-star fields; this may be buggy!
;    stetfile = getenv('MOSAIC_DIR')+'/data/phstd_stetson.tfits'
;    if (file_test(stetfile) eq 0) then begin
;       splog, 'Stetson photometric catalog not found! Check that ${MOSAIC_DIR} is set!'
;    endif else begin
;       splog, 'Reading '+stetfile
;       stet = mrdfits(stetfile,1)
;;      grp = spheregroup(stet.ra,stet.dec,0.5)
;       spherematch, stet.ra, stet.dec, struct.ra, $
;         struct.dec, 0.5, m1, m2
;       if (m2[0] ne -1) then struct[m2].type = 'STD'
;    endelse

; write out a convenient summary file
    if not keyword_set( NOFILE ) then begin
       sel = ['img_root','obj','type','exp',$
         'am','filter','naxis1','naxis2']
       struct_print, struct_trimtags(struct,select=sel), $
         file=repstr(outfil,'.fits','.txt')
    endif

; Write the structure to FITS
    splog, 'Writing '+outfil
    mwrfits, struct, outfil, /create

return
end
