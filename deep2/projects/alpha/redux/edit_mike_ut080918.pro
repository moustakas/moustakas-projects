function edit_mike_ut080918, rawmike
; jm09jan16nyu - edit the MIKE structure for this night 
    
    mike = rawmike
    mike.obj = strtrim(mike.obj,2)
    mike.obj = repstr(mike.obj,' ','_')

    mike.type = strtrim(mike.type,2)
    
; only analyze objects with 1x1 binning and the 0.7" slit; also
; only analyze the red side for now!

    good = where((mike.rowbin eq 1) and (mike.colbin eq 1) and $
      (mike.slit eq 0.7) and (mike.side ne 0),comp=bad,ncomp=nbad)
    mike[good].setup = 1
    if (nbad ne 0L) then mike[bad].flg_anly = 0

    realgood = where(mike.setup eq 1 and mike.flg_anly eq 1)
    struct_print, struct_trimtags(mike[realgood],sel=['frame',$
      'img_root','obj','colbin','rowbin','type','slit','obj_id',$
      'flg_anly','setup'])

; pointing problem & clouds - do not reduce
    crap = where($
      (strtrim(mike.img_root,2) eq 'r0022.fits') or $
      (strtrim(mike.img_root,2) eq 'r0023.fits') or $
      (strtrim(mike.img_root,2) eq 'r0024.fits') or $
      (strtrim(mike.img_root,2) eq 'r0025.fits') or $
      (strtrim(mike.img_root,2) eq 'r0026.fits'))
    mike[crap].flg_anly = 0

; object not in slit; do not reduce; also toss out the first
; lamp
    crap = where($
      (strtrim(mike.img_root,2) eq 'r0027.fits') or $
      (strtrim(mike.img_root,2) eq 'r0028.fits') or $
      (strtrim(mike.img_root,2) eq 'r0029.fits') or $
      (strtrim(mike.img_root,2) eq 'r0032.fits') or $
      (strtrim(mike.img_root,2) eq 'r0033.fits'))
    mike[crap].flg_anly = 0
    
; the internal quartz flats are generally misidentified as milky flats
; (which are the flats observed *with* the diffuser)
    lamp = where((strmatch(mike.obj,'*qz_flats_no_diffuser*',/fold) eq 1))
    mike[lamp].type = 'TFLT'

; all the science spectra are mis-identified as zeros; also fix the
; object names to be compatible with the matching
    fixme = where((strtrim(mike.type,2) eq 'ZRO') or $
      (strmatch(mike.obj,'*obj*',/fold) and $
      ((strmatch(mike.obj,'*th*',/fold) eq 0B) or $
      (strmatch(mike.obj,'*ar*',/fold) eq 0B))),nfixme)
    mike[fixme].type = 'OBJ'
    for jj = 0L, nfixme-1L do begin
       obj = strtrim(mike[fixme[jj]].obj,2)
       len = strlen(obj)
       mike[fixme[jj]].obj = 'obj_'+strmid(obj,3,3)+'_'+strmid(obj,0,/reverse)
    endfor
    
; only keep objects that we need    
    
    realgood = where(mike.setup eq 1 and mike.flg_anly eq 1)
    struct_print, struct_trimtags(mike[realgood],sel=['frame',$
      'img_root','obj','colbin','rowbin','type','slit','obj_id',$
      'flg_anly','setup'])

return, mike[realgood]
end
