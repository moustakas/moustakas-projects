function edit_mike_ut080417, rawmike
; jm09jan07nyu - edit the MIKE structure for this night 
    
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

;   realgood = where(mike.setup eq 1 and mike.flg_anly eq 1)
;   struct_print, struct_trimtags(mike[realgood],sel=['frame',$
;     'img_root','obj','colbin','rowbin','type','slit','obj_id',$
;     'flg_anly','setup'])

; this was supposed to be a lamp (it's identified as a TFLT by
; the code), but it looks like the lamp was never turned on; reject 
    fixme = where(strtrim(mike.img_root,2) eq 'r0023.fits')
    mike[fixme].flg_anly = 0
    
; the internal quartz flats are generally misidentified as milky flats
; (which are the flats observed *with* the diffuser)
    lamp = where((strmatch(mike.obj,'*qz_lamp*',/fold) eq 1) and $
      (strmatch(mike.obj,'*diffuser*',/fold) eq 0))
    mike[lamp].type = 'TFLT'

; all the science spectra are mis-identified as zeros
    fixme = where((strtrim(mike.type,2) eq 'ZRO') or $
      (strmatch(mike.obj,'*obj*',/fold) and $
      (strmatch(mike.obj,'*lamp*',/fold) eq 0B)),nfixme)
    mike[fixme].type = 'OBJ'
    for jj = 0L, nfixme-1L do begin
       obj = strsplit(mike[fixme[jj]].obj,'_',/extract)
       obj[1] = string(obj[1],format='(I3.3)')
       mike[fixme[jj]].obj = strjoin(obj,'_')
    endfor
    
; only keep objects that we need    
    
    realgood = where(mike.setup eq 1 and mike.flg_anly eq 1)
    struct_print, struct_trimtags(mike[realgood],sel=['frame',$
      'img_root','obj','colbin','rowbin','type','slit','obj_id',$
      'flg_anly','setup'])

return, mike[realgood]
end
