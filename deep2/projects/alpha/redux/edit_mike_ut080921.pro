function edit_mike_ut080921, rawmike
; jm09jan16nyu - edit the MIKE structure for this night 
    
    mike = rawmike
    mike.obj = strtrim(mike.obj,2)
    mike.obj = repstr(mike.obj,' ','_')

    mike.type = strtrim(mike.type,2)

; the slit-width was not properly written to the header
    mike.slit = 0.7

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

; the internal quartz flats are generally misidentified as milky flats
; (which are the flats observed *with* the diffuser)
    lamp = where((strmatch(mike.obj,'*qz_flat*',/fold) eq 1) and $
      (strmatch(mike.obj,'*with_diffuser*',/fold) eq 0))
    mike[lamp].type = 'TFLT'

; only keep one of each of the standard-star observations
    toss = where($
      (strtrim(mike.img_root,2) eq 'r0022.fits') or $
      (strtrim(mike.img_root,2) eq 'r0025.fits'))
    mike[toss].flg_anly = 0

; fix the object names individually
;   realgood = where(mike.setup eq 1 and mike.flg_anly eq 1 and $
;     (strmatch(mike.obj,'*obj*') eq 1) and (strmatch(mike.obj,'*thar*') eq 0))
;   struct_print, struct_trimtags(mike[realgood],sel=['frame',$
;     'img_root','obj','colbin','rowbin','type','slit','obj_id',$
;     'flg_anly','setup'])
   mike[where(strtrim(mike.img_root,2) eq 'r0029.fits')].obj = 'obj_077_a'
   mike[where(strtrim(mike.img_root,2) eq 'r0031.fits')].obj = 'obj_077_b'
   mike[where(strtrim(mike.img_root,2) eq 'r0034.fits')].obj = 'obj_152_a'
   mike[where(strtrim(mike.img_root,2) eq 'r0037.fits')].obj = 'obj_160_a'
   mike[where(strtrim(mike.img_root,2) eq 'r0039.fits')].obj = 'obj_160_b'
   mike[where(strtrim(mike.img_root,2) eq 'r0041.fits')].obj = 'obj_160_c'
   mike[where(strtrim(mike.img_root,2) eq 'r0044.fits')].obj = 'obj_138_a'
   mike[where(strtrim(mike.img_root,2) eq 'r0046.fits')].obj = 'obj_138_b'
   mike[where(strtrim(mike.img_root,2) eq 'r0049.fits')].obj = 'obj_142_a'
   mike[where(strtrim(mike.img_root,2) eq 'r0051.fits')].obj = 'obj_142_b'
   mike[where(strtrim(mike.img_root,2) eq 'r0054.fits')].obj = 'obj_129_a'
   mike[where(strtrim(mike.img_root,2) eq 'r0056.fits')].obj = 'obj_129_b'
   mike[where(strtrim(mike.img_root,2) eq 'r0058.fits')].obj = 'obj_129_c'
   mike[where(strtrim(mike.img_root,2) eq 'r0061.fits')].obj = 'obj_162_a'
   mike[where(strtrim(mike.img_root,2) eq 'r0063.fits')].obj = 'obj_162_b'
   mike[where(strtrim(mike.img_root,2) eq 'r0065.fits')].obj = 'obj_162_c'
    
; toss out some crap observations that aren't worth reducing, as
; well as the corresponding lamps
    crap = where($
      (strtrim(mike.img_root,2) eq 'r0048.fits') or $
      (strtrim(mike.img_root,2) eq 'r0049.fits') or $
      (strtrim(mike.img_root,2) eq 'r0050.fits') or $
      (strtrim(mike.img_root,2) eq 'r0051.fits') or $
      (strtrim(mike.img_root,2) eq 'r0052.fits'))
    mike[crap].flg_anly = 0

; some of the science spectra are mis-identified as zeros
    fixme = where((strtrim(mike.type,2) eq 'ZRO'))
    mike[fixme].type = 'OBJ'

; only keep objects that we need    
    
    realgood = where(mike.setup eq 1 and mike.flg_anly eq 1)
    struct_print, struct_trimtags(mike[realgood],sel=['frame',$
      'img_root','obj','colbin','rowbin','type','slit','obj_id',$
      'flg_anly','setup'])

return, mike[realgood]
end
