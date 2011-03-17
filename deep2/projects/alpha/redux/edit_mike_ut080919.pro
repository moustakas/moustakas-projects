function edit_mike_ut080919, rawmike
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
      (strmatch(mike.obj,'*diffuser*',/fold) eq 0))
    mike[lamp].type = 'TFLT'

; only keep one of each of the standard-star observations; also get
; rid of some of the corresponding lamps
    toss = where($
      (strtrim(mike.img_root,2) eq 'r0022.fits') or $
      (strtrim(mike.img_root,2) eq 'r0023.fits') or $
      (strtrim(mike.img_root,2) eq 'r0026.fits') or $
      (strtrim(mike.img_root,2) eq 'r0027.fits'))
    mike[toss].flg_anly = 0
    
; fix the object names individually
   mike[where(strtrim(mike.img_root,2) eq 'r0031.fits')].obj = 'obj_152_a'
   mike[where(strtrim(mike.img_root,2) eq 'r0033.fits')].obj = 'obj_152_b'
   mike[where(strtrim(mike.img_root,2) eq 'r0036.fits')].obj = 'obj_133_a'
   mike[where(strtrim(mike.img_root,2) eq 'r0038.fits')].obj = 'obj_133_b'
   mike[where(strtrim(mike.img_root,2) eq 'r0040.fits')].obj = 'obj_133_c'
   mike[where(strtrim(mike.img_root,2) eq 'r0042.fits')].obj = 'obj_133_d'
   mike[where(strtrim(mike.img_root,2) eq 'r0045.fits')].obj = 'obj_136_a'
   mike[where(strtrim(mike.img_root,2) eq 'r0048.fits')].obj = 'obj_136_b'
   mike[where(strtrim(mike.img_root,2) eq 'r0050.fits')].obj = 'obj_136_c'
   mike[where(strtrim(mike.img_root,2) eq 'r0052.fits')].obj = 'obj_136_d'
   mike[where(strtrim(mike.img_root,2) eq 'r0055.fits')].obj = 'obj_147_a'
   mike[where(strtrim(mike.img_root,2) eq 'r0058.fits')].obj = 'obj_150_a'
   mike[where(strtrim(mike.img_root,2) eq 'r0060.fits')].obj = 'obj_150_b'
   mike[where(strtrim(mike.img_root,2) eq 'r0063.fits')].obj = 'obj_131_a'
   mike[where(strtrim(mike.img_root,2) eq 'r0065.fits')].obj = 'obj_131_b'
    
; toss out some crap observations that aren't worth reducing, as
; well as the corresponding lamps
    crap = where($
      (strtrim(mike.img_root,2) eq 'r0057.fits') or $
      (strtrim(mike.img_root,2) eq 'r0058.fits'))
    mike[crap].flg_anly = 0

; only keep objects that we need    
    realgood = where(mike.setup eq 1 and mike.flg_anly eq 1)
    struct_print, struct_trimtags(mike[realgood],sel=['frame',$
      'img_root','obj','colbin','rowbin','type','slit','obj_id',$
      'flg_anly','setup'])

return, mike[realgood]
end
