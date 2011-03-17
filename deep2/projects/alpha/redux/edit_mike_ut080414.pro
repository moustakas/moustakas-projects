function edit_mike_ut080414, rawmike
; jm08apr13nyu - edit the MIKE structure for this night 

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

; these were water vapor observations for LCO    

    vapor = where(strmatch(mike.obj,'*water*',/fold))
    mike[vapor].flg_anly = 0

; twilight flats mis-identified as objects in the red side    
    twi = where(strmatch(mike.obj,'*twilight*',/fold))
    mike[twi].type = 'TWI' 

; rename the "internal flats" as "trace flats"; also, on the red side
; they were mis-identified as "milky flats"
    lamp = where((strmatch(mike.obj,'*qz_lamp*',/fold) eq 1) and $
      (strmatch(mike.obj,'*diffuser*',/fold) eq 0))
    mike[lamp].type = 'TFLT'

; object mis-identified as a zero
;   fixme = where(strtrim(mike.obj,2) eq '071_a')
;   mike[fixme].type = 'OBJ'

; because of the pre-processing, all the science spectra were
; mis-identified as zeros; also prepend "obj_" to the object name 
    fixme = where(strtrim(mike.type,2) eq 'ZRO')
    mike[fixme].type = 'OBJ'
    mike[fixme].obj = 'obj_'+strtrim(mike[fixme].obj,2)
    
; only keep objects that we need    
    
    realgood = where(mike.setup eq 1 and mike.flg_anly eq 1)
    struct_print, struct_trimtags(mike[realgood],sel=['frame',$
      'img_root','obj','colbin','rowbin','type','slit','obj_id',$
      'flg_anly','setup'])

return, mike[realgood]
end
