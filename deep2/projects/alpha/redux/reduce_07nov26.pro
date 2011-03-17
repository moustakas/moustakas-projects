pro edit_mike_07nov26, rawmike, mike
; jm08apr12nyu - edit the mike structure for the 07nov26 data 

    date = '07nov26'
    if (n_elements(rawmike) eq 0L) then $
      rawmike = mrdfits('rawmike_'+date+'.fits',1)
    
    mike = rawmike
    mike.obj = strtrim(mike.obj,2)
    mike.obj = repstr(mike.obj,' ','_')

; only analyze objects with 1x1 binning    

;   good = where((mike.rowbin eq 1) and (mike.colbin eq 1),comp=bad)
    good = where((mike.rowbin eq 1) and (mike.colbin eq 1) and $
      (mike.slit eq 0.7),comp=bad)
    mike[good].setup = 1
    mike[bad].flg_anly = 0

; these objects (and corresponding arcs) were too faint/diffuse to be
; useful, so don't bother analyzing them; also the last object
; of the night is no good
    mike[102:114].flg_anly = 0
    mike[137:138].flg_anly = 0
    
; reject additional crap
    mike[33:42].flg_anly = 0 ; wrong slit width
    mike[93:94].flg_anly = 0 ; something weird with the power supply
    mike[95].flg_anly = 0    ; don't use external lamp
    mike[96].flg_anly = 0    ; probably ok, but science arcs should do
    mike[104].flg_anly = 0   ; from log: "wrong source"

; could possibly use the sky flats as the milky flats (MFLT)
    mike[101].type = 'TWI' ; 'MFLT'
    mike[144:161].type = 'TWI' ; 'MFLT'
    
; milky flats - wrong/different slit, but the binning is right    
    mike[83:92].flg_anly = 1
    mike[83:92].setup = 1
   
;; mis-identified as MFLT flats; actually internal flats
     mike[53:62].type = 'TFLT'
;
;; mis-identified twilight flats (apparently not used)
;    mike[101].type = 'TWI'
;    mike[144:161].type = 'TWI'

; mis-identified as zero's; actually object images
    mike[[106,108,111,113]].type = 'OBJ'

; don't use darks
    mike[139:143].flg_anly = 0

; we want individual exposures of the same object to be reduced and
; extracted independently; to do this, we have to trick MIKE_SETUP by
; changing the object names (otherwise it will group them together and
; give them common obj_id values)    

    suffix = '_'+['a','b','c','d','e']
    
    obj = where(strtrim(mike.obj,2) eq '003',nobj)
    mike[obj].obj = mike[obj].obj+suffix[0L:nobj-1L]
    obj = where(strtrim(mike.obj,2) eq '004',nobj)
    mike[obj].obj = mike[obj].obj+suffix[0L:nobj-1L]
    obj = where(strtrim(mike.obj,2) eq '028',nobj)
    mike[obj].obj = mike[obj].obj+suffix[0L:nobj-1L]
    obj = where(strtrim(mike.obj,2) eq '040',nobj)
    mike[obj].obj = mike[obj].obj+suffix[0L:nobj-1L]
    
; only write out the objects that we need    
    
;   realgood = where(mike.setup eq 1 and mike.flg_anly eq 1 and 
;     mike.slit eq 0.7)
    realgood = where(mike.setup eq 1 and mike.flg_anly eq 1)
    struct_print, struct_trimtags(mike[realgood],sel=['frame',$
      'obj','colbin','rowbin','type','slit','obj_id','flg_anly','setup'])
    ;; Output
    mwrfits, mike[realgood], 'mike_07nov26.fits', /create

return
end


pro reduce_07nov26

; pre-processing

    date = '07nov26'

    mike_strct, rawmike, /noedit, outfil='rawmike_'+date+'.fits'
    edit_mike_07nov26
    mike = mike_ar()
    mike_setup, mike, outfil='mike_summ_'+date+'.txt'
    mike_wrstrct, mike, outfil='mike_'+date+'.list', fits='mike_'+date+'.fits' ; overwrite!

; process the flats and optionally verify

    mike_allflat, mike, 1, 2, /clobber
;   xatv, 'Flats/Flat_R_01_T.fits'
;   mike_chktrcflat, mike, 1, 2, /nostop, /fit

; process the arcs
    mike_allarc, mike, 1, 2, fits='mike_'+date+'.fits', /clobber, chk=chk ; overwrite!

; construct the slit profile
    mike_slitflat, mike, 1, 2, /clobber, chk=chk

; now extract object spectra!
    obj_id = 1
    mike_allobj, mike, 1, obj_id, 2, /clobber, /procall, $
      /nocr, /docomb, /do1d, checkall=checkall, /noflux, $
      /nohelio, /novac
    
    mike_box, mike, 1, obj_id, 2, ordrs=66, /ochk, /chk, /reschk, $
      base_aper=[0.20,0.20]

return
end

    
