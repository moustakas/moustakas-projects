pro trgb_daomaster, objname, nr, verbose=verbose
;+
; NAME:
;	TRGB_DAOMASTER
;
; PURPOSE:
;	Create a .raw ALLFRAME photometry file by applying the
;	transformation equations.
;
; INPUTS:
;	objname : galaxy name
;	nr	: number of galaxy images (eg, 4 I, 4 V: nr = 8)
;
; KEYWORD PARAMETERS:
;	verbose : verbose (to the screen) DAOMASTER output
;
; OUTPUTS:
;
; PROCEDURES USED:
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 15, UCB
;-

	spawn, ['pwd'], datapath	; current directory
        datapath = datapath[0]

; check for a raw file

        spawn, ['find '+datapath+' -name "'+objname+'.raw" -print'], rawfile
        rawfile = rawfile[0]

; DAOMASTER script to create the raw files
        
        mnumber = nr-2 > 1      ; minimum number
        eframes = nr-1 > 1      ; enough frames
        
        openw, lun1, 'daomaster_script', /get_lun
        printf, lun1, objname+'.mch'
        printf, lun1, strn(mnumber)+',0.5,'+strn(eframes)
        printf, lun1, '5.'	; maximum sigma
        printf, lun1, '6'	; number of transformation coefficients
        printf, lun1, '1'       ; critical matchup radius
        printf, lun1, '1'   
        printf, lun1, '1'   
        printf, lun1, '1'   
        printf, lun1, '1'   
        printf, lun1, '1'   
        printf, lun1, '1'   
        printf, lun1, '1'   
        printf, lun1, '1'   
        printf, lun1, '1'   
        printf, lun1, '1'   
        printf, lun1, '0'       ; exit
        printf, lun1, 'n'	; no: new ids
        printf, lun1, 'n'	; no: mean mags & scatter
        printf, lun1, 'n'	; no: corrected mags & errors
        printf, lun1, 'y'	; yes: raw magnitudes & errors
        printf, lun1, objname+'.raw'
        if rawfile ne '' then printf, lun1, ' ' ; overwrite 
        printf, lun1, 'n'	; no: new transformations
        printf, lun1, 'n'	; no: transfer table
        printf, lun1, 'n'	; no: individual coo files
        printf, lun1, 'n'	; no: transfer star ids
        free_lun, lun1
        
        if keyword_set(verbose) then spawn, ['daomaster < daomaster_script'] else $
          spawn, ['daomaster < daomaster_script > daomaster_log']

return
end
