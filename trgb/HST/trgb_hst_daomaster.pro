; jm00june15ucb

; script to match stars based on the ap files (should be done
; manually, but i'm lazy)

pro trgb_hst_daomaster, objname

        spawn, ['pwd'], datapath
        datapath = datapath[0]

        imnames = rdtxt(datapath+'/image.names')        
        nr = n_elements(imnames)

; DAOMASTER script
        
        mnumber = nr-2 > 1      ; minimum number
        eframes = nr-1 > 1      ; enough frames
        
        openw, lun1, 'daomaster_script', /get_lun
        printf, lun1, objname+'.mch'
        printf, lun1, strn(mnumber)+',0.5,'+strn(eframes)
        printf, lun1, '5.'	; maximum sigma
        printf, lun1, '6'	; number of transformation coefficients
        printf, lun1, '5'       ; critical matchup radius
        printf, lun1, '5'
        printf, lun1, '5'
        printf, lun1, '4'
        printf, lun1, '4'
        printf, lun1, '4'
        printf, lun1, '3'
        printf, lun1, '3'
        printf, lun1, '3'
        printf, lun1, '2'
        printf, lun1, '2'
        printf, lun1, '2'
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
        printf, lun1, 'n'	; no: raw magnitudes & errors
        printf, lun1, 'y'	; yes: new transformations
        printf, lun1, objname+'.mch'
        printf, lun1, ' '	; overwrite
        printf, lun1, 'n'	; no: transfer table
        printf, lun1, 'n'	; no: individual coo files
        printf, lun1, 'n'	; no: transfer star ids
        free_lun, lun1
        
        spawn, ['daomaster < daomaster_script > daomaster_log']

return
end
