pro trgb_hst_threshold, rootname, pc=pc, postscript=postscript
;+
; NAME:
;	TRGB_HST_THRESHOLD
;
; PURPOSE:
;	Generate a plot of threshold versus number of stars found.
;
; MODIFICATION HISTORY:
;	Written, J. Moustakas, 2000 June 10
;-

        npar = n_params()
        if npar ne 1L then begin
            print, 'Syntax: trgb_hst_threshold, rootname, pc=pc, postscript=postscript'
            return
        endif
        
        spawn, ['pwd'], datapath  ; current directory
        datapath = datapath[0]
        
;       chipnum = strmid(datapath,strpos(datapath,'CHIP')+4L) ; chip number
;       objname = rootname+'_'+chipnum+'m' ; median image

        imnames = rdtxt('image.names')
        objname = imnames[0]	; first image

        trgb_hst_makeopt, pc=pc        ; create the .opt files
        
        threshold = findgen(30)+1.
        nstars = threshold-threshold
        
; delete any old files in this directory
            
        spawn, ['\rm '+objname+'.coo*']
            
        openw, lun1, 'find_thresh_script', /get_lun
        printf, lun1, 'OPTIONS'
        printf, lun1, ' '
        printf, lun1, 'THRESH = 1'
        printf, lun1, ' '
        printf, lun1, 'ATTACH '+objname+'.fits'
        printf, lun1, 'FIND'
        printf, lun1, '1,1'
        printf, lun1, objname+'.coo'
        printf, lun1, 'N'       ; answer to "Are you happy with this?"
        for k = 1L, n_elements(threshold)-1L do begin
            printf, lun1, strn(threshold[k])
            printf, lun1, objname+'.coo'
            printf, lun1, ' '
            if k eq n_elements(threshold)-1L then $
              printf, lun1, 'Y' else printf, lun1, 'N'
        endfor
        printf, lun1, 'EXIT'
        free_lun, lun1
            
        spawn, ['daophot < find_thresh_script > find_thresh_log']
            
        spawn, ['grep stars find_thresh_log'], logline
            
        if n_elements(logline) ne n_elements(threshold) then return

        for k = 0L, n_elements(logline)-1L do begin
            stars = strmid(logline[k],0,rstrpos(logline[k],' stars.'))
            nstars[k] = float(stars)
        endfor

        spawn, ['\rm '+objname+'.coo*']
        spawn, ['\rm *jnk.fits']
            
; plot the diagnostic diagram
        
        if keyword_set(postscript) then ps_open, 'thresh_vs_num' else window, 0, xs=450, ys=450
        plot, threshold, nstars, ps=2, xsty=3, ysty=3, $
          xtit='Threshold', ytit='Number of stars', tit=strupcase(objname), $
          thick=2, charthick=1.5, charsize=1.5, xthick=2, ythick=2
        if keyword_set(postscript) then ps_close

return
end







