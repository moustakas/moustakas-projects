; if the mean keyword is set then preliminary starlist stuff (like
; FIND and PHOT) can be run on the mean image, otherwise image.names
; is read in and run on all the images.

pro trgb_hst_pre_starlist, rootname, verbose=verbose, mean=mean, pc=pc
;+
; NAME:
;	TRGB_HST_PRE_STARLIST
;
; MODIFICATION HISTORY:
;	Written, J. Moustakas, 2000 June 7
;-

;	on_error, 2	; return to user
        npar = n_params()

        if npar ne 1L then begin
            print, 'Syntax: trgb_hst_pre_starlist, objname, verbose=verbose, mean=mean, pc=pc'
            return
        endif

;	spawn, ['clear']
	spawn, ['pwd'], datapath	; current directory
        datapath = datapath[0]

        chipnum = strmid(datapath,strpos(datapath,'CHIP')+4L) ; chip number

        imnames = rdtxt(datapath+'/image.names')
        iband = imnames[where(strpos(imnames,'_I') gt 0,ni)] ; I-band frames
        nr = n_elements(imnames)

        find_pars = fltarr(2)

        if keyword_set(mean) then begin

            objname = rootname+'_Im_'+chipnum ; mean I-band image
        
            find_pars[0] = ni ; number of frames averaged
            find_pars[1] = 1.00	  ; number of frames summed	

            nr = 1L

        endif else begin

            find_pars[0] = 1.
            find_pars[1] = 1.

        endelse

        trgb_hst_makeopt, pc=pc        ; create the .opt files

        for k = 0L, nr-1L do begin ; loop on the images

            if not keyword_set(mean) then objname = imnames[k]

; delete any old files in this directory

            spawn, ['\rm '+objname+'.coo*']
            spawn, ['\rm '+objname+'.ap*']

; FIND script
            openw, lun1, 'find_script', /get_lun
            printf, lun1, 'ATTACH '+objname+'.fits'
            printf, lun1, 'FIND'
            printf, lun1, strn(find_pars[0])+','+strn(find_pars[1])
            printf, lun1, objname+'.coo'
            printf, lun1, 'Y'   ; answer to "Are you happy with this?"
            printf, lun1, 'EXIT'
            free_lun, lun1
            
            if keyword_set(verbose) then spawn, ['daophot < find_script'] else $
               spawn, ['daophot < find_script > find_log']

; PHOT script
            openw, lun1, 'phot_script', /get_lun
            printf, lun1, 'ATTACH '+objname+'.fits' ; image
            printf, lun1, 'PHOT'
            printf, lun1, 'photo.opt' ; photometry file
            printf, lun1, ' '
            printf, lun1, objname+'.coo' ; coordinate file
            printf, lun1, objname+'.ap' ; output aperture photometry file
            printf, lun1, 'EXIT'
            free_lun, lun1
            
            if keyword_set(verbose) then spawn, ['daophot < phot_script'] else $
              spawn, ['daophot < phot_script > phot_log']

        endfor

; delete the junk files

        spawn, ['\rm *jnk.fits']
            
return
end


