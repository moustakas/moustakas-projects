; jm00june30ucb

pro multimch, verbose=verbose

	spawn, ['pwd'], datapath	; current directory
        datapath = datapath[0]

        imnames = rdtxt(datapath+'/image.names')
        nr = n_elements(imnames)

        trgb_hst_makeopt

        for k = 0L, nr-1L do begin ; loop on the images

            objname = imnames[k]

            spawn, ['\rm '+objname+'.coo*']
            spawn, ['\rm '+objname+'.ap*']

            openw, lun1, 'find_script', /get_lun
            printf, lun1, 'ATTACH '+objname+'.fits'
            printf, lun1, 'FIND'
            printf, lun1, '1,1'
            printf, lun1, objname+'.coo'
            printf, lun1, 'Y'   ; answer to "Are you happy with this?"
            printf, lun1, 'EXIT'
            free_lun, lun1
            
            if keyword_set(verbose) then spawn, ['daophot < find_script'] else $
               spawn, ['daophot < find_script > find_log']

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

        spawn, ['\rm *jnk.fits']
            
return
end


