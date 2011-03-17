pro trgb_pre_allframe, objname, verbose=verbose
;+
; NAME:
;	TRGB_PRE_ALLFRAME
;
; PURPOSE:
;	Iteratively runs ALLFRAME on a user-selected modified
;	(smaller) starlist to determine the appropriate PROFILE and
;	PERCENT ERROR parameters which give the best CHI
;	distribution.  
;
; INPUTS:
;	objname	: galaxy name (string)
;
; KEYWORD PARAMETERS:
;	verbose : verbose (to the screen) ALLFRAME output
;
; OUTPUTS:
;	Updates the ALLFRAME.OPT parameter file with the best PROFILE
;	and PERCENT values.
;
; COMMON BLOCKS:
;
; RESTRICTIONS:
;	The procedure should be run from the directory where the fits
;	files reside.
;
; EXAMPLE:
;	trgb_pre_allframe, 'SextansB', /verbose
;
; PROCEDURES USED:
;	RDTXT(), TRGB_MAKEMCH, TRGB_MAKEMINI, TRGB_CHECKCHI
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 May 25, UCB
;-

        npar = n_params()
        if npar ne 1 then begin
            print, 'Syntax: trgb_pre_allframe, objname, verbose=verbose'
            return
        endif

	spawn, ['pwd'], datapath	; current directory
        datapath = datapath[0]

        imnames = rdtxt(datapath+'/image.names') ; read in the image names
        nr = n_elements(imnames)
        iband = imnames[where(strpos(imnames,'_I') gt 0)] ; I-band frames
        ni = n_elements(iband)
        
        trgb_makemch, objname, /als ; create the als match files

; iterate on a miniature starlist until the PROFILE and PERCENT ERROR
; give satisfactory CHI values.

        counter = 0L
        happy = ' '
        while strupcase(happy) ne 'Y' do begin

            if counter eq 0L then begin
                print
                read, pro_err, prompt='Profile error (suggested: 3.0) : '
                read, per_err, prompt='Percent error (suggested: 0.55): '
                print
            endif else begin
                print
                print, 'Remember, RAISING (LOWERING) percent error brings chi DOWN (UP).'
                print
                print, 'Old profile error: ', pro_err 
                print, 'Old percent error: ', per_err
                print
                read, pro_err, prompt='New profile error: '
                read, per_err, prompt='New percent error: '
            endelse

            print
            print, 'Profile Error = ', pro_err
            print, 'Percent Error = ', per_err
            print
                        
            trgb_makemini, objname ; mini starlist

; remove any old ALLFRAME files
            spawn, ['\rm '+objname+'.bck']
            spawn, ['\rm '+objname+'*j.fits']
            spawn, ['\rm '+objname+'*k.fits']

; ALLFRAME script on the mini-starlist
            openw, lun1, 'allframe_mini_script', /get_lun
            printf, lun1, 'profile = '+strn(pro_err)
            printf, lun1, 'percent = '+strn(per_err)
            printf, lun1, ' '
            printf, lun1, objname+'.mch'
            printf, lun1, objname+'_mini.mag'
            free_lun, lun1
            
            print
            print, 'Running ALLFRAME on the mini starlist (this '
            print, 'will take about ten minutes for 500 stars).'
            if keyword_set(verbose) then spawn, ['allframe < allframe_mini_script'] else $
              spawn, ['allframe < allframe > allframe_mini_log'] ; this should take about 10 minutes

            trgb_makemch, objname, /alf, /mini	; update the transformation file

; check for a mini raw file

            spawn, ['find '+datapath+' -name "'+objname+'.raw" -print'], miniraw
            miniraw = miniraw[0]

; DAOMASTER script to create the raw files
            mnumber = nr-2 > 1 ; minimum number
            eframes = nr-1 > 1 ; enough frames

            openw, lun1, 'daomaster_script', /get_lun
            printf, lun1, objname+'_mini.mch'
            printf, lun1, strn(mnumber)+',0.5,'+strn(eframes)
            printf, lun1, '5.'	; maximum sigma
            printf, lun1, '6'	; number of transformation coefficients
            printf, lun1, '1'   ; critical matchup radius
            printf, lun1, '1'   
            printf, lun1, '1'   
            printf, lun1, '1'   
            printf, lun1, '1'   
            printf, lun1, '1'   
            printf, lun1, '0'   ; exit
            printf, lun1, 'n'	; no: new ids
            printf, lun1, 'n'	; no: mean mags & scatter
            printf, lun1, 'n'	; no: corrected mags & errors
            printf, lun1, 'y'	; yes: raw magnitudes & errors
            printf, lun1, objname+'_mini.raw'
            if miniraw ne ' ' then printf, lun1, ' ' ; overwrite 
            printf, lun1, 'n'	; no: new transformations
            printf, lun1, 'n'	; no: transfer table
            printf, lun1, 'n'	; no: individual coo files
            printf, lun1, 'n'	; no: transfer star ids
            free_lun, lun1

            if keyword_set(verbose) then spawn, ['daomaster < daomaster_script'] else $
              spawn, ['daomaster < daomaster > daomaster_log']

            ncol = 5+2*ni       ; number of columns in the raw file

            trgb_checkchi, objname, ncol, /mini ; check the chi values

            print
            print, 'Type Y if you are happy with this, '
            print, 'otherwise press any key to change '
            read, happy, prompt='the error parameters: '

            counter = counter+1L 
           
        endwhile
        
; update the allframe.opt file with the best profile and percent error

        allfopt = rdtxt(datapath+'/allframe.opt',nline=10)

        openw, lun2, datapath+'/allframe.opt', /get_lun
        for i = 0L, 7L do printf, lun2, allfopt[i]
        printf, lun2, 'PER = ', strn(per_err)
        printf, lun2, 'PRO = ', strn(pro_err)
        free_lun, lun2

return
end 
