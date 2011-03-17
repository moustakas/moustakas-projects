; jm00may25ucb
pro checkfiles, objname, imnames, hst=hst
; routine to ensure all the correct files exist in this directory

        spawn, ['find '+objname+'.mch'], mchfile ; transformation file
        spawn, ['find '+objname+'.mag'], magfile ; starlist
        if keyword_set(hst) then $
          spawn, ['find *_*_*.psf'], psfiles else $ ; psfs 
          spawn, ['find *.psf'], psfiles
        spawn, ['find *.als'], alsfiles 	 ; als header files

        endit = ' ' 
        print
        print, 'The transformation file: ', mchfile
        print, 'The starlist           : ', magfile
        print
        print, 'The images             : ', imnames
        print, 'The psfs               : ', psfiles
        print, 'The als files          : ', alsfiles
        print
        print, '-------------------------'
        print
        read, endit, prompt='Type Q to quit or ENTER to go on: '
        if strupcase(endit) eq 'Q' then retall
        
return
end

pro trgb_allframe, objname, verbose=verbose, hst=hst, nocheck=nocheck
;+
; NAME:
;	TRGB_ALLFRAME
;
; PURPOSE:
;	Run ALLFRAME.
;
; INPUTS:
;	objname : galaxy name
;
; KEYWORD PARAMETERS:
;	verbose : verbose (to the screen) ALLFRAME output
;	hst	: keyword specifying HST data
;	nocheck	: do not check to make sure all the necessary ALLFRAME
;		  files are present
;
; OUTPUTS:
;
; COMMON BLOCKS:
;
; PROCEDURES USED:
;	RDTXT(), TRGB_MAKEALS, TRGB_MAKEMCH
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 May 26, UCB
;-

        npar = n_params()
        if npar ne 1 then begin
            print, 'Syntax: trgb_allframe, objname, verbose=verbose, hst=hst, nocheck=nocheck'
            return
        endif

	spawn, ['pwd'], datapath	; current directory
        datapath = datapath[0]

        imnames = rdtxt(datapath+'/image.names') ; read in the image names
        nr = n_elements(imnames)
        
        trgb_makeals, hst=hst	    ; create the als header files

; check that the right files are in datapath
        
        if not keyword_set(nocheck) then checkfiles, objname, imnames, hst=hst

; pre-ALLFRAME: checks for the right files and iterates on the profile
; and percent error parameters until the chi values are satisfactory.

        if not keyword_set(hst) then trgb_pre_allframe, objname, verbose=verbose

        trgb_makemch, objname, /als ; update the transformation file

; remove any old ALLFRAME files
        spawn, ['\rm '+objname+'.bck']
        spawn, ['\rm '+objname+'*j.fits']
        spawn, ['\rm '+objname+'*k.fits']

; ALLFRAME script
        openw, lun1, 'allframe_script', /get_lun
        printf, lun1, ' '
        printf, lun1, objname+'.mch'
        printf, lun1, objname+'.mag'
        free_lun, lun1
        
        print
        print, 'Running ALLFRAME.  This will take a while!'
        if keyword_set(verbose) then spawn, ['allframe < allframe_script'] else $
          spawn, ['allframe < allframe_script > allframe_log']

; create a raw file

;       trgb_makemch, objname, /alf ; update the transformation file
;       trgb_daomaster, objname, nr, verbose=verbose

; remove large, unnecessary ALLFRAME output files
;       spawn, ['\rm '+objname+'*j.fits']
        spawn, ['\rm '+objname+'*k.fits']

return
end 
