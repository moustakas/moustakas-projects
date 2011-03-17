pro trgb_starlist, objname, verbose=verbose, hst=hst, pc=pc
;+
; NAME:
;	TRGB_STARLIST
;
; PURPOSE:
;	Iteratively generate a master object starlist with ALLSTAR.
;
; INPUTS:
;	The file master.dir, which refers to the directory where the
;	data are kept, must exist in the current directory.  The 
;	input objname is the name of the object used to create the
;	master starlist (for example, a median image).
;
; OUTPUTS:
;	Creates a file called objname.mag in the master.dir directory:
;	this is the master starlist.  The files created during the
;	iterative process of making the starlist are written to the
;	SUBDIRECTORY of master.dir called STARLIST.
;
; KEYWORDS:
;	verbose : Set to show the DAOPHOT output to the screen
;	       	  (recommended).  
;
; PROCEDURE:
;	First any old files from a previous run of this code are
;	deleted.  An environment variable and alias is created to 
;	the specified image.  ALLSTAR is run, creating the first
;	starlist.  Then the user is prompted to search for the next
;	level of stars.  At this point DAOPHOT is called, then ALLSTAR
;	again, and the procedure iterates until the user is satisfied
;	with the starlist.
;
; PROCEDURES USED:
;	TRGB_MAKEOPT, TRGB_HST_MAKEOPT(), RDTXT()

; WARNINGS:
;	This program must be run from an *empty* subdirectory of
;	master.dir called STARLIST (all capitals).  Also, in
;	master.dir must exist the, .ap, .psf files corresponding
;	to the image, and the image itself.
;
; MODIFICATION HISTORY:
;	Written, J. Moustakas, 2000 March 18-19
;	generalized to include HST data, jm00jun10ucb
;-

        npar = n_params()
        if npar ne 1 then begin
            print, 'Syntax: trgb_starlist, objname, verbose=verbose, hst=hst, pc=pc'
            return
        endif

; delete any old files in this directory

        spawn, ['\rm '+objname+'.als*']
        spawn, ['\rm '+objname+'s.coo*']
        spawn, ['\rm '+objname+'s.ap*']
        spawn, ['\rm '+objname+'.cmb*']
        spawn, ['\rm '+objname+'.grp*']
        spawn, ['\rm '+objname+'s.grp*']
        spawn, ['\rm '+objname+'s*.fits']

;	spawn, ['clear']
	spawn, ['pwd'], datapath	; current directory
        masterdir = rdtxt(datapath+'/master.dir')
        masterdir = masterdir[0]

; create the environment setting and temporary alias to the image

        setenv, 'imdir='+masterdir
        spawn, ['alias imdir cd $imdir']

 ; create the .opt files

        if not keyword_set(hst) then begin ; KECK data

; pop out the fits header for this image and the IRAF psf

            spawn, ['find '+masterdir+' -name "'+objname+'.psf.fits" -print'], irafpsf
            psfhead = headfits(irafpsf[0])
            spawn, ['find '+masterdir+' -name "'+objname+'.pars" -print'], parfile

            trgb_makeopt, parfile, psfhead

        endif else trgb_hst_makeopt, pc=pc ; HST data                        

        allgood = 0L
        i = 1L

        repeat begin

            offset = 20000*i
            
; ALLSTAR script
            openw, lun1, 'allstar_script', /get_lun
            printf, lun1, ' '
            printf, lun1, 'imdir:'+objname+'.fits'    ; input image
            printf, lun1, 'imdir:'+objname+'.psf'     ; psf
            if i eq 1 then begin
                printf, lun1, 'imdir:'+objname+'.ap'  ; photometry file
            endif else begin
                printf, lun1, objname+'.cmb.'+strn(i-1L) ; most recent photometry file
            endelse
            printf, lun1, objname+'.als.'+strn(i)      ; output starlist
            printf, lun1, objname+'s.'+strn(i)+'.fits' ; output subtracted image
            printf, lun1, 'EXIT'
            free_lun, lun1
            
            if keyword_set(verbose) then begin
                spawn, ['allstar < allstar_script'] 
            endif else spawn, ['allstar < allstar_script > allstar_log']
            
            if not keyword_set(hst) then begin
                print & print, 'To search for the next level of stars,'
                read, 'type 0, otherwise type 1 to stop: ', allgood
            endif else if i eq 2L then allgood = 1L

            if long(allgood) eq 1L then goto, FINISH
            
; create the next level of stars if necessary

            if not keyword_set(hst) then begin
                print
                read, threshold, prompt='New threshold: '
                print
            endif else threshold = '5'

            openw, lun1, 'daophot_script', /get_lun
            printf, lun1, 'ATTACH '+objname+'s.'+strn(i)+'.fits'
            printf, lun1, 'OPTIONS'
            printf, lun1, ' '
            printf, lun1, 'thr = '+strn(threshold)
            printf, lun1, ' '

            printf, lun1, 'FIND'
            printf, lun1, '1,1'
            printf, lun1, objname+'s.coo.'+strn(i)  ; output coordinate file
            printf, lun1, 'Y'			  ; answer to "Are you happy with this?"
            printf, lun1, 'PHOT'
            printf, lun1, 'photo.opt'
            printf, lun1, ' '
            printf, lun1, objname+'s.coo.'+strn(i)  ; input coordinate file
            printf, lun1, objname+'s.ap.'+strn(i)   ; output photometry file
            printf, lun1, 'OFFSET' 		  ; offset the ID numbers
            printf, lun1, objname+'s.ap.'+strn(i)
            printf, lun1, strn(offset)+'/'
            printf, lun1, objname+'s.ap.'+strn(i)
            printf, lun1, ' '   		  ; overwrite
            printf, lun1, 'GROUP' 		  ; put the photometry files...
            printf, lun1, objname+'s.ap.'+strn(i)   ; ...in the same format
            printf, lun1, 'imdir:'+objname+'.psf'; psf
            printf, lun1, '3'   		  ; critical overlap (pixels)
            printf, lun1, objname+'s.grp.'+strn(i)  ; output group file
            printf, lun1, 'GROUP'			
            printf, lun1, objname+'.als.'+strn(i)
            printf, lun1, 'imdir:'+objname+'.psf'
            printf, lun1, '3'
            printf, lun1, objname+'.grp.'+strn(i)
            printf, lun1, 'APPEND' 		   ; merge the starlists
            printf, lun1, objname+'s.grp.'+strn(i)
            printf, lun1, objname+'.grp.'+strn(i)
            printf, lun1, objname+'.cmb.'+strn(i)	   ; output starlist
            printf, lun1, 'EXIT'
            free_lun, lun1
            
            print & print, 'Searching for the next level of stars.' & print
            if keyword_set(verbose) then begin
                spawn, ['daophot < daophot_script'] 
            endif else spawn, ['daophot < daophot_script > daophot_log']            
            
            i = i+1

            FINISH:

        endrep until allgood

        rootname = strmid(objname,0,strpos(objname,'_'))	; root object name
            
        print & print, 'Creating '+rootname+'.mag in the master directory . . . '
        spawn, ['\cp '+objname+'.als.'+strn(i)+' '+masterdir+'/'+rootname+'.mag']

; destroy the temporary alias and environment variable

        spawn, ['unalias imdir']
        spawn, ['unsetenv imdir']

        spawn, ['\rm *jnk*.fits']

        print & print, 'Thanks for playing!' & print

return
end


