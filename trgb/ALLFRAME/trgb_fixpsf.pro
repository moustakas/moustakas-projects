function checkfiles, pstlist, coolist, parslist, nsublist, irafpsflist, imnames
; jm00mar16ucb
; function to check with the user that the file names match up.

	print
        print, 'I found the following files: ' & print
        print, transpose(pstlist) & print		; IRAF psf starlist
        print, transpose(coolist) & print		; IRAF coordinate list
        print, transpose(parslist) & print		; IRAF parameter files
        print, transpose(nsublist) & print		; neighbor-subtracted image
        print, transpose(irafpsflist) & print		; IRAF psf
        print, 'And they all need to match with: ' & print
        print, transpose(imnames)			; current fits image
        okay = ' ' & print
        print, 'If this is okay type Y, otherwise press any key to quit: '
        read, okay
        return, okay

end

pro trgb_fixpsf, verbose=verbose
;+
; NAME:
;	TRGB_FIXPSF
;
; PURPOSE:
;	Convert IRAF DAOPHOT point-spread functions to FORTRAN DAOPHOT psfs.
;
; INPUTS:
;	The program assumes the following files exist in the current
;	directory:
; 		master.dir  : Lowest directory name of the current
; 			      object being analyzed.
;		image.names : A list, one per line, of all the frames
;			      of the current object.
;	The user should be in the directory where all the FITS images
;	reside. 
;
; OUTPUTS:
;	The program generates a new .coo and .lst file for every
;	image.  Then DAOPHOT is called to create .ap and .psf files.
;	The log files for PHOT and PSF are written to photlog and
;	psflog.  Files called phot_fixpsf.script and psf_fixpsf.script
;	which contain the DAOPHOT scripts are also created.
;
; PROCEDURE:
;	First checks to make sure the .pars, .coo.1, .pst files, the
;	neighbor-subtracted images, and the IRAF psfs all match up (in
;	order) with the list given in images.names.  Then, for each
;	image, creates the .opt files for that image, a new .coo, a
;	new .lst, a .ap, and a .psf file.  The psf is created from the
;	image that has had all its neighbors subtracted. 
;
; PROCEDURES USED:
;	TRGB_MAKEOPT, IRAF2DAO, TRGB_CONVERTALL, RDTXT(), HEADFITS()
;
; RESTRICTIONS:
;	(**) Any known old .coo.1 or .pst files should be placed in an
;	     OLD/ subdirectory or deleted.
;	(**) Deletes all files created by a previous run of this code
;            (to prevent DAOPHOT from crashing.
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 March 15-19, UCB
;-

        npar = n_params()
        if npar gt 0 then begin
            print, 'Syntax: trgb_fixpsf, verbose=verbose'
            return
        endif

        convert_type = ' ' 
        print & print, 'Type C to convert the data type in the'
        read, convert_type, prompt='header, or press ENTER to continue: '
        if strupcase(convert_type) eq 'C' then trgb_convertall

	spawn, ['pwd'], datapath	; current directory
        datapath = datapath[0]
        masterdir = rdtxt(datapath+'/master.dir')
        masterdir = masterdir[0]

        imnames = rdtxt(datapath+'/image.names') ; read in the image names

; find the iraf parameter files to generate the DAOPHOT header

        spawn, ['find '+masterdir+' -name "*.pars" -print'], parslist
        nold = where(strpos(parslist,'OLD') eq -1,count)
        if count ne 0 then parslist = parslist[nold]

; find the iraf master psf starlist

        spawn, ['find '+masterdir+' -name "*pst" -print'], pstlist
        nold = where(strpos(pstlist,'OLD') eq -1,count)
        if count ne 0 then pstlist = pstlist[nold]

; find the iraf coordinate files

        spawn, ['find '+masterdir+' -name "*coo.1" -print'], coolist
        nold = where(strpos(coolist,'OLD') eq -1,count)
        if count ne 0 then coolist = coolist[nold]

; find the neighbor-subtracted images

        spawn, ['find '+masterdir+' -name "*nsub.fits" -print'], nsublist
        nold = where(strpos(nsublist,'OLD') eq -1,count)
        if count ne 0 then nsublist = nsublist[nold]

; find the IRAF psf

        spawn, ['find '+masterdir+' -name "*psf.fits" -print'], irafpsflist
        nold = where(strpos(irafpsflist,'OLD') eq -1,count)
        if count ne 0 then irafpsflist = irafpsflist[nold]

; make sure all the files match up

        okay = checkfiles(pstlist,coolist,parslist,nsublist,irafpsflist,imnames)
        if strupcase(okay) ne 'Y' then return

        for i = 0L, n_elements(imnames)-1L do begin  ; loop on the images

            image = imnames[i]
            pst = pstlist[i]
            coo = coolist[i]
            nsub = nsublist[i]
            irafpsf = irafpsflist[i]
            pars = parslist[i]

            nname = strmid(nsub,strpos(nsub,image+'.')) ; name of neighbor-subtracted image
            nndir = strmid(nsub,0,strpos(nsub,nname))   ; directory of neighbor-subtracted image

; create the environment setting and temporary alias to the neighbor-subtracted image

            setenv, 'dir='+nndir	; dir is the image directory
            spawn, ['alias dir cd $dir']
            
; delete any old files for this image in this directory

            spawn, ['\rm '+image+'.coo']
            spawn, ['\rm '+image+'.ap']
            spawn, ['\rm '+image+'.psf']
            spawn, ['\rm '+image+'.nei']
;           spawn, ['\rm phot_fixpsf.script']
;           spawn, ['\rm psf_fixpsf.script']

; pop out the fits header for this image and the IRAF psf

            fhead = headfits(datapath+'/'+image+'.fits')
            psfhead = headfits(irafpsf)

; create the .opt files for DAOPHOT appropriate to this image

            trgb_makeopt, pars, psfhead

; create the new .lst and the new .coo files

            lstfile = masterdir+'/'+image+'.lst'
            iraf2dao, pst, pars, fhead, lstfile, /lst

            coofile = masterdir+'/'+image+'.coo'
            iraf2dao, coo, pars, fhead, coofile, /coo            

; PHOT script with the new .coo file
            openw, lun1, 'phot_fixpsf_script', /get_lun
            printf, lun1, 'ATTACH '+image ; image
            printf, lun1, 'PHOT'
            printf, lun1, 'photo.opt' 	  ; photometry file
            printf, lun1, ' '
            printf, lun1, image+'.coo' 	  ; coordinate file
            printf, lun1, image+'.ap' 	  ; output aperture photometry file
            printf, lun1, 'EXIT'
            free_lun, lun1

            print & print, 'Doing aperture photometry on '+image+'.' & print
            if keyword_set(verbose) then begin
                spawn, ['daophot < phot_fixpsf_script']
            endif else spawn, ['daophot < phot_fixpsf_script > phot_log']

; PSF script and make the psf
            openw, lun2, 'psf_fixpsf_script', /get_lun
            printf, lun2, 'ATTACH dir:'+nname ; neighbor-subtracted image
            printf, lun2, 'PSF'		
            printf, lun2, image+'.ap'	  ; aperture photometry file
            printf, lun2, image+'.lst' 	  ; psf starlist
            printf, lun2, image+'.psf'	  ; output psf
            printf, lun2, 'EXIT'
            free_lun, lun2

            print & print, 'Generating the PSF for '+image+'.' & print
            if keyword_set(verbose) then begin
                spawn, ['daophot < psf_fixpsf_script']
            endif else spawn, ['daophot < psf_fixpsf_script > psf_log']

; delete the neighbors output file

            spawn, ['\rm '+image+'.nei']

        endfor

; delete the temporary alias and environment variable

        spawn, ['unsetenv dir']
        spawn, ['unalias dir']

        print & print, 'All done.' & print

return
end




