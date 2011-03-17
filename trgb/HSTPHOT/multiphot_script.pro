pro multiphot_script, objname, log=log, nocheck=nocheck
;+
; NAME:
;	MULTIPHOT_SCRIPT
;
; PURPOSE:
;	Generate a script that runs all pre-MULTIPHOT and MULTIPHOT
;	routines. 
;
; CALLING SEQUENCE:
;	multiphot_script, objname
;
; INPUTS:
;	objname : string name of the object
;
; KEYWORD PARAMETERS:
;	None.
;
; OUTPUTS:
;	Generates a script (text file) called multiphot_script in the
;	object directory.
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;	Dolphin's HSTPHOT photometry package is assumed to be
;	installed in the directory /deep1/ioannis/trgb/hstphot.  The
;	datapath is assumed to be /deepscr1/ioannis/trgb.  
;
; PROCEDURE:
;	The MULTIPHOT script is created based on the recipe provided
;	in the HSTPHOT documentation and includes running the routines
;	crclean, getsky, hotpixels, and multiphot.  The script is
;	closed and sourced, calling the appropriate HSTPHOT routines.
;	WARNING:  HSTPHOT generates lots of files!  The final
;	MULTIPHOT photometry file is called objname_multi.dat.
;
; EXAMPLE:
;	multiphot_script, 'ngc2683'
;
; PROCEDURES USED:
;	TRGB_DATAPATH()
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 August 1, UCB
;-

	if (not keyword_set(objname)) then $
          message, 'Syntax: multiphot_script, objname'

        paths = trgb_datapath()
        hstphot_path = '/deep1/ioannis/trgb/hstphot/'

        pushd, paths[0]+objname
        multidir = 'MULTIPHOT'
        spawn, ['mkdir -m 0755 -p '+multidir]

; copy the data to a temporary MULTIPHOT subdirectory

        if keyword_set(nocheck) then $
          spawn, ['\rm MULTIPHOT/*'] else begin
            check = 'N'
            read, check, prompt='Delete everything in the MULTIPHOT subdirectory (Y/[N])? '
            if strupcase(check) eq 'Y' then $
              spawn, ['\rm MULTIPHOT/*'] else return ; quit
        endelse 
        spawn, ['\cp *.fits '+multidir] ; copy

        flist = findfile(objname+'*.fits')
        flist = flist[where((strpos(flist,'deep') eq -1L) and (strpos(flist,'mosaic') eq -1L))]
        fdata = flist[where(strpos(flist,'_dq') eq -1L)]
        fmask = flist[where(strpos(flist,'_dq') gt 0L)]
        ndata = n_elements(fdata)

        if (n_elements(fmask) ne ndata) then begin
            print & print, 'The data files do not match the data-quality files!' & popd
            return
        endif

; index the I- and the V-band data

        idata = where(strpos(fdata,'_I') gt 0L,icount)
        vdata = where(strpos(fdata,'_V') gt 0L,vcount)

; multiphot output photometry file
        datfile = paths[0]+objname+'/'+objname+'_multi.dat'
        sfile = objname+'_multiphot_script'	; multiphot script name
        root = strarr(ndata)			; rootname string array
        multiarr = strarr(1)			; multiphot argument
        
        openw, lun1, paths[0]+objname+'/'+sfile, /get_lun
        printf, lun1, 'date'
        printf, lun1, 'cd '+multidir
        
        dfile = findfile(objname+'_dither.dat',count=fcount) ; dithering structure
        if fcount eq 1L then dither = sread(dfile[0])

        print & print, 'Generating the MULTIPHOT script . . . '

        if icount ne 0L then begin
            rooticr = ''
            for k = 0L, icount-1L do $
              rooticr = rooticr+' '+strmid(fdata[idata[k]],0,rstrpos(fdata[idata[k]],'.fits'))+'_cr.fits'
        endif
        if vcount ne 0L then begin
            rootvcr = ''
            for k = 0L, vcount-1L do $
              rootvcr = rootvcr+strmid(fdata[vdata[k]],0,rstrpos(fdata[vdata[k]],'.fits'))+'_cr.fits'+' '
        endif

        for k = 0L, ndata-1L do begin

            root[k] = strmid(fdata[k],0,rstrpos(fdata[k],'.fits'))

            if fcount eq 1L then $
              multiarr = multiarr+' '+root[k]+'_cr '+$
;             strn(dither[k].star_shift[0],form='(F6.2)')+' '+strn(dither[k].star_shift[1],form='(F6.2)') else $
              strn(dither[k].ast_shift[0],form='(F6.2)')+' '+strn(dither[k].ast_shift[1],form='(F6.2)') else $
              multiarr = multiarr+' '+root[k]+'_cr 0.0 0.0'             
            
            printf, lun1, hstphot_path+'mask '+fdata[k]+' '+fmask[k]

        endfor

        printf, lun1, 'echo Cleaning cosmic rays . . . '
        for k = 0L, ndata-1L do $	; clean cosmic rays
          printf, lun1, hstphot_path+'crclean '+fdata[k]+' 1 3 0 '+root[k]+'_cr.fits'

        if icount ne 0L then begin
            printf, lun1, 'echo Creating a deep I-band image . . . '
            printf, lun1, hstphot_path+'crclean '+rooticr+' 1 3 0 '+'../'+objname+'_deep_I.fits'
        endif

        if vcount ne 0L then begin
            printf, lun1, 'echo Creating a deep V-band image . . . '
            printf, lun1, hstphot_path+'crclean '+rootvcr+' 1 3 0 '+'../'+objname+'_deep_V.fits'
        endif

        printf, lun1, 'echo Estimating the sky in the deep images . . . '
        if icount ne 0L then printf, lun1, hstphot_path+'getsky '+'../'+objname+'_deep_I'
        if vcount ne 0L then printf, lun1, hstphot_path+'getsky '+'../'+objname+'_deep_V'

        printf, lun1, 'echo Generating a WFPC2 mosaic of the deep images . . . '
        if icount ne 0L then printf, lun1, hstphot_path+'hstmosaic '+'../'+objname+'_deep_I '+$
          '../'+objname+'_mosaic_I.fits'
        if vcount ne 0L then printf, lun1, hstphot_path+'hstmosaic '+'../'+objname+'_deep_V '+$
          '../'+objname+'_mosaic_V.fits'

;        if ((icount ne 0L) and (vcount ne 0L)) then begin
;            printf, lun1, 'echo Generating an RGB mosaic . . . '
;            printf, lun1, hstphot_path+'makergb '+'../'+objname+'_mosaic_V.fits 0.0015 0.052 "" 0.65 0.35 ../'+$
;              objname+'_mosaic_I.fits 0.0016 0.056 ../'+objname+'.bmp 2'
;        endif

        printf, lun1, 'echo Estimating the sky in the remaining images . . . '
        for k = 0L, ndata-1L do $	; get sky estimate
          printf, lun1, hstphot_path+'getsky '+root[k]+'_cr'

        printf, lun1, 'echo Flagging hot pixels . . . '
        for k = 0L, ndata-1L do $	; flag hot pixels
          printf, lun1, hstphot_path+'hotpixels '+root[k]+'_cr'

        printf, lun1, 'echo Running multiphot . . . '
        if vcount ne 0L then $	; use the deep V-band image as the reference image
          printf, lun1, hstphot_path+'multiphot '+datfile+' 3.5 5.0 -1 0 0 0 0 0 '+multiarr+' ../'+$
          objname+'_deep_V' else $
          if icount ne 0L then $ ; use the I-band image as the reference image
          printf, lun1, hstphot_path+'multiphot '+datfile+' 3.5 5.0 -1 0 0 0 0 0 '+multiarr+' ../'+$
          objname+'_deep_I' else begin
            print & print, 'There is a problem!' & return
        endelse
          
        printf, lun1, 'cd ..'
        printf, lun1, 'date'
        printf, lun1, 'echo Done!'

        free_lun, lun1

        print & print, 'Launching MULTIPHOT . . . '
        spawn, ['chmod +x '+objname+'_multiphot_script']
        if keyword_set(log) then $
          spawn, [paths[0]+objname+'/'+objname+'_multiphot_script > '+$
                  paths[0]+objname+'/'+objname+'_log'] else $
          spawn, [paths[0]+objname+'/'+objname+'_multiphot_script']
          
        print, '. . . done!'

        popd
        
return
end








