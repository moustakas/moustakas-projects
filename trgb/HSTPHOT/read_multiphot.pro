pro read_multiphot, objname, ni, nv
;+
; NAME:
;	READ_MULTIPHOT
;
; PURPOSE:
;	Read the photometry file from multiphot and create an
;	objname_starlist.dat file (also see TRGB_HST_COMBINE_STARLISTS).  
;
; INPUTS:
;	objname : string galaxy name
;	ni 	: total number of I-band images
;	nv 	: total number of V-band images
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	Creates a file called objname_starlist.dat in the object's directory.
;
; COMMON BLOCKS:
;
; RESTRICTIONS:
;	If there is no color information (no V-band data) then nv
;	should be set to 0.
;
; PROCEDURES USED:
;	READFAST()
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 August 2, UCB
;-

        npar = n_params()
        if npar ne 3 then begin
            print
            print, 'The total number of I- and V-band images need to be specified.'
            print
            return
        endif

        paths = trgb_datapath()
        pushd, paths[0]+objname

        if (long(ni) eq 0L) and (long(nv) eq 0L) then begin
            
            print
            print, 'There must be at least one I-band image!'
            return

        endif else begin

            pcscale = 0.045 / 60. ; PC plate scale (arcmin/pixel)
            wfscale = 0.1 / 60. ; WFC plate scale (arcmin/pixel)

            if long(nv) eq 0L then begin ; no V data, so calibrate the I data using the TRGB color

                ncol = 8*(ni+1)+7
        
                print & print, 'Reading in the data . . . '

                readfast, objname+'_multi.dat', data, ncols=ncol, skip=0
                ntot = n_elements(data[0,*])

; make cuts based on the following criteria (as suggested by dolphin):
; overall properties:  chi < 5.0, -0.5 < sharp < 0.5, object < 3, mag,
; < 50., S/N > 5. individual properties: chi < 5.0, -0.5 < sharp <
; 0.5, S/N > 3.5

                good = where((data[3,*] lt 5.) and (data[5,*] gt -0.5) and $
                             (data[5,*] lt 0.5) and (data[6,*] lt 3) and $
                             (data[9,*] lt 90.) and (data[4,*] gt 5.0),nstars)
                
                imag = data[9,good] 			 ; instrumental magnitude
                itrue = imag - 0.063*1.42 + 0.025*1.42^2 ; (V-I)_TRGB = 1.42

                starlist = fltarr(4,nstars)

                starlist[0,*] = data[1,good]  ; x-center (0-relative)
                starlist[1,*] = data[2,good]  ; y-center (0-relative)
                starlist[2,*] = itrue  	      ; I magnitude (transformed)
                starlist[3,*] = data[11,good] ; I error

            endif else begin	; I and V data

                ncol = 8*(ni+nv+2)+7

                print & print, 'Reading in the data . . . '

                readfast, objname+'_multi.dat', data, ncols=ncol, skip=0
                ntot = n_elements(data[0,*])

                good = where((data[3,*] lt 5.0) and (data[5,*] gt -0.5) and $
                             (data[5,*] lt 0.5) and (data[6,*] lt 3) and $
                             (data[4,*] gt 5.0) and (data[18,*] lt 30.) and $
                             (data[10,*] lt 30.),nstars)
        
                starlist = fltarr(7,nstars) ; add ID numbers afterwards
                
                starlist[0,*] = data[1,good] ; x-center (0-relative)
                starlist[1,*] = data[2,good] ; y-center (0-relative)
                starlist[2,*] = data[18,good] ; I magnitude (transformed)
                starlist[3,*] = data[19,good] ; I error
                starlist[4,*] = data[10,good] ; V magnitude (transformed)
                starlist[5,*] = data[11,good] ; V error
                starlist[6,*] = data[10,good]-data[18,good] ; V-I color
                
            endelse

; transform the pixel coordinates to a global wfpc2 system (0,0)->(1599,1599)

            chip0 = where(long(data[0,good]) eq 0L)
            chip1 = where(long(data[0,good]) eq 1L)
            chip2 = where(long(data[0,good]) eq 2L)
            chip3 = where(long(data[0,good]) eq 3L)

            starlist[0:1,chip0] =  starlist[0:1,chip0]*pcscale/wfscale + 799.
            arrchip1 = starlist[0:1,chip1]
            starlist[0,chip1]   = -arrchip1[1,*]+799.
            starlist[1,chip1]   =  arrchip1[0,*]+799.   
            starlist[0:1,chip2] = -starlist[0:1,chip2]+799.
            arrchip3 = starlist[0:1,chip3]
            starlist[0,chip3]   =  arrchip3[1,*]+799.
            starlist[1,chip3]   = -arrchip3[0,*]+799.

            id = long(findgen(nstars)) ; sequential ID numbers
        
            print & print, 'Writing '+objname+'_starlist.dat . . .'

            openw, lun1, objname+'_starlist.dat', /get_lun
            if long(nv) ne 0L then $
              printf, lun1, '#  ID   Xcenter   Ycenter    I       Ierr      V       Verr     V-I' else $
              printf, lun1, '#  ID   Xcenter   Ycenter    I       Ierr'
            printf, lun1, ' '
            for j = 0L, nstars-1L do $
              printf, lun1, id[j], starlist[*,j], format = '(1x,I6,7F9.3)'
            free_lun, lun1
        
        endelse

        popd
        
return
end
