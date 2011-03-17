function writefile, daoform, newhead, daofile, frmt
; jm00mar16ucb
; function to create the output file

        print, 'Writing '+daofile+'.'
        openw, lun1, daofile, /get_lun
        printf, lun1, newhead[0] ; insert the daophot header
        printf, lun1, newhead[1] ; insert the header values
        printf, lun1, ' '        ; insert a blank line
        for j = 0L, n_elements(daoform[0,*])-1L do $
          printf, lun1, daoform[*,j], format = frmt
        free_lun, lun1

	return, 'writeit'
end

pro iraf2dao, irafile, parfile, fhead, daofile, lst=lst, coo=coo
;+
; NAME:
;	IRAF2DAO
;
; PURPOSE:
;	Convert IRAF format .coo or .pst files to DAOPHOT format .coo
;	or .lst files (including headers).
;
; INPUTS:
;	irafile : The IRAF file to be converted (e.g. image.coo.1)
;		  (scalar filename of the full directory path).
;	parfile : The parameter file corresponding to the current
;		  image (scalar filename of the full directory path).
;	fhed    : The fits header of the current image.
;	daofile : The name of the output DAOPHOT file (scalar filename
;		  of the full directory path). 
;
; KEYWORD PARAMETERS:
; 	lst : Keyword specifying a .lst file is to be created.
;	coo : Keyword specifying a .coo file is to be created.
;
; OUTPUTS:
;	The program creates a new .coo or .lst file in the current
;	directory. 
;
; PROCEDURE:
;	Reads in the IRAF file, rearranges the columns, creates a new
;	header, and writes out the DAOPHOT file.
;
; WARNINGS:
;	Any known old .coo.1 or .pst.4 files should be placed in an
;	OLD/ subdirectory or deleted.
;
; PROCEDURES USED:
;	READCOL, DAOHEADER

; MODIFICATION HISTORY:
;	John Moustakas, 2000 March 15-16, UCB
;-

        if keyword_set(lst) then begin

; read in the .pst file

            readcol, irafile, id, xcenter, ycenter, mag, msky, /silent, $
              format='I,F,F,F,F'
            daoform = [[id],[xcenter],[ycenter],[mag],[msky]]
            daoform = transpose(daoform)            

            daoheader, parfile, fhead, newhead, /lst

            frmt = '(1x,I5,3F9.3,F9.3)'
            writeit = writefile(daoform, newhead, daofile, frmt)


        endif else begin

            if keyword_set(coo) then begin

; read in the .coo file
                
                readcol, irafile, xcenter, ycenter, mag, sharp, sround, ground, id, /silent, $
                  format='F,F,F,F,F,F,I'
                daoform = [[id],[xcenter],[ycenter],[mag],[sharp],[sround]]
                daoform = transpose(daoform)

                daoheader, parfile, fhead, newhead, /coo

                frmt = '(1x,I5,3F9.3,2F9.3)'
                writeit = writefile(daoform, newhead, daofile, frmt)

            endif else begin

                print, 'You need to specify a file to convert with the keyword!'
                return

            endelse

        endelse
                           
return        
end
