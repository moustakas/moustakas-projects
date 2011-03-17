;+
; NAME:
;	READ_DAOPHOT()
;
; PURPOSE:
;	Read a two-line DAOPHOT photometry file into a data array.
;	(The daophot photometry files have two lines of information
;	per star, followed by a blank space.  Arg!)
;
; INPUTS:
;	filename : complete path and name of the file to read
;
; KEYWORD PARAMETERS:
;	naps : number of PHOT apertures used to create the file.
;
; OUTPUTS:
;	data - returns a data array in the following format
;		ID x-position y-position mag (ap1) mag (ap2) . . . magerr (ap1) magerr (ap2) . . . 
;
; COMMON BLOCKS:
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 28, UCB
;       jm05may30uofa - converted to a function 
;-

function read_daophot, filename, naps=naps

    openr, lun1, filename, /get_lun
    
    spawn, ['wc -l '+filename], wcount
    reads, wcount, nrows
    nrows = long(nrows) - 3L

    ndata = nrows/3L            ; 2 data lines and 1 blank space compressed to 1 line
    ncols = 2L*naps+3L          ; mag and magerr for each aperture, plus 3 for ID and xy
    
    temp1 = fltarr(naps+3L)	; for each star
    temp2 = fltarr(naps)
    junk = strarr(1)
    data = fltarr(ncols,ndata)
    
    head = strarr(4)            ; the header is always 4 lines
    readf, lun1, head           ; read the header

    k = 0L

    while not eof(lun1) do begin

       readf, lun1, temp1, format='(1x,I5,'+string(2+naps,format='(I0)')+'F9.3)'
       readf, lun1, temp2, format='(25x,'+string(naps,format='(I0)')+'F9.3)'
       if k ne ndata-1L then readf, lun1, junk
       
       data[0:3+naps-1L,k] = temp1[0:3+naps-1L] ; ID, xy positions and mag
       data[3+naps:ncols-1L,k] = temp2[0:naps-1L] ; mag error

       k = k+1L

    endwhile

    free_lun, lun1

return, data
end
