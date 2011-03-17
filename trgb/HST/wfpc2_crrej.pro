pro wfpc2_crrej, objname
;+
; NAME:
;	WFPC2_CRREJ
;
; PURPOSE:
;	This program takes raw WFPC2 data and data-quality masks and
;	uses the IDL routine cr_reject to flag the cosmic rays in each
;	image.  Bad pixels (cosmic rays AND those pixels flagged as
;	bad in the data-quality mask) are set to a saturation value of
;	1.0E6.  The fits images are written back out to the working
;	directory with mwrfits.
;
; CALLING SEQUENCE:
;	wfpc2_crrej, objname
;
; INPUTS:
;	objname : string name of the galaxy
;
; OPTIONAL INPUTS:
;	None.
;
; KEYWORD PARAMETERS:
;	None.
;
; OUTPUTS:
;	The data fits files are written to the object directory, but
;	see SIDE EFFECTS below.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	The fits file written out by MWRFITS is larger than the
;	original fits file by 2.88 kb.  This may be because of a bug
;	in MWRFITS somewhere.  The WFPC2 image, however, does not seem
;	to be any different, and all the header information is still
;	present. 
;
; PROCEDURE:
;	Symbolic links to the raw data should be created using
;	wfpc2_links.  This routine then reads in the data and the
;	data-quality mask and calls IDL's cr_reject to flag bad pixels
;	using a noise model (the routine works fairly well even when
;	only two images exist).  The data are written back out to the
;	working directory.  
;
; EXAMPLE:
;	wfpc2_crrej, 'ngc2903'
;
; PROCEDURES USED:
;	MRDFITS, MWRFITS, CR_REJECT
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 August 3, UCB
;
;-
	
	npar = n_params()
        if npar ne 1 then begin
            print, 'Syntax: wfpc2_crrej, objname'
            return
        endif

	flist = findfile(objname+'*.fits')
        fdata = flist[where(strpos(flist,'_dq') eq -1L)]
        fmask = flist[where(strpos(flist,'_dq') gt 0L)]
        ndata = n_elements(fdata)
        
        if n_elements(fmask) ne ndata then begin
            print
            print, 'The data files do not match the data-quality files!'
            return
        endif

; read in the data, the headers, and the data-quality masks

        imdata = fltarr(800,800,4,ndata)
        imst = ptrarr(ndata)
        imask = fltarr(800,800,4,ndata)
        imhead = ptrarr(ndata)
        imsthead = ptrarr(ndata)

        print & print, 'Reading in the data . . . '

        for k = 0L, ndata-1L do begin

            imdata[*,*,*,k] = mrdfits(fdata[k],0,head,/silent)
            imhead[k] = ptr_new(head)
            imask[*,*,*,k] = readfits(fmask[k],/silent)

            temp = mrdfits(fdata[k],1,shead,/silent)
            imst[k] = ptr_new(temp)
            imsthead[k] = ptr_new(shead)

        endfor
        
; the readnoise and gain are WFPC2-specific values

        rdnoise = 0.71
        gain = 7.00

        chip = fltarr(800,800,ndata)
        mask = bytarr(800,800,ndata)

        for j = 0L, 3L do begin
            
            print & print, 'Flagging cosmic rays in chip '+strn(j)+' . . .'

            chip[*,*,*] = imdata[*,*,j,*]
            mask[*,*,*] = byte(imask[*,*,j,*]) eq 0B ; 0B is bad, 1B is good

            cr_reject, chip, rdnoise, 0, gain, 0.1, combined_image, $
              combined_noise, combined_npix, $ ; exptime=exp, $
              nsig = [8,6,4], input_mask=mask, dilation=1, dfactor=0.5, $
              /noskyadjust, /noclearmask, mask_cube=cr_mask;, /verbose

            chip[where(cr_mask eq 0B)] = 1.0E6	; set bad pixels to a big number
            imdata[*,*,j,*] = chip

        endfor
            
; write the data back out

        for i = 0L, ndata-1L do begin

            print & print, 'Writing '+fdata[i]+'.'
            spawn, ['\rm '+fdata[i]]
            mwrfits, imdata[*,*,*,i], fdata[i], *imhead[i], /create
            mwrfits, *imst[i], fdata[i], *imsthead[i], /ascii

        endfor

        ptr_free, imhead
        ptr_free, imst
        ptr_free, imsthead

        print & print, 'Done.' & print

return
end




