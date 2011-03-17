function make_astr, st
; jm00aug2ucb
; create an astrometry structure using the header information

	cd = [[st.cd1_1,st.cd1_2],[st.cd2_1,st.cd2_2]]	; degrees/pixel
        cdelt = [1,1]
        crpix = [st.crpix1,st.crpix2]   ; reference pixel coordinates (pixels)
        crval = [st.crval1,st.crval2]	; reference pixel coordinates (RA,DEC)
        ctype = [st.ctype1,st.ctype2]	; coordinate type (RA--TAN,DEC--TAN)

        longpole = 180.
        projp1 = -1.
        projp2 = -2.

	astr = {cd: double(cd), cdelt: double(cdelt), $
                    crpix: float(crpix), crval:double(crval), $
                    ctype: string(ctype), longpole: float(longpole[0]),  $
                    projp1: float(projp1[0]), projp2: float(projp2[0])}

	return, astr

end

pro hst_align, objname, shift
;+
; NAME:
;	HST_ALIGN
;
; PURPOSE:
;	Determine the amount of dithering between any number of images
;	relative to the first image for HST/WFPC2 data using the
;	astrometry in the header.
;
; INPUTS:
;	objname : object name (string)
;
; KEYWORD PARAMETERS:
;	None.
;
; OUTPUTS:
;	shift : two-dimensional array containing the x and y shifts in
;		pixels of each image relative to the first image
;
; COMMON BLOCKS:
;;
; SIDE EFFECTS:
;;
; PROCEDURE:
;	Calls mrdfits to read in the ascii extension for a WFPC2
;	image, creates an astrometry structure for the first image,
;	and calculates the shift in pixels of each subsequent image
;	using ad2xy which converts the RA and DEC of a reference pixel
;	to pixels. 
;
; EXAMPLE:
;	hst_align, 'ngc1705'
;
; MODIFICATION HISTORY:
;	John Moustakas, 1 August 2000, UCB
;
;-

	path = '/deepscr1/ioannis/trgb/'
        pushd, path+objname

; generate the file list

        flist = findfile(objname+'*.fits')
        fdata = flist[where(strpos(flist,'_dq') eq -1L)]
        ndata = n_elements(fdata)

; solve for the shifts using chip 0

        im = mrdfits(fdata[0],0,imhead,/silent)
        st = mrdfits(fdata[0],1,head,/silent) ; reference (first) image
        astr = make_astr(st[0]) ; create an astrometry structure

        raref = st[0].crval1    ; RA of reference pixel
        decref = st[0].crval2   ; DEC of reference pixel
        paref = sxpar(imhead,'ORIENTAT') ; position angle of the pointing
        
        ad2xy, raref, decref, astr, xref, yref
        
        shift = fltarr(2,ndata) ; dithering wrt the first image
        shift[0,0] = [0.,0.]
            
        for k = 1L, ndata-1L do begin

            im = mrdfits(fdata[0],0,imhead,/silent)
            st = mrdfits(fdata[k],1,/silent)
            
            ra = st[0].crval1  
            dec = st[0].crval2 
            pa = sxpar(imhead,'ORIENTAT')
                
            ad2xy, ra, dec, astr, x, y
            
            shift[*,k] = [x-xref,y-yref]
            print, k, paref, pa, x-xref, y-yref
            
        endfor
        
; ----------------------------------------------------------------------

        

        popd

stop            
return
end
