;pro check_solution, refim, im, xstars, ystars, xst, yst, dither
;
;        window, 2, xs=300, ys=300
;        plot, xstars, ystars, ps=3
;        oplot, xst, yst, ps=2, syms=0.2
;
;stop
;return
;end

function make_astr, st
; jm00aug2ucb
; create an astrometry structure using the header information.  this
; function replaces EXTAST because of the way WFPC2 fits files are
; written. 

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

pro wfpc2_dithering, objname
;+
; NAME:
;	WFPC2_DITHERING
;
; PURPOSE:
;	Determine the amount of dithering between any number of images
;	relative to the a first image for HST/WFPC2 data using two
;	methods:  the astrometry in the header (expected pointing);
;	and starlist matching (empirical solution).
;
; INPUTS:
;	objname : object name (string)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	Writes an IDL save set structure in the object directory with
;	the dithering solutions.
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; PROCEDURE:
;	Calls MRDFITS to read in the ascii extension for a WFPC2
;	image, creates an astrometry structure for the first image,
;	and calculates the shift in pixels of each subsequent image
;	using AD2XY which converts the RA and DEC of a reference pixel
;	to pixels.  Also generates a starlist for each image using
;	FIND and solves for the offsets using OFFSET_FROM_PAIRS.  All
;	the offset are with respect to the first image.
;
; COMMENTS:
;	It is strongly recommended that WFPC2_CRREJ be run on all the
;	images before this procedure because FIND is easily confused
;	by cosmic rays.
;
; EXAMPLE:
;	wfpc2_dithering, 'ngc1705'
;
; PROCEDURES USED:
;	MRDFITS, SXPAR(), AD2XY, FIND, OFFSET_FROM_PAIRS(), SWRITE()
;
; MODIFICATION HISTORY:
;	John Moustakas, 1 August 2000, UCB
;
;-

	if (not keyword_set(objname)) then $
            message, 'Syntax: wfpc2_dithering, objname'

        paths = trgb_datapath()
        pushd, paths[0]+objname

; generate the file list

        flist = findfile(objname+'*.fits')
        fdata = flist[where(strpos(flist,'_dq') eq -1L)]
        ndata = n_elements(fdata)

	template = {name:		' ', $
                    image_name:		strarr(1), $ ; image name
                    pa:			fltarr(1), $ ; position angle of pointing
                    ra_dec_ref:		dblarr(2), $ ; RA & DEC of PC reference pixel
                    ast_shift:		dblarr(2), $ ; astrometric x & y shifts (pixels)
                    star_shift:		dblarr(2)}   ; star-matching x & y shifts (pixels)
                    
        dither = replicate(template,ndata)

; solve for the shifts using chip 0, relative to the first image in the file list

        refim = mrdfits(fdata[0],0,imhead,/silent)
        st = mrdfits(fdata[0],1,head,/silent)
        astr = make_astr(st[0]) 		; create an astrometry structure for CHIP 0

        ad2xy, st[0].crval1, st[0].crval2, astr, xref, yref
        
        dither[0].image_name = fdata[0]
        dither[0].pa = sxpar(imhead,'ORIENTAT') ; position angle of the pointinp
        dither[0].ra_dec_ref = [st[0].crval1,st[0].crval2]
        dither[0].ast_shift = [0.,0.]
        dither[0].star_shift = [0.,0.]

        jm_find, refim[*,*,0], xstars, ystars, flux, sharp, round, 4.0, 1.4, $ ; generate a starlist
          [-1.0,1.0], [0.2,1.0], /silent

        fsort = sort(flux)			; sort by brightness
        nmax = long(0.5*n_elements(xstars))	; use 50% of the starlist
        xstars = xstars[fsort[0:nmax-1L]]
        ystars = ystars[fsort[0:nmax-1L]]
        
        for k = 1L, ndata-1L do begin

            im = mrdfits(fdata[k],0,imhead,/silent)
            st = mrdfits(fdata[k],1,/silent)

            ad2xy, st[0].crval1, st[0].crval2, astr, x, y ; astrometric solution

; ----------------------------------------------------------------------            
            im[*,*,0] = shift(im[*,*,0],xref-x,yref-y)
; ----------------------------------------------------------------------            
            
            jm_find, im[*,*,0], xst, yst, flux, sharp, round, 4.0, 1.4, $ ; match starlists
              [-1.0,1.0], [0.2,1.0], /silent

            fsort = sort(flux)
            nmax = long(0.5*n_elements(xst))
            xst = xst[fsort[0:nmax-1L]]
            yst = yst[fsort[0:nmax-1L]]

            print
            print, 'Matching '+strn(long(0.5*n_elements(xst)))+' stars . . .'
            xyshift = offset_from_pairs(xstars,ystars,xst,yst,dmax=8,$
                                       minpeak=4.,binsz=0.5)

; the sense of the shift is correct in the astrometric solution

;            if (long(x-xref) gt 0L) then xyshift[0] = abs(xyshift[0]) else $
;              xyshift[0] = -abs(xyshift[0])
;            if (long(y-yref) gt 0L) then xyshift[1] = abs(xyshift[1]) else $
;              xyshift[1] = -abs(xyshift[1])

            dither[k].image_name = fdata[k]
            dither[k].pa = sxpar(imhead,'ORIENTAT')
            dither[k].ra_dec_ref = [st[0].crval1,st[0].crval2]
            dither[k].ast_shift = [xref-x,yref-y]
            dither[k].star_shift = [xref-x+xyshift[0],yref-y+xyshift[1]]

            print
            print, '       Astrometric     Empirical   '
            print, 'Image    x     y        x     y    '
            print, '-----------------------------------'
            print, ' '+strn(k,form='(I3)')+'   '+strn(xref-x,form='(F8.3)')+'   '+$
              strn(yref-y,form='(F8.3)')+'  '+strn(xref-x+xyshift[0],form='(F8.3)')+$
              '  '+strn(yref-y+xyshift[1],form='(F8.3)')

;           check_solution, refim, im, xstars, ystars, xst, yst, dither

        endfor        

        okay = 'N' & print
        read, okay, prompt='Write the solutions to a structure? (Y/[N])? '
        if strupcase(okay) eq 'Y' then begin
            print
            print, 'Writing '+objname+'_dither.dat.'
            swrite, dither, objname+'_dither.dat'
        endif

        popd

        while !d.window ne -1L do wdelete, !d.window ; delete all windows
        
return
end
