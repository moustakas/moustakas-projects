pro daoheader, paramfile, fhead, daohead, lst=lst, coo=coo, als=als
;+
; NAME:
;	DAOHEADER
;
; PURPOSE:
;	Create a two-line DAOPHOT header.
;
; INPUTS:
;	The program requires a single line text file with the names of
;	each DAOPHOT header value to exist, for example:
;		/deep1/ioannis/trgb/daohead.sample
; 	The other inputs are:
;		paramfile : The parameter file corresponding to the current
;		   	    image (scalar filename of the full directory path).
;		fhead     : The fits header of the current image.
;
; KEYWORD PARAMETERS:
; 	lst : Keyword specifying a .lst header is to be created.
;	coo : Keyword specifying a .coo header is to be created.
;	als : Keyword specifying a .als header is to be created.
;
; OUTPUTS:
;	daohead : A two-element string array of the new DAOPHOT header.
;
; PROCEDURE:
;	Extracts the necessary header information from the fits header
;	and from the parameter file.
;
; PROCEDURES USED:
;	RDTXT(), SXPAR(), IRAFHVAL()
; MODIFICATION HISTORY:
;	John Moustakas, 2000 March 16, UCB
;-

; read in the parameters file

        pars = rdtxt(paramfile)

; read in the DAOPHOT sample header

	daohead = strarr(2)

        dum = ' '
        if keyword_set(als) then begin
            file = '/deep1/ioannis/trgb/daohead.sample_als' 
            frmt = '(2x,I1,1x,I4,1x,I4,1x,F7.1,1x,F7.1,2x,F6.1,3x,f5.2,2x,f5.2,3x,f5.2,3x,f5.2)'
        endif else begin
            file = '/deep1/ioannis/trgb/daohead.sample'          
            frmt = '(2x,I1,1x,I4,1x,I4,1x,F7.1,1x,F7.1,2x,F6.1,3x,f5.2,2x,f5.2,3x,f5.2)'
        endelse

	openr, lun1, file, /get_lun
        readf, lun1, dum 
        daohead[0] = dum
        free_lun, lun1

; pull out information from the fits header

        nx = sxpar(fhead,'NAXIS1')
        ny = sxpar(fhead,'NAXIS2')
        gain = sxpar(fhead,'GAIN')
        rdnoise = sxpar(fhead,'RDNOISE')

; now pull out information from the IRAF headers

        highbad = irafhval(pars,'datapars.datamax')
        lowbad = irafhval(pars,'datapars.datamin')
        tsigmas = irafhval(pars,'findpars.threshold')
        skysig = irafhval(pars,'datapars.sigma')
        fwhm = irafhval(pars,'datapars.fwhmpsf')
        tcounts = skysig*tsigmas

        if keyword_set(lst) then $
          daohead[1] = string([1],[nx],[ny],[lowbad],[highbad],$
                              [tcounts],[fwhm],[gain],[rdnoise],format=frmt)
        
        if keyword_set(coo) then $
          daohead[1] = string([1],[nx],[ny],[lowbad],[highbad],$
                              [tcounts],[0.0],[gain],[rdnoise],format=frmt)

        if keyword_set(als) then $
          daohead[1] = string([1],[nx],[ny],[lowbad],[highbad],$
                              [tcounts],[fwhm],[gain],[rdnoise],[fwhm],format=frmt)

return        
end




