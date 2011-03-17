; ###########################################################################
; RELEGATED BY SINGS_NED_WEBGET!!!
; ###########################################################################

stop

;+
; NAME:
;       SINGS_DIAMETERS
;
; PURPOSE:
;       Collect size and position angle information from NED for the
;       SINGS. 
;
; CALLING SEQUENCE:
;       sings_diameters, diameters, /write
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       write - write out
;
; OUTPUTS:
;       diameters - output results from NED_WEBGET_DIAMETERS 
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;       NED_WEBGET_DIAMETERS
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2005 Jul 25, U of A
;-

pro sings_diameters, diameters, write=write
    
    analysis_path = sings_path(/analysis)
    sings = mrdfits(analysis_path+'sings_ned.fits.gz',1,/silent)

    galaxy = strtrim(sings.nedgalaxy,2)
    ngalaxy = n_elements(galaxy)

    outfile = 'sings_diameters.fits'
    
    ned_webget_diameters, galaxy, diameters, outfile=outfile, $
      outpath=analysis_path, write=write

return
end
    
