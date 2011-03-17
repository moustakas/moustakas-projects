;+
; NAME:
;       KENN92_DIAMETERS
;
; PURPOSE:
;       Collect size and position angle information from NED for the
;       Kennicutt (1992) sample.
;
; CALLING SEQUENCE:
;       kenn92_diameters, diameters, /write
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
;       J. Moustakas, 2005 Aug 02, U of A
;-

pro kenn92_diameters, diameters, write=write
    
    analysis_path = kenn92_path(/analysis)
    kenn92 = mrdfits(analysis_path+'kenn92_ned.fits.gz',1,/silent)

    galaxy = strtrim(kenn92.nedgalaxy,2)
    ngalaxy = n_elements(galaxy)

    outfile = 'kenn92_diameters.fits'
    
    ned_webget_diameters, galaxy, diameters, outfile=outfile, $
      outpath=analysis_path, write=write

return
end
    
