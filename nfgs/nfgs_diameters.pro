;+
; NAME:
;       NFGS_DIAMETERS
;
; PURPOSE:
;       Collect size and position angle information from NED for the
;       NFGS. 
;
; CALLING SEQUENCE:
;       nfgs_diameters, diameters, /write
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
;       J. Moustakas, 2005 Jul 22, U of A
;-

pro nfgs_diameters, diameters, write=write
    
    analysis_path = nfgs_path(/analysis)
    nfgs = mrdfits(analysis_path+'nfgs_ned.fits.gz',1,/silent)

    galaxy = strtrim(nfgs.nedgalaxy,2)
    ngalaxy = n_elements(galaxy)

    outfile = 'nfgs_diameters.fits'
    
    ned_webget_diameters, galaxy, diameters, outfile=outfile, $
      outpath=analysis_path, write=write

return
end
    
