; ###########################################################################
; RELEGATED BY ATLAS_NED_WEBGET!!!
; ###########################################################################

stop

;+
; NAME:
;       ATLAS_DIAMETERS
;
; PURPOSE:
;       Collect size and position angle information from NED for the
;       spectral atlas.
;
; CALLING SEQUENCE:
;       atlas_diameters, diameters, /write
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

pro atlas_diameters, diameters, write=write
    
    analysis_path = atlas_path(/analysis)
    atlas = mrdfits(analysis_path+'atlas_ned.fits.gz',1,/silent)

    galaxy = strtrim(atlas.nedgalaxy,2)
    ngalaxy = n_elements(galaxy)

    outfile = 'atlas_diameters.fits'
    
    ned_webget_diameters, galaxy, diameters, outfile=outfile, $
      outpath=analysis_path, write=write

return
end
    
