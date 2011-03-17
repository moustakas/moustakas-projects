;+
; NAME:
;       KENNICUTT92_DIAMETERS
;
; PURPOSE:
;       Collect size and position angle information from NED for the
;       KENNICUTT92. 
;
; CALLING SEQUENCE:
;       kennicutt92_diameters, diameters, /write
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
;       J. Moustakas, 2005 Oct 20, U of A
;-

pro kennicutt92_diameters, diameters, write=write
    
    root = '92kennicutt'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'

    kennicutt92 = mrdfits(path+'92kennicutt_ned.fits.gz',1,/silent)

    galaxy = strtrim(kennicutt92.nedgalaxy,2)
    ngalaxy = n_elements(galaxy)

    outfile = '92kennicutt_diameters.fits'
    
    ned_webget_diameters, galaxy, diameters, outfile=outfile, $
      outpath=path, write=write

return
end
    
