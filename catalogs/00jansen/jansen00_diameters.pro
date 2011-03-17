;+
; NAME:
;       JANSEN00_DIAMETERS
;
; PURPOSE:
;       Collect size and position angle information from NED for the
;       JANSEN00. 
;
; CALLING SEQUENCE:
;       jansen00_diameters, diameters, /write
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

pro jansen00_diameters, diameters, write=write
    
    root = '00jansen'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'

    jansen00 = mrdfits(path+'00jansen_ned.fits.gz',1,/silent)

    galaxy = strtrim(jansen00.nedgalaxy,2)
    ngalaxy = n_elements(galaxy)

    outfile = '00jansen_diameters.fits'
    
    ned_webget_diameters, galaxy, diameters, outfile=outfile, $
      outpath=path, write=write

return
end
    
