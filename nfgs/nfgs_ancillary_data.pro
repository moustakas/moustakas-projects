;+
; NAME:
;       NFGS_ANCILLARY_DATA
;
; PURPOSE:
;       Compile ancillary data for the NFGS.
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; COMMENTS:
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 May 14 - originally written
;       jm02-04 - many updates
;       jm05jul24uofa - more updates
;-

pro nfgs_ancillary_data, table, write=write

    datapath = nfgs_path(/analysis)
    outpath = datapath

    outname = 'nfgs_ancillary_data.fits'
    
    basicname = 'nfgs_ned.fits.gz'
    photoname = 'nfgs_ned_photo.fits.gz'
    distname = 'nfgs_distances.fits.gz'
    diamname = 'nfgs_diameters.fits.gz'

    leda = nfgs_read_leda()
    
    write_ancillary_data, table, datapath=datapath, outpath=outpath, $
      basicname=basicname, photoname=photoname, distname=distname, $
      diamname=diamname, leda=leda, outname=outname, write=write

return
end
