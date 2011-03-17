;+
; NAME:
;       KENN92_ANCILLARY_DATA
;
; PURPOSE:
;       Compile ancillary data for the Kennicutt (1992) spectral
;       atlas. 
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
;       J. Moustakas, 2005 Aug 02
;-

pro kenn92_ancillary_data, table, write=write

    datapath = kenn92_path(/analysis)
    outpath = datapath

    outname = 'kenn92_ancillary_data.fits'
    
    basicname = 'kenn92_ned.fits.gz'
    photoname = 'kenn92_ned_photo.fits.gz'
    distname = 'kenn92_distances.fits.gz'
    diamname = 'kenn92_diameters.fits.gz'

    write_ancillary_data, table, datapath=datapath, outpath=outpath, $
      basicname=basicname, photoname=photoname, distname=distname, $
      diamname=diamname, outname=outname, write=write

return
end
