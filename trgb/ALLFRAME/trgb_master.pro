pro trgb_master, objname, verbose=verbose
;+
; NAME:
;	TRGB_MASTER
;
; PURPOSE:
;	Master ALLFRAME photometric analysis of the Keck TRGB data:
;	(1) convert IRAF psf's to DAOPHOT psf's; (2) generate a master
;	starlist; (3) run ALLFRAME; (4) determine the aperture
;	corrections; (5) calibrate the photometry to a standard
;	system. 
;
; INPUTS:
;	objname : galaxy name
;
; KEYWORD PARAMETERS:
;	verbose : verbose (to the screen) ALLFRAME output
;
; OUTPUTS:
;
; PROCEDURES USED:
;	TRGB_FIXPSF, TRGB_STARLIST, TRGB_ALLFRAME, TRGB_APCORRECT,
;	TRGB_CALIBRATE 
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 May 26, UCB
;-

	trgb_fixpsf, verbose=verbose	; convert IRAF psf to DAOPHOT psf

; check to make sure the psf's were created.  if there were any bad
; pixels, make the psf's by hand using trgb_manualpsf.

        trgb_starlist, verbose=verbose	

        trgb_allframe, verbose=verbose

        trgb_apcorrect, verbose=verbose

        trgb_calibrate, objname, dataout
        
return
end
