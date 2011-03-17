function trgb_robust, mags, sigs
;+
; NAME:
;	TRGB_ROBUST
;
; PURPOSE:
;	Calculate the weighted average of HST photometric magnitudes. 
;
; INPUTS:
;	mags : magnitude array to average
;	sigs : magnitude error array corresponding to mags
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 10, UCB, based on J. Newman's ROBUST 
;-

; flag saturated stars

	bad = where(mags gt 90,nbad)

        if nbad eq n_elements(mags) then begin

            newmag = 99.9999
            return, newmag

        endif else begin

            if total(bad) gt -1L then sigs[bad] = 10.E10

; convert to data numbers (ADU)

            flx = 10D^((double(25) - mags)/2.5D)
            sigflx = (alog(10)/2.5D)*flx*sigs

            wghts = 1D/sigflx^2 ; compute the weighted average        
            flxavg = total(flx*wghts)/total(wghts)
        
            newmag =  25D - 2.5D*alog10(flxavg)

            return, newmag

        endelse

return
end



