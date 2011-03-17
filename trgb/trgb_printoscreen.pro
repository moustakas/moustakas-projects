pro trgb_printoscreen, trgb, trgberr, infobase, distance, vpec, log=log
;+
; NAME:
;	TRGB_PRINTOSCREEN
;
; PURPOSE:
;	Print some basic results of a TRGB edge-detection.
;
; INPUTS:
;	trgb	 : magnitude of the TRGB
;	trgberr  : magnitude error of the TRGB
;	infobase : object informational structure (TRGB_OBJECT_DATA)
;
; KEYWORD PARAMETERS:
;	log	: set if the TRGB was detected from a logarithmic
;		  luminosity function
;
; OUTPUTS:
;	distance : galaxy distance
;	vpec	 : galaxy peculiar velocity
;
; COMMON BLOCKS:
;
; RESTRICTIONS:
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 26, UCB
;-

        h_0 = 69.               ; Hubble's constant
        light = 2.99E5          ; speed of light
        itrgb = -4.06           ; absolute magnitude of the TRGB

        distance = 10^((trgb-itrgb-25.)/5.) ; Mpc

        dfar = 10^((trgb-itrgb-trgberr-25.)/5.)
        dclose = 10^((trgb-itrgb+trgberr-25.)/5.)
        derror = abs(dfar-dclose)

        vpec = light*infobase.z_helio - h_0*distance

        vpecfar = light*infobase.z_helio - h_0*dfar
        vpecclose = light*infobase.z_helio - h_0*dclose
        vpecerror = abs(vpecfar-vpecclose)
        
        if keyword_set(log) then begin
            print
            print, 'Logarithmic response: '
            print, '--------------------------------------------------'
        endif else begin
            print
            print, 'Linear response: '
            print, '--------------------------------------------------'
        endelse

        print, 'TRGB-mean (mag)     = '+strn(trgb,form='(F5.2)')+' \pm '+strn(trgberr,form='(F5.2)')
;       print, 'TRGB-median (mag)   = '+strn(trgbmed,form='(F5.2)')+' \pm '+strn(trgberr,form='(F5.2)')
        print, 'Distance (Mpc)      = '+strn(distance,form='(F5.2)')+' \pm '+strn(derror,form='(F5.2)')
        print, 'Peculiar vel (km/s) = '+strn(vpec,form='(F7.1)')+ ' \pm '+strn(vpecerror,form='(F5.2)')
        print

return
end

