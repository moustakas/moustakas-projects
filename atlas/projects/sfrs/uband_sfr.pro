function uband_sfr, loglb, logluband
; jm05aug26uofa - given the B-band luminosity and an observed U-band
;                 luminosity return the SFR and the 1-sigma
;                 uncertainty

    nloglb = n_elements(loglb)
    nuband = n_elements(logluband)
    
    if (nloglb eq 0L) or (nuband eq 0L) then begin
       print, 'sfr = uband_sfr(loglb,logluband)'
       return, -1L
    endif

    if (nloglb ne nuband) then begin
       splog, 'LOGLB and LOGLUBAND must have the same number of elements!'
       return, -1L
    endif

; restore the binary FITS table    
    
    uband_sfrpath = atlas_path(/projects)+'sfrs/'
    uband_sfrfile = 'uband_sfr.fits'

    uband_sfr = mrdfits(uband_sfrpath+uband_sfrfile,1,/silent)

;   logsfr = interpol(uband_sfr.mean,uband_sfr.loglb,loglb) + logluband - 41.0
    logsfr = interpol(uband_sfr.p50,uband_sfr.loglb,loglb) + logluband - 42.0
    sfr = 10.0^logsfr
    
return, sfr
end
    
