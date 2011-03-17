function read_lf
;+
; NAME:
;	READ_LF()
;
; PURPOSE:
;	Read in a luminosity function.
;
; INPUTS:
;	None.
; OUTPUTS
;	lf - luminosity function structure
;    
; COMMENTS:
;
; PROCEDURES USED:
;	READCOL
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 October 25, U of A
;-

    path = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='lib')
    readfast, path+'fir_lf.dat', lfdata, skip=4, /double
    lum_lf = reform(lfdata[0,*]) & phi_lf = reform(lfdata[1,*])

    nlf = n_elements(lum_lf)
    logbinsz = 0.4D ; mag
    
    lum = 10D^(lum_lf) * 3.826D26 * 60D / 2.99793e14 ; L_nu (W/Hz) at 60 mu
    phi = 10D^(phi_lf) * 2.5D / (lum * alog(10))     ; Phi (#/Mpc^3/L_nu)

; pack everything into a structure
    
    lf = {name: 'Luminosity Function', lum_lf: double(lum_lf), $
          phi_lf: double(phi_lf), lum: double(lum), phi: double(phi), $
          logbinsz: double(logbinsz), nlf: n_elements(lum)}

return, lf
end

