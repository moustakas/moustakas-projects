; ----------------------------------------------------------------------
;   read_fit.pro
;
; Created to utilized the flux ratio method to compute the dust 
; attenuation of (starburst) galaxies at different wavelengths.  The
; details of this method are contained in Gordon K. D., Clayton,
; G. C., Witt, A. N., & Misselt, K. A. 2000, ApJ, 533, 236.  
;
; v 1.0 - Karl D. Gordon (Steward Obs., 9 June 2000)
; ----------------------------------------------------------------------

PRO read_fit,filename,fit_info

openr,unit1,filename,/get_lun
readf,unit1,n_bands,n_ages

bands = strarr(n_bands)
band_waves = fltarr(n_bands)
breaks = fltarr(2,n_ages,n_bands)
ages = fltarr(n_ages)
n_coeff1 = 4
n_coeff2 = 3
fit_coeff1 = fltarr(n_coeff1,n_ages,n_bands)
fit_coeff2 = fltarr(n_coeff2,n_ages,n_bands)

tstr = ''
t1 = fltarr(n_coeff1)
t2 = fltarr(n_coeff2)
t3 = fltarr(2)
tflt = 0.0
tband = ''
FOR i = 0,(n_bands-1) DO BEGIN 
    readf,unit1,tstr
    readf,unit1,tflt,tband,format='(F9.1,2x,A)'
    bands[i] = tband
    band_waves[i] = tflt
    FOR j = 0,(n_ages-1) DO BEGIN
        readf,unit1,tflt
        ages[j] = tflt
        readf,unit1,t3
        breaks[*,j,i] = t3
        readf,unit1,t1
        fit_coeff1[*,j,i] = t1
        readf,unit1,t2
        fit_coeff2[*,j,i] = t2
    ENDFOR 
ENDFOR

free_lun,unit1

; make record to hold details

fit_info = {bands : bands, $
            band_waves : band_waves, $
            ages : ages, $
            breaks : breaks, $
            fit_coeff1 : fit_coeff1, $
            fit_coeff2 : fit_coeff2, $
            n_bands : n_bands, $
            n_ages : n_ages, $
            n_coeff1 : n_coeff1, $
            n_coeff2 : n_coeff2}

END 

; ----------------------------------------------------------------------
