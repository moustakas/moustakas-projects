; ----------------------------------------------------------------------
;   fr2atten.pro
;
; Created to utilized the flux ratio method to compute the dust 
; attenuation of (starburst) galaxies at different wavelengths.  The
; details of this method are contained in Gordon K. D., Clayton,
; G. C., Witt, A. N., & Misselt, K. A. 2000, ApJ, 533, 236.  
;
; v 1.0 - Karl D. Gordon (Steward Obs., 9 June 2000)
; v 1.1 - Karl D. Gordon (Steward Obs., 3 Oct 2000)
;    fixed the error in that fit is done to fr + 1 and so targ_fr
;    needs to have 1 added before being sent to the fit
; ----------------------------------------------------------------------

FUNCTION fr2atten_calc,fr,fit_info,age_k,met_k,ad_k,band_k

; add one to fr as the fit was done to fr+1

in_fr = fr
fr = fr + 1.0

IF (fr LE fit_info.breaks[1,age_k,met_k,ad_k,band_k]) THEN BEGIN
    atten = total(fit_info.fit_coeff1[*,age_k,met_k,ad_k,band_k]* $
                  fr^indgen(fit_info.n_coeff1+1)) 
ENDIF ELSE IF (fr LE fit_info.breaks[0,age_k,met_k,ad_k,band_k]) THEN BEGIN 
    width = fit_info.breaks[0,age_k,met_k,ad_k,band_k] - $
      fit_info.breaks[1,age_k,met_k,ad_k,band_k]
    a1_weight = (fit_info.breaks[0,age_k,met_k,ad_k,band_k] - fr)/width
    a2_weight = 1.0 - a1_weight
    atten = a1_weight*total(fit_info.fit_coeff1[*,age_k,met_k,ad_k,band_k]* $
                            fr^indgen(fit_info.n_coeff1+1)) + $
      a2_weight*total(fit_info.fit_coeff2[*,age_k,met_k,ad_k,band_k]* $
                      alog10(fr)^indgen(fit_info.n_coeff2+1))
ENDIF ELSE BEGIN 
    atten = total(fit_info.fit_coeff2[*,age_k,met_k,ad_k,band_k]* $
                  alog10(fr)^indgen(fit_info.n_coeff2+1))
ENDELSE

; reset fr to the input fr

fr = in_fr

return,atten

END

; ----------------------------------------------------------------------

FUNCTION fr2atten, targ_fr, targ_band, targ_wave, targ_age, targ_met, $
  targ_ad, fit_info, silent=silent

npar = n_params(0)
if (npar eq 0) then begin
    print,"attenuation = fr2atten(flux_ratio,band,wavelength,age," + $
      "metallicity,a_d,"
    print,"                       fit_info,/silent)"
    print,""
    print,"inputs: flux_ratio [float]  = measured flux ratio as defined in "
    print,"           Gordon et al. 2000, ApJ, 533, 236"
    print,"        band [string]       = name of band " + $
      "[mainly for use with H-alpha = Ha]"
    print,"        wavelength [float]  = wavelength of band in Angstroms"
    print,"        age [float]         = age of SES model [in Myrs]"
    print,"        metallicity [float] = metallicity of SES model " + $
      "[solar = 0.0]"
    print,"        a_d [float]         = fraction of Lyman " + $
      "continuum (ionizing) "
    print,"                photons absorbed by the dust in H II regions"
    print,"        fit_info [struct]   = structure returned by the IDL" + $
      " procedure "
    print,"                              read_many_fit.pro"
    print,""
    print,"keywords: silent = set to supress messages"
retall & end


; determine the two fits with ages braketting the target age

IF (targ_age LT min(fit_info.ages)) THEN targ_age = min(fit_info.ages)
IF (targ_age GE max(fit_info.ages)) THEN targ_age = max(fit_info.ages)*0.99

age_k = 0
WHILE (fit_info.ages[age_k] LE targ_age) DO age_k = age_k + 1
age2_weight = (targ_age - fit_info.ages[age_k-1])/ $
  (fit_info.ages[age_k] - fit_info.ages[age_k-1])
age1_weight = 1.0 - age2_weight

; determine the two metallicities breaketting the target metallicity

IF (targ_met LT min(fit_info.mets)) THEN targ_met = min(fit_info.mets)
IF (targ_met GE max(fit_info.mets)) THEN targ_met = max(fit_info.mets)*0.99

met_k = 0
WHILE (fit_info.mets[met_k] LE targ_met) DO met_k = met_k + 1
met2_weight = (targ_met - fit_info.mets[met_k-1])/ $
  (fit_info.mets[met_k] - fit_info.mets[met_k-1])
met1_weight = 1.0 - met2_weight

; determine the two ad values braketting the target ad

IF (targ_ad LT min(fit_info.ads)) THEN targ_ad = min(fit_info.ads)
IF (targ_ad GE max(fit_info.ads)) THEN targ_ad = max(fit_info.ads)*0.99
ad_k = 0
WHILE (fit_info.ads[ad_k] LE targ_ad) DO ad_k = ad_k + 1
ad2_weight = (targ_ad - fit_info.ads[ad_k-1])/ $
  (fit_info.ads[ad_k] - fit_info.ads[ad_k-1])
ad1_weight = 1.0 - ad2_weight

; if the band is in the list, use it

indxs = where(fit_info.bands EQ targ_band,n_indxs)
IF (n_indxs GT 0) THEN BEGIN

    band_k = indxs[0]
    atten1 = age1_weight*fr2atten_calc(targ_fr,fit_info,age_k-1,met_k-1, $
                                       ad_k-1,band_k) + $
      age2_weight*fr2atten_calc(targ_fr,fit_info,age_k,met_k-1,ad_k-1,band_k)
    atten2 = age1_weight*fr2atten_calc(targ_fr,fit_info,age_k-1,met_k, $
                                       ad_k-1,band_k) + $
      age2_weight*fr2atten_calc(targ_fr,fit_info,age_k,met_k,ad_k-1,band_k)

    atten_ad1 = met1_weight*atten1 + met2_weight*atten2

    atten1 = age1_weight*fr2atten_calc(targ_fr,fit_info,age_k-1,met_k-1, $
                                       ad_k,band_k) + $
      age2_weight*fr2atten_calc(targ_fr,fit_info,age_k,met_k-1,ad_k,band_k)
    atten2 = age1_weight*fr2atten_calc(targ_fr,fit_info,age_k-1,met_k, $
                                       ad_k,band_k) + $
      age2_weight*fr2atten_calc(targ_fr,fit_info,age_k,met_k,ad_k,band_k)

    atten_ad2 = met1_weight*atten1 + met2_weight*atten2

    atten = ad1_weight*atten_ad1 + ad2_weight*atten_ad2

; otherwise use the nearest two bands

ENDIF ELSE BEGIN

    k = 0
    while (targ_wave GT fit_info.band_waves[k]) do begin
        k = k + 1
    endwhile
    k1 = k - 1
    k2 = k
    wave2_weight = (targ_wave - fit_info.band_waves[k1])/ $
      (fit_info.band_waves[k2] - fit_info.band_waves[k1])
    wave1_weight = 1.0 - wave2_weight

    if (not keyword_set(silent)) then begin
        print,'interpolating between bands (wavelengths)'
        print,'1 = ' + fit_info.bands[k1] + ' (' + $
          strtrim(string(fit_info.band_waves[k1]),2) + ')'
        print,'2 = ' + fit_info.bands[k2] + ' (' + $
          strtrim(string(fit_info.band_waves[k2]),2) + ')'
    endif

; get the attten for the first band

    band_k = k1
    atten1 = age1_weight*fr2atten_calc(targ_fr,fit_info,age_k-1,met_k-1, $
                                       ad_k-1,band_k) + $
      age2_weight*fr2atten_calc(targ_fr,fit_info,age_k,met_k-1,ad_k-1,band_k)
    atten2 = age1_weight*fr2atten_calc(targ_fr,fit_info,age_k-1,met_k, $
                                       ad_k-1,band_k) + $
      age2_weight*fr2atten_calc(targ_fr,fit_info,age_k,met_k,ad_k-1,band_k)

    atten_ad1 = met1_weight*atten1 + met2_weight*atten2

    atten1 = age1_weight*fr2atten_calc(targ_fr,fit_info,age_k-1,met_k-1, $
                                       ad_k,band_k) + $
      age2_weight*fr2atten_calc(targ_fr,fit_info,age_k,met_k-1,ad_k,band_k)
    atten2 = age1_weight*fr2atten_calc(targ_fr,fit_info,age_k-1,met_k, $
                                       ad_k,band_k) + $
      age2_weight*fr2atten_calc(targ_fr,fit_info,age_k,met_k,ad_k,band_k)

    atten_ad2 = met1_weight*atten1 + met2_weight*atten2

    atten_band1 = ad1_weight*atten_ad1 + ad2_weight*atten_ad2

; get the attten for the second band

    band_k = k2
    atten1 = age1_weight*fr2atten_calc(targ_fr,fit_info,age_k-1,met_k-1, $
                                       ad_k-1,band_k) + $
      age2_weight*fr2atten_calc(targ_fr,fit_info,age_k,met_k-1,ad_k-1,band_k)
    atten2 = age1_weight*fr2atten_calc(targ_fr,fit_info,age_k-1,met_k, $
                                       ad_k-1,band_k) + $
      age2_weight*fr2atten_calc(targ_fr,fit_info,age_k,met_k,ad_k-1,band_k)

    atten_ad1 = met1_weight*atten1 + met2_weight*atten2

    atten1 = age1_weight*fr2atten_calc(targ_fr,fit_info,age_k-1,met_k-1, $
                                       ad_k,band_k) + $
      age2_weight*fr2atten_calc(targ_fr,fit_info,age_k,met_k-1,ad_k,band_k)
    atten2 = age1_weight*fr2atten_calc(targ_fr,fit_info,age_k-1,met_k, $
                                       ad_k,band_k) + $
      age2_weight*fr2atten_calc(targ_fr,fit_info,age_k,met_k,ad_k,band_k)

    atten_ad2 = met1_weight*atten1 + met2_weight*atten2

    atten_band2 = ad1_weight*atten_ad1 + ad2_weight*atten_ad2

; conbine the different band attenuations

    atten = wave1_weight*atten_band1 + wave2_weight*atten_band2

ENDELSE 

return,atten

END 
