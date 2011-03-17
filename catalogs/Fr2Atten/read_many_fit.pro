PRO read_many_fit,fit_info,sf_type=sf_type, path=path
; jm05feb21uofa - added variable datapath
    
IF (NOT keyword_set(sf_type)) THEN sf_type = 'c'

if (n_elements(path) eq 0L) then path = filepath('',root=getenv('CATALOGS_DIR'),subdirectory='Fr2Atten')
;if (n_elements(path) eq 0L) then path = './'
;path = '~/Pro/Fr2Atten/'
mod_val = ['m23','m17','m07','m04','p00','p04','p07']
n_mod_vals = n_elements(mod_val)
mets = [-2.3,-1.7,-0.7,-0.4,0.0,0.4,0.7]
n_mets = n_mod_vals
ad_val = [0.01,0.25,0.50,0.75,0.95]
n_ad_vals = n_elements(ad_val)
FOR i = 0,(n_mod_vals-1) DO BEGIN
    FOR j = 0,(n_ad_vals-1) DO BEGIN 
        read_fit,path+'fr_A_' + mod_val[i]+'_' + sf_type + '_' + $
          strtrim(string(ad_val[j],format='(F4.2)'),2)+'.dat',ofit_info
        IF ((i EQ 0) AND (j EQ 0)) THEN BEGIN
            fit_coeff1 = fltarr(ofit_info.n_coeff1,ofit_info.n_ages,n_mets, $
                                n_ad_vals,ofit_info.n_bands)
            fit_coeff2 = fltarr(ofit_info.n_coeff2,ofit_info.n_ages,n_mets, $
                                n_ad_vals,ofit_info.n_bands)
            breaks = fltarr(2,ofit_info.n_ages,n_mets,n_ad_vals, $
                            ofit_info.n_bands)
        ENDIF
        fit_coeff1[*,*,i,j,*] = ofit_info.fit_coeff1
        fit_coeff2[*,*,i,j,*] = ofit_info.fit_coeff2
        breaks[*,*,i,j,*] = ofit_info.breaks
    ENDFOR 
ENDFOR

; make record to hold details

fit_info = {bands : ofit_info.bands, $
            band_waves : ofit_info.band_waves, $
            ages : ofit_info.ages, $
            mets : mets, $
            ads : ad_val, $
            breaks : breaks, $
            fit_coeff1 : fit_coeff1, $
            fit_coeff2 : fit_coeff2, $
            n_bands : ofit_info.n_bands, $
            n_ages : ofit_info.n_ages, $
            n_mets : n_mets, $
            n_ads : n_ad_vals, $
            n_coeff1 : ofit_info.n_coeff1, $
            n_coeff2 : ofit_info.n_coeff2}

END 
