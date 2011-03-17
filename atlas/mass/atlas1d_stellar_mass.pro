;+
; NAME:
;       ATLAS1D_STELLAR_MASS
;
; PURPOSE:
;       Estimate the stellar masses of the K/M Atlas galaxies using
;       SED fitting. 
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
;       J. Moustakas, 2005 Mar 03, U of A
;-

pro atlas1d_stellar_mass, atlas, _extra=extra

    if (n_elements(atlas) eq 0L) then atlas = read_integrated()
    ngalaxy = n_elements(atlas)

    datapath = atlas_path(/analysis)+'mass/'

    prefix = 'atlas1d_UBVJHKs'
    filterlist = ['bessell_U','bessell_B','bessell_V','twomass_J','twomass_H','twomass_Ks']+'.par'
    filtinfo = im_filterspecs(filterlist=filterlist)
    nfilter = n_elements(filterlist)

    H0 = 70.0                   ; [km/s/Mpc]
    light = 2.99792458D5        ; speed of light [km/s]
    
; generate the model magnitudes for BwRIK and BwRI

;;    isedfit_models, filterlist, datapath=datapath, prefix=prefix, $
;;      minredshift=0.01, maxredshift=0.08, deltaredshift=0.01, $
;;      minebv=0.0, maxebv=1.0, deltaebv=0.1, minwave=2900.0, $
;;      _extra=extra, /write
;;
;;stop
    
; convert the photometry in each band to MAGGIES; set the inverse
; variance in objects with no photometry to a large number 

    vegamags = transpose([$
      [atlas.U],$
      [atlas.B],$
      [atlas.V],$
      [atlas.twomass_J],$
      [atlas.twomass_H],$
      [atlas.twomass_Ks]])

    vegamags_err = transpose([$
      [atlas.U_err],$
      [atlas.B_err],$
      [atlas.V_err],$
      [atlas.twomass_J_err],$
      [atlas.twomass_H_err],$
      [atlas.twomass_Ks_err]])

    ABmags = vegamags + rebin(filtinfo.vega2ab,nfilter,ngalaxy)
    ABmags_err = vegamags_err

    maggies        = 10.0^(-0.4D*(ABmags>0.0))
    maggies_err    = 0.4D*alog(10.0)*maggies*(ABmags_err>0.0)

    for ifilter = 0L, nfilter-1L do begin
       bad = where(vegamags[ifilter,*] lt -900.0,nbad)
       if (nbad ne 0L) then begin
          maggies[ifilter,bad] = 0.0D
          maggies_err[ifilter,bad] = 1D16
       endif
    endfor

    maggies_invvar = 1.0/maggies_err^2

    zobj = atlas.distance*H0/light
    galaxy = strlowcase(strtrim(atlas.galaxy,2))

; do the fit and generate QA plots

    isedfit, float(maggies[*,0:1]), float(maggies_invvar[*,0:1]), zobj[0:1], result, $
      filterlist=filterlist, galaxy=galaxy[0:1], datapath=datapath, prefix=prefix, $
      _extra=extra
    
;   isedfit, float(maggies), float(maggies_invvar), zobj, result, $
;     galaxy=galaxy, datapath=datapath, prefix=prefix, _extra=extra

stop    
    
    atlas1d_isedfit_qaplot, atlas, datapath=datapath, prefix=prefix


    
stop    
    
return
end
    
