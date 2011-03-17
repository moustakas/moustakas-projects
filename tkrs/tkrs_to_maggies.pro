;+
; NAME:
;   tkrs_to_maggies
; PURPOSE:
;   convert TKRS catalog input to Galactic-extcintion corrected AB maggies 
; CALLING SEQUENCE:
;   tkrs_to_maggies,tkrs,maggies,ivar
; INPUTS:
;   tkrs - [N] TKRS style input
;               .RA (J2000 degrees)
;               .DEC (J2000 degrees)
;               .BMAG_MAGAUTO
;               .BMAGERR_MAGAUTO
;               .VMAG_MAGAUTO
;               .VMAGERR_MAGAUTO
;               .IMAG_MAGAUTO
;               .IMAGERR_MAGAUTO
;               .ZMAG_MAGAUTO
;               .ZMAGERR_MAGAUTO
;               .JMAG_MAGAUTO
;               .JMAGERR_MAGAUTO
;               .HMAG_MAGAUTO
;               .HMAGERR_MAGAUTO
;               .KMAG_MAGAUTO
;               .KMAGERR_MAGAUTO
; OUTPUTS:
;   maggies - [7, N] output in AB maggies in BVizJHK
;   ivar - [7, N] inverse variance of maggies
; COMMENTS:
;   It ALWAYS applies a minimum error of 0.02 in all
;   bandpasses. Except in H band it ALWAYS sets the inverse variance
;   to zero since I don't believe the H band calibration.
;
;   Requires you to have the dust maps so that dust_getval can find
;   them. (If somebody wants me to set "default" columns in the
;   tkrs structure that this code looks for, let me know).
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
; jm08apr24nyu - identical to GOODS_TO_MAGGIES, but converts the
;                photometry to AB first
;-
;------------------------------------------------------------------------------
pro tkrs_to_maggies, tkrs, maggies, ivar, useh=useh

minerrors=replicate(0.02, 7)
dfactors=[4.32, 3.32, 2.00, 1.54, 0.90, 0.58, 0.37]
names=['B','V','I','Z','J','H','K']

filterlist = 'goods_'+['acs_'+['f435w','f606w','f775w','f850lp'],$
  ['J','H','Ks']+'_isaac_etc']+'.par'
vega2ab=k_vega2ab(filterlist=filterlist, /kurucz,/silent)*0.0 ; already in AB

;; get Galactic extinction
euler,tkrs.ra,tkrs.dec,ll,bb,1
ebv=dust_getval(ll, bb, /interp, /noloop)*0.0 ; already extinction corrected

maggies=fltarr(7,n_elements(tkrs))
ivar=fltarr(7,n_elements(tkrs))
for iband=0L, 6L do begin
    itag_m=tag_indx(tkrs[0], names[iband]+'mag_magauto')
    itag_msig=tag_indx(tkrs[0], names[iband]+'magerr_magauto')
    indx=where(tkrs.(itag_m) gt 0. AND $
               tkrs.(itag_msig) ge 0. AND $
               tkrs.(itag_m) ne 99. and $
               tkrs.(itag_msig) ne 99. and $
               tkrs.(itag_m) ne -99. and $
               tkrs.(itag_msig) ne -99., count)
    if(count gt 0) then begin
        maggies[iband, indx]=10.^(-0.4*(tkrs[indx].(itag_m) $
                                        -ebv[indx]*dfactors[iband])+vega2ab[iband])
        sig=tkrs[indx].(itag_msig)> 0.001
        ivar[iband, indx]= $
          1./(0.4*alog(10.)*(maggies[iband,indx]*sig))^2
    endif
endfor

k_minerror, maggies, ivar, minerrors

if(NOT keyword_set(useh)) then $
  ivar[5,*]=0.

end
