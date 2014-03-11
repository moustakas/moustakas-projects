function deep2_get_oiiflux, ppxf, kised, oiierr=oiierr, $
  oiilimit=oiilimit, snrcut=snrcut
; given the ppxf and isedfit/k-correct data structures, compute the
; total [OII] flux and upper limits

    if n_elements(snrcut) eq 0 then snrcut = 1.5
    
    ngal = n_elements(ppxf)
    oii = fltarr(ngal)-1.0
    oiierr = fltarr(ngal)-1.0
    oiilimit = intarr(ngal)
    
    both = where(ppxf.oii_3727_1[1] gt 0 and ppxf.oii_3727_2[1] gt 0 and $
      ppxf.oii_3727_1[0]/ppxf.oii_3727_1[1] gt snrcut and $
      ppxf.oii_3727_2[0]/ppxf.oii_3727_2[1] gt snrcut,nboth)
    oii[both] = (ppxf[both].oii_3727_1_ew[0]+ppxf[both].oii_3727_2_ew[0])*$
      kised[both].cflux_3727
    oiierr[both] = sqrt(ppxf[both].oii_3727_1_ew[1]^2+ppxf[both].oii_3727_2_ew[1]^2)*$
      kised[both].cflux_3727

    one = where(ppxf.oii_3727_1[1] gt 0 and ppxf.oii_3727_2[1] eq -1 and $
      ppxf.oii_3727_1[0]/ppxf.oii_3727_1[1] gt snrcut,none)
    oii[one] = ppxf[one].oii_3727_1_ew[0]*kised[one].cflux_3727
    oiierr[one] = ppxf[one].oii_3727_1_ew[1]*kised[one].cflux_3727

    two = where(ppxf.oii_3727_1[1] eq -1 and ppxf.oii_3727_2[1] gt 0 and $
      ppxf.oii_3727_2[0]/ppxf.oii_3727_2[1] gt snrcut,ntwo)
    oii[two] = ppxf[two].oii_3727_2_ew[0]*kised[two].cflux_3727
    oiierr[two] = ppxf[two].oii_3727_2_ew[1]*kised[two].cflux_3727

    neither = where(ppxf.oii_3727_1[1] eq -2 or ppxf.oii_3727_2[1] eq -2,nneither)
    oii[neither] = -2.0
    oiierr[neither] = -2.0

; upper limits  
    lim = where(oii eq -1.0,nlim)
    oiilimit[lim] = 1
    oii[lim] = (ppxf[lim].oii_3727_1_ew_limit+ppxf[lim].oii_3727_2_ew_limit)*kised[lim].cflux_3727
    
    splog, ngal, nboth, none, ntwo, nneither, nboth+none+ntwo+nneither
    
return, oii
end

