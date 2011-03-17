pro write_flamingos_ahead_global

    xdist_coeff = [$
      -3.13385180345,$
      1.004072924,$
      -0.00183808434959,$
      1.90899554759D-06,$
      -7.02600909987D-07,$
      1.31772038336D-06,$
      -1.46253633834D-08,$
      8.74646975868D-10,$
      -1.27477572409D-08,$
      3.84044359946D-10,$
      -9.23163130406D-13,$
      -3.02699688831D-13,$
      -2.17804944819D-12,$
      7.13140243591D-13,$
      -4.22992633579D-13,$
      7.50779676661D-15,$
      -1.05067728167D-16,$
      1.32865593246D-14,$
      -1.39547603335D-15,$
      4.79680525157D-15,$
      5.27901629365D-16,$
      -5.28209885681D-19,$
      -2.15829018477D-19,$
      -1.69082032126D-19,$
      2.98401155856D-19,$
      1.42989656823D-19,$
      -2.34978528457D-19,$
      -1.4139870516D-19]

    ydist_coeff = [$
      -1.59317700345,$
      -7.58792414804D-05,$
      1.00403420094,$
      -2.87917603493D-08,$
      1.4431725501D-08,$
      -4.00634996454D-06,$
      -2.46081186836D-10,$
      -1.53044126738D-08,$
      8.91239390236D-10,$
      -1.53282151607D-08,$
      -3.25515890748D-14,$
      -3.40198566384D-13,$
      2.21451114579D-12,$
      -8.13094365249D-13,$
      3.39460425228D-12,$
      2.9734493965D-16,$
      7.34769704289D-15,$
      2.74279596123D-16,$
      1.33306619088D-14,$
      -1.57274385218D-15,$
      7.21586951996D-15,$
      5.0088021963D-19,$
      3.53436803775D-19,$
      -4.51355862217D-19,$
      -1.471392533D-19,$
      2.80336330843D-19,$
      -2.1021589111D-19,$
      -3.01518328277D-19]

    mkhdr, hdr, 0

    extast, headfits('ks.fits'), astr
    putast, hdr, astr
    sxdelpar, hdr, 'COMMENT'
    sxdelpar, hdr, 'DATE'
    sxdelpar, hdr, 'HISTORY'
    sxdelpar, hdr, 'SIMPLE'
    sxdelpar, hdr, 'BITPIX'
    sxdelpar, hdr, 'NAXIS'
    sxdelpar, hdr, 'EXTEND'
    sxdelpar, hdr, 'CTYPE1'
    sxdelpar, hdr, 'CTYPE2'
    sxdelpar, hdr, 'CRVAL1'
    sxdelpar, hdr, 'CRVAL2'
    sxdelpar, hdr, 'LONPOLE'
    sxdelpar, hdr, 'LATPOLE'
    sxdelpar, hdr, 'PV2_1'
    sxdelpar, hdr, 'PV2_2'
    sxdelpar, hdr, 'EQUINOX'

    for ix = 0L, n_elements(xdist_coeff)-1L do sxaddpar, hdr, 'PV1_'+$
      string(ix,format='(I0)'), xdist_coeff[ix]
    for iy = 0L, n_elements(xdist_coeff)-1L do sxaddpar, hdr, 'PV2_'+$
      string(iy,format='(I0)'), ydist_coeff[iy]
    
;   niceprint, hdr
    openw, lun, 'flamingos.ahead', /get_lun
    for ii = 0L, n_elements(hdr)-1L do printf, lun, hdr[ii]
    free_lun, lun
    
stop    

return
end
