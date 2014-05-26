function deep2_colorcut, magb, magr, magi, mlim=mlim
    if n_elements(mlim) eq 0 then mlim = 24.1
    return, magr lt mlim and $
      ((magb-magr lt 2.45*(magr-magi)-0.5) or $
      (magr-magi gt 1.1) or (magb-magr lt 0.5))
;   return, magr lt mlim and $
;     ((magb-magr lt 2.45*(magr-magi)-0.311) or $
;     (magr-magi gt 1.211) or (magb-magr lt 0.389))
end

pro build_deep2_targeting_weights
; jm14may21siena - compute the targeting weights for DEEP2/DR4
; for details see: Newman+13 and
; http://deep.ps.uci.edu/dr4/completeness.html 

    catpath = deep2_path(/catalogs)
    winpath = deep2_path(/window)

; read the full redshift catalog    
    zcat = mrdfits(catpath+'zcat.deep2.dr4.uniq.fits.gz',1)
    ngal = n_elements(zcat)

; obj_weight - object weight
; targ_weight - probability of being selected for observation (final
;   targeting weight) 
    
    out = struct_addtags(struct_trimtags(zcat,select=$
      ['objno','objname','ra','dec','pgal','magr']),replicate({$
      obj_weight: -1.0, targ_weight: -1.0, $
      deep2_nmatch: 0, deep2_psel: -1.0, deep2_probcut: -1.0, $
      deep2_weight: -1.0},ngal))
    
; ####################
; Field 1/EGS targeting weights: W = Wsg * WR; note we
; select on OBJNO here (not OBJNAME) because there are ~6 objects with
; a pointing prefactor of "6" that we would otherwise miss
    field1 = where(strmid(strtrim(zcat.objno,2),0,1) eq '1',nobj)
    zcat1 = zcat[field1]
    magb = zcat1.magb
    magr = zcat1.magr
    magi = zcat1.magi

    wsg = fltarr(nobj)          ; star-galaxy weight    
    wsg[where(zcat1.pgal ge 2)] = 1.0
    ww = where(zcat1.pgal lt 2 and zcat1.pgal gt 0.2)
    wsg[ww] = zcat1[ww].pgal<1.0

    wc = deep2_colorcut(magb,magr,magi)        ; color weight
    wr = (0.75*10^(-0.4*(magr-24.1)))<1        ; magnitude weight

    wr_final = fltarr(nobj)
    w1 = where((magr le 21.5) or (wc eq 1),n1)
    w2 = where((magr gt 21.5) and (wc eq 0),n1)
    wr_final[w1] = wr[w1]
    wr_final[w2] = (0.1*10^(-0.4*(magr[w2]-24.1)))<1

; object weight and final targeting weight
    out[field1].obj_weight = wsg*wr_final
    out[field1].targ_weight = 0.33398 + 0.42687*out[field1].obj_weight
;   djs_plot, magr, out[field1].targ_weight, psym=3, xsty=3, ysty=3

; ####################
; Fields 2-4 targeting weights: W = Wsg * Wc * WR * Wbc
    fields24 = where(strmid(strtrim(zcat.objno,2),0,1) ne '1',nobj)
    zcat24 = zcat[fields24]
    magb = zcat24.magb
    magr = zcat24.magr
    magi = zcat24.magi

    wsg = fltarr(nobj)          ; star-galaxy weight    
    wsg[where(zcat24.pgal ge 2)] = 1.0
    ww = where(zcat24.pgal lt 2 and zcat24.pgal gt 0.2)
    wsg[ww] = zcat24[ww].pgal<1.0

    wc = deep2_colorcut(magb,magr,magi)        ; color weight
    wr = (0.75*10^(-0.4*(magr-24.1)))<1        ; magnitude weight
    wbc = (poly(magr-magi,[0.0,2.2222])>0.1)<1 ; "blue color" weight

;   djs_plot, magr-magi, magb-magr, psym=3
;   ww = where(wc eq 0)
;   djs_oplot, magr[ww]-magi[ww], magb[ww]-magr[ww], psym=3, color='red'

; final weight; don't apply the color weight because it has already
; been incorporated into the redshift catalog for these fields
    out[fields24].obj_weight = wsg*wr*wbc
;   out[fields24].obj_weight = wsg*wc*wr*wbc
    out[fields24].targ_weight = 0.27976 + 0.44717*out[fields24].obj_weight-$
      0.09137*out[fields24].obj_weight^2

;   djs_plot, magr, out[fields24].targ_weight, psym=3, xsty=3, ysty=3

;;; check!    
;;    deep2 = rsex(catpath+'deep2_selection.dat')
;;;   match, out.objno, deep2.objno, m1, m2
;;;   out[m1].deep2_weight = deep2[m2].weight
;;    mult = 0
;;    for ii = 0, ngal-1 do begin
;;       match = where(out[ii].objno eq deep2.objno,nmatch)
;;       out[ii].deep2_nmatch = nmatch
;;;      if nmatch eq 0 then splog, 'No match for ', out[ii].objno
;;       if nmatch eq 1 then begin
;;          out[ii].deep2_psel = deep2[match].psel
;;          out[ii].deep2_probcut = deep2[match].prob_cut
;;          out[ii].deep2_weight = deep2[match].weight
;;       endif
;;       if nmatch gt 1 then begin
;;          splog, 'Multiple matches for ', out[ii].objno
;;          mult++
;;          struct_print, deep2[match]
;;          this = where(deep2[match].weight gt 0)
;;          out[ii].deep2_psel = deep2[match[this[0]]].psel
;;          out[ii].deep2_probcut = deep2[match[this[0]]].prob_cut
;;          out[ii].deep2_weight = deep2[match[this[0]]].weight
;;       endif
;;    endfor
    
; write out
    im_mwrfits, out, catpath+'weight.zcat.deep2.dr4.uniq.fits', /clobber

; also write out files that are line-matched to my catalogs which have
; been processed through deep2_check_spec1d_dr4:
    zcat = read_deep2_zcat(/all)
    match, zcat.objno, out.objno, m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    help, zcat, m2

    im_mwrfits, out[m2], catpath+'weight.zcat.dr4.goodspec1d.fits', /clobber
    
    zcat = read_deep2_zcat()
    match, zcat.objno, out.objno, m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    help, zcat, m2

    im_mwrfits, out[m2], catpath+'weight.zcat.dr4.goodspec1d.Q34.fits', /clobber

return
end
    
