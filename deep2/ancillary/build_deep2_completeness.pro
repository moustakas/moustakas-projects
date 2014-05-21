function deep2_colorcut, magb, magr, magi, mlim=mlim
    if n_elements(mlim) eq 0 then mlim = 24.1
    return, magr lt mlim and $
      ((magb-magr lt 2.45*(magr-magi)-0.5) or $
      (magr-magi gt 1.1) or (magb-magr lt 0.5))
;   return, magr lt mlim and $
;     ((magb-magr lt 2.45*(magr-magi)-0.311) or $
;     (magr-magi gt 1.211) or (magb-magr lt 0.389))
end

pro build_deep2_completeness
; jm14may21siena - compute the targeting weights for DEEP2/DR4
; for details see: Newman+13 and
; http://deep.ps.uci.edu/dr4/completeness.html 

    catpath = deep2_path(/catalogs)
    winpath = deep2_path(/window)

; read the full redshift catalog    
    zcat = mrdfits(catpath+'zcat.deep2.dr4.uniq.fits.gz',1)
    ngal = n_elements(zcat)

    out = struct_addtags(struct_trimtags(zcat,select=$
      ['objno','objname','ra','dec']),replicate({targ_weight: -1.0},ngal))
    
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

; final weight
    out[field1].targ_weight = wsg*wr_final
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
    out[fields24].targ_weight = wsg*wr*wbc
;   out[fields24].targ_weight = wsg*wc*wr*wbc
;   djs_plot, magr, out[fields24].targ_weight, psym=3, xsty=3, ysty=3

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
    
