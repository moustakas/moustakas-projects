pro build_deep2_window_completeness
; jm14mar04siena - compute the statistical completeness for DEEP2/DR4
; using the distributed window functions; for details see:
;   http://deep.ps.uci.edu/dr4/completeness.html
;   http://deep.ps.uci.edu/dr4/skycoverage.html

; this code uses the "weights" maps, which it turns out are useless
; because they also include redshift success

    catpath = deep2_path(/catalogs)
    winpath = deep2_path(/window)

; read the full redshift catalog    
    zcat = mrdfits(catpath+'zcat.deep2.dr4.uniq.fits.gz',1)
    ngal = n_elements(zcat)

    out = struct_addtags(struct_trimtags(zcat,select=$
      ['objno','objname','ra','dec']),replicate({weight: -1.0},ngal))
;   weight = rsex(catpath+'deep2_target_weights.dat')
    
; Field 1/EGS; note we select on OBJNO here (not OBJNAME) because
; there are ~6 objects with a pointing prefactor of "6" that we would
; otherwise miss 
    win = mrdfits(winpath+'windowf.egs.fits.gz',0,hdr)
    extast, hdr, astr

    these = where(strmid(strtrim(zcat.objno,2),0,1) eq '1')
    tot += n_elements(these)
    ad2xy, zcat[these].ra, zcat[these].dec, astr, xgal, ygal

;   xim = findgen(astr.naxis[0])#(fltarr(astr.naxis[1])+1)
;   yim = (fltarr(astr.naxis[0])+1)#reverse(findgen(astr.naxis[1]))
    out[these].weight = interpolate(win,xgal,ygal,missing=0.0)

; Fields 2-4; note: no statistical weights for pointing 43, which is
; severely incomplete; also note that we select on OBJNAME due to the
; intricacies of how the masks were designed (see the footnote at this
; site: http://deep.ps.uci.edu/dr4/ztags.html)
    point = ['21','22','31','32','33','41','42']
    for ii = 0, n_elements(point)-1 do begin
       win = mrdfits(winpath+'windowf.'+point[ii]+'.fits.gz',0,hdr)
       extast, hdr, astr

       these = where(strmid(strtrim(zcat.objname,2),0,2) eq point[ii])
       tot += n_elements(these)
       ad2xy, zcat[these].ra, zcat[these].dec, astr, xgal, ygal
       out[these].weight = interpolate(win,xgal,ygal,missing=0.0)
    endfor

; check to make sure all the missing weights are from pointing 43
    ww = where(out.weight eq -1.0)
;   splog, ngal, tot
;   print, zcat[ww].objno
    out[ww].weight = 0.0 ; replace

    im_mwrfits, out, catpath+'weight.zcat.deep2.dr4.uniq.fits', /clobber
    
; also write out files that are line-matched to my catalogs which have
; been processed through deep2_check_spec1d_dr4:
    zcat = read_deep2_zcat(/all)
    match, zcat.objno, out.objno, m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    help, m2

    im_mwrfits, out[m2], catpath+'weight.zcat.dr4.goodspec1d.fits', /clobber
    
    zcat = read_deep2_zcat()
    match, zcat.objno, out.objno, m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    help, m2

    im_mwrfits, out[m2], catpath+'weight.zcat.dr4.goodspec1d.Q34.fits', /clobber

return
end
    
