function read_bcgmstar_sample, zsort=zsort
; jm13oct19siena - read the sample

    propath = bcgmstar_path(/propath)
    
;   sample = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    sample = rsex(bcgmstar_path(/propath)+'bcgmstar_sample.sex')
    ncl = n_elements(sample)

; add the M(500) masses from Merten+14 and from Donahue+14; note that
; Merten+14 does not give masses for Abell1423 and MACS2129; side
; note: for the A1423 virial mass Lemze+14 gets 0.385E15+/-0.078E15
    sample = struct_addtags(sample,replicate({m500_merten: 0.0, m500_merten_err: 0.0, $
      m500_xray: 0.0, m500_xray_err: 0.0, m500: 0.0, m500_err: 0.0},ncl))
    
    merten = rsex(propath+'merten_lensing_mass.sex')
    match, merten.cluster, strupcase(sample.shortname), m1, m2
    sample[m2].m500_merten = alog10(merten[m1].m_500*1D15)
    sample[m2].m500_merten_err = merten[m1].delta_m_500/merten[m1].m_500/alog(10)

; Donahue+14: use the vmax and radii from Table 4 and equation 5 to
; get M(500)
    xray = mrdfits(propath+'xray_mass_parsed.fits.gz',1)
    match, strtrim(xray.cluster,2), strtrim(strupcase(sample.shortname),2), m1, m2
    sample[m2].m500_xray = alog10(xray[m1].m500)
    sample[m2].m500_xray_err = xray[m1].m500_err/xray[m1].m500/alog(10)

; use the X-ray M(500) where Merten's values are missing
    sample.m500 = sample.m500_merten
    sample.m500_err = sample.m500_merten_err
    zero = where(sample.m500 eq 0.0)
    sample[zero].m500 = sample[zero].m500_xray
    sample[zero].m500_err = sample[zero].m500_xray_err
    
;   splog, 'HACK!!!!!'
;   sample = sample[14]
    if keyword_set(zsort) then sample = sample[sort(sample.z)]
;   struct_print, sample

return, sample
end
