function init_cat_linefit, ngalaxy=ngalaxy
; jm05jan01uofa - initialize the linefit data structure for the
;                 various compilations of high-redshift data from the
;                 literature
; jm08apr11nyu - modified

    if (n_elements(ngalaxy) eq 0L) then ngalaxy = 1L
    if (ngalaxy eq 0L) then ngalaxy = 1L

    moretags = {id: 0L, galaxy: '', alt_galaxy: '', ra: 0.0D, dec: 0.0D, $
      z: -999.0, m_b: -999.0, m_b_err: -999.0, m_g: -999.0, m_g_err: -999.0, $
      m_r: -999.0, m_r_err: -999.0, mass: -999.0, r23branch: '', $
      lit_log12oh: -999.0, lit_log12oh_err: -999.0}
    linename = [$
      'oii_3727','neiii_3869','h_delta','h_gamma','oiii_4363','h_beta',$
      'oiii_4959','oiii_5007','oi_6300','nii_6548','h_alpha','nii_6584',$
      'sii_6716','sii_6731']
    linewave = [$
      3727.4235,3869.06,4101.734,4340.464,4363.21,4861.325,$
      4958.911,5006.843,6300.304,6548.04,6562.8,6583.46,$
      6716.14,6730.81]
    nline = n_elements(linename)
;   niceprint, linename, linewave

    data = struct_addtags(replicate(moretags,ngalaxy),$
      replicate({linename: linename, linewave: linewave},ngalaxy))
    data = struct_addtags(data,mrd_struct(linename,replicate('[0.0,-2.0]',nline),ngalaxy))
    data = struct_addtags(data,mrd_struct(linename+'_wave',strtrim(linewave,2),ngalaxy))
    data = struct_addtags(data,mrd_struct(linename+'_ew',replicate('[0.0,-2.0]',nline),ngalaxy))
    data = struct_addtags(data,mrd_struct(linename+'_continuum',replicate('[0.0,-2.0]',nline),ngalaxy))

return, data
end
