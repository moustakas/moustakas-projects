pro sings_electronic_tables
; jm10jul27ucsd - make electronic versions of various tables for
; distribution
    
    version = sings_log12oh_version()
    outpath = sings_path(/projects)+'log12oh/'

    final1 = mrdfits(outpath+'sings_log12oh_final_'+version+'.fits.gz',1)
    result1 = mrdfits(outpath+'sings_log12oh_'+version+'.fits.gz',1)
    hii1 = mrdfits(outpath+'sings_log12oh_hiiregions_'+version+'.fits.gz',1)

; ###########################################################################
; TABLE 8 - nuclear, circumnuclear, and radial-strip abundances
    kk04select = ['nice_galaxy','r23_branch','hii_kk04_slope','hii_kk04_log12oh_nuclear',$
      'hii_kk04_log12oh_char','hii_kk04_log12oh_avg','hii_kk04_nhii_used']
    pt05select = ['nice_galaxy','r23_branch','hii_pt05_slope','hii_pt05_log12oh_nuclear',$
      'hii_pt05_log12oh_char','hii_pt05_log12oh_avg','hii_pt05_nhii_used']
    kk04 = im_struct_trimtags(hii1[where(hii1.hii_kk04_slope[1] gt -900.0)],$
      select=kk04select,newtags=repstr(repstr(repstr($
      kk04select,'hii_',''),'_nused','_nhii_used'),'kk04_',''))
    pt05 = im_struct_trimtags(hii1[where(hii1.hii_pt05_slope[1] gt -900.0)],$
      select=pt05select,newtags=repstr(repstr(repstr($
      pt05select,'hii_',''),'_nused','_nhii_used'),'pt05_',''))
    nkk04 = n_elements(kk04)
    npt05 = n_elements(pt05)
    result = [kk04,pt05]

    result = struct_addtags(result,replicate({calibration: ''},nkk04+npt05))
    result[0:nkk04-1].calibration = 'KK04'
    result[nkk04:nkk04+npt05-1].calibration = 'PT05'
    result.nice_galaxy = strcompress(result.nice_galaxy,/remove)

    outfile = outpath+'sings_moustakas_table8.dat'
    splog, 'Writing '+outfile
    openw, lun, outfile, /get_lun
    printf, lun, '## John Moustakas, john.moustakas@gmail.com'
    printf, lun, '## '+strmid(im_today(),0,11)
    printf, lun, '## '
    printf, lun, '## This file contains Table 8 in "Optical Spectroscopy and'
    printf, lun, '## Nebular Oxygen Abundances of the Spitzer/SINGS Galaxies",'
    printf, lun, '## ApJS, 2010, in press (astro-ph/1007.4547) by J. Moustakas, '
    printf, lun, '## R.C. Kennicutt, C.A. Tremonti, D.A. Dale, J.D.T. Smith, & '
    printf, lun, '## D. Calzetti.  Publications that make use of these abundances '
    printf, lun, '## must cite this paper.'
    printf, lun, '## '
    printf, lun, '## Abundances have been computed using both the empirical'
    printf, lun, '## Pilyugin & Thuan (2005, PT05) and theoretical Kobulnicky &'
    printf, lun, '## Kewley (2004, KK04) R23 calibrations. There are numerous'
    printf, lun, '## caveats and potential systematic uncertainties affecting'
    printf, lun, '## both calibrations, as discussed extensively in the associated'
    printf, lun, '## paper, and these issues must be kept in mind in any analysis'
    printf, lun, '## making use of these abundances.'
    printf, lun, '## '
    printf, lun, '## Null measurements are indicated by -999.  See the paper '
    printf, lun, '## for a complete description of the meaning of each column.'
    printf, lun, '## '
    printf, lun, '#  1 Galaxy'
    printf, lun, '#  2 R23_Branch'
    printf, lun, '#  3 Slope'
    printf, lun, '#  4 Slope_Err'
    printf, lun, '#  5 log12oh_0'
    printf, lun, '#  6 log12oh_0_Err'
    printf, lun, '#  7 log12oh_Char'
    printf, lun, '#  8 log12oh_Char_Err'
    printf, lun, '#  9 log12oh_Avg'
    printf, lun, '# 10 log12oh_Avg_Err'
    printf, lun, '# 11 NHII'
    printf, lun, '# 12 Calibration'

    struct_print, result, lun=lun, /no_head
    free_lun, lun

; ###########################################################################
; TABLE 7 - nuclear, circumnuclear, and radial-strip abundances
    nucselect = ['nice_galaxy','r23_branch','nuclear_r23','nuclear_o32',$
      'nuclear_p','nuclear_logu_kk04','nuclear_log12oh_kk04','nuclear_log12oh_pt05']
    d20select = ['nice_galaxy','r23_branch','drift20_r23','drift20_o32',$
      'drift20_p','drift20_logu_kk04','drift20_log12oh_kk04','drift20_log12oh_pt05']
    d56select = ['nice_galaxy','r23_branch','drift56_r23','drift56_o32',$
      'drift56_p','drift56_logu_kk04','drift56_log12oh_kk04','drift56_log12oh_pt05']
    nuc = im_struct_trimtags(result1[where(result1.nuclear_log12oh_kk04[0] gt 0.0)],$
      select=nucselect,newtags=repstr(nucselect,'nuclear_',''))
    d20 = im_struct_trimtags(result1[where(result1.drift20_log12oh_kk04[0] gt 0.0)],$
      select=d20select,newtags=repstr(d20select,'drift20_',''))
    d56 = im_struct_trimtags(result1[where(result1.drift56_log12oh_kk04[0] gt 0.0)],$
      select=d56select,newtags=repstr(d56select,'drift56_',''))
    nnuc = n_elements(nuc)
    nd20 = n_elements(d20)
    nd56 = n_elements(d56)
    result = [nuc,d20,d56]
    result = struct_addtags(result,replicate({spectrum: ''},nnuc+nd20+nd56))
    result[0:nnuc-1].spectrum = 'nuclear'
    result[nnuc:nnuc+nd20-1].spectrum = 'circum'
    result[nnuc+nd20:nnuc+nd20+nd56-1].spectrum = 'strip'
    result.nice_galaxy = strcompress(result.nice_galaxy,/remove)

    outfile = outpath+'sings_moustakas_table7.dat'
    splog, 'Writing '+outfile
    openw, lun, outfile, /get_lun
    printf, lun, '## John Moustakas, john.moustakas@gmail.com'
    printf, lun, '## '+strmid(im_today(),0,11)
    printf, lun, '## '
    printf, lun, '## This file contains Table 7 in "Optical Spectroscopy and'
    printf, lun, '## Nebular Oxygen Abundances of the Spitzer/SINGS Galaxies",'
    printf, lun, '## ApJS, 2010, in press (astro-ph/1007.4547) by J. Moustakas, '
    printf, lun, '## R.C. Kennicutt, C.A. Tremonti, D.A. Dale, J.D.T. Smith, & '
    printf, lun, '## D. Calzetti.  Publications that make use of these abundances '
    printf, lun, '## must cite this paper.'
    printf, lun, '## '
    printf, lun, '## Abundances have been computed using both the empirical'
    printf, lun, '## Pilyugin & Thuan (2005, PT05) and theoretical Kobulnicky &'
    printf, lun, '## Kewley (2004, KK04) R23 calibrations. There are numerous'
    printf, lun, '## caveats and potential systematic uncertainties affecting'
    printf, lun, '## both calibrations, as discussed extensively in the associated'
    printf, lun, '## paper, and these issues must be kept in mind in any analysis'
    printf, lun, '## making use of these abundances.'
    printf, lun, '## '
    printf, lun, '## Null measurements are indicated by -999.  See the paper '
    printf, lun, '## for a complete description of the meaning of each column.'
    printf, lun, '## '
    printf, lun, '#  1 Galaxy'
    printf, lun, '#  2 R23_Branch'
    printf, lun, '#  3 R23'
    printf, lun, '#  4 R23_Err'
    printf, lun, '#  5 O32'
    printf, lun, '#  6 O32_Err'
    printf, lun, '#  7 P'
    printf, lun, '#  8 P_Err'
    printf, lun, '#  9 LogU'
    printf, lun, '# 10 LogU_Err'
    printf, lun, '# 11 log12oh_KK04'
    printf, lun, '# 12 log12oh_KK04_Err'
    printf, lun, '# 13 log12oh_PT05'
    printf, lun, '# 14 log12oh_PT05_Err'
    printf, lun, '# 15 Spectrum'

    struct_print, result, lun=lun, /no_head
    free_lun, lun

; ###########################################################################
; TABLE 9 - final abundances    
    select = ['GALAXY',$
      'LOG12OH_KK04_CENTRAL','LOG12OH_KK04_CHAR','LOG12OH_KK04_LZ',$
      'LOG12OH_PT05_CENTRAL','LOG12OH_PT05_CHAR','LOG12OH_PT05_LZ']
    final = im_struct_trimtags(final1,select=select)
    final.galaxy = strcompress(final1.nice_galaxy,/remove)

    outfile = outpath+'sings_moustakas_table9.dat'
    splog, 'Writing '+outfile
    openw, lun, outfile, /get_lun
    printf, lun, '## John Moustakas, john.moustakas@gmail.com'
    printf, lun, '## '+strmid(im_today(),0,11)
    printf, lun, '## '
    printf, lun, '## This file contains Table 9 in "Optical Spectroscopy and'
    printf, lun, '## Nebular Oxygen Abundances of the Spitzer/SINGS Galaxies",'
    printf, lun, '## ApJS, 2010, in press (astro-ph/1007.4547) by J. Moustakas, '
    printf, lun, '## R.C. Kennicutt, C.A. Tremonti, D.A. Dale, J.D.T. Smith, & '
    printf, lun, '## D. Calzetti.  Publications that make use of these abundances '
    printf, lun, '## must cite this paper.'
    printf, lun, '## '
    printf, lun, '## Abundances have been computed using both the empirical'
    printf, lun, '## Pilyugin & Thuan (2005, PT05) and theoretical Kobulnicky &'
    printf, lun, '## Kewley (2004, KK04) R23 calibrations. There are numerous'
    printf, lun, '## caveats and potential systematic uncertainties affecting'
    printf, lun, '## both calibrations, as discussed extensively in the associated'
    printf, lun, '## paper, and these issues must be kept in mind in any analysis'
    printf, lun, '## making use of these abundances. The table also includes abundances'
    printf, lun, '## inferred using the B-band luminosity-metallicity (LZ) relation,'
    printf, lun, '## which may be subject to additional systematic errors, as '
    printf, lun, '## described in the paper.'
    printf, lun, '## '
    printf, lun, '## Null measurements are indicated by -999.  See the paper '
    printf, lun, '## for a complete description of the meaning of each column.'
    printf, lun, '## '
    printf, lun, '#  1 Galaxy'
    printf, lun, '#  2 log12oh_KK04_Cent'
    printf, lun, '#  3 log12oh_KK04_Cent_Err'
    printf, lun, '#  4 log12oh_KK04_Char'
    printf, lun, '#  5 log12oh_KK04_Char_Err'
    printf, lun, '#  6 log12oh_KK04_LZ'
    printf, lun, '#  7 log12oh_KK04_LZ_Err'
    printf, lun, '#  8 log12oh_PT05_Cent'
    printf, lun, '#  9 log12oh_PT05_Cent_Err'
    printf, lun, '# 10 log12oh_PT05_Char'
    printf, lun, '# 11 log12oh_PT05_Char_Err'
    printf, lun, '# 12 log12oh_PT05_LZ'
    printf, lun, '# 13 log12oh_PT05_LZ_Err'

    struct_print, final, lun=lun, /no_head
    free_lun, lun

return
end
    
