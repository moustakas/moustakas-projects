function photoztemplates_filterlist
; jm09oc30ucsd
    return, ['galex_'+['FUV','NUV'],'sdss_'+['u0','g0','r0','i0','z0'],$
      'twomass_'+['J','H','Ks'],'spitzer_irac_'+['ch1','ch2','ch3','ch4'],$
      'spitzer_mips_'+['24','70','160'],$
;     'bessell_'+['B','V','R','I'],$
      'wise_'+['w1','w2','w3','w4'],$
      ['S1','S2']+'_custom'$
      ]+'.par'
end
