function bootes_filterlist
; jm10jan05ucsd
    return, ['lbc_blue_ufilter',$
      'ndwfs_'+['Bw','R','I'],$
      'bok_90prime_z',$
      'newfirm_'+['J','H','Ks'],$
      'spitzer_irac_'+['ch1','ch2','ch3','ch4']$
      ]+'.par'
end
