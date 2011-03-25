function bootes_filterlist, noirac=noirac
; jm10jan05ucsd
    filters = ['lbc_blue_ufilter',$
      'ndwfs_'+['Bw','R','I'],$
      'bok_90prime_z',$
      'newfirm_'+['J','H','Ks'],$
      'spitzer_irac_'+['ch1','ch2','ch3','ch4']$
      ]+'.par'
    if keyword_set(noirac) then filters = filters[where(strmatch(filters,'*irac*',/fold) eq 0)]
return, filters    
end
