function bootes_filterlist, noirac=noirac, bands=bands
; jm10jan05ucsd
    filters = ['lbc_blue_ufilter',$
      'ndwfs_'+['Bw','R','I'],$
      'bok_90prime_z',$
      'lbc_red_yfilter',$
      'newfirm_'+['J','H','Ks'],$
      'spitzer_irac_'+['ch1','ch2','ch3','ch4']$
      ]+'.par'
    bands = ['u','Bw','R','I','z','y','J','H','Ks','ch1','ch2','ch3','ch4']
    if keyword_set(noirac) then begin
       keep = where(strmatch(filters,'*irac*',/fold) eq 0)
       filters = filters[keep]
       bands = bands[keep]
    endif
return, filters    
end
