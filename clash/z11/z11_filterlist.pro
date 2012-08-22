function z11_filterlist, short_filter=short_filter, nice_filter=nice_filter, weff=weff, fwhm=fwhm
; jm12aug14siena

    filterlist = [$
      'clash_wfc3_f225w.par',$
      'clash_wfc3_f275w.par',$
      'clash_wfc3_f336w.par',$
      'clash_wfc3_f390w.par',$
      'clash_acs_f435w.par',$
      'clash_acs_f475w.par',$
      'clash_acs_f555w.par',$ ; archival band!
      'clash_acs_f606w.par',$
      'clash_acs_f625w.par',$
      'clash_acs_f775w.par',$
      'clash_acs_f814w.par',$
      'clash_acs_f850lp.par',$
      'clash_wfc3_f105w.par',$
      'clash_wfc3_f110w.par',$
      'clash_wfc3_f125w.par',$
      'clash_wfc3_f140w.par',$
      'clash_wfc3_f160w.par']
    nice_filter = [$
      'WFC3/UVIS-F225W',$
      'WFC3/UVIS-F275W',$
      'WFC3/UVIS-F336W',$
      'WFC3/UVIS-F390W',$
      'ACS-F435w',$
      'ACS-F475w',$
      'ACS-F555W',$
      'ACS-F606w',$
      'ACS-F625w',$
      'ACS-F775w',$
      'ACS-F814w',$
      'ACS-F850LP',$
      'WFC3/IR-F105W',$
      'WFC3/IR-F110W',$
      'WFC3/IR-F125W',$
      'WFC3/IR-F140W',$
      'WFC3/IR-F160W']
    short_filter = [$
      'f225w',$ 
      'f275w',$
      'f336w',$
      'f390w',$
      'f435w',$
      'f475w',$
      'f555w',$
      'f606w',$
      'f625w',$
      'f775w',$
      'f814w',$
      'f850lp',$
      'f105w',$
      'f110w',$
      'f125w',$
      'f140w',$
      'f160w']

; add IRAC    
    filterlist = [filterlist,(irac_filterlist())[0:1]]
    nice_filter = [nice_filter,'IRAC-[ch1]','IRAC-[ch2]']
    short_filter = [short_filter,'ch1','ch2']

    if arg_present(weff) then weff = k_lambda_eff(filterlist=filterlist)

return, filterlist
end
