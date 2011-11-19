function clash_filterlist, short_filter=short_filter, nice_filter=nice_filter, $
  zpt=zpt, alam=alam, useirac=useirac, weff=weff
; jm11apr24ucsd 

    filterlist = [$
      'clash_wfc3_f225w.par',$
      'clash_wfc3_f275w.par',$
      'clash_wfc3_f336w.par',$
      'clash_wfc3_f390w.par',$
      'clash_acs_f435w.par',$
      'clash_acs_f475w.par',$
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

; zeropoints (as of 2011-Nov-09)    

; # MACS1206 mosdriz 20110815 2011-08-15
; # Zeropoint (AB mag) for each filter
; # with and without galactic extinction included.
; # Extinctions derived using value from Schlegel dust maps:
; # E(B-V) = 0.06283
;    zpt = [$
;      24.09656D,$
;      24.17417,$
;      24.64529,$
;      25.37141,$
;      25.65779,$
;      26.05926,$
;      26.49116,$
;      25.90669,$
;      25.66506,$
;      25.94335,$
;      24.84247,$
;      26.27068,$
;      26.82514,$
;      26.24736,$
;      26.46450,$
;      25.95590]
    zpt = [$
      23.64640D,$
      23.80435,$
      24.33872,$
      25.09953,$
      25.40982,$
      25.83493,$
      26.31478,$
      25.74580,$
      25.54353,$
      25.83358,$
      24.75373,$
      26.20955,$
      26.77237,$
      26.20177,$
      26.42782,$
      25.92759]    
    
; add IRAC    
    if keyword_set(useirac) then begin
       filterlist = [filterlist,(irac_filterlist())[0:1]]
       nice_filter = [nice_filter,'IRAC-[ch1]','IRAC-[ch2]']
       short_filter = [short_filter,'ch1','ch2']
       zpt = [zpt,0.0,0.0]
    endif

    if arg_present(weff) then weff = k_lambda_eff(filterlist=filterlist)
    
return, filterlist
end
