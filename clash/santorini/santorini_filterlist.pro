function santorini_filterlist, short_filter=short_filter, nice_filter=nice_filter, $
  zpt=zpt, alam=alam, weff=weff, fwhm=fwhm, pivotwave=pivotwave, width=width
; jm11apr24ucsd 

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
      'ACS-F435W',$
      'ACS-F475W',$
      'ACS-F555W',$
      'ACS-F606W',$
      'ACS-F625W',$
      'ACS-F775W',$
      'ACS-F814W',$
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

; http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c06_uvis06.html#387761 - WFC3/UVIS
; http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c07_ir06.html - WFC3/IR
; http://etc.stsci.edu/etcstatic/users_guide/appendix_b_acs.html - ACS/WFC
; http://iopscience.iop.org/1538-3881/141/5/173/pdf/aj_141_5_173.pdf - IRAC

    pivotwave = [$
; WFC3/UVIS
      235.9,$
      270.4,$
      335.5,$
      392.1,$
; ACS/WFC
      431.70,$ ; f435w
      474.44,$
      535.964,$ ; f555w
      591.77,$
      631.05,$
      769.30,$
      805.98,$
      905.48,$
; WFC3/IR
      1055.2,$
      1153.4,$
      1248.6,$
      1392.3,$
      1536.9]*10 ; Angstrom

    width = [$
; WFC3/UVIS
      46.7,$
      39.8,$
      51.1,$
      89.6,$
; ACS/WFC ; FWHM from 
      69.108,$
      98.927,$
      84.779,$ ; f555w
      158.32,$
      97.832,$
      102.34,$
      154.16,$
      127.03,$
; WFC3/IR
      265.0,$
      443.0,$
      284.5,$
      384.0,$
      268.3]*10                 ; Angstrom
    
; add IRAC    
    filterlist = [filterlist,(irac_filterlist())[0:1]]
    nice_filter = [nice_filter,'IRAC-[ch1]','IRAC-[ch2]']
    short_filter = [short_filter,'ch1','ch2']
    pivotwave = [pivotwave,1D4*[3.551,4.496,5.724,7.884]]
;   width = [width,[0.7,1.0,1.5,3.0]*1D4] ; approximate for IRAC!
    width = [width,[4313.0,5712,8197,16717]*2] ; approximate for IRAC!

    if arg_present(weff) then weff = k_lambda_eff(filterlist=filterlist)

return, filterlist
end
