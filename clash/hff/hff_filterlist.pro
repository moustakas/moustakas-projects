function hff_filterlist, short_filter=short_filter, nice_filter=nice_filter, $
  zpt=zpt, alam=alam, useirac=useirac, weff=weff, fwhm=fwhm, pivotwave=pivotwave, $
  width=width, instr=instr, usehawki=usehawki, usemoircs=usemoircs
; jm11apr24ucsd 

    filterlist = [$
      'clash_acs_f435w.par',$
      'clash_acs_f606w.par',$
      'clash_acs_f814w.par',$
      'clash_wfc3_f105w.par',$
      'clash_wfc3_f125w.par',$
      'clash_wfc3_f140w.par',$
      'clash_wfc3_f160w.par']
    nice_filter = [$
      'ACS-F435w',$
      'ACS-F606w',$
      'ACS-F814w',$
      'WFC3/IR-F105W',$
      'WFC3/IR-F125W',$
      'WFC3/IR-F140W',$
      'WFC3/IR-F160W']
    short_filter = [$
      'f435w',$
      'f606w',$
      'f814w',$
      'f105w',$
      'f125w',$
      'f140w',$
      'f160w']

    instr = [$
      'acs',$
      'acs',$
      'acs',$
      'wfc3ir',$
      'wfc3ir',$
      'wfc3ir',$
      'wfc3ir']

; Zeropoint (AB mag) for each filter =
;   photzpt = sxpar(hdr,'PHOTZPT')
;   photflam = sxpar(hdr,'PHOTFLAM')
;   photplam = sxpar(hdr,'PHOTPLAM')
;   zpt = photzpt - 2.5*alog10(photflam) - 5.0*alog10(photplam/5475.4)
    zpt = [$
; ACS
      25.65779,$
      26.49116,$
      25.94335,$
; WFC3/IR
      26.27068,$
      26.24736,$
      26.46450,$
      25.95590]

; http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c06_uvis06.html#387761 - WFC3/UVIS
; http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c07_ir06.html - WFC3/IR
; http://etc.stsci.edu/etcstatic/users_guide/appendix_b_acs.html - ACS/WFC
; http://iopscience.iop.org/1538-3881/141/5/173/pdf/aj_141_5_173.pdf - IRAC

    pivotwave = [$
; ACS/WFC
      431.70,$
      591.77,$
      805.98,$
; WFC3/IR
      1055.2,$
      1248.6,$
      1392.3,$
      1536.9]*10 ; Angstrom

    width = [$
; ACS/WFC
      69.108,$
      158.32,$
      154.16,$
; WFC3/IR
      265.0,$
      284.5,$
      384.0,$
      268.3]*10                 ; Angstrom
    
; add IRAC    
    if keyword_set(useirac) then begin
       filterlist = [filterlist,irac_filterlist(/warm)]
       nice_filter = [nice_filter,'IRAC-[ch1]','IRAC-[ch2]']
       short_filter = [short_filter,'irac_ch1','irac_ch2']
       zpt = [zpt,18.803,18.318]
       pivotwave = [pivotwave,1D4*[3.551,4.496]]
;      width = [width,[0.7,1.0,1.5,3.0]*1D4] ; approximate for IRAC!
       width = [width,[4313.0,5712]*2] ; approximate for IRAC!
    endif

; add HAWK-I/Ks
    if keyword_set(usehawki) then begin
       filterlist = [filterlist,'hawki_Ks1.par']
       nice_filter = [nice_filter,'HAWKI-Ks']
       short_filter = [short_filter,'HAWKI_KS']
       zpt = [zpt,0.0]
       pivotwave = [pivotwave,21420.8]
       width = [width,1922.24]
    endif

; add MOIRCS/K
    if keyword_set(usemoircs) then begin
       filterlist = [filterlist,'moircs_K.par']
       nice_filter = [nice_filter,'MOIRCS-K']
       short_filter = [short_filter,'MOIRCS_K']
       zpt = [zpt,0.0]
       pivotwave = [pivotwave,21420.8]
       width = [width,1922.24]
    endif

; get the effective wavelength    
    if arg_present(weff) then weff = k_lambda_eff(filterlist=filterlist)

return, filterlist
end
