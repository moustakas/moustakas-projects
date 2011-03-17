pro mz_withagn
; jm09feb06nyu - 
    
    kcorr1 = read_sdss_vagc_mpa(/kcorr,sample='dr7',$
      letter=letter,poststr=poststr)
    ispec1 = read_sdss_vagc_mpa(/ispec,sample='dr7',$
      letter=letter,poststr=poststr)

    snrcut = 5.0
    keep = where($
      (ispec1.oii_3727[0]/ispec1.oii_3727[1] gt snrcut) and $
      (ispec1.nii_6584[0]/ispec1.nii_6584[1] gt snrcut) and $
      (ispec1.oiii_5007[0]/ispec1.oiii_5007[1] gt snrcut) and $
      (ispec1.h_beta[0]/ispec1.h_beta[1] gt snrcut) and $
      (ispec1.h_alpha[0]/ispec1.h_alpha[1] gt snrcut),nkeep)

    iclass = iclassification(ispec1,snrcut_class=snrcut,$
      ratios=iratios,silent=silent)
    ispecnodust = iunred_linedust(ispec1,snrcut=snrcut,silent=silent)
    iabund = im_abundance(ispecnodust,nmonte=0,snrcut_abund=snrcut,$
      /justflux,/nodensity)

; make the plots

    sf = where(strmatch(strtrim(iratios.final_class,2),'*SF*') and $
      (iabund.zstrong_12oh_kk04_upper gt -900.0))
    agn = where(strmatch(strtrim(iratios.final_class,2),'*AGN*') and $
      (iabund.zstrong_12oh_kk04_upper gt -900.0))

    dfpsplot, 'mz_withagn.ps', /square, /color
    hogg_scatterplot, kcorr1[sf].mass, iabund[sf].zstrong_12oh_kk04_upper, $
      xrange=[8.3,11.7], yrange=[8.4,9.3], /xsty, /ysty
    djs_oplot, kcorr1[agn].mass, iabund[agn].zstrong_12oh_kk04_upper, $
      psym=6, symsize=0.2, color='dark green'
    dfpsclose

; plots to make
; * residuals from MZ relation vs D(BPT): (1) no AGN culling; (2)
;   Kewley AGN culling; (3) Kauffmann AGN culling
    
stop    
    
return
end
    
