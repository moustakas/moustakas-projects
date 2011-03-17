pro hcn_lit_classify
; jm08oct02nyu 

    linedust = {galaxy: '', linename: ['H_ALPHA','H_BETA','OIII_5007','NII_6584'], $
      h_alpha: [0.0,0.0], h_beta: [0.0,0.0], oiii_5007: [0.0,0.0], nii_6584: [0.0,0.0], $
      h_alpha_wave: 6562.8, h_beta_wave: 4861.33, oiii_5007_wave: 5006.84, nii_6584_wave: 6583.46}
    linedust = replicate(linedust,2)

; Veilleux et al. 1999    
    linedust[0].galaxy = 'IRAS12112+0305'
    linedust[0].h_alpha = 0.8E-14*[1.0,0.1] 
    linedust[0].h_beta = 0.13*linedust[0].h_alpha
    linedust[0].oiii_5007 = 0.27*linedust[0].h_alpha
    linedust[0].nii_6584 = 0.64*linedust[0].h_alpha

; Wu et al. 1998
    wu = mrdfits('~/catalogs/98wu/98wu_table1.fits',1)
    indx = where(strmatch(wu.iras,'*05084*'))
    linedust[1].galaxy = 'VIIZw31' ; = IRAS05084+7936
    linedust[1].h_beta = [1.0,0.01]
    linedust[1].h_alpha = [wu[indx].ha_hb,0.01]
    linedust[1].oiii_5007 = [10.0^wu[indx].log__oiii__hb_,0.01]
    linedust[1].nii_6584 = [10.0^wu[indx].log__nii__ha_*wu[indx].ha_hb,0.01]

    line = iunred_linedust(linedust)
    class = iclassification(linedust,/doplot,ratios=rr)
    struct_print, rr

    mwrfits, line, 'hcn_lit_classify.fits', /create
    mwrfits, linedust, 'hcn_lit_classify.fits'
    
stop    
    
return
end
    
