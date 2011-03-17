pro write_97ho, out, out_nodust, _extra=extra
; jm03sep15uofa - read the basic data from NED (ho97.ned) and the
;                 classification table from Ho et al (1997) from
;                 Vizier, cross-match, and write a final binary FITS
;                 file that contains the galaxy name, coordinates, and
;                 nuclear spectral classification
; jm05may17uofa - redo the classification using the Kauffmann et 
;                 al. curve     
; jm06feb08uofa - corrected the lines for Galactic and internal
;                 reddening; derive strong-line nebular abundances 
    
    oratio = 2.984
    nratio = 3.054

;   parse_ned_byname, '97ho.ned', ned, inputnedfile='galaxy_list.txt', outfile='97ho_ned.fits'

    ned = mrdfits('97ho_ned.fits.gz',1,/silent)
    class = mrdfits('97ho_table4.fits.gz',1,/silent)
    ratios = mrdfits('97ho_table2.fits.gz',1,/silent)
    ancil = mrdfits('97ho_table10_matched.fits.gz',1,/silent)
    phot = mrdfits('97ho_table11_matched.fits.gz',1,/silent)
    dist = mrdfits('97ho_distances.fits.gz',1,/silent)

;   phot = mrdfits('97ho_table11.fits.gz',1,/silent)
;   match = match_string(ratios.name+ratios.m_name,phot.name+phot.m_name,/exact,findex=indx)
;   mwrfits, phot[indx], '97ho_table11_matched.fits', /create
;   spawn, ['gzip -f '+'97ho_table11_matched.fits'], /sh
;   niceprint, phot[indx].name+phot[indx].m_name, ratios.name+ratios.m_name
;   
;   ancil = mrdfits('97ho_table10.fits.gz',1,/silent)
;   match = match_string(ratios.name+ratios.m_name,ancil.name+ancil.m_name,/exact,findex=indx)
;   mwrfits, ancil[indx], '97ho_table10_matched.fits', /create
;   spawn, ['gzip -f '+'97ho_table10_matched.fits'], /sh
;   niceprint, ancil[indx].name+ancil[indx].m_name, ratios.name+ratios.m_name
    
    ned = struct_trimtags(ned,select=['GALAXY','NED_GALAXY','RA','DEC'])

    ned = struct_addtags(ned,struct_trimtags(dist,select=['DISTANCE*']))
;   class = struct_trimtags(class,select=['CLASS'])

    glactc, 15.0D*im_hms2dec(ned.ra), im_hms2dec(ned.dec), 2000.0, gl, gb, 1, /degree
    ebv_mw = dust_getval(gl,gb,/interp)

    ngalaxy = n_elements(class)
    
; use ICLASSIFICATION to classify the galaxies in this sample; the
; ratios are *observed*, not reddening-corrected!

    line = {$
;     galaxy:           '',        $
      linename:         ['H_BETA','OIII_4959','OIII_5007','OI_6300',$
      'NII_6548','H_ALPHA','NII_6584','SII_6716','SII_6731'], $
      h_beta_wave:      4861.33, $
      oiii_4959_wave:   4958.91, $
      oiii_5007_wave:   5006.84, $
      oi_6300_wave:     6300.304,$
      nii_6548_wave:    6548.04, $
      h_alpha_wave:     6562.80, $
      nii_6584_wave:    6583.46, $
      sii_6716_wave:    6716.14, $
      sii_6731_wave:    6730.81, $
      h_beta:        [0.0,-2.0], $
      oiii_4959:     [0.0,-2.0], $
      oiii_5007:     [0.0,-2.0], $
      oi_6300:       [0.0,-2.0], $
      nii_6548:      [0.0,-2.0], $
      h_alpha:       [0.0,-2.0], $
      nii_6584:      [0.0,-2.0], $
      sii_6716:      [0.0,-2.0], $
      sii_6731:      [0.0,-2.0], $
      photflag:             'Y', $
      ho_class:              ''}
    line = replicate(line,ngalaxy)

    line.ho_class = class.class

    k_hb        = k_lambda(line[0].h_beta_wave,/odonnell)
    k_ha        = k_lambda(line[0].h_alpha_wave,/odonnell)
    k_oiii_4959 = k_lambda(line[0].oiii_4959_wave,/odonnell)
    k_oiii_5007 = k_lambda(line[0].oiii_5007_wave,/odonnell)
    k_oi_6300   = k_lambda(line[0].oi_6300_wave,/odonnell)  
    k_nii_6548  = k_lambda(line[0].nii_6548_wave,/odonnell) 
    k_nii_6584  = k_lambda(line[0].nii_6584_wave,/odonnell) 
    k_sii_6716  = k_lambda(line[0].sii_6716_wave,/odonnell) 
    k_sii_6731  = k_lambda(line[0].sii_6731_wave,/odonnell) 
    
;   k_hb_ha        = k_lambda(line[0].h_beta_wave,/odonnell)   -k_lambda(line[0].h_alpha_wave,/odonnell)
;   k_oiii_4959_ha = k_lambda(line[0].oiii_4959_wave,/odonnell)-k_lambda(line[0].h_alpha_wave,/odonnell)
;   k_oiii_5007_ha = k_lambda(line[0].oiii_5007_wave,/odonnell)-k_lambda(line[0].h_alpha_wave,/odonnell)
;   k_oi_6300_ha   = k_lambda(line[0].oi_6300_wave,/odonnell)  -k_lambda(line[0].h_alpha_wave,/odonnell)
;   k_nii_6548_ha  = k_lambda(line[0].nii_6548_wave,/odonnell) -k_lambda(line[0].h_alpha_wave,/odonnell)
;   k_nii_6584_ha  = k_lambda(line[0].nii_6584_wave,/odonnell) -k_lambda(line[0].h_alpha_wave,/odonnell)
    
    nonphot = where((strmatch(ratios.n_logf_ha_,'*L*') eq 1B),nindx)
    line[nonphot].photflag = 'N'

    indx = where((ratios.logf_ha_ ne 0.0),nindx)
;   indx = where((ratios.logf_ha_ ne 0.0) and (strmatch(ratios.n_logf_ha_,'*u*') eq 0B) and $
;     (strmatch(ratios.n_logf_ha_,'*L*') eq 0B) and (strmatch(ratios.n_logf_ha_,'*s*') eq 0B),nindx)
    line[indx].h_alpha[0] = 10^ratios[indx].logf_ha_
    line[indx].h_alpha[1] = line[indx].h_alpha[0]*0.1 ; 10% error
    hierr = where(ratios[indx].n_logf_ha_ eq 'b')
    line[indx[hierr]].h_alpha[1] = line[indx[hierr]].h_alpha[0]*0.4 ; 40% error
    hierr = where(ratios[indx].n_logf_ha_ eq 'c')
    line[indx[hierr]].h_alpha[1] = line[indx[hierr]].h_alpha[0]*1.0 ; 100% error
    
    line[indx].h_alpha = line[indx].h_alpha*rebin(rebin(10^(0.4*ebv_mw[indx]*k_ha),1,nindx),2,nindx)

    ul = where((strmatch(ratios[indx].n_logf_ha_,'*u*') eq 1B),nul) ; upper limit
    line[indx[ul]].h_alpha[1] = -3.0

; H-beta    
    
    indx = where((ratios.hbeta gt 0.0),nindx)
    line[indx].h_beta[0] = ratios[indx].hbeta*line[indx].h_alpha[0]
    line[indx].h_beta[1] = line[indx].h_beta[0]*0.1 ; 10% error
    hierr = where(ratios[indx].n_hbeta eq 'b')
    line[indx[hierr]].h_beta[1] = line[indx[hierr]].h_beta[0]*0.4 ; 40% error
    hierr = where(ratios[indx].n_hbeta eq 'c')
    line[indx[hierr]].h_beta[1] = line[indx[hierr]].h_beta[0]*1.0 ; 100% error

    line[indx].h_beta = line[indx].h_beta*rebin(rebin(10^(0.4*ebv_mw[indx]*k_hb),1,nindx),2,nindx)
    
    ul = where((strmatch(ratios[indx].n_hbeta,'*u*') eq 1B),nul) ; upper limit
    line[indx[ul]].h_beta[1] = -3.0

; [O I] 6300
    
    indx = where((ratios._oi_ gt 0.0),nindx)
    line[indx].oi_6300[0] = ratios[indx]._oi_*line[indx].h_alpha[0]
    line[indx].oi_6300[1] = line[indx].oi_6300[0]*0.1 ; 10% error
    hierr = where(ratios[indx].n__oi_ eq 'b')
    line[indx[hierr]].oi_6300[1] = line[indx[hierr]].oi_6300[0]*0.4 ; 40% error
    hierr = where(ratios[indx].n__oi_ eq 'c')
    line[indx[hierr]].oi_6300[1] = line[indx[hierr]].oi_6300[0]*1.0 ; 100% error

    line[indx].oi_6300 = line[indx].oi_6300*rebin(rebin(10^(0.4*ebv_mw[indx]*k_oi_6300),1,nindx),2,nindx)
    
    ul = where((strmatch(ratios[indx].n__oi_,'*u*') eq 1B),nul)
    line[indx[ul]].oi_6300[1] = -3.0

; [O III] 5007 and 4959
    
    indx = where((ratios._oiii_ gt 0.0),nindx)
    line[indx].oiii_5007[0] = ratios[indx]._oiii_*line[indx].h_alpha[0]
    line[indx].oiii_5007[1] = line[indx].oiii_5007[0]*0.1 ; 10% error
    hierr = where(ratios[indx].n__oiii_ eq 'b')
    line[indx[hierr]].oiii_5007[1] = line[indx[hierr]].oiii_5007[0]*0.4 ; 40% error
    hierr = where(ratios[indx].n__oiii_ eq 'c')
    line[indx[hierr]].oiii_5007[1] = line[indx[hierr]].oiii_5007[0]*1.0 ; 100% error

    line[indx].oiii_4959 = line[indx].oiii_5007/oratio
    line[indx].oiii_4959 = line[indx].oiii_4959*rebin(rebin(10^(0.4*ebv_mw[indx]*k_oiii_4959),1,nindx),2,nindx)
    line[indx].oiii_5007 = line[indx].oiii_5007*rebin(rebin(10^(0.4*ebv_mw[indx]*k_oiii_5007),1,nindx),2,nindx)

    ul = where((strmatch(ratios[indx].n__oiii_,'*u*') eq 1B),nul)
    line[indx[ul]].oiii_4959[1] = -3.0
    line[indx[ul]].oiii_5007[1] = -3.0

; [N II] 6584 and 6548
    
    indx = where((ratios._nii_ gt 0.0),nindx)
    line[indx].nii_6584[0] = ratios[indx]._nii_*line[indx].h_alpha[0]
    line[indx].nii_6584[1] = line[indx].nii_6584[0]*0.1 ; 10% error
    hierr = where(ratios[indx].n__nii_ eq 'b')
    line[indx[hierr]].nii_6584[1] = line[indx[hierr]].nii_6584[0]*0.4 ; 40% error
    hierr = where(ratios[indx].n__nii_ eq 'c')
    line[indx[hierr]].nii_6584[1] = line[indx[hierr]].nii_6584[0]*1.0 ; 100% error

    line[indx].nii_6548 = line[indx].nii_6584/nratio
    line[indx].nii_6548 = line[indx].nii_6548*rebin(rebin(10^(0.4*ebv_mw[indx]*k_nii_6548),1,nindx),2,nindx)
    line[indx].nii_6584 = line[indx].nii_6584*rebin(rebin(10^(0.4*ebv_mw[indx]*k_nii_6584),1,nindx),2,nindx)

    ul = where((strmatch(ratios[indx].n__nii_,'*u*') eq 1B),nul)
    line[indx[ul]].nii_6548[1] = -3.0
    line[indx[ul]].nii_6584[1] = -3.0

; [S II] 6716
    
    indx = where((ratios._sii_a gt 0.0),nindx)
    line[indx].sii_6716[0] = ratios[indx]._sii_a*line[indx].h_alpha[0]
    line[indx].sii_6716[1] = line[indx].sii_6716[0]*0.1 ; 10% error
    hierr = where(ratios[indx].n__sii_a eq 'b')
    line[indx[hierr]].sii_6716[1] = line[indx[hierr]].sii_6716[0]*0.4 ; 40% error
    hierr = where(ratios[indx].n__sii_a eq 'c')
    line[indx[hierr]].sii_6716[1] = line[indx[hierr]].sii_6716[0]*1.0 ; 100% error

    line[indx].sii_6716 = line[indx].sii_6716*rebin(rebin(10^(0.4*ebv_mw[indx]*k_sii_6716),1,nindx),2,nindx)
    
    ul = where((strmatch(ratios[indx].n__sii_a,'*u*') eq 1B),nul)
    line[indx[ul]].sii_6716[1] = -3.0

; [S II] 6731
    
    indx = where((ratios._sii_b gt 0.0),nindx)
    line[indx].sii_6731[0] = ratios[indx]._sii_b*line[indx].h_alpha[0]
    line[indx].sii_6731[1] = line[indx].sii_6731[0]*0.1 ; 10% error
    hierr = where(ratios[indx].n__sii_b eq 'b')
    line[indx[hierr]].sii_6731[1] = line[indx[hierr]].sii_6731[0]*0.4 ; 40% error
    hierr = where(ratios[indx].n__sii_b eq 'c')
    line[indx[hierr]].sii_6731[1] = line[indx[hierr]].sii_6731[0]*1.0 ; 100% error

    line[indx].sii_6731 = line[indx].sii_6731*rebin(rebin(10^(0.4*ebv_mw[indx]*k_sii_6731),1,nindx),2,nindx)
    
    ul = where((strmatch(ratios[indx].n__sii_b,'*u*') eq 1B),nul)
    line[indx[ul]].sii_6731[1] = -3.0

; compute the emission-line luminosities

    lums = compute_linelums(line,ancillary=ned,select_lines=['OIII_5007'])
    line = struct_addtags(line,lums)
    
; deredden the line fluxes

    line_nodust = iunred_linedust(line,/silent)
;   niceprint, line_nodust.hahb, class.ha_hb
    
; classify using the observed and the de-reddened line fluxes; some of
; the objects may change classification if they are on the boundary,
; so adopt the class based on the de-reddened fluxes    
    
    iclass = iclassification(line,ratios=ratios,_extra=extra,snrcut=1.0)

;   iclass_nodust = iclassification(line_nodust,ratios=ratios_nodust,_extra=extra)
;   niceprint, ned.galaxy, iclass.bpt_nii_mixture_class, iclass_nodust.bpt_nii_mixture_class
;   w = where((iclass.bpt_nii_mixture_class ne iclass_nodust.bpt_nii_mixture_class) and $
;     (iclass.bpt_nii_mixture_class ne 'Unknown') and (iclass_nodust.bpt_nii_mixture_class) ne 'Unknown')
;   niceprint, ned[w].galaxy, iclass[w].bpt_nii_mixture_class, iclass_nodust[w].bpt_nii_mixture_class
;   iclass[w].bpt_nii_mixture_class = iclass_nodust[w].bpt_nii_mixture_class

;   plot_kewley_grids, /hii
;   plotsym, 0, 0.5, /fill
;   djs_oplot, ratios.nii_ha, ratios.oiii_hb, ps=8, color='yellow'
;   djs_oplot, ratios_nodust.nii_ha, ratios_nodust.oiii_hb, ps=8, color='red'
    
; no [O II]!
    
    abund = im_abundance(line,snrcut=0.0)
    abund_nodust = im_abundance(line_nodust,snrcut=0.0)

    out = struct_addtags(struct_addtags(ned,struct_addtags(line,$
      im_struct_trimtags(iclass,select='BPT_NII_MIXTURE_CLASS',newtags='CLASS'))),abund)
    out_nodust = struct_addtags(struct_addtags(ned,struct_addtags(line_nodust,$
      im_struct_trimtags(iclass,select='BPT_NII_MIXTURE_CLASS',newtags='CLASS'))),abund_nodust)

; write out    
    
    srt = sort(im_hms2dec(out.ra))
    out = out[srt]
    out_nodust = out_nodust[srt]
    
    mwrfits, out, '97ho.fits', /create
    spawn, ['gzip -f 97ho.fits'], /sh

    mwrfits, out_nodust, '97ho_nodust.fits', /create
    spawn, ['gzip -f 97ho_nodust.fits'], /sh

; simplify the spectral classification    
    
;   oldclass = strcompress(class.class,/remove)
;   newclass = oldclass
;
;   newclass[where(newclass eq '')] = 'U'
;   newclass[where(strmatch(newclass,'H*') eq 1B)] = 'H'
;   newclass[where(strmatch(newclass,'*S2*') eq 1B)] = 'S2'
;   newclass[where(strmatch(newclass,'*T*') eq 1B)] = 'L'
;   newclass[where(strmatch(newclass,'*L*') eq 1B)] = 'L'
;
;   newclass = repstr(newclass,':','')    
;   
;   niceprint, ned.galaxy, oldclass, newclass
;
;   class.class = newclass
;   out = struct_addtags(ned,class)

return
end    
