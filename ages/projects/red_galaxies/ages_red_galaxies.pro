pro ages_red_galaxies
; jm08jul03nyu - look into Kate's red-galaxy selection

    ii1 = read_ages(/ispec)
    aa1 = read_ages(/ancillary)
    kate1 = mrdfits('obj_cat_10jul08.fits',1)
    nkate = n_elements(kate1)

; select the parent sample    
    
    badplates = [106,110,209,310,311]
    good = lonarr(n_elements(aa1))+1L
    good[where_array(badplates,aa1.pass)] = 0L
    keep = where((aa1.ubvri_absmag[2] le -19.4) and (aa1.z ge 0.15) and $
      (aa1.z lt 0.3) and (aa1.main_flag) and (good eq 1),nkeep)

    aa = aa1[keep] & ii = ii1[keep]
    match, aa.galaxy, kate1.galaxy, m1, m2
;   niceprint, aa[m1].galaxy, kate1[m2].galaxy

; resort KATE to match my structures
    
    kate = kate1[m2]
    moretags = replicate({final_class: 'unclas'},nkate)
    out = struct_addtags(kate1,moretags)

; now identify red (quiescent and infrared-excess), blue, and green
; galaxies 

    z = aa.z
    mv = aa.ubvri_absmag[2]
    uv = aa.ubvri_absmag[0]-aa.ubvri_absmag[2]
    f24 = aa.phot_mips24

    allred = where(uv gt (1.55 - 0.25 - 0.04*(mv+20.0) - 0.42*(z-0.05) + 0.07*(z-0.05)^2))
    ir_red = where(uv gt (1.55 - 0.25 - 0.04*(mv+20.0) - 0.42*(z-0.05) + 0.07*(z-0.05)^2) and (f24 gt 0.3))
    red = where(uv gt (1.55 - 0.25 - 0.04*(mv+20.0) - 0.42*(z-0.05) + 0.07*(z-0.05)^2) and (f24 lt 0.3))
    green = where((uv le (1.55 - 0.25 - 0.04*(mv+20.0) - 0.42*(z-0.05) + 0.07*(z-0.05)^2)) and $
      (uv gt (1.35 - 0.25 - 0.04*(mv+20.0) - 0.42*(z-0.05) + 0.07*(z-0.05)^2)))
    blue = where(uv le (1.35 - 0.25 - 0.04*(mv+20.0) - 0.42*(z-0.05) + 0.07*(z-0.05)^2))
    help, allred, ir_red, red, green, blue

    im_plothist, uv[allred], line=2, bin=0.1, xr=[-0.5,2.0]
    im_plothist, uv[blue], /over                           
    im_plothist, uv[green], /over, line=5, thick=3

; now classify    

    dfpsplot, 'ir_red.bpt.ps', /square, /color
    bpt = iclassification(ii,snrcut_class=3.0,niihacut=-0.3,/doplot,/kauffmann)
    dfpsclose
    nbpt = n_elements(bpt)
    
    limit = where((strtrim(bpt.bpt_pure_nii_class,2) eq 'Unknown') and $
      (ii.h_alpha[0]/ii.h_alpha[1] gt 3.0) and $
;     (ii.oii_3727[0]/ii.oii_3727[1] gt 3.0) and $
;     (ii.h_beta[0]/ii.h_beta[1] gt 3.0) and $
      (ii.nii_6584[0]/ii.nii_6584[1] lt 3.0))
    sf_limit = where(alog10(ii[limit].nii_6584_limit/ii[limit].h_alpha[0]) lt -0.3)
;   im_plothist, alog10(ii[limit].nii_6584_limit/ii[limit].h_alpha[0]), bin=0.1
    bpt[limit[sf_limit]].bpt_pure_nii_class = 'HII'
    
    unk = where(strtrim(bpt.bpt_pure_nii_class,2) eq 'Unknown',nunk)
    sf = where(strtrim(bpt.bpt_pure_nii_class,2) eq 'HII',nsf)
    agn = where(strtrim(bpt.bpt_pure_nii_class,2) eq 'AGN',nagn)
    help, bpt, sf, agn, unk & print
    splog, 'All galaxies', ' sf: ', nsf/float(nbpt), ' agn: ', $
      nagn/float(nbpt), ' unk: ', nunk/float(nbpt)

    out[sf].final_class = 'sf'
    out[agn].final_class = 'agn'

    irred = where((strtrim(out.gal_type,2) eq 'ir_red'),nirred)
    red = where((strtrim(out.gal_type,2) eq 'red'),nred)
    blu = where((strtrim(out.gal_type,2) eq 'blu'),nblu)
    gre = where((strtrim(out.gal_type,2) eq 'gre'),ngre)
    splog, 'N(IR Red) = ', nirred, ', ', 'N(red) = ', nred, ', ', $
      'N(blue) = ', nblu, ', ', 'N(green) = ', ngre

    sfirred = where((strtrim(out.gal_type,2) eq 'ir_red') and $
      (strtrim(out.final_class,2) eq 'sf'),nsfirred)
    agnirred = where((strtrim(out.gal_type,2) eq 'ir_red') and $
      (strtrim(out.final_class,2) eq 'agn'),nagnirred)
    unkirred = where((strtrim(out.gal_type,2) eq 'ir_red') and $
      (strtrim(out.final_class,2) eq 'unclas'),nunkirred)
    splog, 'IR Red: ', ' sf: ', nsfirred, nsfirred/float(nirred), ' agn: ', $
      nagnirred, nagnirred/float(nirred), ' unk: ', nunkirred, $
      nunkirred/float(nirred)

    sfred = where((strtrim(out.gal_type,2) eq 'red') and $
      (strtrim(out.final_class,2) eq 'sf'),nsfred)
    agnred = where((strtrim(out.gal_type,2) eq 'red') and $
      (strtrim(out.final_class,2) eq 'agn'),nagnred)
    unkred = where((strtrim(out.gal_type,2) eq 'red') and $
      (strtrim(out.final_class,2) eq 'unclas'),nunkred)
    splog, 'Red: ', ' sf: ', nsfred, nsfred/float(nred), ' agn: ', $
      nagnred, nagnred/float(nred), ' unk: ', nunkred, $
      nunkred/float(nred)

    sfblu = where((strtrim(out.gal_type,2) eq 'blu') and $
      (strtrim(out.final_class,2) eq 'sf'),nsfblu)
    agnblu = where((strtrim(out.gal_type,2) eq 'blu') and $
      (strtrim(out.final_class,2) eq 'agn'),nagnblu)
    unkblu = where((strtrim(out.gal_type,2) eq 'blu') and $
      (strtrim(out.final_class,2) eq 'unclas'),nunkblu)
    splog, 'Blu: ', ' sf: ', nsfblu, nsfblu/float(nblu), ' agn: ', $
      nagnblu, nagnblu/float(nblu), ' unk: ', nunkblu, $
      nunkblu/float(nblu)

    sfgre = where((strtrim(out.gal_type,2) eq 'gre') and $
      (strtrim(out.final_class,2) eq 'sf'),nsfgre)
    agngre = where((strtrim(out.gal_type,2) eq 'gre') and $
      (strtrim(out.final_class,2) eq 'agn'),nagngre)
    unkgre = where((strtrim(out.gal_type,2) eq 'gre') and $
      (strtrim(out.final_class,2) eq 'unclas'),nunkgre)
    splog, 'Gre: ', ' sf: ', nsfgre, nsfgre/float(ngre), ' agn: ', $
      nagngre, nagngre/float(ngre), ' unk: ', nunkgre, $
      nunkgre/float(ngre)

    mwrfits, out, 'moustakas_obj_cat_11jul08.fits', /create
    
stop    
    
; write out    
    


; ---------------------------------------------------------------------------    
    
    
    kred = kate[where(strmatch(kate.gal_type,'red*',/fold))] ; ir-quiet
    kred = kate[where(strmatch(kate.gal_type,'red*',/fold))] ; ir-quiet
    kered = kate[where(strmatch(kate.gal_type,'ir_red*',/fold))] ; ir-excess
    kblu = kate[where(strmatch(kate.gal_type,'*blu*',/fold))]
    kgre = kate[where(strmatch(kate.gal_type,'*gre*',/fold))]
    
    spherematch, aa1.ra, aa1.dec, kate.ra, kate.dec, $
      1.0/3600.0, all, m2 & ii = ii1[all] & aa = aa1[all]
    
    spherematch, aa1.ra, aa1.dec, kred.ra, kred.dec, $
      1.0/3600.0, ared, m2 & red = ii1[ared] & kred = kred[m2]
    spherematch, aa1.ra, aa1.dec, kered.ra, kered.dec, $
      1.0/3600.0, aered, m2 & ered = ii1[aered] & kered = kered[m2]
    spherematch, aa1.ra, aa1.dec, kblu.ra, kblu.dec, $
      1.0/3600.0, ablu, m2 & blu = ii1[ablu] & kblu = kblu[m2]
    spherematch, aa1.ra, aa1.dec, kgre.ra, kgre.dec, $
      1.0/3600.0, agre, m2 & gre = ii1[agre] & kgre = kgre[m2]
    
    class = iclassification(ii,snrcut_class=3.0,/doplot,/kauffmann)
    unk = where(strtrim(class.bpt_class,2) eq 'Unknown')

    
    
    ired = iclassification(red,snrcut_class=3.0,/doplot)
    iered = iclassification(ered,snrcut_class=3.0,/doplot)
    iblu = iclassification(blu,snrcut_class=3.0,/doplot)
    igre = iclassification(gre,snrcut_class=3.0,/doplot)

    snrcut = 3.0

    these = where((ii[unk].nii_6584[0]/ii[unk].nii_6584[1] lt snrcut) and $
      (ii[unk].h_alpha[0]/ii[unk].h_alpha[1] gt snrcut) and $
      (ii[unk].oiii_5007[0]/ii[unk].oiii_5007[1] gt snrcut) and $
      (ii[unk].h_beta[0]/ii[unk].h_beta[1] gt snrcut))
    
    these = where($
      ((ii[unk].h_alpha[0]/ii[unk].h_alpha[1] gt snrcut) and $
      (ii[unk].h_beta[0]/ii[unk].h_beta[1] gt snrcut)) or $
      (ii[unk].oiii_5007[0]/ii[unk].oiii_5007[1] gt snrcut) or $
      (ii[unk].nii_6584[0]/ii[unk].nii_6584[1] gt snrcut))

    these = where($
      (ii.h_alpha[0]/ii.h_alpha[1] gt snrcut) and $
;     (ii.h_beta[0]/ii.h_beta[1] gt snrcut)) and $
;     (ii.oiii_5007[0]/ii.oiii_5007[1] lt snrcut) and $
      (ii.nii_6584[0]/ii.nii_6584[1] gt snrcut))

    these = where($
      ((ii[unk].h_alpha[0]/ii[unk].h_alpha[1] gt snrcut) and $
      (ii[unk].h_beta[0]/ii[unk].h_beta[1] lt snrcut)) and $
;     (ii[unk].oiii_5007[0]/ii[unk].oiii_5007[1] lt snrcut) and $
      (ii[unk].nii_6584[0]/ii[unk].nii_6584[1] gt snrcut))

    plot, alog10(ii[these].nii_6584_limit/ii[these].h_alpha[0]), $
      alog10(ii.oiii_5007[0]/ii.h_beta[0]), ps=4, $
      xrange=[-2,1.0], yrange=[-1.5,1.5], xsty=1, ysty=1
    kew = kewley_bpt_lines()
    oplot, kew.x_nii, kew.y_nii, line=2, thick=2.0

;   kred = kate[where(strmatch(kate.gal_type,'red*',/fold))] ; ir-quiet
;   kred = kate[where(strmatch(kate.gal_type,'red*',/fold))] ; ir-quiet
;   kered = kate[where(strmatch(kate.gal_type,'ir_red*',/fold))] ; ir-excess
;   kblu = kate[where(strmatch(kate.gal_type,'*blu*',/fold))]
;   kgre = kate[where(strmatch(kate.gal_type,'*gre*',/fold))]
;   help, kred, kered, kblu, kgre

; ---------------------------------------------------------------------------    
    
    badplates = [106,110,209,310,311]
    good = lonarr(n_elements(aa1))+1L
    good[where_array(badplates,aa1.pass)] = 0L
    
    keep = where((aa1.ubvri_absmag[2] le -19.4) and (aa1.z ge 0.15) and $
      (aa1.z lt 0.3) and (aa1.main_flag) and (good eq 1),nkeep)

    ii = ii1[keep]
    aa = aa1[keep]

    z = aa.z
    mv = aa.ubvri_absmag[1]
    uv = aa.ubvri_absmag[0]-aa.ubvri_absmag[2]
    f24 = aa.phot_mips24

    red = where(uv gt (1.55 - 0.25 - 0.04*(mv+20.0) - 0.42*(z-0.05) + 0.07*(z-0.05)^2),nred)
    green = where((uv le (1.55 - 0.25 - 0.04*(mv+20.0) - 0.42*(z-0.05) + 0.07*(z-0.05)^2)) and $
      (uv gt (1.35 - 0.25 - 0.04*(mv+20.0) - 0.42*(z-0.05) + 0.07*(z-0.05)^2)),ngreen)
    blue = where(uv le (1.35 - 0.25 - 0.04*(mv+20.0) - 0.42*(z-0.05) + 0.07*(z-0.05)^2),nblue)
    
    qq = where((f24 ge 0.3) and (uv gt (1.55 - 0.25 - 0.04*(mv+20.0) - 0.42*(z-0.05) + 0.07*(z-0.05)^2)),nqq)
    help, red, green, blue, qq

; cut

    snrcut = 3.0
    lred = where((ii[red].nii_6584[0]/ii[red].nii_6584[1] gt snrcut) and $
      (ii[red].h_alpha[0]/ii[red].h_alpha[1] gt snrcut) and $
      (ii[red].h_beta[0]/ii[red].h_beta[1] gt snrcut),nlred)
    lblue = where((ii[blue].nii_6584[0]/ii[blue].nii_6584[1] gt snrcut) and $
      (ii[blue].h_alpha[0]/ii[blue].h_alpha[1] gt snrcut) and $
      (ii[blue].h_beta[0]/ii[blue].h_beta[1] gt snrcut),nlblue)
    lgreen = where((ii[green].nii_6584[0]/ii[green].nii_6584[1] gt snrcut) and $
      (ii[green].h_alpha[0]/ii[green].h_alpha[1] gt snrcut) and $
      (ii[green].h_beta[0]/ii[green].h_beta[1] gt snrcut),nlgreen)
    help, lred, red, lgreen, green, lblue, blue

    ired = iclassification(ii[red[lred]],snrcut_class=0.1,/doplot)
    iblue = iclassification(ii[blue[lblue]],snrcut_class=0.1,/doplot)
    igreen = iclassification(ii[green[lgreen]],snrcut_class=0.1,/doplot)
    
    
; ---------------------------------------------------------------------------
; CMD
; ---------------------------------------------------------------------------

    postthick1 = 2.0
    postthick2 = 2.0
    
    mv = aa.ubvri_absmag[1]
    uv = aa.ubvri_absmag[0]-aa.ubvri_absmag[2]
    
    xrange = [-18.8,-22.8]
    yrange = [-0.3,1.9]

    xtitle = 'M_{V}'
    ytitle = 'U - V'

    dfpsplot, 'cmd.ps', /square, /color
    hogg_scatterplot, mv, uv, xsty=1, ysty=1, levels=errorf(0.5*[1.0,2.0,3.0]), $
      /outliers, /internal_weight, outsymsize=0.5, $
      xrange=xrange, yrange=yrange, xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=1.9, charthick=postthick2, xthick=postthick1, ythick=postthick1;, position=pos[*,0], $
    dfpsclose
    
    
    ages_zed=ages_cat2.z
    mu=ages_cat.synth_M_U
    mv=ages_cat.synth_M_V
    mb=ages_cat.synth_M_B
    class=ages_cat.class
    f24=ages_cat.phot_mips24
    f3p6=ages_cat.phot_ch1
    f4p5=ages_cat.phot_ch2
    f5p8=ages_cat.phot_ch3
    f8=ages_cat.phot_ch4

    red = where(mu-mv gt 1.55 - 0.25 - 0.04*(mv+20.0) - 0.42*(ages_zed-0.05) + 0.07*(ages_zed-0.05)^2,n_red)
    print,n_red,'red sequence galaxies'
    gre=where(mu-mv le 1.55 - 0.25 - 0.04*(mv+20.0) - 0.42*(ages_zed-0.05) + 0.07*(ages_zed-0.05)^2 and mu-mv gt 1.35 - 0.25 - 0.04*(mv+20.0) - 0.42*(ages_zed-0.05) + 0.07*(ages_zed-0.05)^2,n_gre)
    print,n_gre,'green valley galaxies'
    blu=where(mu-mv le 1.35 - 0.25 - 0.04*(mv+20.0) - 0.42*(ages_zed-0.05) + 0.07*(ages_zed-0.05)^2,n_blu)
    print,n_blu,'blue sequence galaxies'
    qq=where(f24 ge 0.3 and (mu-mv gt 1.55 - 0.25 - 0.04*(mv+20.0) - 0.42*(ages_zed-0.05) + 0.07*(ages_zed-0.05)^2),n_qq)
    print,n_qq,'infrared-bright red galaxies'



    
return
end
    
