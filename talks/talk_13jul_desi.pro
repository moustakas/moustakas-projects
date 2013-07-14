pro talk_13jul_desi
; jm13jul14siena - make some plots for the talk

    path = getenv('IM_RESEARCH_DIR')+'/talks/2013/13jul_desi/'

; ---------------------------------------------------------------------------
; plot the EW([OII]) distribution from DEEP2

    zz = read_deep2_zcat()
    ss = read_deep2(/ispec)
    ww = where(zz.z gt 0.7 and zz.z lt 1.4 and $
      ss.oii_3727_1[0]/ss.oii_3727_1[1] gt 5.0 and $
      ss.oii_3727_2[0]/ss.oii_3727_2[1] gt 5.0 and $
      ss.oii_3727_1_ew[0] gt 0.0 and $
      ss.oii_3727_2_ew[0] gt 0.0,ngal)
    zz = zz[ww]
    ss = ss[ww]

    psfile = path+'ewoii_hist.ps'
    im_plotconfig, 0, pos, psfile=psfile

    
    

    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
; ---------------------------------------------------------------------------
; see ISEDFIT_CALIBRATE_EMLINES
    
    datapath = getenv('IM_PROJECTS_DIR')+'/isedfit/'
;   path = getenv('IM_PROJECTS_DIR')+'/desi/'

    allatlas = mrdfits(datapath+'integrated_atlas_speclinefit_v3.1_nodust.fits.gz',1)
;   allatlas = mrdfits(datapath+'integrated_atlas_speclinefit_v3.1.fits.gz',1)
    allsdss = read_vagc_garching(/ispecline)
    allhiiinfo = read_hiiregions(linefit=allhii)

    snrcut = 5.0
    
; remove AGN from the galaxy samples and also require S/N>5 on all
; emission lines of interest 
    sclass = iclassification(allsdss,snrcut=snrcut,ratios=sratios)
    keep = where($
      strmatch(sratios.final_class,'*AGN*') eq 0 and $
;     strmatch(sratios.final_class,'*SF*') and $
      allsdss.oiii_5007[0]/allsdss.oiii_5007[1] gt snrcut and $
      allsdss.oii_3727[0]/allsdss.oii_3727[1] gt snrcut and $
      allsdss.h_beta[0]/allsdss.h_beta[1] gt snrcut and $
      allsdss.h_alpha[0]/allsdss.h_alpha[1] gt snrcut and $
      allsdss.nii_6584[0]/allsdss.nii_6584[1] gt snrcut and $
      allsdss.sii_6716[0]/allsdss.sii_6716[1] gt snrcut and $
      allsdss.sii_6731[0]/allsdss.sii_6731[1] gt snrcut,nsdss)
    splog, 'SDSS ', nsdss
    sdss = allsdss[keep]
    
    aclass = iclassification(allatlas,snrcut=snrcut,ratios=aratios)
    keep = where($
      strmatch(aratios.final_class,'*AGN*') eq 0 and $
;     strmatch(aratios.final_class,'*SF*') and $
      allatlas.oiii_5007[0]/allatlas.oiii_5007[1] gt snrcut and $
      allatlas.oii_3727[0]/allatlas.oii_3727[1] gt snrcut and $
      allatlas.h_beta[0]/allatlas.h_beta[1] gt snrcut and $
      allatlas.h_alpha[0]/allatlas.h_alpha[1] gt snrcut and $
      allatlas.nii_6584[0]/allatlas.nii_6584[1] gt snrcut and $
      allatlas.sii_6716[0]/allatlas.sii_6716[1] gt snrcut and $
      allatlas.sii_6731[0]/allatlas.sii_6731[1] gt snrcut,natlas)
    splog, 'ATLAS ', natlas
    atlas = allatlas[keep]

    keep = where($
      allhii.oiii_5007[0]/allhii.oiii_5007[1] gt snrcut and $
      allhii.oii_3727[0]/allhii.oii_3727[1] gt snrcut and $
      allhii.h_beta[0]/allhii.h_beta[1] gt snrcut and $
      allhii.h_alpha[0]/allhii.h_alpha[1] gt snrcut and $
      allhii.nii_6584[0]/allhii.nii_6584[1] gt snrcut and $
      allhii.sii_6716[0]/allhii.sii_6716[1] gt snrcut and $
      allhii.sii_6731[0]/allhii.sii_6731[1] gt snrcut,nhii)
    splog, 'HII ', nhii
    hii = allhii[keep]
    hiiinfo = allhiiinfo[keep]
    
; build the line-ratios
    lineratio, sdss, 'OIII_5007', 'H_BETA', 'OII_3727', $
      'H_BETA', soiiihb, soiiihb_err, soiihb, soiihb_err, snrcut=snrcut
    lineratio, atlas, 'OIII_5007', 'H_BETA', 'OII_3727', $
      'H_BETA', aoiiihb, aoiiihb_err, aoiihb, aoiihb_err, snrcut=snrcut
    lineratio, hii, 'OIII_5007', 'H_BETA', 'OII_3727', $
      'H_BETA', hoiiihb, hoiiihb_err, hoiihb, hoiihb_err, snrcut=snrcut

    lineratio, sdss, 'OIII_5007', 'H_BETA', 'NII_6584', $
      'H_ALPHA', soiiihb, soiiihb_err, sniiha, sniiha_err, snrcut=snrcut
    lineratio, atlas, 'OIII_5007', 'H_BETA', 'NII_6584', $
      'H_ALPHA', aoiiihb, aoiiihb_err, aniiha, aniiha_err, snrcut=snrcut
    lineratio, hii, 'OIII_5007', 'H_BETA', 'NII_6584', $
      'H_ALPHA', hoiiihb, hoiiihb_err, hniiha, hniiha_err, snrcut=snrcut

    lineratio, sdss, 'OIII_5007', 'H_BETA', ['SII_6716','SII_6731'], $
      'H_ALPHA', soiiihb, soiiihb_err, ssiiha, ssiiha_err, snrcut=snrcut
    lineratio, atlas, 'OIII_5007', 'H_BETA', ['SII_6716','SII_6731'], $
      'H_ALPHA', aoiiihb, aoiiihb_err, asiiha, asiiha_err, snrcut=snrcut
    lineratio, hii, 'OIII_5007', 'H_BETA', ['SII_6716','SII_6731'], $
      'H_ALPHA', hoiiihb, hoiiihb_err, hsiiha, hsiiha_err, snrcut=snrcut

; fit the median relations, giving equal weight to every sample
    xx = [soiiihb,aoiiihb,hoiiihb]
    weight = [soiiihb*0.0+1.0/nsdss,aoiiihb*0.0+1.0/natlas,hoiiihb*0.0+1.0/nhii]

    binsz = 0.1
    oiihbmed = im_medxbin(xx,[soiihb,aoiihb,hoiihb],binsz,weight=weight)
    niihamed = im_medxbin(xx,[sniiha,aniiha,hniiha],binsz,weight=weight)
    siihamed = im_medxbin(xx,[ssiiha,asiiha,hsiiha],binsz,weight=weight)

    oiihbcoeff = poly_fit(oiihbmed.medx,oiihbmed.medy,3)
    niihacoeff = poly_fit(niihamed.medx,niihamed.medy,3)
    siihacoeff = poly_fit(siihamed.medx,siihamed.medy,3)

    splog, '[OII]/Hb coeff ', oiihbcoeff
    splog, '[NII]/Ha coeff ', niihacoeff
    splog, '[SII]/Ha coeff ', siihacoeff
    
; make the plots
    nogrey = 0
    oiiihbrange = [-1.0,1.0]
    oiiihbaxis = range(oiiihbrange[0]+0.1,oiiihbrange[1]-0.1,200)
    
    psfile = path+'oiiihb_oiihb.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.7, charsize=1.9, height=5.0
    
; [OII]/Hb vs [OIII]/Hb
    im_hogg_scatterplot, soiiihb, soiihb, position=pos[*,0], nogrey=nogrey, $
      xsty=1, ysty=1, xrange=oiiihbrange, yrange=[-1,1.2], $
      ytitle=textoidl('log ([OII] \lambda3727 / H\beta)'), $
      xtitle=textoidl('log ([O III] \lambda5007 / H\beta)')
    djs_oplot, hoiiihb, hoiihb, psym=symcat(6,thick=4), color='orange', symsize=0.5
    djs_oplot, aoiiihb, aoiihb, psym=symcat(16), color='blue'
;   oploterror, oiihbmed.medx, oiihbmed.medy, oiihbmed.sigy, $
;     psym=symcat(9,thick=8), symsize=3.0, errthick=6.0

;   djs_oplot, oiihbmed.medx, oiihbmed.medy, line=0, thick=10
;   djs_oplot, oiihbmed.medx, oiihbmed.medy+0.2, line=5, thick=10
;   djs_oplot, oiihbmed.medx, oiihbmed.medy-0.2, line=5, thick=10

;   djs_oplot, oiihbmed.medx, oiihbmed.quant75, line=5, thick=10
;   djs_oplot, oiihbmed.medx, oiihbmed.quant25, line=5, thick=10
    djs_oplot, oiiihbaxis, poly(oiiihbaxis,oiihbcoeff), line=0, thick=10
    djs_oplot, oiiihbaxis, poly(oiiihbaxis,oiihbcoeff)+0.2, line=5, thick=10
    djs_oplot, oiiihbaxis, poly(oiiihbaxis,oiihbcoeff)-0.2, line=5, thick=10
    
    im_legend, ['SDSS Galaxies','Integrated Spectra','HII Regions'], $
      psym=[15,16,6], color=['grey','blue','orange'], /left, /top, $
      box=0, margin=0, symthick=10, charsize=1.8
    im_legend, 'Median Fit', /right, /bottom, box=0, line=0, $
      thick=10, pspacing=1.9, charsize=2.1, margin=0
    
    im_plotconfig, /psclose, psfile=psfile, /pdf
    
stop    
    
return
end
    
