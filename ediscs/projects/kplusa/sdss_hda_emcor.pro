pro sdss_hda_emcor, nodust, postscript=postscript

    sdss = read_sdss_main(/ispec)
    if (n_elements(nodust) eq 0L) then nodust = $
      iunred_linedust(sdss,snrcut=5.0,/silent)

; plotting variables    
    
    postthick1 = 2.0
    postthick2 = 2.0
    textcolor1 = 'white'

;   pspath = getenv('PAPERSPATH')+'/literature/'
    pspath = './'
    if keyword_set(postscript) then begin
       postthick1 = 4.0
       postthick2 = 3.0
       textcolor1 = 'black'
    endif

    charsize_0 = 1.0
    charsize_1 = 1.1
    charsize_2 = 1.2
    charsize_3 = 1.3
    charsize_4 = 1.4
    charsize_5 = 1.5
    charsize_6 = 1.6
    charsize_7 = 1.7
    charsize_8 = 1.8
    charsize_9 = 1.9
    singlecharsize_0 = 2.0
    singlecharsize_1 = 2.1
    singlecharsize_2 = 2.2
    singlecharsize_3 = 2.3
    singlecharsize_4 = 2.4
    singlecharsize_5 = 2.5
    charsize_30 = 3.0

    armexten = '.ps'
    
; ---------------------------------------------------------------------------    
; selection plot + BPT diagram
; ---------------------------------------------------------------------------    

    psname = 'pb_agn'
    if keyword_set(postscript) or keyword_set(pdf) then begin
       psfile = pspath+psname+armexten
       if keyword_set(postscript) then splog, 'Writing '+psfile
    endif else delvarx, psfile

    xmargin = [1.1,0.3] & ymargin = [0.3,1.0]
    xspace = 1.0 & yspace = 0.0
    width = 4.5*[1,1] & height = 4.5
    xpage = total(height)+total(ymargin)+total(yspace)
    ypage = total(width)+total(xmargin)+total(xspace)

    xspace1 = xspace & yspace1 = yspace
    xmargin1 = xmargin & ymargin1 = ymargin
    width1 = width & height1 = height
    xpage1 = xpage & ypage1 = ypage
    
    arm_plotconfig, landscape=1, nx=2, ny=1, xmargin=xmargin1, $
      ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
      height=height1, coord=pos, xpage=xpage1, ypage=ypage1, $
      psfile=psfile, /writeover, bw=0l;, /show
    cleanplot, /silent

    xtitle1 = 'log EW([O II]) (\AA)'
    ytitle1 = 'H\delta_{A} (\AA, corrected)'
    xrange1 = alog10([0.1,500])
    yrange1 = [-3,9]

    xtitle2 = 'log ([N II] \lambda6584/H\alpha)'
    ytitle2 = 'log ([O III] \lambda5007/H\beta)'
    xrange2 = [-1.7,0.9]
    yrange2 = [-1.1,1.2]

    snrcut = 3.0
    allindx = where((sdss.oii_3727_ew[1] gt 0.0) and $
      (sdss.nii_6584[0]/sdss.nii_6584[1] gt snrcut) and $
      (sdss.oiii_5007[0]/sdss.oiii_5007[1] gt snrcut) and $
      (sdss.h_alpha[0]/sdss.h_alpha[1] gt snrcut) and $
      (sdss.h_beta[0]/sdss.h_beta[1] gt snrcut),nallindx)

    allniiha = alog10(sdss[allindx].nii_6584[0]/sdss[allindx].h_alpha[0])
    alloiiihb = alog10(sdss[allindx].oiii_5007[0]/sdss[allindx].h_beta[0])
    allewoii = alog10(sdss[allindx].oii_3727_ew[0])
    allhda = sdss[allindx].lick_hd_a_model[0]

    indx_bpt = where((sdss.oii_3727_ew[1] gt 0.0) and $
      (sdss.nii_6584[0]/sdss.nii_6584[1] gt snrcut) and $
      (sdss.oiii_5007[0]/sdss.oiii_5007[1] gt snrcut) and $
      (sdss.h_alpha[0]/sdss.h_alpha[1] gt snrcut) and $
      (sdss.h_beta[0]/sdss.h_beta[1] gt snrcut) and $
      (sdss.oii_3727_ew[0] lt 0.9) and $
      (sdss.lick_hd_a_model[0] gt 3.0),nindx_bpt)
    indx = where((sdss.oii_3727_ew[1] gt 0.0) and $
      (sdss.oii_3727_ew[0] lt 0.9) and $
      (sdss.lick_hd_a_model[0] gt 3.0),nindx)

    cc = iclassification(sdss[indx],$
      snrcut_class=snrcut,ratios=rr,/kauffmann)
    sf = where(cc.bpt_nii_kauffmann_class eq 4)
    agn = where(cc.bpt_nii_kauffmann_class eq 2)
    unk = where(cc.bpt_nii_kauffmann_class eq 0)

    niiha = rr.nii_ha
    niiha_err = rr.nii_ha_err
    oiiihb = rr.oiii_hb
    oiiihb_err = rr.oiii_hb_err
;   niiha = alog10(sdss[indx].nii_6584[0]/sdss[indx].h_alpha[0])
;   oiiihb = alog10(sdss[indx].oiii_5007[0]/sdss[indx].h_beta[0])
    ewoii = alog10(sdss[indx].oii_3727_ew[0])
    hda = sdss[indx].lick_hd_a_model[0]

    im_symbols, 108, psize=0.15, fill=1, color=fsc_color(textcolor1,100)
    im_hogg_scatterplot, allewoii, allhda, outliers=1, xsty=1, ysty=1, $
      label=0, outpsym=8, outsymsize=1.0, outcolor=fsc_color(textcolor1,100), $
      levels=errorf(0.5*[1.0,2.0,3.0]), $
      xrange=xrange1, yrange=yrange1, xtitle=textoidl(xtitle1), ytitle=textoidl(ytitle1), $
      charsize=charsize_4, charthick=postthick2, xthick=postthick1, ythick=postthick1, position=pos[*,0], $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(textcolor1,100), $
      contour_color=fsc_color(textcolor1,100)

    im_symbols, 108, psize=0.8, fill=1, color=fsc_color('forest green',101)
    djs_oplot, ewoii[agn], hda[agn], ps=8
    im_symbols, 108, psize=0.8, fill=1, color=fsc_color('dodger blue',102)
    djs_oplot, ewoii[sf], hda[sf], ps=8
    im_symbols, 108, psize=0.5, fill=1, color=fsc_color('firebrick',103)
    djs_oplot, ewoii[unk], hda[unk], ps=8

    djs_oplot, [xrange1[0],alog10(0.9)], 3.0*[1,1], line=5, $
      thick=postthick1
    djs_oplot, alog10(0.9)*[1,1], [3.0,yrange1[1]], line=5, $
      thick=postthick1
    
    im_symbols, 108, psize=0.15, fill=1, color=fsc_color(textcolor1,100)
    im_hogg_scatterplot, allniiha, alloiiihb, /noerase, outliers=1, xsty=1, ysty=1, $
      label=0, outpsym=8, outsymsize=1.0, outcolor=fsc_color(textcolor1,100), $
      levels=errorf(0.5*[1.0,2.0,3.0]), $
      xrange=xrange2, yrange=yrange2, xtitle=textoidl(xtitle2), ytitle=textoidl(ytitle2), $
      charsize=charsize_4, charthick=postthick2, xthick=postthick1, ythick=postthick1, position=pos[*,1], $
      color=fsc_color(textcolor1,100), axis_color=fsc_color(textcolor1,100), $
      contour_color=fsc_color(textcolor1,100)

    models = kewley_bpt_lines(/kauffmann,_extra=extra)
    oplot, models.x_nii, models.y_nii, line=0, thick=postthick1
    models = kewley_bpt_lines(_extra=extra)
    oplot, models.x_nii, models.y_nii, line=2, thick=postthick1
    
    im_symbols, 108, psize=1.0, fill=1, color=fsc_color('forest green',101)
    oploterror, niiha[agn], oiiihb[agn], niiha_err[agn], oiiihb_err[agn], ps=8, $
      errthick=postthick1, errcolor=fsc_color('forest green',101)
    im_symbols, 108, psize=1.0, fill=1, color=fsc_color('dodger blue',102)
    oploterror, niiha[sf], oiiihb[sf], niiha_err[sf], oiiihb_err[sf], ps=8, $
      errthick=postthick1, errcolor=fsc_color('dodger blue',101)
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    


return
end
