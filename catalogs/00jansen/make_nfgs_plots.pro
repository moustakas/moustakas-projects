pro make_nfgs_plots, postscript=postscript
; jm02may28uofa
; jm04dec10uofa - updated
; generate plots from Jansen et al. 2001, ApJ, 551, 825
    
    root = '00jansen'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'

    nfgs = read_00jansen(datanodust=nfgs_cor)

    if keyword_set(postscript) then begin
       dfpsplot, path+'nfgs_plots.ps', /square
       postthick = 8.0
    endif else begin
       window, 0, xs=550, ys=550
       postthick = 2.0
    endelse

    plotsym, 0, 1, /fill

; ------------------------------------------------------------
; EW([N II]) versus [N II] ratio
; ------------------------------------------------------------

    indx = where((nfgs.nii_6548[1] gt 0.0) and (nfgs.nii_6584[1] gt 0.0) and $
      (nfgs.nii_6584_ew[1] gt 0.0),nindx)

    x = nfgs[indx].nii_6584_ew[0]
    y = nfgs[indx].nii_6584[0]/nfgs[indx].nii_6548[0]
    
    xtitle = 'EW([N II] \lambda6584)'
    ytitle = '[N II] \lambda6584 / [N II] \lambda6548'

    xrange = minmax(x)
    yrange = [0,6]
    
    djs_plot, x, y, ps=8, xsty=3, ysty=3, charsize=2.0, charthick=postthick, xthick=postthick, $
      ythick=postthick, xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, $
      thick=postthick, title=title, /xlog
    djs_oplot, 10^!x.crange, 3*[1,1], line=0, thick=postthick
    
    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; ------------------------------------------------------------
; EW([O III]) versus [O III] ratio
; ------------------------------------------------------------

    indx = where((nfgs.oiii_4959[1] gt 0.0) and (nfgs.oiii_5007[1] gt 0.0) and $
      (nfgs.oiii_5007_ew[1] gt 0.0),nindx)

    x = nfgs[indx].oiii_5007_ew[0]
    y = nfgs[indx].oiii_5007[0]/nfgs[indx].oiii_4959[0]
    
    xtitle = 'EW([O III] \lambda5007)'
    ytitle = '[O III] \lambda5007 / [O III] \lambda4959'

    xrange = minmax(x)
    yrange = [0,6]
    
    djs_plot, x, y, ps=8, xsty=3, ysty=3, charsize=2.0, charthick=postthick, xthick=postthick, $
      ythick=postthick, xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, $
      thick=postthick, title=title, /xlog
    djs_oplot, 10^!x.crange, 3*[1,1], line=0, thick=postthick
    
    if (not keyword_set(postscript)) then cc = get_kbrd(1)
    
    if keyword_set(postscript) then dfpsclose

stop
    
; --------
; FIGURE 1
; --------
; M_B versus log OII/Ha with EW(Ha)>10

    good = where((nfgs.oii_3727[1] gt 0.0) and (nfgs.h_alpha[1] gt 0.0) and (nfgs.H_ALPHA_EW[0] ge 10.0),ngood)

    x = nfgs[good].M_B
    y = alog10(nfgs[good].OII_3727[0]/nfgs[good].h_alpha[0])
    
    xrange = [-13,-23] & yrange = [-1.4,0.8]
    xtitle = 'M_{B}'
    ytitle = 'log ([OII]/H\alpha)_{obs}'
    title='Figure 1'

    plotsym, 0, 1, /fill
    djs_plot, [0], [0], xsty=1, ysty=1, charsize=2.0, charthick=postthick, xthick=postthick, $
      ythick=postthick, xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, $
      thick=postthick, title=title
    djs_oplot, x, y, ps=8
    coeff = [1.37,0.087] & fit = poly(x,coeff)
    djs_oplot, x, fit, line=0, thick=postthick

    splog, 'Points: '+string(ngood,format='(I0)')
    
    if (not keyword_set(postscript)) then cc = get_kbrd(1)
    
; ---------
; FIGURE 2a
; ---------
; M_B versus log OII/Ha with EW(Ha)>10 and EW(Hb)>5

    good = where((nfgs.oii_3727[1] gt 0.0) and (nfgs.h_alpha[1] gt 0.0) and $
      (nfgs.H_ALPHA_EW[0] ge 10.0) and (nfgs.H_BETA_EW[0] ge 5.0),ngood)

    x = nfgs[good].M_B
    y = alog10(nfgs[good].OII_3727[0]/nfgs[good].h_alpha[0])
    
    xrange = [-13,-23] & yrange = [-1.4,0.8]
    xtitle = 'M_{B}'
    ytitle = 'log ([OII]/H\alpha)_{obs}'
    title='Figure 2a'

    plotsym, 0, 1, /fill
    djs_plot, [0], [0], xsty=1, ysty=1, charsize=2.0, charthick=postthick, xthick=postthick, $
      ythick=postthick, xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, $
      thick=postthick, title=title
    djs_oplot, x, y, ps=8
    coeff = [1.40,0.087] & fit = poly(x,coeff)
    djs_oplot, x, fit, line=0, thick=postthick

    splog, 'Points: '+string(ngood,format='(I0)')

    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; ---------
; FIGURE 2b
; ---------
; M_B versus log OII/Ha corrected for dust

    good = where((nfgs_cor.OII_3727[1] gt 0.0) and (nfgs_cor.h_alpha[1] gt 0.0) and $
      (nfgs_cor.ebv_hahb_err gt 0.0) and (nfgs.H_ALPHA_EW[0] ge 10.0) and (nfgs.H_BETA_EW[0] ge 5.0),ngood)

    x = nfgs_cor[good].M_B
    y = alog10(nfgs_cor[good].OII_3727[0]/nfgs_cor[good].h_alpha[0])
    
    xrange = [-13,-23] & yrange = [-1.4,0.8]
    xtitle = 'M_B'
    ytitle = 'log ([OII]/H\alpha)_{cor}'
    title='Figure 2b'

    plotsym, 0, 1, /fill
    djs_plot, [0], [0], xsty=1, ysty=1, charsize=2.0, charthick=postthick, xthick=postthick, $
      ythick=postthick, xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, $
      thick=postthick, title=title
    djs_oplot, x, y, ps=8
    coeff = [0.65,0.035] & fit = poly(x,coeff)
    djs_oplot, x, fit, line=0, thick=postthick

    splog, 'Points: '+string(ngood,format='(I0)')
    
    if (not keyword_set(postscript)) then cc = get_kbrd(1)

; ---------
; FIGURE 2c
; ---------
; E(B-V) versus log OII/Ha_obs

    good = where((nfgs.OII_3727[1] gt 0.0) and (nfgs.h_alpha[1] gt 0.0) and (nfgs.H_ALPHA_EW[0] ge 10.0) and (nfgs.H_BETA_EW[0] ge 5.0) and $
      (nfgs.H_ALPHA_EW[0] ge 10.0) and (nfgs.H_BETA_EW[0] ge 5.0) and (nfgs_cor.ebv_hahb_err gt 0.0),ngood)

    x = nfgs_cor[good].ebv_hahb
    y = alog10(nfgs[good].OII_3727[0]/nfgs[good].h_alpha[0])
    
    xrange = [-0.1,0.9] & yrange = [-1.4,0.8]
    xtitle = 'E(B-V)'
    ytitle = 'log ([OII]/H\alpha)_{obs}'
    title='Figure 2c'

    plotsym, 0, 1, /fill
    djs_plot, [0], [0], xsty=1, ysty=1, charsize=2.0, charthick=postthick, xthick=postthick, $
      ythick=postthick, xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, $
      thick=postthick, title=title
    djs_oplot, x, y, ps=8

; overplot the reddening curve

    ratio = 0.0

    xebv = findgen((1.1-(-0.5))/0.1)*0.1+(-0.5)
    yseaton =  -0.4*xebv*(k_lambda(3727.0,/seaton)-k_lambda(6563.0,/seaton)) + ratio
    djs_oplot, xebv, yseaton, line=0, thick=postthick
    
    splog, 'Points: '+string(ngood,format='(I0)')

    if (not keyword_set(postscript)) then cc = get_kbrd(1)

    if keyword_set(postscript) then dfpsclose

    stop

return
end    
