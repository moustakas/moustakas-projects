; ---------------------------------------------------------------------------    
; JD's elliptical
; ---------------------------------------------------------------------------    

    ispec, 'fdwscra.1123_2.fits', /tracespec, aperture=3.0, traceorder=2, $
      /wfits, /noplot, outname='ngc3962_nuclear'
    specfile = 'ngc3962_nuclear_003.ms.fits'
    head = headfits(specfile)
    sxaddpar, head, 'Z', float(0.006054), 'NED redshift', before='HISTORY'
    modfits, specfile, 0, head

    specdata = ispeclinefit(specfile,specres=8.0,snrcut=1.0,dustmodel=0,$
      datapath=datapath,linepath=atlas_path(/specfit),suffix=suffix,Zmulti=Zmulti,$
      /charlot,/zcrosscor,/postscript,/write,vmaxshift=500.0,/nologfile,$
      starvdisp=100.0)

    hii = read_hii_regions(/nosdss,/samplerefs)
    models_kauffmann = kewley_bpt_lines(/kauffmann)
    models = kewley_bpt_lines()

    dfpsplot, 'ngc3962_bpt.ps', /square, /color
;   im_window, 0, xratio=0.5, /square

; NII/Ha vs OIII/Hb    
    
    djs_plot, models.x_nii, models.y_nii, xsty=3, ysty=3, thick=5.0, $
      xrange=[-2.3,1.0], yrange=[-1.5,1.5], xthick=5.0, ythick=5.0, $
      charsize=1.8, charthick=5.0, xtitle='log ([N II]/H\alpha)', $
      ytitle='log ([O III]/H\beta)'
    oplot, models_kauffmann.x_nii, models_kauffmann.y_nii, line=2, thick=5.0

    plotsym, 0, 0.5, /fill
    indx = where((hii.nii_6584_h_alpha gt -900) and (hii.oiii_5007_h_beta gt -900))
    djs_oplot, hii[indx].nii_6584_h_alpha, hii[indx].oiii_5007_h_beta, ps=8, color='blue'

    plotsym, 8, 2.0, /fill
    plots, alog10(specdata.nii_6584[0]/specdata.h_alpha[0]), alog10(specdata.oiii_5007[0]/specdata.h_beta[0]), $
      psym=8, color=djs_icolor('red')

; SII/Ha vs OIII/Hb    
    
    djs_plot, models.x_sii, models.y_sii, xsty=3, ysty=3, thick=5.0, $
      xrange=[-2.3,1.0], yrange=[-1.5,1.5], xthick=5.0, ythick=5.0, $
      charsize=1.8, charthick=5.0, xtitle='log ([S II]/H\alpha)', $
      ytitle='log ([O III]/H\beta)'
    oplot, models_kauffmann.x_sii, models_kauffmann.y_sii, line=2, thick=5.0

    plotsym, 0, 0.5, /fill
    indx = where((hii.sii_h_alpha gt -900) and (hii.oiii_5007_h_beta gt -900))
    djs_oplot, hii[indx].sii_h_alpha, hii[indx].oiii_5007_h_beta, ps=8, color='blue'

    plotsym, 8, 2.0, /fill
    plots, alog10((specdata.sii_6716[0]+specdata.sii_6731[0])/specdata.h_alpha[0]), $
      alog10(specdata.oiii_5007[0]/specdata.h_beta[0]), psym=8, color=djs_icolor('red')

; OI/Ha vs OIII/Hb    
    
    djs_plot, models.x_oi, models.y_oi, xsty=3, ysty=3, thick=5.0, $
      xrange=[-3.0,0.0], yrange=[-1.5,1.5], xthick=5.0, ythick=5.0, $
      charsize=1.8, charthick=5.0, xtitle='log ([O I]/H\alpha)', $
      ytitle='log ([O III]/H\beta)'
    oplot, models_kauffmann.x_oi, models_kauffmann.y_oi, line=2, thick=5.0

;   plotsym, 0, 0.5, /fill
;   indx = where((hii.oi_6300_h_alpha gt -900) and (hii.oiii_5007_h_beta gt -900))
;   djs_oplot, hii[indx].oi_6300_h_alpha, hii[indx].oiii_5007_h_beta, ps=8, color='blue'

    plotsym, 8, 2.0, /fill
    plots, alog10(specdata.oi_6300[0]/specdata.h_alpha[0]), alog10(specdata.oiii_5007[0]/specdata.h_beta[0]), $
      psym=8, color=djs_icolor('red')

    dfpsclose
    
