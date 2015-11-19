pro nsf_2014_nov
; jm14sep17siena - build plots for my 2014/Nov NSF proposal

    common com_redmapper, bcgs

    nsfpath = getenv('IM_PAPERS_DIR')+'/grants/nsf-2014-nov/'
    redmapperpath = redmapper_path(ver=ver)
    massprofpath = bcgmstar_path(/massprofiles)

; read the catalogs
    if n_elements(bcgs) eq 0L then bcgs = mrdfits(redmapperpath+$
      'redmapper_isedfit_'+ver+'_centrals.fits.gz',1)
    zbins = redbaryons_zbins(nzbins)

; fit a power-law to the sample of BCGs
    these = where(bcgs.z gt 0.1 and bcgs.z lt 0.3,nbcgs)
    mbcg_redmapper = bcgs[these].mstar_50
    m500_redmapper = rykoff_mass(bcgs[these].lambda_chisq)
    redmapper_coeff = robust_linefit(m500_redmapper,mbcg_redmapper,yfit,redmapper_sig)

    redmapper_med = im_medxbin(m500_redmapper,mbcg_redmapper,0.05,/ver)

    big = where(bcgs[these].lambda_chisq gt 100)
    
; --------------------------------------------------    
; BCG stellar mass vs virial mass using the CLASH and redmapper
; samples 
    clash = read_bcgmstar_sample(/getmbcg)
    ncl = n_elements(clash)

; make the plot    
    xrange = [13.5,15.3]
    yrange = [10.5,13.0]
    
    psfile = nsfpath+'mbcg_m500.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, xmargin=[1.3,0.4], width=6.8

    djs_plot, [0], [0], /nodata, position=pos, xsty=5, ysty=5, $
      xrange=xrange, yrange=yrange

; overplot the power-law fit from Kravtsov (Table 2)    
    scat = 0.17
    xx = range(13.6,15.2,50)
;   xx = range(xrange[0],xrange[1],50)
    yy = poly(xx-14.5,[12.24,0.33])
;   polyfill, [xx,reverse(xx)], [yy+scat,reverse(yy-scat)], $
;     /fill, color=cgcolor('light grey'), noclip=0
    djs_oplot, xx, yy, line=3, thick=8;, color=cgcolor('firebrick')

; overplot the redmapper sample
    xx = range(13.5,15.1,50)
    yy = poly(xx,redmapper_coeff)
;   djs_oplot, xx, yy, line=0, color=cgcolor('navy'), thick=8
;   djs_oplot, xx, yy+redmapper_sig, line=2, color=cgcolor('navy'), thick=8
;   djs_oplot, xx, yy-redmapper_sig, line=2, color=cgcolor('navy'), thick=8
;   polyfill, [xx,reverse(xx)], [yy+redmapper_sig,reverse(yy-redmapper_sig)], $
;     /fill, color=cgcolor('powder blue'), noclip=0
    djs_oplot, redmapper_med.medx, redmapper_med.medy, $
      line=0, thick=8, color=cgcolor('dodger blue')
    djs_oplot, redmapper_med.medx, redmapper_med.medy+redmapper_med.sigy, line=5, $
      thick=8, color=cgcolor('dodger blue')
    djs_oplot, redmapper_med.medx, redmapper_med.medy-redmapper_med.sigy, line=5, $
      thick=8, color=cgcolor('dodger blue')
;   djs_oplot, redmapper_med.medx, redmapper_med.quant75, line=5, thick=6
;   djs_oplot, redmapper_med.medx, redmapper_med.quant25, line=5, thick=6
    djs_oplot, m500_redmapper[big], mbcg_redmapper[big], psym=symcat(17), $
      color=cgcolor('dodger blue'), symsize=1
    djs_oplot, m500_redmapper[big], mbcg_redmapper[big], psym=symcat(5), $
      color=cgcolor('navy'), symsize=1.1
      
; overlay Kravtsov (h=0.7, Chabrier)
    mbcg_krav = [3.12,4.14,3.06,1.47,0.79,1.26,1.09,0.91,1.38]*1D12
    mbcg_krav_err = [0.36,0.3,0.3,0.13,0.05,0.11,0.06,0.05,0.14]*1D12
    m500_krav = [15.6,10.3,7,5.34,2.35,1.86,1.34,0.46,0.47]*1D14

    oploterror, alog10(m500_krav), alog10(mbcg_krav), mbcg_krav_err/mbcg_krav/alog(10), $
      psym=symcat(15), color=cgcolor('tomato'), errcolor=cgcolor('tomato'), $
      symsize=1.5

; overplot Gonzalez+13 (Table 3)
    mbcg_gonz = [0.84,0.87,0.33,0.57,0.85,0.60,0.86,0.93,0.71,0.81,0.70,0.57]*1D12*2.65
    mbcg_gonz_err = [0.03,0.09,0.01,0.01,0.14,0.03,0.03,0.05,0.07,0.12,0.02,0.01]*1D12*2.65
;     [2.57,3.11,1.30,2.81,2.13,1.33,1.72,3.33,2.42,2.65,2.46,1.30]
;   mbcg_gonz = [0.68,0.82,0.35,0.74,0.56,0.35,0.46,0.88,0.64,0.7,0.65,0.35]*1D12*mbcg_frac
;   mbcg_gonz_err = [0.04,0.06,0.03,0.06,0.05,0.02,0.02,0.06,0.05,0.06,0.04,0.02]*1D13
    m500_gonz = [2.26,5.15,0.95,3.46,3.59,0.99,0.95,3.23,2.26,2.41,2.37,1.45]*1D14
    m500_gonz_err = [0.19,0.42,0.1,0.32,0.28,0.11,0.1,0.19,0.23,0.18,0.24,0.21]*1D14
    
    oploterror, alog10(m500_gonz), alog10(mbcg_gonz), mbcg_gonz_err/mbcg_gonz/alog(10), $
      psym=symcat(14), color=cgcolor('forest green'), errcolor=cgcolor('forest green'), $
      symsize=1.8

; overplot CLASH    
    oploterror, clash.m500, clash.mbcg-0.25, clash.m500_err, clash.mbcg_err, $
      psym=symcat(16), symsize=1.6, color=cgcolor('grey')
    oploterror, clash.m500, clash.mbcg-0.25, clash.m500_err, clash.mbcg_err, $
      psym=symcat(9), symsize=1.7

; make a legend
    im_legend, ['CLASH','Kravtsov+14','Gonzalez+13','redMaPPer (\lambda>100)'], /right, /bottom, box=0, $
      psym=[16,15,14,17], color=['grey','tomato','forest green','dodger blue'], spacing=2.2, $
      charsize=1.6, symsize=[1.0,1.0,1.3,1.0]*1.5
    im_legend, ['CLASH','Kravtsov+14','Gonzalez+13','redMaPPer (\lambda>100)'], /right, /bottom, box=0, $
      psym=[9,15,14,5], color=['black','tomato','forest green','navy'], spacing=2.2, $
      charsize=1.6, symsize=[1.0,1.0,1.3,1.0]*1.5, symthick=5

;   im_legend, ['Kravtsov+14'], /right, /bottom, box=0, $
;     psym=[16,15], color=['black','tomato'], spacing=2.5
    xyouts, 13.85, 11.06, 'redMaPPer', align=0.5, /data, $
      charsize=1.6, color=cgcolor('dodger blue')
    xyouts, 13.85, 10.93, '(0.1<z<0.3)', align=0.5, /data, $
      charsize=1.3, color=cgcolor('dodger blue')

;   xyouts, 14.3, 12.45, textoidl('\propto'+'M_{500}^{0.33}'), align=0.0, $
;     /data, charsize=2.0, orientation=10
    im_legend, 'M_{*,BCG}\propto'+'M_{500}^{0.33}', /left, /top, box=0, line=3, $
      pspacing=1.7, thick=8
    
;   im_hogg_scatterplot, rykoff_mass(bcgs.lambda_chisq), bcgs.mstar_50, $
;     /overplot, /nogrey, xrange=xrange, yrange=yrange, outliers=0, $
;     /internal, levels=[0.1,0.25,0.5,0.75,0.9,0.95], xnpix=101, ynpix=101, $
;     /noerase, contour_color=cgcolor('navy')
    
; redraw the axes
    djs_plot, [0], [0], /nodata, /noerase, position=pos, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, $
      xtitle='log_{10} (M_{500} / M_{\odot})', $
      ytitle='log_{10} (M_{*,BCG} / M_{\odot})'
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop
    
    
return
end
    
