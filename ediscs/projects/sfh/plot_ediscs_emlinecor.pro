pro plot_ediscs_emlinecor, ps=ps
; jm10may03ucsd -  develop the emission-line corrections

    ps = 1
    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

    sfhpath = ediscs_path(/projects)+'sfh/'
    paperpath = sfhpath
;   cluster = read_ediscs_sfh_sample(/cluster)

    info = read_ediscs(/spec1d)
    ppxf = read_ediscs(/ppxf)
    ppxf = ppxf[where(strmatch(info.memberflag,'*1*'))]

; cross-match against Patricia's emission-line free red galaxies 
    readcol, sfhpath+'greg/table3.dat', gal, format='A', $
      comment='#', /silent
    match, strtrim(ppxf.galaxy,2), strtrim(gal,2), isred, m2
    isblue = cmset_op(lindgen(n_elements(ppxf)),'and',/not2,isred)

; QAplot of the spectra that nominally have H-delta well-detected
   hd = where(ppxf.h_delta[1] gt 0.0)
   ww = cmset_op(isred,'and',hd)
   qaplot_ediscs_gandalf_specfit, ppxf[hd[ww]], psfile=sfhpath+'qaplots/qaplot_hdelta.ps'
stop    
; --------------------------------------------------
; make a bunch of plots
    psfile = paperpath+'emline_cor'+suffix
    im_plotconfig, 6, pos, psfile=psfile, yspace=1.0, xmargin=[1.3,0.2]

; ###############
; D(4000) BC03 vs raw
    xtitle1 = 'D_{n}(4000)_{raw}'
    ytitle1 = 'D_{n}(4000)_{BC03}'
    xtitle2 = 'EW(H\beta) (\AA)'
    ytitle2 = 'D_{n}(4000)_{raw}-D_{n}(4000)_{BC03}'

    xrange1 = [0.8,2.4]
    yrange1 = xrange1
    xrange2 = [0.1,100]
    yrange2 = [-0.4,0.4]
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      xrange=xrange1, yrange=yrange1, xtitle=xtitle1, ytitle=ytitle1
    djs_oplot, ppxf[isblue].d4000_narrow_raw[0], ppxf[isblue].d4000_narrow_model[0], $
      psym=symcat(9,thick=6), color=fsc_color('dodger blue',101)
    djs_oplot, ppxf[isred].d4000_narrow_raw[0], ppxf[isred].d4000_narrow_model[0], $
      psym=symcat(16), color=fsc_color('firebrick',101)
    djs_oplot, !x.crange, !y.crange, line=0, thick=4

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      xrange=xrange2, yrange=yrange2, xtitle=xtitle2, ytitle=ytitle2, /xlog
    resid = ppxf.d4000_narrow_model[0]-ppxf.d4000_narrow_raw[0]
    djs_oplot, ppxf[isblue].h_beta_ew[0], resid[isblue], $
      psym=symcat(9,thick=6), color=fsc_color('dodger blue',101)
    djs_oplot, ppxf[isred].h_beta_ew[0], resid[isred], $
      psym=symcat(16), color=fsc_color('firebrick',101)
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=4

; ###############
; Hd_A, BC03 vs raw
    xtitle1 = 'H\delta_{A,raw}'
    ytitle1 = 'H\delta_{A,BC03}'
    xtitle2 = 'EW(H\beta) (\AA)'
    ytitle2 = 'H\delta_{A,raw}-H\delta_{BC03}'

    xrange1 = [-6,12]
    yrange1 = xrange1
    xrange2 = [0.1,100]
    yrange2 = [-8,8]
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      xrange=xrange1, yrange=yrange1, xtitle=xtitle1, ytitle=ytitle1
    djs_oplot, ppxf[isblue].lick_hd_a_raw[0], ppxf[isblue].lick_hd_a_model[0], $
      psym=symcat(9,thick=6), color=fsc_color('dodger blue',101)
    djs_oplot, ppxf[isred].lick_hd_a_raw[0], ppxf[isred].lick_hd_a_model[0], $
      psym=symcat(16), color=fsc_color('firebrick',101)
    djs_oplot, !x.crange, !y.crange, line=0, thick=4

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      xrange=xrange2, yrange=yrange2, xtitle=xtitle2, ytitle=ytitle2, /xlog
    resid = ppxf.lick_hd_a_model[0]-ppxf.lick_hd_a_raw[0]
    djs_oplot, ppxf[isblue].h_beta_ew[0], resid[isblue], $
      psym=symcat(9,thick=6), color=fsc_color('dodger blue',101)
    djs_oplot, ppxf[isred].h_beta_ew[0], resid[isred], $
      psym=symcat(16), color=fsc_color('firebrick',101)
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=4

; ###############
; Hd_A, BC03 vs cor
    xtitle1 = 'H\delta_{A,cor}'
    ytitle1 = 'H\delta_{A,BC03}'
    xtitle2 = 'EW(H\beta) (\AA)'
    ytitle2 = 'H\delta_{A,cor}-H\delta_{BC03}'

    xrange1 = [-6,12]
    yrange1 = xrange1
    xrange2 = [0.1,100]
    yrange2 = [-8,8]
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      xrange=xrange1, yrange=yrange1, xtitle=xtitle1, ytitle=ytitle1
    djs_oplot, ppxf[isblue].lick_hd_a_cor[0], ppxf[isblue].lick_hd_a_model[0], $
      psym=symcat(9,thick=6), color=fsc_color('dodger blue',101)
    djs_oplot, ppxf[isred].lick_hd_a_cor[0], ppxf[isred].lick_hd_a_model[0], $
      psym=symcat(16), color=fsc_color('firebrick',101)
    djs_oplot, !x.crange, !y.crange, line=0, thick=4

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      xrange=xrange2, yrange=yrange2, xtitle=xtitle2, ytitle=ytitle2, /xlog
    resid = ppxf.lick_hd_a_model[0]-ppxf.lick_hd_a_cor[0]
    djs_oplot, ppxf[isblue].h_beta_ew[0], resid[isblue], $
      psym=symcat(9,thick=6), color=fsc_color('dodger blue',101)
    djs_oplot, ppxf[isred].h_beta_ew[0], resid[isred], $
      psym=symcat(16), color=fsc_color('firebrick',101)
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=4

; ###############
; Hd_A, raw vs cor
    xtitle1 = 'H\delta_{A,cor}'
    ytitle1 = 'H\delta_{A,raw}'
    xtitle2 = 'EW(H\beta) (\AA)'
    ytitle2 = 'H\delta_{A,cor}-H\delta_{raw}'

    xrange1 = [-6,12]
    yrange1 = xrange1
    xrange2 = [0.1,100]
    yrange2 = [-8,8]
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      xrange=xrange1, yrange=yrange1, xtitle=xtitle1, ytitle=ytitle1
    djs_oplot, ppxf[isblue].lick_hd_a_cor[0], ppxf[isblue].lick_hd_a_raw[0], $
      psym=symcat(9,thick=6), color=fsc_color('dodger blue',101)
    djs_oplot, ppxf[isred].lick_hd_a_cor[0], ppxf[isred].lick_hd_a_raw[0], $
      psym=symcat(16), color=fsc_color('firebrick',101)
    djs_oplot, !x.crange, !y.crange, line=0, thick=4

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      xrange=xrange2, yrange=yrange2, xtitle=xtitle2, ytitle=ytitle2, /xlog
    resid = ppxf.lick_hd_a_raw[0]-ppxf.lick_hd_a_cor[0]
    djs_oplot, ppxf[isblue].h_beta_ew[0], resid[isblue], $
      psym=symcat(9,thick=6), color=fsc_color('dodger blue',101)
    djs_oplot, ppxf[isred].h_beta_ew[0], resid[isred], $
      psym=symcat(16), color=fsc_color('firebrick',101)
    djs_oplot, 10^!x.crange, [0,0], line=0, thick=4

    im_plotconfig, /psclose
    
stop    

return
end
