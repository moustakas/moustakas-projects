pro oplot_imf, mass, gamma=gamma, mcut=mcut, _extra=extra
    timf = imf(mass,gamma,mcut)
;   timf = 1E3*timf/total(timf,/double)
    timf = timf/interpol(timf,mass,1.0)
    good = where(timf ne 0)
    djs_oplot, alog10(mass[good]), alog10(timf[good]), $
      _extra=extra
return
end    

pro plot_cosmicimf_imfs, ps=ps, keynote=keynote
; jm10mar29ucsd - plot some example IMFs

    cosmicimfpath = ages_path(/projects)+'cosmicimf/'
    paperpath = ages_path(/papers)+'cosmicimf/'
    if keyword_set(keynote) then paperpath = $
      getenv('RESEARCHPATH')+'/meetings/10apr_florida/keynote/'

    peg = read_cosmicimf_sample(/pegase)
    dmsalp = mrdfits(cosmicimfpath+'dmsalp.fits.gz',1)
    nimf = n_elements(peg)

    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

; --------------------------------------------------
; show the dependence of various parameters on the high-mass slope 
    psfile = paperpath+'props_vs_gamma'+suffix
    im_plotconfig, 4, pos, psfile=psfile, xmargin=[1.3,0.2], $
      keynote=keynote
    if keyword_set(keynote) then begin
       keycolor = djs_icolor('white')
       hacolor = 'green'
       l24color = 'dodger blue'
       ircolor = 'firebrick'
       symthick2 = 10
       salpcolor = keycolor
    endif else begin
       hacolor = 'forest green'
       l24color = 'navy'
       ircolor = 'firebrick'
       symthick2 = 8.0
       salpcolor = 'grey'
    endelse

    psym = 6
    symsize1 = 1.8
    symsize2 = 3.0
    symthick1 = 5.0
    gamma = peg.slope-1
    xrange = minmax(gamma)

; delta mass with respect to Salpeter    
    djs_plot, [0], [0], /nodata, position=pos[*,0], $
      xsty=3, ysty=1, xrange=xrange, yrange=[-0.2,0.6], $
      xtitle='', ytitle='\Delta'+'log(m_{salp})', $
      xtickname=replicate(' ',10)
    polyfill, [gamma[1:nimf-1],reverse(gamma[1:nimf-1])], $
      [dmsalp[1:nimf-1].dmsalp_sfhgrid03-dmsalp[1:nimf-1].dmsalp_sfhgrid03_err/2,$
      reverse(dmsalp[1:nimf-1].dmsalp_sfhgrid03+dmsalp[1:nimf-1].dmsalp_sfhgrid03_err/2)], $
      /data, /fill, color=fsc_color('tan',100), noclip=0
;   djs_oplot, gamma[1:nimf-1], dmsalp[1:nimf-1].dmsalp_sfhgrid03, thick=symthick2
;   djs_oplot, gamma[1:nimf-1], dmsalp[1:nimf-1].dmsalp, $
;     psym=-symcat(psym,thick=symthick1), symsize=symsize1, $
;     thick=symthick1
    djs_oplot, !x.crange, [0,0], line=1, color=salpcolor, thick=symthick2

; SFR conversion factors    
    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], $
      xsty=3, ysty=3, xrange=xrange, yrange=[-0.99,1.99], $
      xtitle='', ytitle=' \Delta'+'log(\psi/L)_{salp}', $
      xtickname=replicate(' ',10)
    djs_oplot, gamma[1:nimf-1], alog10(peg[1:nimf-1].c_ha/peg[0].c_ha), $
;     psym=-symcat(psym,thick=symthick1), symsize=symsize1, $
      thick=symthick2, color=fsc_color(hacolor,101), line=5
    djs_oplot, gamma[1:nimf-1], alog10(peg[1:nimf-1].c_ir/peg[0].c_ir), $
;     psym=-symcat(psym,thick=symthick1), symsize=symsize1, $
      thick=symthick2, color=fsc_color(ircolor,101), line=3
    djs_oplot, gamma[1:nimf-1], alog10(peg[1:nimf-1].c_24_lo/peg[0].c_24_lo), $
;     psym=-symcat(psym,thick=symthick1), symsize=symsize1, $
      thick=symthick2, color=fsc_color(l24color,101), line=0
    djs_oplot, !x.crange, [0,0], line=1, color=salpcolor, thick=symthick2
    im_legend, ['L(H\alpha)','L(IR)','L(24\mu'+'m)'], /left, /top, box=0, $
      charsize=1.6, color=[hacolor,ircolor,l24color], $
      line=[5,3,0], thick=symthick2, pspacing=1.5, textcolor=keycolor
    
; return fraction    
    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,2], $
      xsty=3, ysty=1, xrange=xrange, yrange=[0,1], $
      xtitle='IMF Slope \Gamma (0.5-120 M_{\odot})', $
      ytitle='Return Fraction'

    polyfill, [gamma[1:nimf-1],reverse(gamma[1:nimf-1])], $
      [peg[1:nimf-1].rfrac_ssp,$
      reverse(peg[1:nimf-1].rfrac+peg[1:nimf-1].rfrac_err/2)],$
      /data, /fill, color=fsc_color('tan',100), noclip=0
;   polyfill, [gamma[1:nimf-1],reverse(gamma[1:nimf-1])], $
;     [peg[1:nimf-1].rfrac-peg[1:nimf-1].rfrac_err/2,$
;     reverse(peg[1:nimf-1].rfrac+peg[1:nimf-1].rfrac_err/2)],$
;     /data, /fill, color=fsc_color('tan',100), noclip=0
;   djs_oplot, gamma[1:nimf-1], peg[1:nimf-1].rfrac_ssp, $
;     thick=symthick2
;   djs_oplot, gamma[1:nimf-1], peg[1:nimf-1].rfrac_ssp, $
;     psym=-symcat(psym,thick=symthick1), symsize=symsize1, $
;     thick=symthick1
    djs_oplot, !x.crange, peg[0].rfrac_ssp*[1,1], $
      line=1, color=salpcolor, thick=symthick2
;   plots, peg[0].slope-1, peg[0].rfrac_ssp, psym=symcat(15,thick=symthick1), $
;     color=fsc_color('firebrick',101), symsize=symsize2

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,suffix,'.pdf'), /sh
       rmfile, psfile
    endif

; --------------------------------------------------
; plot some example IMFs
    psfile = paperpath+'example_imfs'+suffix
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote

    mass1 = range(0.1D,100D,500,/log)
    mass2 = range(0.1D,120D,500,/log)

    djs_plot, [0], [0], /nodata, position=pos, $
      xsty=3, ysty=1, xrange=alog10([0.1,120]), yrange=[-3.2,1.4], $
      xtitle='log (m/M_{\odot})', $
      ytitle='log (\xi_{log m}) (dex^{-1})', $
      charsize=2.2
; xi = Number per unit logarithmic mass interval

; Salpeter
    salpcolor = 'red'
    oplot_imf, mass1, gamma=-1.35, mcut=[0.1,100], line=0, $
      color=fsc_color(salpcolor,101), thick=10
    gamma = [0.5,1.5,2.5]
    line = [3,4,5]
    if keyword_set(keynote) then $
      color = ['cyan','wheat','orange'] else $
        color = ['blue','dark green','orange']
    for ii = 0, n_elements(gamma)-1 do oplot_imf, mass2, gamma=[0.5,-gamma[ii]], $
      mcut=[0.1,0.5,120], line=line[ii], color=fsc_color(color[ii],100), thick=10
    im_legend, ['Salpeter','\Gamma='+string(gamma,format='(F3.1)')], $
      line=[0,line], color=[salpcolor,color], thick=10, /left, $
      /bottom, box=0, charsize=1.6, pspacing=1.8, textcolor=keycolor

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,suffix,'.pdf'), /sh
       rmfile, psfile
    endif

return
end
