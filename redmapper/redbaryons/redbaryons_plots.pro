function rykoff_mass, lambda
; relation between richness and M200 taken from equation B4 in
; Rykoff+12 
    return, alog10(1D14*exp(1.72+1.08*alog(lambda/60.0)))
end

pro redbaryons_plots
; jm13aug27siena - build some plots for the redmapper/baryons paper

    common com_redmapper, centrals, satellites

    paperpath = redmapper_path(/redbaryons,/paper)
    massprofpath = bcgmstar_path(/massprofiles)

; read the catalogs
    if n_elements(centrals) eq 0L then read_redmapper, $
      centrals=centrals, satellites=satellites
    zbins = redbaryons_zbins(nzbins)

; --------------------------------------------------    
; plot total stellar mass vs virial mass using the CLASH and redmapper
; samples 
    sample = read_bcgmstar_sample(/getmbcg)
    ncl = n_elements(sample)

; make the plot    
    xrange = alog10([3D13,3D15])
    yrange = [11.0,13.2]
    
    psfile = paperpath+'mbcg_vs_m500.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, xmargin=[1.3,0.4], $
      width=6.8

    djs_plot, [0], [0], /nodata, position=pos, xsty=5, ysty=5, $
      xrange=xrange, yrange=yrange
      
; overplot the power-law fit from Kravtsov (Table 2)    
    scat = 0.17
    xx = range(xrange[0],xrange[1],50)
    yy = poly(xx-14.5,[12.24,0.33])
    polyfill, [xx,reverse(xx)], [yy+scat,reverse(yy-scat)], $
      /fill, color=cgcolor('light grey'), noclip=0
;   djs_oplot, xx, yy, thick=4

; overlay Kravtsov
    mbcg_krav = [3.12,4.14,3.06,1.47,0.79,1.26,1.09,0.91,1.38]*1D12
    mbcg_krav_err = [0.36,0.3,0.3,0.13,0.05,0.11,0.06,0.05,0.14]*1D12
    m500_krav = [15.6,10.3,7,5.34,2.35,1.86,1.34,0.46,0.47]*1D14

    oploterror, alog10(m500_krav), alog10(mbcg_krav), mbcg_krav_err/mbcg_krav/alog(10), $
      psym=symcat(15), color=cgcolor('tomato'), errcolor=cgcolor('tomato'), $
      symsize=1.5

; overplot Gonzalez+13
    mbcg_gonz = [0.68,0.82,0.35,0.74,0.56,0.35,0.46,0.88,0.64,0.7,0.65,0.35]*1D13
    mbcg_gonz_err = [0.04,0.06,0.03,0.06,0.05,0.02,0.02,0.06,0.05,0.06,0.04,0.02]*1D13
    m500_gonz = [2.26,5.15,0.95,3.46,3.59,0.99,0.95,3.23,2.26,2.41,2.37,1.45]*1D14
    m500_gonz_err = [0.19,0.42,0.1,0.32,0.28,0.11,0.1,0.19,0.23,0.18,0.24,0.21]*1D14
    
;    oploterror, alog10(m500_gonz), alog10(mbcg_gonz), mbcg_gonz_err/mbcg_gonz/alog(10), $
;      psym=symcat(14), color=cgcolor('forest green'), errcolor=cgcolor('forest green'), $
;      symsize=1.8

; overplot CLASH    
    oploterror, sample.m500, sample.mbcg, sample.m500_err, sample.mbcg_err, psym=symcat(16), symsize=1.6

; make a legend
    im_legend, ['CLASH','Kravtsov+14'], /right, /bottom, box=0, $
      psym=[16,15], color=['black','tomato'], spacing=2.5

; redraw the axes
    djs_plot, [0], [0], /nodata, /noerase, position=pos, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, $
      xtitle='log_{10} (M_{500} / M_{\odot})', $
      ytitle='log_{10} (M_{*,BCG} / M_{\odot})'
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop
    
; ---------------------------------------------------------------------------
; stellar mass vs redshift
    psfile = paperpath+'mstar_vs_z'+suffix
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, $
      ymargin=[0.4,1.1], xmargin=[1.3,0.4], width=6.8

    im_hogg_scatterplot, satellites.z, satellites.mstar_avg, position=pos, $
      xsty=1, ysty=1, /internal, outlier=0, xrange=[0.05,0.6], yrange=[9.0,12.5], $
      xtitle='Redshift', ytitle=redbaryons_masstitle(), $
      outcolor=im_color('dodger blue'), levels=[0.25,0.5,0.75,0.95,0.99], $
      xnpix=101, ynpix=101, /nogrey, contour_color=im_color('dodger blue'), $
      cthick=8, ytickinterval=1
    im_hogg_scatterplot, centrals.z, centrals.mstar_avg, /noerase, /nogrey, $
      /overplot, /internal, outlier=0, outcolor=im_color('firebrick'), $
      contour_color=im_color('firebrick'), levels=[0.25,0.5,0.75,0.95,0.99], $
      cthick=8
    im_plotconfig, psfile=psfile, /psclose, pdf=pdf
    
; ---------------------------------------------------------------------------
; CSMF for satellites & centrals in six redshift bins
;   xrange = [9.8,12.5]
;   yrange = [1,5.99]
;   yrange = [-8.0,-2.5]

    xrange = [10.8,12.2]
    yrange = [-2.5,1]
;   yrange = [1,4.5]

    mingal = 3
    cen_psym = 16 & cen_color = 'firebrick'
    sat_psym = 15 & sat_color = 'dodger blue'
    
    psfile = paperpath+'csmf_censat'+suffix
    im_plotconfig, 17, pos, psfile=psfile, charsize=1.6, $
      xmargin=[1.1,0.3], width=3.2*[1,1,1], height=2.6*[1,1]

    csmf_sat = mrdfits(datapath+'csmf_sat.fits.gz',1)
    csmf_cen = mrdfits(datapath+'csmf_cen.fits.gz',1)
    for iz = 0, nzbins-1 do begin
       if (iz le 2) then begin
          xtickname = replicate(' ',10)
          xtitle = ''
       endif else begin
          delvarx, xtickname
          xtitle = redbaryons_masstitle()
       endelse
       if (iz eq 0) or (iz eq 3) then delvarx, ytickname else ytickname = replicate(' ',10)
;      if odd(iz) then delvarx, ytickname else ytickname = replicate(' ',10)
       
       plot, [0], [0], /nodata, noerase=iz gt 0, position=pos[*,iz], xsty=1, ysty=1, $
         yrange=yrange, xrange=xrange, ytickname=ytickname, xtitle=xtitle, $
         ytitle='', xtickinterval=0.5, xminor=5, xtickname=xtickname, yminor=5, $
         ytickinterval=1.0
       im_legend, 'z='+strtrim(string(zbins[iz].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[iz].zup,format='(F12.2)'),2), $
         /right, /top, box=0, margin=-0.1, charsize=1.3

       if (iz gt 0) then begin
          good = where(csmf_cen[0].number ge mingal and csmf_cen[0].limit eq 1)
          djs_oplot, csmf_cen[0].meanmass[good], alog10(csmf_cen[0].phi[good]), $
            line=0, color=im_color('grey')
       endif

;; satellites       
;       good = where(csmf_sat[iz].number ge mingal and csmf_sat[iz].limit eq 1)
;       oploterror, csmf_sat[iz].meanmass[good], alog10(csmf_sat[iz].phi[good]), $
;         csmf_sat[iz].phierr[good]/csmf_sat[iz].phi[good]/alog(10), $
;         psym=symcat(sat_psym,thick=3), symsize=0.75, color=im_color(sat_color), $
;         errcolor=im_color(sat_color)

; centrals
       good = where(csmf_cen[iz].number ge mingal and csmf_cen[iz].limit eq 1)
       oploterror, csmf_cen[iz].meanmass[good], alog10(csmf_cen[iz].phi[good]), $
         csmf_cen[iz].phierr[good]/csmf_cen[iz].phi[good]/alog(10), $
         psym=symcat(cen_psym,thick=3), symsize=0.75, color=im_color(cen_color), $
         errcolor=im_color(cen_color)
;      print, alog10(csmf_cen[iz].phi[good])
       
;      djs_oplot, 11.6*[1,1], !y.crange, line=0
    endfor 

; ytitle    
    xyouts, pos[0,0]-0.06, pos[1,0], textoidl('log (\Phi / Mpc^{-3} dex^{-1})'), $
      align=0.5, /normal, orientation=90
;   xyouts, pos[0,0]-0.06, pos[1,0], textoidl('log (\Phi / Mpc^{-3} dex^{-1})'), $
;     align=0.5, /normal, orientation=90
    
    im_plotconfig, psfile=psfile, /psclose, pdf=pdf

; ---------------------------------------------------------------------------
; central galaxy stellar mass vs cluster richness
    good = where(centrals.z gt 0.1 and centrals.z lt 0.35,ngal)
    xx = centrals[good].mstar_avg
    yy = centrals[good].lambda_chisq
    mm = im_medxbin(xx,yy,0.1,minpts=100,/verbose)

    splog, 'Code by redshift!!'
    richrange = [20,200]
    
;   xtitle = textoidl('log_{10}(Stellar Mass) (M_{'+sunsymbol()+'})')
    
    psfile = paperpath+'mstarcen_vs_rich'+suffix
    im_plotconfig, 0, pos, psfile=psfile, height=4.5, $
      ymargin=[0.5,1.1], xmargin=[1.3,1.3], width=5.9
;   djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;     xrange=[9,12], yrange=[20,500], /ylog

    im_hogg_scatterplot, xx, yy, position=pos, xsty=1, ysty=9, $
      /internal, /outlier, xrange=[10.8,12.2], yrange=richrange, $
      xtickinterval=0.5, ytitle=textoidl('Cluster Richness \lambda'), $
      xtitle=redbaryons_masstitle(), $
      ytick_get=yticks, outcolor=im_color('dodger blue')
    axis, /yaxis, ysty=1, yrange=rykoff_mass(richrange), $
;     ytickv=rykoff_mass([richrange,yticks]), yticks=n_elements(yticks)+1, $
      ytitle=textoidl('log_{10} (M_{200m} / h_{70}^{-1} M_{'+sunsymbol()+'})')
;     ytickinterval=0.25
;   djs_oplot, !x.crange, 20*[1,1], line=0
;   djs_oplot, mm.medx, mm.medy, line=0
;   djs_oplot, mm.medx, mm.quant25, line=5
;   djs_oplot, mm.medx, mm.quant75, line=5
    im_plotconfig, psfile=psfile, /psclose, pdf=pdf
    
stop    
    
return
end
    
