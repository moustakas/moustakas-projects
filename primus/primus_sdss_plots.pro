pro primus_sdss_plots, atlas, sdss, sdssancillary, atlasclass, sdssclass, ps=ps
; jm07sep05nyu - 

    if (n_elements(atlas) eq 0L) then begin
       int = read_integrated()
       nfgs = read_nfgs()
       atlas = struct_append(int,nfgs)
    endif
    if (n_elements(atlasclass) eq 0L) then atlasclass = iclassification(atlas,snrcut=5.0,/silent,/kauffmann)
    if (n_elements(sdss) eq 0L) then sdss = mrdfits(sdss_path()+'sdss_specdata_main_dr4.fits.gz',1,/silent)
    if (n_elements(sdssancillary) eq 0L) then sdssancillary = mrdfits(sdss_path()+'sdss_ancillary_main_dr4.fits.gz',1,/silent)
    if (n_elements(sdssclass) eq 0L) then sdssclass = iclassification(sdss,snrcut=5.0,/silent,/kauffmann)

    snrcut = 5.0
    
    if keyword_set(ps) then begin
       postthick1 = 4.0
       postthick2 = 8.0
    endif else begin
       im_window, 0, xratio=0.4, /square
       postthick1 = 2.0
       postthick2 = 2.0
    endelse

; ---------------------------------------------------------------------------    
; Dn(4000) vs EW([O II])
; ---------------------------------------------------------------------------    

    if keyword_set(ps) then dfpsplot, 'dn4000_vs_ewoii.ps', /square, /color

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

; SDSS
    
    indx_sdss = where((sdss.d4000_narrow[0]/sdss.d4000_narrow[1] gt snrcut) and $
      (sdss.oii_3727_ew[0]/sdss.oii_3727_ew[1] gt snrcut))
    x_sdss = sdss[indx_sdss].d4000_narrow_model[0]
;   y_sdss = sdss[indx_sdss].oii_3727_ew[0]
    y_sdss = alog10(sdss[indx_sdss].oii_3727_ew[0])

    sf_sdss = where(strtrim(sdssclass[indx_sdss].bpt_class,2) eq 'HII')
    agn_sdss = where(strtrim(sdssclass[indx_sdss].bpt_class,2) eq 'AGN')
    
; Atlas
    
    indx_atlas = where((atlas.d4000_narrow[0]/atlas.d4000_narrow[1] gt snrcut) and $
      (atlas.oii_3727_ew[0]/atlas.oii_3727_ew[1] gt snrcut))
    x_atlas = atlas[indx_atlas].d4000_narrow_model[0]
;   y_atlas = atlas[indx_atlas].oii_3727_ew[0]
    y_atlas = alog10(atlas[indx_atlas].oii_3727_ew[0])
    
    sf_atlas = where(strtrim(atlasclass[indx_atlas].bpt_class,2) eq 'HII')
    agn_atlas = where(strtrim(atlasclass[indx_atlas].bpt_class,2) eq 'AGN')
    
    xrange = [0.7,2.3]
    yrange = [0.0,2.5]
;   yrange = [0.0,100]

    xtitle = 'D_{n} (4000)'
;   ytitle = 'EW([O II] \lambda3727) (\AA)'
    ytitle = 'log EW([O II] \lambda3727) (\AA)'

    im_symbols, 106, psize=0.2, thick=postthick1, fill=1, color='dark green'
    hogg_scatterplot, x_sdss[sf_sdss], y_sdss[sf_sdss], /outliers, label=0, $
      outpsym=8, outsymsize=1.0, outcolor='dark green', xrange=xrange, yrange=yrange, $
      xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), charsize=1.8, charthick=postthick1, $
      xthick=postthick1, ythick=postthick1, position=pos
;   im_symbols, 106, psize=0.6, thick=postthick1, fill=1, color='blue'
;   djs_oplot, x_atlas[sf_atlas], y_atlas[sf_atlas], ps=8
    legend, 'Star-Forming', /right, /top, box=0, charsize=1.8, charthick=postthick1

    if (not keyword_set(ps)) then cc = get_kbrd(1)

    im_symbols, 106, psize=0.2, thick=postthick1, fill=1, color='dark green'
    hogg_scatterplot, x_sdss[agn_sdss], y_sdss[agn_sdss], /outliers, label=0, $
      outpsym=8, outsymsize=1.0, outcolor='dark green', xrange=xrange, yrange=yrange, $
      xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), charsize=1.8, charthick=postthick1, $
      xthick=postthick1, ythick=postthick1, position=pos
;   im_symbols, 106, psize=0.6, thick=postthick1, fill=1, color='blue'
;   djs_oplot, x_atlas[agn_atlas], y_atlas[agn_atlas], ps=8
    legend, 'AGN', /right, /top, box=0, charsize=1.8, charthick=postthick1

    if keyword_set(ps) then dfpsclose

; ---------------------------------------------------------------------------    
; distribution of Dn(4000) for galaxies with poorly measured EW([O II])
; ---------------------------------------------------------------------------    

    if keyword_set(ps) then dfpsplot, 'dn4000_hist.ps', /square, /color

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

; SDSS

    junk = where((sdss.d4000_narrow[0]/sdss.d4000_narrow[1] gt snrcut) and $
      (sdssancillary.z gt 0.033),normfactor)
    
    indx_sdss_lo = where((sdss.d4000_narrow[0]/sdss.d4000_narrow[1] gt snrcut) and $
      (sdss.oii_3727_ew[0]/sdss.oii_3727_ew[1] lt snrcut) and (sdssancillary.z gt 0.033))
    x_sdss_lo = sdss[indx_sdss_lo].d4000_narrow_model[0]
    indx_sdss_hi = where((sdss.d4000_narrow[0]/sdss.d4000_narrow[1] gt snrcut) and $
      (sdss.oii_3727_ew[0]/sdss.oii_3727_ew[1] gt snrcut) and (sdssancillary.z gt 0.033))
    x_sdss_hi = sdss[indx_sdss_hi].d4000_narrow_model[0]

    binsize = 0.03
    im_plothist, x_sdss_lo, bin=binsize, xbin_lo, ybin_lo, /noplot, norm=normfactor
    im_plothist, x_sdss_hi, bin=binsize, xbin_hi, ybin_hi, /noplot, norm=normfactor

    xrange = [0.7,2.3]
;   xrange = [min(xbin_lo)<min(xbin_hi),max(ybin_lo)>max(ybin_hi)]
    yrange = [0.0,max(ybin_lo)>max(ybin_hi)]*1.05
    
    xtitle = 'D_{n} (4000)'
    ytitle = 'Fraction'

    djs_plot, [0], [0], /nodata, yrange=yrange, xrange=xrange, $
      position=pos[*,0], xsty=1, ysty=1, xtitle=xtitle, ytitle=ytitle, $
      xthick=postthick1, ythick=postthick1, charsize=1.8, charthick=postthick1
    im_plothist, x_sdss_lo, bin=binsize, /overplot, line=0, thick=postthick1, norm=normfactor
    im_plothist, x_sdss_hi, bin=binsize, /overplot, line=2, thick=postthick2, norm=normfactor

    legend, ['EW([O II]) S/N < 5','EW([O II]) S/N > 5'], /left, /top, box=0, charsize=1.8, $
      charthick=postthick1, line=[0,2], thick=postthick2
    
    if keyword_set(ps) then dfpsclose

stop    
    
return
end

    
