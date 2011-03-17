pro atlas_distributions, atlas, encapsulated=encapsulated, postscript=postscript
; jm05sep07uofa

    if (n_elements(atlas) eq 0L) then atlas = read_integrated()

    k92 = read_92kennicutt()
;   k922 = rsex('/home/ioannis/catalogs/92kennicutt/92kennicutt.dat')
    nfgs = read_00jansen()
    
    pspath = atlas_path(/papers)+'atlas/FIG_ATLAS/'
    if keyword_set(postscript) then begin
       postthick = 8.0
       dfpsplot, pspath+'atlas_distributions.eps', encapsulated=encapsulated, xsize=8.5, ysize=9.25
    endif else begin
       im_window, 0, xratio=0.5, /square
       postthick = 2.0
    endelse

    pagemaker, nx=2, ny=2, xspace=1.0, yspace=1.0, width=3.1*[1,1], height=3.1*[1,1], $
      xmargin=[1.0,0.3], ymargin=[0.8,1.1], xpage=8.5, ypage=9.25, position=pos, /normal

; #######
; Panel 1    
; #######
    
    indx = where(atlas.rc3_b_lum gt -900.0)
    x = atlas[indx].rc3_b_lum
    xabs = atlas[indx].rc3_m_b
    stats = im_stats(x,/verbose)
    
    xtitle = 'log L(B) [L(B)'+sunsymbol()+']'
    xrange = [6.9,11.6]
    binsize = 0.25
;   xtitle = 'M_{B} [mag]'
;   xrange = [-12,-23]
;   binsize = 0.75

    im_plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    yrange = minmax(ybin)*[1.0,1.1]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=1.7, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=11, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,0]
    axis, /xaxis, xrange=interpol(xabs,x,!x.crange), xthick=postthick, xstyle=1, $
      charsize=1.7, charthick=postthick
    im_xyouts_title, xtitle='M_{B} [mag]', charsize=1.7, charthick=postthick

; overplot K92 and NFGS

    g = where(k92.rc3_b_lum gt -900)
    im_plothist, k92[g].rc3_b_lum, bin=binsize, /overplot, line=1, thick=postthick, $
      /halfbin, /fill, /fline, forientation=135, fcolor=djs_icolor('grey'), $
      color=djs_icolor('grey'), fspacing=0.05
;     /halfbin, /fill, fcolor=djs_icolor('black'), color=djs_icolor('black')

    g = where(nfgs.rc3_b_lum gt -900)
    im_plothist, nfgs[g].rc3_b_lum, bin=binsize, /overplot, line=2, thick=postthick, $
      /halfbin, /fill, /fline, forientation=45, fcolor=djs_icolor('grey'), $
      color=djs_icolor('grey'), fspacing=0.15

; now overplot my data    
    
    im_plothist, x, bin=binsize, /overplot, thick=postthick, /halfbin

    legend, '(a)', /left, /top, charthick=postthick, charsize=1.5, box=0, $
      clear=keyword_set(postscript), spacing=0

; #######
; Panel 2
; ####### 
    
    indx = where((atlas.rc3_b gt -900) and (atlas.rc3_v gt -900))
    x = atlas[indx].rc3_b-atlas[indx].rc3_v
    stats = im_stats(x,/verbose)

    xtitle = 'B - V [mag]'
    xrange = [0.0,1.1]
    binsize = 0.1
    
    im_plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    yrange = minmax(ybin)*[1.0,1.1]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=1.7, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=3, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,1], /noerase, $
      yminor=3

; overplot K92 and NFGS

    g = where((k92.rc3_b gt -900) and (k92.rc3_v gt -900))
    im_plothist, k92[g].rc3_b-k92[g].rc3_v, bin=binsize, /overplot, line=1, thick=postthick, $
      /halfbin, /fill, /fline, forientation=135, fcolor=djs_icolor('grey'), $
      color=djs_icolor('grey'), fspacing=0.05
;     /halfbin, /fill, fcolor=djs_icolor('black'), color=djs_icolor('black')

    g = where((nfgs.b gt -900) and (nfgs.v gt -900))
    im_plothist, nfgs[g].b-nfgs[g].v, bin=binsize, /overplot, line=2, thick=postthick, $
      /halfbin, /fill, /fline, forientation=45, fcolor=djs_icolor('grey'), $
      color=djs_icolor('grey'), fspacing=0.15

;   p = rsex('/home/ioannis/catalogs/00jansen/nfgs_synth_photometry.dat')
;   im_plothist, p.bvspec, bin=binsize, /overplot, line=2, thick=postthick, $
;     /halfbin, /fill, /fline, forientation=45, fcolor=djs_icolor('green'), $
;     color=djs_icolor('green'), fspacing=0.15
    
; now overplot my data    
    
    im_plothist, x, bin=binsize, /overplot, thick=postthick, /halfbin

    legend, '(b)', /left, /top, charthick=postthick, charsize=1.5, box=0, clear=postscript, spacing=0

; #######
; Panel 3
; #######

    indx = where(atlas.l_fir_l_b gt -900)
    x = alog10(atlas[indx].l_fir_l_b)

    stats = im_stats(x,/verbose)

    xtitle = 'log L(FIR)/L(B)'
    xrange = [-2.2,2.8]
    binsize = 0.25
    
    im_plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    yrange = minmax(ybin)*[1.0,1.1]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=1.7, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,2], /noerase

; overplot K92 and NFGS

    g = where((k92.l_fir_l_b gt -900))
    im_plothist, alog10(k92[g].l_fir_l_b), bin=binsize, /overplot, line=1, thick=postthick, $
      /halfbin, /fill, /fline, forientation=135, fcolor=djs_icolor('grey'), $
      color=djs_icolor('grey'), fspacing=0.05
;     /halfbin, /fill, fcolor=djs_icolor('black'),
;     color=djs_icolor('black')
    print, 'K92 ', median(alog10(k92[g].l_fir_l_b)), 10^djsig(alog10(k92[g].l_fir_l_b),sigrej=4.0)

    g = where(nfgs.l_fir_l_b gt -900)
    im_plothist, alog10(nfgs[g].l_fir_l_b), bin=binsize, /overplot, line=2, thick=postthick, $
      /halfbin, /fill, /fline, forientation=45, fcolor=djs_icolor('grey'), $
      color=djs_icolor('grey'), fspacing=0.15
    print, 'NFGS ', median(alog10(nfgs[g].l_fir_l_b)), 10^djsig(alog10(nfgs[g].l_fir_l_b),sigrej=4.0)

; now overplot my data    
    
    im_plothist, x, bin=binsize, /overplot, thick=postthick, /halfbin

    legend, '(c)', /left, /top, charthick=postthick, charsize=1.5, box=0, clear=postscript, spacing=0

; #######
; Panel 4
; #######

    indx = where(atlas.h_alpha_ew[1] gt 0.0)
    x = alog10(atlas[indx].h_alpha_ew[0])

    stats = im_stats(10^x,/verbose)
    help, x, where(10^x gt 50.0)

    xtitle = 'log EW(H\alpha) ['+angstrom()+']'
    xrange = [-0.5,3.1]
    binsize = 0.15
    
    im_plothist, x, bin=binsize, xbin, ybin, /noplot, /halfbin
    yrange = minmax(ybin)*[1.0,1.1]

    djs_plot, [0], [0], /nodata, xthick=postthick, ythick=postthick, charsize=1.7, $
      charthick=postthick, thick=postthick, xtitle=textoidl(xtitle), ytitle='Number', $
      xstyle=1, ystyle=1, xrange=xrange, yrange=yrange, position=pos[*,3], /noerase

; overplot K92 and NFGS

    g = where((k92.h_alpha_ew[1] gt 0.0))
    print, 'K92 ', 10^median(alog10(k92[g].h_alpha_ew[0]))
    help, g, where(k92[g].h_alpha_ew[0] gt 50.0)
    im_plothist, alog10(k92[g].h_alpha_ew[0]), bin=binsize, /overplot, line=1, thick=postthick, $
      /halfbin, /fill, /fline, forientation=135, fcolor=djs_icolor('grey'), $
      color=djs_icolor('grey'), fspacing=0.05
;     /halfbin, /fill, fcolor=djs_icolor('black'), color=djs_icolor('black')

    g = where(nfgs.h_alpha_ew[1] gt 0.0)
    print, 'NFGS ', 10^median(alog10(nfgs[g].h_alpha_ew[0]))
    help, g, where(nfgs[g].h_alpha_ew[0] gt 50.0)
    im_plothist, alog10(nfgs[g].h_alpha_ew[0]), bin=binsize, /overplot, line=2, thick=postthick, $
      /halfbin, /fill, /fline, forientation=45, fcolor=djs_icolor('grey'), $
      color=djs_icolor('grey'), fspacing=0.15

; now overplot my data    
    
    im_plothist, x, bin=binsize, /overplot, thick=postthick, /halfbin

    legend, '(d)', /left, /top, charthick=postthick, charsize=1.5, box=0, clear=postscript, spacing=0

    if keyword_set(postscript) then dfpsclose

stop
    
return
end
    
