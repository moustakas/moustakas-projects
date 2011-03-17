;+
; NAME:
;       COMPARE_SYNTHMAGS_RC3
;
; PURPOSE:
;       Assess the spectrophotometric accuracy of the atlas by
;       comparing synthetic (spectroscopic) apparent magnitudes and
;       colors with the RC3.
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES USED:
;
; TODO:
;       Flag truly integrated galaxies from partially integrated
;       ones. 
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 Sep 10, U of A, re-written and updated from
;         an earlier version of the code
;       jm04feb05uofa - major re-write
;       jm04mar11uofa - renamed PHOTOSYNTH; include 2MASS photometry  
;       jm05aug01uofa - renamed COMPARE_SYNTHMAGS_RC3; cleaned up 
;-

pro compare_synthmags_rc3, atlas, result, paper=paper, postscript=postscript

    pspath = atlas_path(/papers)+'atlas/FIG_ATLAS/'

    if (n_elements(atlas) eq 0L) then atlas = read_integrated()
    if keyword_set(paper) then postscript = 1

    lcharsize = 1.2

    if (n_elements(xfrac) eq 0L) then xfrac = 2.0
    if (n_elements(yfrac) eq 0L) then yfrac = 2.0

; ---------------------------------------------------------------------------
; Figure: B, V [RC3] vs B, V [Synthetic]
; ---------------------------------------------------------------------------

; initialize some plotting variables

    xmargin = [1.2,1.2]
    ymargin = [0.3,1.2]

    width = [3.0,3.0]
    height = [3.0,2.5]

    xspace = 0.0
    yspace = 1.0
    
    xpage = total(xmargin)+total(xspace)+total(width)
    ypage = total(ymargin)+total(yspace)+total(height)

    pagemaker, nx=2, ny=2, position=pos, /normal, xspace=xspace, width=width, $
      height=height, xmargin=xmargin, ymargin=ymargin, yspace=yspace, $
      xpage=xpage, ypage=ypage

    magrange = [6.5,18.5]
    d25range = [-1.0,1.4]
    
    if keyword_set(postscript) then begin
       dfpsplot, pspath+'compare_synthmags_rc3_bvmags.eps', xsize=xpage, ysize=ypage, /encapsulated
       postthick = 5.0
    endif else begin
       im_window, 0, xratio=0.5, /square
       postthick = 2.0
    endelse

; ---------------------------------------------------------------------------    
; B-band    
; ---------------------------------------------------------------------------    

    photo = where((atlas.rc3_b gt -900.0) and (atlas.synth_b_obs gt -900.0) and $
      (atlas.drift_photflag eq 'Y') and (atlas.d25_maj gt -900.0),nphoto)
    nophoto = where((atlas.rc3_b gt -900.0) and (atlas.synth_b_obs gt -900.0) and $
      (atlas.drift_photflag eq 'N') and (atlas.d25_maj gt -900.0),nnophoto)

    xphoto = atlas[photo].rc3_b 
    yphoto = atlas[photo].synth_b_obs
    d25photo = alog10(atlas[photo].d25_maj)
    
    xnophoto = atlas[nophoto].rc3_b 
    ynophoto = atlas[nophoto].synth_b_obs
    d25nophoto = alog10(atlas[nophoto].d25_maj)

    resphoto = yphoto-xphoto
    resnophoto = ynophoto-xnophoto
    
    im_symbols, 108, psize=0.8, /fill
    djs_plot, xphoto, yphoto, ps=8, xsty=3, ysty=3, charsize=1.7, $
      charthick=postthick, xthick=postthick, ythick=postthick, xrange=magrange, $
      yrange=magrange, position=pos[*,0], $
      ytitle=textoidl('m_{B} [Synthesized]'), xtitle=textoidl('m_{B} [Literature]')
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    im_symbols, 106, psize=psize, fill=0, thick=postthick
    djs_oplot, xnophoto, ynophoto, ps=8
    
;   im_legend, ['Photometric','Non-Photometric'], psym=[108,106], $
;     fill=[1,0], /left, /top, box=0, charsize=1.2, charthick=postthick, $
;     spacing=1.5

; residuals

    resrange = (max(abs(resphoto))>max(abs(resnophoto)))*[-1.1,1.1]
    
    im_symbols, 108, psize=0.8, /fill
    djs_plot, d25photo, resphoto, xsty=3, ysty=3, ps=8, xrange=d25range, $
      position=pos[*,2], /noerase, yrange=resrange, $
      charsize=1.7, charthick=postthick, xthick=postthick, ythick=postthick, $
      xtitle=textoidl('log D_{25} [arcmin]'), yminor=3, /nohat, ytitle='Residuals [mag]'
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    im_symbols, 106, psize=psize, fill=0, thick=postthick
    djs_oplot, d25nophoto, resnophoto, ps=8

;   djs_plot, xphoto, resphoto, xsty=3, ysty=3, ps=8, xrange=magrange, $
;     position=pos[*,2], /noerase, yrange=resrange, $
;     charsize=1.7, charthick=postthick, xthick=postthick, ythick=postthick, $
;     xtitle=textoidl('m_{B} [Literature]'), yminor=3, /nohat, ytitle='Residuals [mag]'
;   djs_oplot, xnophoto, resnophoto, ps=8

    splog, 'B Residual Statistics [Photometric, Non-Photometric, Combined]'
    stats = im_stats(resphoto,/verbose,sigrej=3.0)
    junk = im_stats(resnophoto,/verbose,sigrej=3.0,/no_head)
    junk = im_stats([resphoto,resnophoto],/verbose,sigrej=3.0,/no_head)

    running = im_medxbin([d25photo,d25nophoto],[resphoto,resnophoto],0.5,$
      minpts=10,minx=-0.6,/verbose)

    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   legend, textoidl(xstr), /left, /bottom, box=0, charsize=lcharsize, $
;     charthick=postthick, clear=keyword_set(postscript)

; ---------------------------------------------------------------------------    
; V-band    
; ---------------------------------------------------------------------------    

    photo = where((atlas.rc3_v gt -900.0) and (atlas.synth_v_obs gt -900.0) and $
      (atlas.drift_photflag eq 'Y') and (atlas.d25_maj gt -900),nphoto)
    nophoto = where((atlas.rc3_v gt -900.0) and (atlas.synth_v_obs gt -900.0) and $
      (atlas.drift_photflag eq 'N') and (atlas.d25_maj gt -900),nnophoto)

    xphoto = atlas[photo].rc3_v 
    yphoto = atlas[photo].synth_v_obs
    d25photo = alog10(atlas[photo].d25_maj)
    
    xnophoto = atlas[nophoto].rc3_v 
    ynophoto = atlas[nophoto].synth_v_obs
    d25nophoto = alog10(atlas[nophoto].d25_maj)

    resphoto = yphoto-xphoto
    resnophoto = ynophoto-xnophoto
    
    im_symbols, 108, psize=0.8, /fill
    djs_plot, xphoto, yphoto, ps=8, xsty=3, ysty=11, charsize=1.7, $
      charthick=postthick, xthick=postthick, ythick=postthick, xrange=magrange, $
      yrange=magrange, position=pos[*,1], xtitle=textoidl('m_{V} [Literature]'), $ ;xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), ytitle='', /noerase
    axis, /yaxis, yrange=magrange, ysty=3, ytitle=textoidl('m_{V} [Synthesized]'), $
      yminor=3, ythick=postthick, charsize=1.7, charthick=postthick
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    im_symbols, 106, psize=psize, fill=0, thick=postthick
    djs_oplot, xnophoto, ynophoto, ps=8
    
;   im_legend, ['Photometric','Non-Photometric'], psym=[108,106], $
;     fill=[1,0], /left, /top, box=0, charsize=1.2, charthick=postthick, $
;     spacing=1.5

; residuals

    resrange = (max(abs(resphoto))>max(abs(resnophoto)))*[-1.1,1.1]
    
    im_symbols, 108, psize=0.8, /fill
    djs_plot, d25photo, resphoto, xsty=3, ysty=11, ps=8, $
      position=pos[*,3], /noerase, yrange=resrange, xrange=d25range, $
      charsize=1.7, charthick=postthick, xthick=postthick, ythick=postthick, $
      xtitle=textoidl('log D_{25} [arcmin]'), yminor=3, ytickname=replicate(' ',10)
    axis, /yaxis, yrange=resrange, ysty=3, ytitle='Residuals [mag]', $
      yminor=3, ythick=postthick, charsize=1.7, charthick=postthick
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    im_symbols, 106, psize=psize, fill=0, thick=postthick
    djs_oplot, d25nophoto, resnophoto, ps=8

;   djs_plot, xphoto, resphoto, xsty=3, ysty=11, ps=8, $
;     position=pos[*,3], /noerase, yrange=resrange, xrange=magrange, $
;     charsize=1.7, charthick=postthick, xthick=postthick, ythick=postthick, $
;     xtitle=textoidl('m_{V} [Literature]'), yminor=3, ytickname=replicate(' ',10)
;   axis, /yaxis, yrange=resrange, ysty=3, ytitle='Residuals [mag]', $
;     yminor=3, ythick=postthick, charsize=1.7, charthick=postthick
;   djs_oplot, xnophoto, resnophoto, ps=8

    splog, 'V Residual Statistics [Photometric, Non-Photometric, Combined]'
    stats = im_stats(resphoto,/verbose,sigrej=3.0)
    junk = im_stats(resnophoto,/verbose,sigrej=3.0,/no_head)
    junk = im_stats([resphoto,resnophoto],/verbose,sigrej=3.0,/no_head)

    running = im_medxbin([d25photo,d25nophoto],[resphoto,resnophoto],0.5,$
      minpts=10,minx=-0.6,/verbose)

    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   legend, textoidl(xstr), /left, /bottom, box=0, charsize=lcharsize, $
;     charthick=postthick, clear=keyword_set(postscript)

    if keyword_set(postscript) then dfpsclose else cc = get_kbrd(1)

; ---------------------------------------------------------------------------
; Figure: (B-V) [RC3] vs (B-V) [Synthetic]
; ---------------------------------------------------------------------------

; initialize some plotting variables

    xmargin = [1.2,0.3]
    ymargin = [0.3,1.2]

    width = 6.0
    height = [3.0,2.5]

    xpage = total(xmargin)+total(width)
    ypage = total(ymargin)+total(height)

    pagemaker, nx=1, ny=2, position=pos, /normal, xspace=0.0, width=width, $
      height=height, xmargin=xmargin, ymargin=ymargin, yspace=0.0, $
      xpage=xpage, ypage=ypage

    magrange = [6.5,18.5]
    colorrange = [-0.3,1.2]
    
    if keyword_set(postscript) then begin
       dfpsplot, pspath+'compare_synthmags_rc3_bvcolor.eps', xsize=xpage, ysize=ypage, /encapsulated
       postthick = 8.0
    endif else begin
       im_window, 0, xratio=0.5, /square
       postthick = 2.0
    endelse

    photo = where((atlas.rc3_b gt -900.0) and (atlas.synth_b_obs gt -900.0) and $
      (atlas.rc3_v gt -900.0) and (atlas.synth_v_obs gt -900.0) and $
      (atlas.drift_photflag eq 'Y'),nphoto)
    nophoto = where((atlas.rc3_b gt -900.0) and (atlas.synth_b_obs gt -900.0) and $
      (atlas.rc3_v gt -900.0) and (atlas.synth_v_obs gt -900.0) and $
      (atlas.drift_photflag eq 'N'),nnophoto)

    xphoto = atlas[photo].rc3_b-atlas[photo].rc3_v
    xerrphoto = sqrt(atlas[photo].rc3_b_err^2+atlas[photo].rc3_v_err^2)
    yphoto = atlas[photo].synth_b_obs-atlas[photo].synth_v_obs
    yerrphoto = sqrt(atlas[photo].synth_b_obs_err^2+atlas[photo].synth_v_obs_err^2)
    
    xnophoto = atlas[nophoto].rc3_b-atlas[nophoto].rc3_v
    xerrnophoto = sqrt(atlas[nophoto].rc3_b_err^2+atlas[nophoto].rc3_v_err^2)
    ynophoto = atlas[nophoto].synth_b_obs-atlas[nophoto].synth_v_obs
    yerrnophoto = sqrt(atlas[nophoto].synth_b_obs_err^2+atlas[nophoto].synth_v_obs_err^2)

    resphoto = yphoto-xphoto
    resnophoto = ynophoto-xnophoto

    xerravg = djs_median([xerrphoto,xerrnophoto])
    yerravg = djs_median([yerrphoto,yerrnophoto])

    im_symbols, 108, psize=0.8, /fill
    djs_plot, xphoto, yphoto, ps=8, xsty=3, ysty=3, charsize=1.7, $
      charthick=postthick, xthick=postthick, ythick=postthick, xrange=colorrange, $
      yrange=colorrange, position=pos[*,0], xtickname=replicate(' ',10), $
      ytitle=textoidl('(B-V) [Synthesized]'), /nohat
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    im_symbols, 106, psize=psize, fill=0, thick=postthick
    djs_oplot, xnophoto, ynophoto, ps=8

    xoff = (xfrac*xerravg) > (0.075*(!x.crange[1]-!x.crange[0]))
    yoff = (yfrac*yerravg) > (0.075*(!y.crange[1]-!y.crange[0]))
    xc = !x.crange[1]-xoff & yc = !y.crange[0]+yoff
;   oploterror, xc, yc, xerravg, yerravg, ps=3, errthick=errthick;, /nohat

    xyouts, 1.1, 0.05, 'MRK 0475', align=0.5, /data, charsize=1.0, $
      charthick=postthick

; residuals

    resrange = (max(abs(resphoto))>max(abs(resnophoto)))*[-1.1,1.1]
    
    im_symbols, 108, psize=0.8, /fill
    djs_plot, xphoto, resphoto, xsty=3, ysty=3, ps=8, $
      position=pos[*,1], /noerase, yrange=resrange, xrange=colorrange, $
      charsize=1.7, charthick=postthick, xthick=postthick, ythick=postthick, $
      xtitle=textoidl('(B-V) [Literature]'), yminor=3, ytitle='Residuals [mag]'
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    im_symbols, 106, psize=psize, fill=0, thick=postthick
    djs_oplot, xnophoto, resnophoto, ps=8

    xyouts, 1.1, -0.925, 'MRK 0475', align=0.5, /data, charsize=1.0, $
      charthick=postthick

    splog, '(B-V) Residual Statistics [Photometric, Non-Photometric, Combined]'
    stats = im_stats(resphoto,sigrej=3.0,/verbose)
    junk = im_stats(resnophoto,sigrej=3.0,/verbose,/no_head)
    junk = im_stats([resphoto,resnophoto],sigrej=3.0,/verbose,/no_head)

    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
;   legend, textoidl(xstr), /right, /top, box=0, charsize=lcharsize, $
;     charthick=postthick, clear=keyword_set(postscript)
    
    if keyword_set(postscript) then dfpsclose else cc = get_kbrd(1)

stop
    
return
end
