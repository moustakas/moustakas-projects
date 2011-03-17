pro deep2_zaritsky, kcorr, oiii, ediscs, ediscs_info, ages, ages_info, ps=ps
; jm07feb22nyu 
    
    analysis_path = deep2_path(/analysis)

    if (n_elements(kcorr) eq 0L) then kcorr = mrdfits(analysis_path+'deep2.dr2.kcorr.fits.gz',1,/silent)
    if (n_elements(oiii) eq 0L) then oiii = mrdfits(analysis_path+'deep2.dr2.oiii.fits.gz',1,/silent)
    if (n_elements(ediscs) eq 0L) then ediscs = read_ediscs(ancillary=ediscs_info)
    if (n_elements(ages) eq 0L) then ages_info = read_ages(specdata=ages)

    snrcut = 3.0

; ---------------------------------------------------------------------------    
; DEEP-2
; ---------------------------------------------------------------------------    

    indx = where((oiii.oiii_5007_ew[0] gt 0.0) and (oiii.oiii_4959_ew[0] gt 0.0) and $
      (oiii.oiii_5007_ew[0]/oiii.oiii_5007_ew[1] gt snrcut) and (oiii.dec lt 10.0),nindx)
    sample = struct_addtags(oiii[indx],kcorr[indx])

    sample.oiii_4959 = sample.oiii_4959_ew*rebin(reform(sample.cflux_4959,1,nindx),2,nindx)
    sample.oiii_5007 = sample.oiii_5007_ew*rebin(reform(sample.cflux_5007,1,nindx),2,nindx)

    cut = where(sample.oiii_5007[0] ge 2.6D-17,ncut)

    z = sample[cut].z
    flux = 1D17*sample[cut].oiii_5007[0]

    psname = 'deep2_zhist_alpha.ps'
    if keyword_set(ps) then dfpsplot, analysis_path+psname, /square, /color
    
    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    postthick = 4.0
    
    bin = 0.01
    plothist, z, xbin, ybin, bin=bin, /noplot, /halfbin
    yrange = [0,max(ybin)]*[1.0,1.05]
    xrange = minmax(xbin)*[0.9,1.1]

    xtitle = 'Redshift'
    ytitle = 'Number'

    djs_plot, [0], [0], xsty=3, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle='', ytitle='', xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      position=pos[*,0], title='DEEP-2'
    plothist, z, bin=bin, line=0, thick=postthick, /overplot, $
      /halfbin, color=djs_icolor('red'), fcolor=djs_icolor('red'), /fill
    djs_plot, [0], [0], /noerase, xsty=3, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,0]

    legend, textoidl(['S/N EW([OIII] 5007) > 3',$
      'F_{[O III] 5007} > 2.6x10^{-17} cgs','LCO/Fall Accessible']), $
      /left, /top, box=0, charsize=1.5, charthick=postthick, spacing=2.2, margin=0
    legend, textoidl(['N_{gal} = '+strtrim(ncut,2)]), /right, /top, box=0, charsize=1.8, charthick=postthick

    if keyword_set(ps) then begin
       dfpsclose
       spawn, 'convert '+analysis_path+psname+' '+analysis_path+repstr(psname,'.ps','.png'), /sh
    endif
    
; ---------------------------------------------------------------------------    
; EDISCS
; ---------------------------------------------------------------------------    

    cut = where((ediscs.oiii_5007_ew[0] gt 0.0) and (ediscs.oiii_4959_ew[0] gt 0.0) and $
      (ediscs.oiii_5007_ew[0]/ediscs.oiii_5007_ew[1] gt snrcut),ncut)
    sample = ediscs[cut]

    z = sample.z_obj

    psname = 'ediscs_zhist_alpha.ps'
    if keyword_set(ps) then dfpsplot, analysis_path+psname, /square, /color
    
    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    postthick = 4.0
    
    bin = 0.01
    plothist, z, xbin, ybin, bin=bin, /noplot, /halfbin
    yrange = [0,max(ybin)]*[1.0,1.05]
    xrange = minmax(xbin)*[0.9,1.1]

    xtitle = 'Redshift'
    ytitle = 'Number'

    djs_plot, [0], [0], xsty=3, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle='', ytitle='', xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      position=pos[*,0], title='EDisCS'
    plothist, z, bin=bin, line=0, thick=postthick, /overplot, $
      /halfbin, color=djs_icolor('red'), fcolor=djs_icolor('red'), /fill
    djs_plot, [0], [0], /noerase, xsty=3, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle=xtitle, ytitle=ytitle, charsize=2.0, charthick=postthick, $
      xthick=postthick, ythick=postthick, position=pos[*,0]

    legend, textoidl(['S/N EW([OIII] 5007) > 3',$;'EW([OIII 5007]) > 0',$
      'F_{[O III] 5007} > 2.6x10^{-17} cgs']), $
      /left, /top, box=0, charsize=1.5, charthick=postthick, spacing=2.0, margin=0
;   legend, textoidl(['N_{gal} = '+strtrim(ncut,2),'F_{[O III] 5007} > 2.6x10^{-17} cgs']), $
;     /left, /top, box=0, charsize=1.8, charthick=postthick, spacing=2.0, margin=0
    legend, textoidl(['N_{gal} = '+strtrim(ncut,2)]), /right, /top, box=0, charsize=1.8, charthick=postthick

    if keyword_set(ps) then begin
       dfpsclose
       spawn, 'convert '+analysis_path+psname+' '+analysis_path+repstr(psname,'.ps','.png'), /sh
    endif
    
stop

return
end
    
