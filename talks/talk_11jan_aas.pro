pro talk_11jan_aas, keynote=keynote
; jm10jan07ucsd - miscellaneous plots for my 2011/Jan AAS talk 

    common com_talk, sdss_parent
    
    mfpath = mf_path()
    isedpath = mf_path(/isedfit)
    
    rootpath = getenv('IM_RESEARCH_DIR')+'/talks/2011/11jan_aas/'
    if keyword_set(keynote) then talkpath = rootpath+'keynote/' else $
      talkpath = rootpath+'figures/'

; --------------------------------------------------
; remake Fig 8 in Ken's paper
    psfile = talkpath+'11wong_fig7.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, $
      height=5.2, xmargin=[1.7,0.2], width=6.6, charsize=2

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='Redshift', xrange=[0.27,0.71], $
      yrange=[0.09,-0.4], ytitle=textoidl('\Delta<FUV-r>'), $
      color=keycolor, ytickinterval=0.1
    oplot, !x.crange, [0,0], line=1, thick=4, color=keycolor
    oploterror, [mean([0.25,0.5]),mean([0.5,0.75])], $  ; 30 kpc
      [-0.232,-0.109], [0.080,0.112], psym=-symcat(6,thick=8), line=0, $
      color=fsc_color('tan',101), errcolor=fsc_color('tan',101), $
      thick=8, errthick=8, symsize=3
    oploterror, [mean([0.25,0.5]),mean([0.5,0.75])]-0.03, $ ; 50 kpc
      [-0.191,-0.034], [0.055,0.076], psym=-symcat(15), line=5, $
      color=fsc_color('dodger blue',101), errcolor=fsc_color('dodger blue',101), $
      thick=8, errthick=8, symsize=3
    im_legend, ['r_{p}<30 h^{-1} kpc','r_{p}<50 h^{-1} kpc'], $
      /right, /top, box=0, charsize=1.7, $
      pspacing=1.9, line=[0,5], color=['tan','dodger blue'], $
      textcolor=keycolor, psym=[6,15], symsize=2, thick=8, symthick=8
    
    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

stop    


; --------------------------------------------------
; zvz - read the data
    zerod = primus_read_zerod(field='zspec',zmerge=zmerge)
    zcut = where((zmerge.z_spec gt 0.0) and (zmerge.z_spec lt 1.2) and $
      (zmerge.good_spec eq 1) and (zerod.zprimus_zconf ge 3) and $
      (strtrim(zerod.zprimus_class,2) eq 'GALAXY') and $
      ((zerod.fieldnum lt 3 and zerod.targ_mag le 22.5) or $
      (zerod.fieldnum eq 4 and zerod.targ_mag le 22.8) or $
      (zerod.fieldnum eq 5 and zerod.targ_mag le 22.8) or $
      (zerod.fieldnum eq 6 and zerod.targ_mag le 22.5)),nzconf3)

    allfield = strtrim(zmerge[zcut].field,2)
    allzconf = zerod[zcut].zprimus_zconf
    allzprimus = zerod[zcut].zprimus
    allzspec = zmerge[zcut].z_spec

; use just a fraction of the deep2 sources
    isdeep = where(strmatch(allfield,'*deep*'),nisdeep,comp=notdeep)
    nkeep = long(nisdeep*0.2)
    keep = shuffle_indx(nisdeep,num_sub=nkeep,seed=seed)

    zconf = [allzconf[notdeep],allzconf[isdeep[keep]]]
    zprimus = [allzprimus[notdeep],allzprimus[isdeep[keep]]]
    zspec = [allzspec[notdeep],allzspec[isdeep[keep]]]
    field = [allfield[notdeep],allfield[isdeep[keep]]]

; shuffle the indices so that they don't overlap on the plot
    npts = n_elements(zprimus)
    indx = shuffle_indx(npts)
    zconf = zconf[indx]
    zprimus = zprimus[indx]
    zspec = zspec[indx]
    field = field[indx]

; make the plot    
    psfile = talkpath+'zvz.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, charsize=1.7, $
      width=7.0, height=5.0

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')
    
; plotting preferences    
    nfield = 4
    prefs = replicate({field: '', nicefield: '', color: '', $
      psym: 0, symsize: 1.0},nfield)
    prefs.field = ['cosmos','calib','deep2_02hr','deep2_23hr']
    prefs.nicefield = ['zCOSMOS','VVDS','DEEP2 02^{hr}','DEEP2 23^{hr}']
    if keyword_set(keynote) then $
      prefs.color = ['white','medium orchid','dodger blue','tan'] else $
        prefs.color = ['black','medium orchid','dodger blue','tan']
    prefs.psym = [14,9,15,4]
    prefs.symsize = [1.1,0.9,0.9,1.1]*0.7

    zrange = [-0.06,1.27]
    xtitle1 = 'Redshift from High-Resolution Spectroscopy'
    ytitle1 = 'PRIMUS Redshift'

    zaxis = findgen((2.0-(-1.0))/0.01)*0.01+(-1.0)
    dzcut_03 = 0.03*(1.0+zaxis)
    dzcut_10 = 0.10*(1.0+zaxis)

; Q=4
    plot, [0], [0], /nodata, /xsty, /ysty, position=pos[*,0], $
      xrange=zrange, yrange=zrange, xtitle=xtitle1, ytitle=ytitle1, $
      color=keycolor
    best = where(zconf eq 4,nbest)
    for ii = 0, nbest-1 do begin
       this = where(field[best[ii]] eq strtrim(prefs.field,2))
       plots, zspec[best[ii]], zprimus[best[ii]], psym=symcat(prefs[this].psym,thick=6), $
         color=fsc_color(prefs[this].color,101), symsize=prefs[this].symsize
    endfor
    oplot, zaxis, zaxis, line=0, thick=5, color=keycolor
    oplot, zaxis, zaxis-dzcut_03, line=5, thick=5, color=keycolor
    oplot, zaxis, zaxis+dzcut_03, line=5, thick=5, color=keycolor
    oplot, zaxis, zaxis-dzcut_10, line=1, thick=5, color=keycolor
    oplot, zaxis, zaxis+dzcut_10, line=1, thick=5, color=keycolor
    
    stats = primus_zvz_stats(zspec[best],zprimus[best])
    label = [$
      'F(>0.03)='+strtrim(string(100.0*stats.frac_bad_03,format='(F12.1)'),2)+'%',$
      'F(>0.10)='+strtrim(string(100.0*stats.frac_bad_10,format='(F12.1)'),2)+'%',$
      '\sigma_{\Delta'+'z/(1+z)}='+strtrim(string(stats.mad,format='(F12.4)'),2)]
    legend, textoidl(label), /left, /top, box=0, charsize=1.4, $
      margin=0, textcolor=keycolor
;   xyouts, 0.43, 0.82, 'Q=4', align=0.5, /data, color=keycolor
    
; legend
    im_legend, prefs.nicefield, textcolor=keycolor, $
      psym=prefs.psym, color=strtrim(prefs.color,2), $
      /right, /bottom, box=0, margin=0, charsize=1.4, $
      symsize=prefs.symsize*2.3, symthick=6
    
    im_plotconfig, psfile=psfile, keynote=keynote, /psclose

return
end
