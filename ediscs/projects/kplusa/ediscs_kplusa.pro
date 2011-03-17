function form_model, model1, zform=zform, zobserve=zobserve, $
  theseages=theseages, zaxis=zaxis, lookback=lookback
; apply a formation redshift to a given BC03 model

    if (n_elements(zform) eq 0L) then zform = 3.0
    if (n_elements(zobserve) eq 0L) then zobserve = 0.0
    if (n_elements(theseages) eq 0L) then $
      ages = model1.age/1E9 else $ ; [Gyr]
      ages = theseages
    tform = getage(zform)
    finalage = getage(zobserve)-getage(zform)

    get_element, model1.age/1E9, ages, these
    model = model1[these]
    
    zaxis = getredshift(model.age/1E9+tform)
    lookback = getage(zobserve)-getage(zaxis)

    good = where(zaxis gt zobserve,ngood)
    if (ngood eq 0L) then message, 'Go directly to jail - do not pass go.'
    zaxis = zaxis[good]
    lookback = lookback[good]
    
return, model[good]
end
    

pro ediscs_kplusa, postscript=postscript, encapsulated=encapsulated, pdf=pdf
; jm08may18nyu - 

; read the data

    cluster1 = read_ediscs_kplusa_sample(/cluster1)
    field1 = read_ediscs_kplusa_sample(/field1)

; constants

    red, h100=0.7, omega0=0.3, omega_lambda=0.7
    h100 = redh100()
    
    timelabel1 = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0] ; [Gyr]

; initialize some path names and plotting variables

    kplusapath = ediscs_path(/projects)+'kplusa/'
    paperpath = ediscs_path(/papers)+'kplusa/'
    pspath = paperpath+'FIG_KPLUSA/'
    bc03path = getenv('BC03_SFHGRID_DIR')+'/measure/' ; SFH grid path

    postthick1 = 2.0
    postthick2 = 2.0
    postthick3 = 2.0
    postthick4 = 2.0
    postthick5 = 2.0
    postthick6 = 1.0
    postthick7 = 2.0
    
    textcolor1 = 'white'
    axis_color = 'white'

    cluster1psize = 0.9 & cluster1sym = 106 & cluster1color = 'dark grey'
    field1psize = 0.9 & field1sym = 106 & field1color = 'dodger blue'

    if (not keyword_set(postscript)) and (not keyword_set(pdf)) then $
      im_window, 0, xratio=0.5, /square

    if keyword_set(pdf) then begin

       postscript = 0L
       encapsulated = 0L
       pspath = paperpath+'/keynote/'

       textcolor1 = 'white'
       axis_color = 'black'

       postthick1 = 6.0 
       postthick2 = 6.0 
       postthick3 = 8.0 
       postthick4 = 4.0 
       postthick5 = 8.0 
       postthick6 = 1.0 
       postthick7 = 10.0 

    endif
    
    if keyword_set(postscript) then begin
    
       textcolor1 = 'black'
       axis_color = 'black'
 
       postthick1 = 4.0 
       postthick2 = 3.0 
       postthick3 = 8.0 
       postthick4 = 4.0 
       postthick5 = 6.0 
       postthick6 = 1.0 
       postthick7 = 10.0 

    endif

    if keyword_set(encapsulated) then armexten = '.eps' else armexten = '.ps'

    charsize_0 = 1.0
    charsize_1 = 1.1
    charsize_2 = 1.2
    charsize_3 = 1.3
    charsize_4 = 1.4
    charsize_5 = 1.5
    charsize_6 = 1.6
    charsize_7 = 1.7
    charsize_8 = 1.8
    charsize_9 = 1.9
    singlecharsize_0 = 2.0
    singlecharsize_1 = 2.1
    singlecharsize_2 = 2.2
    singlecharsize_3 = 2.3
    singlecharsize_4 = 2.4
    singlecharsize_5 = 2.5
    charsize_30 = 3.0

; read some BC03 models

    labelages = [0.1,0.5,1.0,1.5,3.0,5.0,6.0,10.0]

    const1psize = 1.6 & const1sym = 108 & const1color = 'firebrick'
    const2psize = 1.6 & const2sym = 108 & const2color = 'cyan'
    tau0psize   = 1.6 & tau0sym   = 108 & tau0color   = 'forest green'
    tau1psize   = 1.6 & tau1sym   = 108 & tau1color   = 'navy'
    tau6psize   = 1.6 & tau6sym   = 108 & tau6color   = 'orchid'
    burst1psize = 1.6 & burst1sym = 108 & burst1color = 'orange'

    const1 = mrdfits(bc03path+'salp_m62_const_20.0Gyr.info.fits',1)
    const2 = mrdfits(bc03path+'salp_m62_const_1.0Gyr.info.fits',1)
    tau0 = mrdfits(bc03path+'salp_m62_tau_00.0Gyr.info.fits',1)
    tau1 = mrdfits(bc03path+'salp_m62_tau_01.0Gyr.info.fits',1)
    tau6 = mrdfits(bc03path+'salp_m62_tau_02.0Gyr.info.fits',1)
    burst1 = mrdfits(bc03path+'salp_m62_tau_01.0Gyr_tb_05.0Gyr_dtb_0.10Gyr_fb_0.20.info.fits',1)

    zform = 3.0    ; formation redshift
    zobserve = 0.4 ; redshift of observation
    
    const1 = form_model(const1,zform=zform,zobserve=zobserve,$
      theseages=theseages,zaxis=const1_zaxis)
    const2 = form_model(const2,zform=zform,zobserve=zobserve,$
      theseages=theseages,zaxis=const2_zaxis)
    tau0 = form_model(tau0,zform=zform,zobserve=zobserve,$
      theseages=theseages,zaxis=tau0_zaxis)
    tau1 = form_model(tau1,zform=zform,zobserve=zobserve,$
      theseages=theseages,zaxis=tau1_zaxis)
    tau6 = form_model(tau6,zform=zform,zobserve=zobserve,$
      theseages=theseages,zaxis=tau6_zaxis)
    burst1 = form_model(burst1,zform=zform,zobserve=zobserve,$
      theseages=theseages,zaxis=burst1_zaxis)

; ------------------------------------------------------------
; D(4000) versus H-delta_A - all cluster/field galaxies
; ------------------------------------------------------------
    
    xpage = 8.5 & ypage = 8.5

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=7.0, height=7.0, $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=xpage, ypage=ypage, $
      position=pos, /normal

    xrange = [0.85,2.1]
    yrange = [-2,12]

    xtitle = 'D_{n}(4000)'
    ytitle = 'H\delta_{A} (\AA)'

    clindx = lindgen(n_elements(cluster1))
    xcluster1 = cluster1[clindx].d4000_narrow[0] & xcluster1err = cluster1[clindx].d4000_narrow[1]
    ycluster1 = cluster1[clindx].lick_hd_a[0] & ycluster1err = cluster1[clindx].lick_hd_a[1]
    xcluster2 = cluster1[clindx].d4000_narrow_model[0] & xcluster2err = cluster1[clindx].d4000_narrow_model[1]
    ycluster2 = cluster1[clindx].lick_hd_a_model[0] & ycluster2err = cluster1[clindx].lick_hd_a_model[1]

; for the field sample we want: which field galaxies are within DZ of
; each *cluster* redshift (not cluster *galaxy* redshift)

    cl_z = cluster1[clindx[uniq(cluster1[clindx].cluster_z,$
      sort(cluster1[clindx].cluster_z))]].cluster_z

    dz = 0.05
    for ii = 0L, n_elements(cl_z)-1L do begin
       these = where((abs(cl_z[ii]-field1.z) lt dz),nthese)
       if (ii eq 0L) then fldindx = these else fldindx = [fldindx,these]
    endfor
    fldindx = fldindx[uniq(fldindx,sort(fldindx))]
    
    fldindx = lindgen(n_elements(field1))
    xfield1 = field1[fldindx].d4000_narrow_model[0]
    yfield1 = field1[fldindx].lick_hd_a_model[0]
    
; make the plot - cluster galaxies
    
    psname = 'd4000_vs_hda_cluster'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

    plot, [0], [0], /nodata, ysty=1, xsty=1, xrange=xrange, yrange=yrange, $
      xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=singlecharsize_0, charthick=postthick2, xthick=postthick1, $
      ythick=postthick1, position=pos[*,0], color=fsc_color(textcolor1,100)

    im_symbols, cluster1sym, psize=cluster1psize, thick=postthick1, $
      fill=1, color=fsc_color(cluster1color,1)
    djs_oplot, xcluster1, ycluster1, ps=8, color=fsc_color(cluster1color,50)
    im_symbols, cluster1sym, psize=cluster1psize, thick=postthick1, $
      fill=1, color=fsc_color('forest green',50)
    djs_oplot, xcluster2, ycluster2, ps=8, sym=0.8

; overlay some BC03 models

    djs_oplot, const1.d4000_narrow, const1.lick_hd_a, line=0, $
      thick=postthick3, color=fsc_color(const1color,90)
    im_symbols, const1sym, psize=const1psize, thick=postthick1, $
      fill=1, color=fsc_color(const1color,90)
    for iage = 0L, n_elements(labelages)-1L do begin
       if (labelages[iage] le max(const1.age/1E9)) and $
         (labelages[iage] ge min(const1.age/1E9)) then begin
          lxx = interpol(const1.d4000_narrow,const1.age/1E9,labelages[iage])
          lyy = interpol(const1.lick_hd_a,const1.age/1E9,labelages[iage])
          plots, lxx, lyy, ps=8
          xyouts, lxx, lyy, string(labelages[iage],format='(G0)'), $
            align=0.0, charsize=1.2, charthick=postthick2
       endif
    endfor
    
    djs_oplot, tau0.d4000_narrow, tau0.lick_hd_a, line=0, $
      thick=postthick3, color=fsc_color(tau0color,90)
    im_symbols, tau0sym, psize=tau0psize, thick=postthick1, $
      fill=1, color=fsc_color(tau0color,90)
    for iage = 0L, n_elements(labelages)-1L do begin
       if (labelages[iage] le max(tau0.age/1E9)) and $
         (labelages[iage] ge min(tau0.age/1E9)) then begin
          lxx = interpol(tau0.d4000_narrow,tau0.age/1E9,labelages[iage])
          lyy = interpol(tau0.lick_hd_a,tau0.age/1E9,labelages[iage])
          plots, lxx, lyy, ps=8
          xyouts, lxx, lyy, string(labelages[iage],format='(G0)'), $
            align=0.0, charsize=1.2, charthick=postthick2
       endif
    endfor

    djs_oplot, tau1.d4000_narrow, tau1.lick_hd_a, line=0, $
      thick=postthick3, color=fsc_color(tau1color,90)
    im_symbols, tau1sym, psize=tau1psize, thick=postthick1, $
      fill=1, color=fsc_color(tau1color,90)
    for iage = 0L, n_elements(labelages)-1L do begin
       if (labelages[iage] le max(tau1.age/1E9)) and $
         (labelages[iage] ge min(tau1.age/1E9)) then begin
          lxx = interpol(tau1.d4000_narrow,tau1.age/1E9,labelages[iage])
          lyy = interpol(tau1.lick_hd_a,tau1.age/1E9,labelages[iage])
          plots, lxx, lyy, ps=8
          xyouts, lxx, lyy, string(labelages[iage],format='(G0)'), $
            align=0.0, charsize=1.2, charthick=postthick2
       endif
    endfor

    djs_oplot, tau6.d4000_narrow, tau6.lick_hd_a, line=0, $
      thick=postthick3, color=fsc_color(tau6color,90)
    im_symbols, tau6sym, psize=tau6psize, thick=postthick1, $
      fill=1, color=fsc_color(tau6color,90)
    for iage = 0L, n_elements(labelages)-1L do begin
       if (labelages[iage] le max(tau6.age/1E9)) and $
         (labelages[iage] ge min(tau6.age/1E9)) then begin
          lxx = interpol(tau6.d4000_narrow,tau6.age/1E9,labelages[iage])
          lyy = interpol(tau6.lick_hd_a,tau6.age/1E9,labelages[iage])
          plots, lxx, lyy, ps=8
          xyouts, lxx, lyy, string(labelages[iage],format='(G0)'), $
            align=0.0, charsize=1.2, charthick=postthick2
       endif
    endfor

;   djs_oplot, burst1.d4000_narrow, burst1.lick_hd_a, line=0, $
;     thick=postthick3, color=fsc_color(burst1color,90)
;   im_symbols, burst1sym, psize=burst1psize, thick=postthick1, $
;     fill=1, color=fsc_color(burst1color,90)
;   for iage = 0L, n_elements(labelages)-1L do begin
;      if (labelages[iage] le max(burst1.age/1E9)) and $
;        (labelages[iage] ge min(burst1.age/1E9)) then begin
;         lxx = interpol(burst1.d4000_narrow,burst1.age/1E9,labelages[iage])
;         lyy = interpol(burst1.lick_hd_a,burst1.age/1E9,labelages[iage])
;         plots, lxx, lyy, ps=8
;         xyouts, lxx, lyy, string(labelages[iage],format='(G0)'), $
;           align=0.0, charsize=1.2, charthick=postthick2
;      endif
;   endfor
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; make the plot - field galaxies
    
    psname = 'd4000_vs_hda_field'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

    plot, [0], [0], /nodata, ysty=1, xsty=1, xrange=xrange, yrange=yrange, $
      xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=singlecharsize_0, charthick=postthick2, xthick=postthick1, $
      ythick=postthick1, position=pos[*,0], color=fsc_color(textcolor1,100)

    im_symbols, field1sym, psize=field1psize, thick=postthick1, $
      fill=1, color=fsc_color(field1color,1)
    djs_oplot, xfield1, yfield1, ps=8, color=fsc_color(field1color,50)

; overlay some BC03 models

    djs_oplot, const1.d4000_narrow, const1.lick_hd_a, line=0, $
      thick=postthick3, color=fsc_color(const1color,90)
    im_symbols, const1sym, psize=const1psize, thick=postthick1, $
      fill=1, color=fsc_color(const1color,90)
    for iage = 0L, n_elements(labelages)-1L do begin
       if (labelages[iage] le max(const1.age/1E9)) and $
         (labelages[iage] ge min(const1.age/1E9)) then begin
          lxx = interpol(const1.d4000_narrow,const1.age/1E9,labelages[iage])
          lyy = interpol(const1.lick_hd_a,const1.age/1E9,labelages[iage])
          plots, lxx, lyy, ps=8
          xyouts, lxx, lyy, string(labelages[iage],format='(G0)'), $
            align=0.0, charsize=1.2, charthick=postthick2
       endif
    endfor
    
    djs_oplot, tau0.d4000_narrow, tau0.lick_hd_a, line=0, $
      thick=postthick3, color=fsc_color(tau0color,90)
    im_symbols, tau0sym, psize=tau0psize, thick=postthick1, $
      fill=1, color=fsc_color(tau0color,90)
    for iage = 0L, n_elements(labelages)-1L do begin
       if (labelages[iage] le max(tau0.age/1E9)) and $
         (labelages[iage] ge min(tau0.age/1E9)) then begin
          lxx = interpol(tau0.d4000_narrow,tau0.age/1E9,labelages[iage])
          lyy = interpol(tau0.lick_hd_a,tau0.age/1E9,labelages[iage])
          plots, lxx, lyy, ps=8
          xyouts, lxx, lyy, string(labelages[iage],format='(G0)'), $
            align=0.0, charsize=1.2, charthick=postthick2
       endif
    endfor

    djs_oplot, burst1.d4000_narrow, burst1.lick_hd_a, line=0, $
      thick=postthick3, color=fsc_color(burst1color,90)
    im_symbols, burst1sym, psize=burst1psize, thick=postthick1, $
      fill=1, color=fsc_color(burst1color,90)
    for iage = 0L, n_elements(labelages)-1L do begin
       if (labelages[iage] le max(burst1.age/1E9)) and $
         (labelages[iage] ge min(burst1.age/1E9)) then begin
          lxx = interpol(burst1.d4000_narrow,burst1.age/1E9,labelages[iage])
          lyy = interpol(burst1.lick_hd_a,burst1.age/1E9,labelages[iage])
          plots, lxx, lyy, ps=8
          xyouts, lxx, lyy, string(labelages[iage],format='(G0)'), $
            align=0.0, charsize=1.2, charthick=postthick2
       endif
    endfor
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    

; ------------------------------------------------------------
; D(4000) versus H-delta_A - 4-panel plot; cluster + field
; ------------------------------------------------------------
    
    xpage = 8.5 & ypage = 8.5

    pagemaker, nx=2, ny=2, xspace=0, yspace=0.0, width=3.5*[1,1], height=3.5*[1,1], $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=xpage, ypage=ypage, $
      position=pos, /normal

    xrange = [0.85,2.1]
    yrange = [-1.9,11]

    xtitle = 'D_{n}(4000)'
    ytitle = 'H\delta_{A} (\AA)'

    cl1 = where((cluster1.cluster_z lt 0.6) and (cluster1.cluster_sigma lt 500.0),ncl1)
    cl2 = where((cluster1.cluster_z gt 0.6) and (cluster1.cluster_sigma lt 500.0),ncl1)
    cl3 = where((cluster1.cluster_z lt 0.6) and (cluster1.cluster_sigma gt 500.0),ncl2)
    cl4 = where((cluster1.cluster_z gt 0.6) and (cluster1.cluster_sigma gt 500.0),ncl2)

    cl1_z = cluster1[cl1[uniq(cluster1[cl1].cluster_z,sort(cluster1[cl1].cluster_z))]].cluster_z
    cl2_z = cluster1[cl2[uniq(cluster1[cl2].cluster_z,sort(cluster1[cl2].cluster_z))]].cluster_z
    cl3_z = cluster1[cl3[uniq(cluster1[cl3].cluster_z,sort(cluster1[cl3].cluster_z))]].cluster_z
    cl4_z = cluster1[cl4[uniq(cluster1[cl4].cluster_z,sort(cluster1[cl4].cluster_z))]].cluster_z
    
; for the field sample we want to ask: which field galaxies are within
; dz=0.1 of each cluster in each of the above bins of redshift and
; velocity dispersion

    dz = 0.05

    for ii = 0L, n_elements(cl1_z)-1L do begin
       these = where((abs(cl1_z[ii]-field1.z) lt dz),nthese)
       if (ii eq 0L) then fld1 = these else fld1 = [fld1,these]
    endfor
    fld1 = fld1[uniq(fld1,sort(fld1))]

    for ii = 0L, n_elements(cl2_z)-1L do begin
       these = where((abs(cl2_z[ii]-field1.z) lt dz),nthese)
       if (ii eq 0L) then fld2 = these else fld2 = [fld2,these]
    endfor
    fld2 = fld2[uniq(fld2,sort(fld2))]

    for ii = 0L, n_elements(cl3_z)-1L do begin
       these = where((abs(cl3_z[ii]-field1.z) lt dz),nthese)
       if (ii eq 0L) then fld3 = these else fld3 = [fld3,these]
    endfor
    fld3 = fld3[uniq(fld3,sort(fld3))]

    for ii = 0L, n_elements(cl4_z)-1L do begin
       these = where((abs(cl4_z[ii]-field1.z) lt dz),nthese)
       if (ii eq 0L) then fld4 = these else fld4 = [fld4,these]
    endfor
    fld4 = fld4[uniq(fld4,sort(fld4))]
    
    xcl1 = cluster1[cl1].d4000_narrow_model[0] & ycl1 = cluster1[cl1].lick_hd_a_model[0]
    xcl2 = cluster1[cl2].d4000_narrow_model[0] & ycl2 = cluster1[cl2].lick_hd_a_model[0]
    xcl3 = cluster1[cl3].d4000_narrow_model[0] & ycl3 = cluster1[cl3].lick_hd_a_model[0]
    xcl4 = cluster1[cl4].d4000_narrow_model[0] & ycl4 = cluster1[cl4].lick_hd_a_model[0]

    xfld1 = field1[fld1].d4000_narrow_model[0] & yfld1 = field1[fld1].lick_hd_a_model[0]
    xfld2 = field1[fld2].d4000_narrow_model[0] & yfld2 = field1[fld2].lick_hd_a_model[0]
    xfld3 = field1[fld3].d4000_narrow_model[0] & yfld3 = field1[fld3].lick_hd_a_model[0]
    xfld4 = field1[fld4].d4000_narrow_model[0] & yfld4 = field1[fld4].lick_hd_a_model[0]
    
; cluster galaxies
    
    psname = 'd4000_vs_hda_cluster_4panel'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

    im_symbols, cluster1sym, psize=cluster1psize, thick=postthick1, $
      fill=1, color=fsc_color(cluster1color,1)

    plot, [0], [0], /nodata, ysty=1, xsty=1, xrange=xrange, yrange=yrange, $
      xtitle='', ytitle=textoidl(ytitle), $
      xtickname=replicate(' ',10), $
      charsize=charsize_8, charthick=postthick2, xthick=postthick1, $
      ythick=postthick1, position=pos[*,0], color=fsc_color(textcolor1,100)
    djs_oplot, xcl1, ycl1, ps=8, color=fsc_color(cluster1color,50)
    legend, textoidl(['z < 0.6','\sigma < 500 km s^{-1}']), /left, /bottom, $
      box=0, charsize=charsize_5, charthick=postthick2

    plot, [0], [0], /nodata, /noerase, ysty=1, xsty=1, xrange=xrange, yrange=yrange, $
      xtitle='', ytitle='', $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      charsize=charsize_8, charthick=postthick2, xthick=postthick1, $
      ythick=postthick1, position=pos[*,1], color=fsc_color(textcolor1,100)
    djs_oplot, xcl2, ycl2, ps=8, color=fsc_color(cluster1color,50)
    legend, textoidl(['z > 0.6','\sigma < 500 km s^{-1}']), /left, /bottom, $
      box=0, charsize=charsize_5, charthick=postthick2

    plot, [0], [0], /nodata, /noerase, ysty=1, xsty=1, xrange=xrange, yrange=yrange, $
      xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=charsize_8, charthick=postthick2, xthick=postthick1, $
      ythick=postthick1, position=pos[*,2], color=fsc_color(textcolor1,100)
    djs_oplot, xcl3, ycl3, ps=8, color=fsc_color(cluster1color,50)
    legend, textoidl(['z < 0.6','\sigma > 500 km s^{-1}']), /left, /bottom, $
      box=0, charsize=charsize_5, charthick=postthick2

    plot, [0], [0], /nodata, /noerase, ysty=1, xsty=1, xrange=xrange, yrange=yrange, $
      xtitle=textoidl(xtitle), ytitle='', $
      ytickname=replicate(' ',10), $
      charsize=charsize_8, charthick=postthick2, xthick=postthick1, $
      ythick=postthick1, position=pos[*,3], color=fsc_color(textcolor1,100)
    djs_oplot, xcl4, ycl4, ps=8, color=fsc_color(cluster1color,50)
    legend, textoidl(['z > 0.6','\sigma > 500 km s^{-1}']), /left, /bottom, $
      box=0, charsize=charsize_5, charthick=postthick2

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; field galaxies
    
    psname = 'd4000_vs_hda_field_4panel'
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

    im_symbols, field1sym, psize=field1psize, thick=postthick1, $
      fill=1, color=fsc_color(field1color,1)

    plot, [0], [0], /nodata, ysty=1, xsty=1, xrange=xrange, yrange=yrange, $
      xtitle='', ytitle=textoidl(ytitle), $
      xtickname=replicate(' ',10), $
      charsize=charsize_8, charthick=postthick2, xthick=postthick1, $
      ythick=postthick1, position=pos[*,0], color=fsc_color(textcolor1,100)
    djs_oplot, xfld1, yfld1, ps=8, color=fsc_color(field1color,50)
    legend, textoidl(['z < 0.6','\sigma < 500 km s^{-1}']), /left, /bottom, $
      box=0, charsize=charsize_5, charthick=postthick2

    plot, [0], [0], /nodata, /noerase, ysty=1, xsty=1, xrange=xrange, yrange=yrange, $
      xtitle='', ytitle='', $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      charsize=charsize_8, charthick=postthick2, xthick=postthick1, $
      ythick=postthick1, position=pos[*,1], color=fsc_color(textcolor1,100)
    djs_oplot, xfld2, yfld2, ps=8, color=fsc_color(field1color,50)
    legend, textoidl(['z > 0.6','\sigma < 500 km s^{-1}']), /left, /bottom, $
      box=0, charsize=charsize_5, charthick=postthick2

    plot, [0], [0], /nodata, /noerase, ysty=1, xsty=1, xrange=xrange, yrange=yrange, $
      xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=charsize_8, charthick=postthick2, xthick=postthick1, $
      ythick=postthick1, position=pos[*,2], color=fsc_color(textcolor1,100)
    djs_oplot, xfld3, yfld3, ps=8, color=fsc_color(field1color,50)
    legend, textoidl(['z < 0.6','\sigma > 500 km s^{-1}']), /left, /bottom, $
      box=0, charsize=charsize_5, charthick=postthick2

    plot, [0], [0], /nodata, /noerase, ysty=1, xsty=1, xrange=xrange, yrange=yrange, $
      xtitle=textoidl(xtitle), ytitle='', $
      ytickname=replicate(' ',10), $
      charsize=charsize_8, charthick=postthick2, xthick=postthick1, $
      ythick=postthick1, position=pos[*,3], color=fsc_color(textcolor1,100)
    djs_oplot, xfld4, yfld4, ps=8, color=fsc_color(field1color,50)
    legend, textoidl(['z > 0.6','\sigma > 500 km s^{-1}']), /left, /bottom, $
      box=0, charsize=charsize_5, charthick=postthick2

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; ------------------------------------------------------------
; M_B vs B-V 
; ------------------------------------------------------------
    
    psname = 'mb_vs_bv'
    xpage = 8.5 & ypage = 8.5
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=7.0, height=7.0, $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=xpage, ypage=ypage, $
      position=pos, /normal

    xrange = [-17.0,-23.5]
    yrange = [-0.1,1.15]

    xtitle = 'M_{B}'
    ytitle = 'B - V'

    clindx = where((cluster1.ubvrijhk_absmag[1] gt -900.0) and $
      (cluster1.ubvrijhk_absmag[2] gt -900.0),nclindx)
    xcluster1 = cluster1[clindx].ubvrijhk_absmag[1]
    ycluster1 = cluster1[clindx].ubvrijhk_absmag[1]-cluster1.ubvrijhk_absmag[2]
    
    fldindx = where((field1.ubvrijhk_absmag[1] gt -900.0) and $
      (field1.ubvrijhk_absmag[2] gt -900.0),nfldindx)
    xfield1 = field1[fldindx].ubvrijhk_absmag[1]
    yfield1 = field1[fldindx].ubvrijhk_absmag[1]-field1.ubvrijhk_absmag[2]
    
; make the plot    
    
    plot, [0], [0], /nodata, ysty=1, xsty=1, xrange=xrange, yrange=yrange, $
      xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=singlecharsize_0, charthick=postthick2, xthick=postthick1, $
      ythick=postthick1, position=pos[*,0], color=fsc_color(textcolor1,100)

; cluster galaxies    
    im_symbols, cluster1sym, psize=cluster1psize, thick=postthick1, $
      fill=1, color=fsc_color(cluster1color,1)
    djs_oplot, xcluster1, ycluster1, ps=8, color=fsc_color(cluster1color,50)

; field galaxies    
    im_symbols, field1sym, psize=field1psize, thick=postthick1, $
      fill=1, color=fsc_color(field1color,1)
    djs_oplot, xfield1, yfield1, ps=8, color=fsc_color(field1color,50)

; overlay some BC03 models

    djs_oplot, const1.ubvrijhk[1], const1.ubvrijhk[1]-const1.ubvrijhk[2], line=0, $
      thick=postthick3, color=fsc_color(const1color,90)
    im_symbols, const1sym, psize=const1psize, thick=postthick1, $
      fill=1, color=fsc_color(const1color,90)
    for iage = 0L, n_elements(labelages)-1L do begin
       if (labelages[iage] le max(const1.age/1E9)) and $
         (labelages[iage] ge min(const1.age/1E9)) then begin
          lxx = interpol(const1.ubvrijhk[1],const1.age/1E9,labelages[iage])
          lyy = interpol(const1.ubvrijhk[1]-const1.ubvrijhk[2],const1.age/1E9,labelages[iage])
          plots, lxx, lyy, ps=8
          xyouts, lxx, lyy, string(labelages[iage],format='(G0)'), $
            align=0.0, charsize=1.2, charthick=postthick2
       endif
    endfor
    
    djs_oplot, tau0.ubvrijhk[1], tau0.ubvrijhk[1]-tau0.ubvrijhk[2], line=0, $
      thick=postthick3, color=fsc_color(tau0color,90)
    im_symbols, tau0sym, psize=tau0psize, thick=postthick1, $
      fill=1, color=fsc_color(tau0color,90)
    for iage = 0L, n_elements(labelages)-1L do begin
       if (labelages[iage] le max(tau0.age/1E9)) and $
         (labelages[iage] ge min(tau0.age/1E9)) then begin
          lxx = interpol(tau0.ubvrijhk[1],tau0.age/1E9,labelages[iage])
          lyy = interpol(tau0.ubvrijhk[1]-tau0.ubvrijhk[2],tau0.age/1E9,labelages[iage])
          plots, lxx, lyy, ps=8
          xyouts, lxx, lyy, string(labelages[iage],format='(G0)'), $
            align=0.0, charsize=1.2, charthick=postthick2
       endif
    endfor

    djs_oplot, burst1.ubvrijhk[1], burst1.ubvrijhk[1]-burst1.ubvrijhk[2], line=0, $
      thick=postthick3, color=fsc_color(burst1color,90)
    im_symbols, burst1sym, psize=burst1psize, thick=postthick1, $
      fill=1, color=fsc_color(burst1color,90)
    for iage = 0L, n_elements(labelages)-1L do begin
       if (labelages[iage] le max(burst1.age/1E9)) and $
         (labelages[iage] ge min(burst1.age/1E9)) then begin
          lxx = interpol(burst1.ubvrijhk[1],burst1.age/1E9,labelages[iage])
          lyy = interpol(burst1.ubvrijhk[1]-burst1.ubvrijhk[2],burst1.age/1E9,labelages[iage])
          plots, lxx, lyy, ps=8
          xyouts, lxx, lyy, string(labelages[iage],format='(G0)'), $
            align=0.0, charsize=1.2, charthick=postthick2
       endif
    endfor
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    

; ------------------------------------------------------------
; redshift versus M_B - sample 1
; ------------------------------------------------------------
    
    psname = 'redshift_vs_mb'
    xpage = 8.5 & ypage = 8.5
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.5, height=6.5, $
      xmargin=[1.6,0.4], ymargin=[1.0,1.0], xpage=xpage, ypage=ypage, $
      position=pos, /normal

    xrange = [0.35,1.1]
    yrange = [-17.0,-23.5]

    xtitle = 'Redshift'
    ytitle = 'M_{B}'

    clindx = where((cluster1.z gt -900.0) and $
      (cluster1.ubvrijhk_absmag[1] gt -900.0),nclindx)
    xcluster1 = cluster1[clindx].z
    ycluster1 = cluster1[clindx].ubvrijhk_absmag[1]
    
; make the plot    
    
    plot, [0], [0], /nodata, ysty=1, xsty=9, xrange=xrange, yrange=yrange, $
      xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), $
      charsize=singlecharsize_0, charthick=postthick2, xthick=postthick1, $
      ythick=postthick1, position=pos[*,0], color=fsc_color(textcolor1,100)
    axis, /xaxis, xthick=postthick1, xsty=1, xrange=xrange, charsize=singlecharsize_0, $
      charthick=postthick2, xtitle='Lookback Time (Gyr)', color=fsc_color(textcolor1,100), $
      xtickv=getredshift(getage(0.0)-timelabel1), xticks=n_elements(timelabel1)-1L, $
      xtickname=string(timelabel1,format='(I0)')

    im_symbols, cluster1sym, psize=cluster1psize, thick=postthick1, $
      fill=1, color=fsc_color(cluster1color,1)
    djs_oplot, xcluster1, ycluster1, ps=8, color=fsc_color(cluster1color,50)
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    

return
end
