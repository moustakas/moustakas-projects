pro oplot_hda_d4000, hda, d4000, hdaerr, d4000err, ewoii
; code by EW([OII])

    psym1 = [16,14,15,14]
    color1 = ['firebrick','navy','forest green','orange']
    symsize1 = [1.0,1.3,1.0,1.0]
    locut = [0.0,5.0,10.0,25.0]
    hicut = [5.0,10.0,25.0,500.0]
    ncut = n_elements(locut)
    for ii = 0, ncut-1 do begin
       indx1 = where((ewoii ge locut[ii]) and (ewoii lt hicut[ii]),nindx1)
       if (nindx1 ne 0) then oploterror, d4000[indx1], hda[indx1], $
         d4000err[indx1], hdaerr[indx1], psym=symcat(psym1[ii]), $
         symsize=symsize1[ii], color=fsc_color(color1[ii],101), $
         errcolor=fsc_color(color1[ii],101), errthick=!p.thick
    endfor

    label = string(locut,format='(I0)')+'<EW([O II])<'+string(hicut,format='(I0)')+' \AA'
    label[0] = 'EW([O II])<'+string(hicut[0],format='(I0)')+' \AA'
    label[ncut-1] = 'EW([O II])>'+string(locut[ncut-1],format='(I0)')+' \AA'
    im_legend, label, color=color1, psym=psym1, /left, $
      /bottom, box=0, charsize=1.6, symsize=symsize1*1.5
    
return
end
    
pro plot_ediscs_sfh
; jm10may03ucsd - build plots for the EDisCS/SFH project with Greg 

    sfhpath = ediscs_path()+'sfh/'
    paperpath = sfhpath
    cluster = read_ediscs_sfh_sample(/cluster)
    field = read_ediscs_sfh_sample(/field)
    allcluster = read_ediscs_sfh_sample(/cluster,/all)
    allfield = read_ediscs_sfh_sample(/field,/all)

; --------------------------------------------------
; 4-panel plot for the paper showing the spectra and the spectral fits
; for four objects
;    
;    EDCSNJ1301302-1138187 - old with emission:  EW(OII) = 6.30 +/- 0.80
;    EDCSNJ1018467-1211527
;    EDCSNJ1018454-1212235 or EDCSNJ1040337-1157231 - middle age bin 
;    EDCSNJ1018445-1208545 - young age bin

    gal = [$
      'EDCSNJ1018445-1208545',$ ; young
      'EDCSNJ1018454-1212235',$ ; or EDCSNJ1040337-1157231 intermediate-age
      'EDCSNJ1301302-1138187',$ ; old, with emission
      'EDCSNJ1018467-1211527']  ; old, no emission
    match, strtrim(cluster.galaxy,2), gal, m1, m2
    srt = sort(m2)
    m1 = m1[srt] & m2 = m2[srt]
    nobj = n_elements(gal)

    spec = read_ediscs_gandalf_specfit(cluster[m1])

    psfile = paperpath+'sfh_examples.ps'
    im_plotconfig, 5, pos, psfile=psfile, xmargin=[0.8,0.2], $
;   im_plotconfig, 2, pos, psfile=psfile, xmargin=[1.1,0.2], $
      height=[2.6,2.6], xspace=0.05, yspace=0.05

    label = ['Young','Intermediate','Old (with [OII])','Old (no [OII])']
    scale = 1D
    for ii = 0, nobj-1 do begin
       if odd(ii) then ytickname = replicate(' ',10) else delvarx, ytickname
       if ii le 1 then xtickname = replicate(' ',10) else delvarx, xtickname
       djs_plot, [0], [0], /nodata, xrange=[3420,4900], yrange=[-0.2,2], $
         xsty=1, ysty=1, position=pos[*,ii], noerase=ii gt 0, $
         ytickname=ytickname, xtickname=xtickname, xtickinterval=500, $
         ytickinterval=1.0
       im_legend, gal[ii], /right, /bottom, box=0, margin=0, charsize=1.0, $
         charthick=3.0
       im_legend, label[ii], /right, /top, box=0, margin=0, charsize=1.0, $
         charthick=3.0, position=[pos[2,ii]-0.005,pos[3,ii]-0.03], /norm
       good = where(spec[ii].wave ne 0,npix)
       norm = interpol(spec[ii].flux[good]*scale,exp(spec[ii].wave[good]),4500)
       splog, norm, minmax(exp(spec[ii].wave[good]))
       djs_oplot, exp(spec[ii].wave[good]), spec[ii].flux[good]*scale/norm, $
         color=cgcolor('grey'), psym=10, thick=3
       djs_oplot, exp(spec[ii].wave[good]), (spec[ii].continuum[good]+$
         spec[ii].smooth_continuum[good]+spec[ii].linefit[good])*scale/norm, $
         color=cgcolor('dodger blue'), thick=3, psym=10
       djs_oplot, exp(spec[ii].wave[good]), (spec[ii].continuum[good]+$
         spec[ii].smooth_continuum[good])*scale/norm, $
         color=cgcolor('firebrick'), thick=3, psym=10

; make an inset focused on
       pos2 = [pos[0,ii]+0.01,pos[1,ii]+0.25,pos[0,ii]+0.25,pos[3,ii]-0.01]
       djs_plot, [0], [0], /nodata, /noerase, xsty=1, ysty=1, position=pos2, $
         charsize=1.0, /norm, ytickname=replicate(' ',10), xrange=[3700,4150], $
         yrange=[0,1.2], xtickinterval=200
       djs_oplot, exp(spec[ii].wave[good]), spec[ii].flux[good]*scale/norm, $
         color=cgcolor('grey'), psym=10, thick=3
       djs_oplot, exp(spec[ii].wave[good]), (spec[ii].continuum[good]+$
         spec[ii].smooth_continuum[good]+spec[ii].linefit[good])*scale/norm, $
         color=cgcolor('dodger blue'), thick=3, psym=10
       djs_oplot, exp(spec[ii].wave[good]), (spec[ii].continuum[good]+$
         spec[ii].smooth_continuum[good])*scale/norm, $
         color=cgcolor('firebrick'), thick=3, psym=10
    endfor

    xyouts, pos[0,0]-0.08, pos[1,0], 'Relative Flux', align=0.5, orientation=90, /normal
    xyouts, pos[0,3], pos[1,3]-0.12, textoidl('Rest Wavelength (\AA)'), $
      align=0.5, /normal
    
    im_plotconfig, psfile=psfile, /psclose, /pdf
    
stop    
    
; --------------------------------------------------
; S/N vs magnitude
    psfile = paperpath+'snr_vs_mag.eps'
    im_plotconfig, 6, pos, psfile=psfile

    xrange = [0.2,100]
    yrange1 = [17.5,25]
    yrange2 = [-16,-24.5]

    xtitle = 'Continuum S/N (pixel^{-1})'
    ytitle1 = 'I_{tot} (AB mag)'
    ytitle2 = 'M_{V} (AB mag)'

; vs I-mag    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      xtickname=replicate(' ',10), xtitle='', ytitle=ytitle1, $
      xrange=xrange, yrange=yrange1, /xlog
    djs_oplot, allcluster.continuum_snr, -2.5*alog10(allcluster.maggies[3]), $
      psym=symcat(16), color=fsc_color('firebrick',101), symsize=0.7
    djs_oplot, allfield.continuum_snr, -2.5*alog10(allfield.maggies[3]), $
      psym=symcat(15), color=fsc_color('dodger blue',101), symsize=0.7
    djs_oplot, 3.0*[1,1], !y.crange, line=0, thick=4
    djs_oplot, 10^!x.crange, 23*[1,1], line=5, thick=6

; vs M_V
    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      xtitle=xtitle, ytitle=ytitle2, xrange=xrange, yrange=yrange2, /xlog
    djs_oplot, allcluster.continuum_snr, allcluster.ubvrijhk_absmag_00[2], $
      psym=symcat(16), color=fsc_color('firebrick',101), symsize=0.7
    djs_oplot, allfield.continuum_snr, allfield.ubvrijhk_absmag_00[2], $
      psym=symcat(15), color=fsc_color('dodger blue',101), symsize=0.7
    djs_oplot, 3.0*[1,1], !y.crange, line=0, thick=4
    djs_oplot, 10^!x.crange, -19.0*[1,1], line=5, thick=6

; plot absolute magnitude vs I-band magnitude
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
      xtitle=ytitle1, ytitle=ytitle2, xrange=yrange1, yrange=yrange2
    djs_oplot, -2.5*alog10(allcluster.maggies[3]), allcluster.ubvrijhk_absmag_00[2], $
      psym=symcat(16), color=fsc_color('firebrick',101), symsize=0.7
    djs_oplot, -2.5*alog10(allfield.maggies[3]), allfield.ubvrijhk_absmag_00[2], $
      psym=symcat(15), color=fsc_color('dodger blue',101), symsize=0.7
    djs_oplot, !x.crange, -19.0*[1,1], line=0, thick=4
    djs_oplot, 23*[1,1], !y.crange, line=0, thick=4

    im_plotconfig, /psclose, /pdf, psfile=psfile

stop    
    
; --------------------------------------------------
; S/N vs I-band magnitude and M_V
    psfile = paperpath+'snr_vs_mag.ps'
    im_plotconfig, 6, pos, psfile=psfile

    xrange = [0.2,100]
    yrange1 = [17.5,25]
    yrange2 = [-16,-24.5]

    xtitle = 'Continuum S/N (pixel^{-1})'
    ytitle1 = 'I_{tot} (AB mag)'
    ytitle2 = 'M_{V} (AB mag)'

; vs I-mag    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      xtickname=replicate(' ',10), xtitle='', ytitle=ytitle1, $
      xrange=xrange, yrange=yrange1, /xlog
    djs_oplot, allcluster.continuum_snr, -2.5*alog10(allcluster.maggies[3]), $
      psym=symcat(16), color=fsc_color('firebrick',101), symsize=0.7
    djs_oplot, allfield.continuum_snr, -2.5*alog10(allfield.maggies[3]), $
      psym=symcat(15), color=fsc_color('dodger blue',101), symsize=0.7
    djs_oplot, 3.0*[1,1], !y.crange, line=0, thick=4
    djs_oplot, 10^!x.crange, 23*[1,1], line=5, thick=6

; vs M_V
    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      xtitle=xtitle, ytitle=ytitle2, xrange=xrange, yrange=yrange2, /xlog
    djs_oplot, allcluster.continuum_snr, allcluster.ubvrijhk_absmag_00[2], $
      psym=symcat(16), color=fsc_color('firebrick',101), symsize=0.7
    djs_oplot, allfield.continuum_snr, allfield.ubvrijhk_absmag_00[2], $
      psym=symcat(15), color=fsc_color('dodger blue',101), symsize=0.7
    djs_oplot, 3.0*[1,1], !y.crange, line=0, thick=4
    djs_oplot, 10^!x.crange, -19.0*[1,1], line=5, thick=6

; plot absolute magnitude vs I-band magnitude
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
      xtitle=ytitle1, ytitle=ytitle2, xrange=yrange1, yrange=yrange2
    djs_oplot, -2.5*alog10(allcluster.maggies[3]), allcluster.ubvrijhk_absmag_00[2], $
      psym=symcat(16), color=fsc_color('firebrick',101), symsize=0.7
    djs_oplot, -2.5*alog10(allfield.maggies[3]), allfield.ubvrijhk_absmag_00[2], $
      psym=symcat(15), color=fsc_color('dodger blue',101), symsize=0.7
    djs_oplot, !x.crange, -19.0*[1,1], line=0, thick=4
    djs_oplot, 23*[1,1], !y.crange, line=0, thick=4

    im_plotconfig, /psclose, /pdf, psfile=psfile

; --------------------------------------------------
; D(4000) vs Hd_A

    yrange = [-5.5,10]
    xrange = [0.8,2.3]
    ytitle = 'H\delta_{A} (\AA)'
    xtitle = 'D_{n}(4000)'

; #########################
; field
    psfile = paperpath+'d4000_hda_field.ps'
    im_plotconfig, 0, pos, psfile=psfile
; raw
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle=xtitle+ '(raw)', ytitle=ytitle+ '(raw)', xrange=xrange, yrange=yrange
    oplot_hda_d4000, field.lick_hd_a_cor[0], field.d4000_narrow_cor[0], $
      field.lick_hd_a_cor[1], field.d4000_narrow_cor[1], field.oii_3727_ew[0]
    im_legend, 'Field', /right, /top, box=0, charsize=1.8
; model
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle=xtitle+' (BC03)', ytitle=ytitle+ ' (BC03)', xrange=xrange, yrange=yrange
    oplot_hda_d4000, field.lick_hd_a_model[0], field.d4000_narrow_model[0], $
      field.lick_hd_a_model[1]*0.0, field.d4000_narrow_model[1]*0.0, field.oii_3727_ew[0]
    im_legend, 'Field', /right, /top, box=0, charsize=1.8
    im_plotconfig, /psclose, /pdf, psfile=psfile

; #########################
; cluster
    psfile = paperpath+'d4000_hda_cluster.ps'
    im_plotconfig, 0, pos, psfile=psfile
; raw
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle=xtitle+' (raw)', ytitle=ytitle+' (raw)', xrange=xrange, yrange=yrange
    oplot_hda_d4000, cluster.lick_hd_a_cor[0], cluster.d4000_narrow_cor[0], $
      cluster.lick_hd_a_cor[1], cluster.d4000_narrow_cor[1], cluster.oii_3727_ew[0]
    im_legend, 'Cluster', /right, /top, box=0, charsize=1.8
    !p.multi = 0
; model
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle=xtitle+' (BC03)', ytitle=ytitle+' (BC03)', xrange=xrange, yrange=yrange
    oplot_hda_d4000, cluster.lick_hd_a_model[0], cluster.d4000_narrow_model[0], $
      cluster.lick_hd_a_model[1]*0.0, cluster.d4000_narrow_model[1]*0.0, cluster.oii_3727_ew[0]
    im_legend, 'Cluster', /right, /top, box=0, charsize=1.8
    im_plotconfig, /psclose, /pdf, psfile=psfile
    
return
end
    


;;; --------------------------------------------------
;;; test the emission-line corrections
;;
;;    psfile = paperpath+'emline_cor.ps'
;;    im_plotconfig, 6, pos, psfile=psfile, xmargin=[1.3,0.2]
;;
;;; ###############
;;; D(4000)
;;    xtitle1 = 'D_{n}(4000)_{raw}'
;;    ytitle1 = 'D_{n}(4000)_{BC03}'
;;    xtitle2 = 'EW(H\beta) (\AA)'
;;    ytitle2 = 'D_{n}(4000)_{raw}-D_{n}(4000)_{BC03}'
;;
;;    xrange1 = [0.8,2.4]
;;    yrange1 = xrange1
;;    xrange2 = [0.1,100]
;;    yrange2 = [-0.3,0.4]
;;    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
;;      xrange=xrange1, yrange=yrange1, xtitle=xtitle1, ytitle=ytitle1
;;
;;; cluster
;;    resid = cluster.d4000_narrow_model[0]-cluster.d4000_narrow_raw[0]
;;    resid_err = cluster.d4000_narrow_raw[1]
;;    djs_oplot, cluster.h_beta_ew[0], resid, psym=symcat(16), $
;;      color=fsc_color('firebrick',101), symsize=1.4
;;;   oploterror, cluster.h_beta_ew[0], resid, cluster.h_beta_ew[1], $
;;;     resid_err, psym=symcat(16), color=fsc_color('firebrick',101), $
;;;     errcolor=fsc_color('firebrick',101), errthick=!p.thick, $
;;;     symsize=1.4
;;; field
;;    resid = field.d4000_narrow_model[0]-field.d4000_narrow_raw[0]
;;    resid_err = field.d4000_narrow_raw[1]
;;    djs_oplot, field.h_beta_ew[0], resid, psym=symcat(14), $
;;      color=fsc_color('dodger blue',101), symsize=1.4
;;;   oploterror, field.h_beta_ew[0], resid, field.h_beta_ew[1], $
;;;     resid_err, psym=symcat(14), color=fsc_color('dodger blue',101), $
;;;     errcolor=fsc_color('dodger blue',101), errthick=!p.thick, $
;;;     symsize=1.4
;;    im_legend, ['Cluster','Field'], /top, /left, box=0, $
;;      charsize=1.6, color=['firebrick','dodger blue'], psym=[16,14], $
;;      symsize=[1.6,1.9]
;;
;;; ###############
;;; Hd_{A}
;;    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;;      xtitle='EW(H\beta) (\AA)', ytitle='H\delta_{A, BC03} - H\delta_{A, raw}', $
;;      xrange=[0.1,100], yrange=[-3,8], /xlog
;;; cluster
;;    resid = cluster.lick_hd_a_model[0]-cluster.lick_hd_a_raw[0]
;;    resid_err = cluster.lick_hd_a_raw[1]
;;    djs_oplot, cluster.h_beta_ew[0], resid, $
;;      psym=symcat(16), color=fsc_color('firebrick',101), $
;;      symsize=1.4
;;;   oploterror, cluster.h_beta_ew[0], resid, cluster.h_beta_ew[1], $
;;;     resid_err, psym=symcat(16), color=fsc_color('firebrick',101), $
;;;     errcolor=fsc_color('firebrick',101), errthick=!p.thick, $
;;;     symsize=1.4
;;; field
;;    resid = field.lick_hd_a_model[0]-field.lick_hd_a_raw[0]
;;    resid_err = field.lick_hd_a_raw[1]
;;    djs_oplot, field.h_beta_ew[0], resid, $
;;      psym=symcat(14), color=fsc_color('dodger blue',101), $
;;      symsize=1.4
;;;   oploterror, field.h_beta_ew[0], resid, field.h_beta_ew[1], $
;;;     resid_err, psym=symcat(14), color=fsc_color('dodger blue',101), $
;;;     errcolor=fsc_color('dodger blue',101), errthick=!p.thick, $
;;;     symsize=1.4
;;    im_legend, ['Cluster','Field'], /top, /left, box=0, $
;;      charsize=1.6, color=['firebrick','dodger blue'], psym=[16,14], $
;;      symsize=[1.6,1.9]
;;
;;    im_plotconfig, /psclose
;;
;;stop    
;;    
