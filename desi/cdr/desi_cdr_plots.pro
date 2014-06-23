pro desi_cdr_plots
; jm14jun16siena - build plots for the DESI/CDR section on ELGs

    catpath = deep2_path(/cat)
    cdrpath = getenv('IM_SVNREPOS')+'/desi/cdr/4_Targets/plots/'
    targpath = getenv('IM_PROJECTS_DIR')+'/desi/targeting/'

    magcut1 = 23.0
    magcut2 = 23.5

    allphot = mrdfits(catpath+'deep2.pcat_ext.fits.gz',1)
    keep = where(allphot.zquality ge 3 and allphot.g gt 0 and $
      allphot.r gt 0 and allphot.z gt 0 and allphot.gerr lt 1 and $
      allphot.rerr lt 1 and allphot.zerr lt 1 and $
      allphot.badflag eq 0 and allphot.pgal ge 1)
    allphot = deep2_get_ugriz(allphot[keep])

    phot = mrdfits(targpath+'deep2egs-photparent.fits.gz',1)
;   zcat = mrdfits(targpath+'deep2egs_parent_zcat.fits.gz',1)
    zcat = mrdfits(targpath+'deep2egs-zcatparent.Q34.fits.gz',1)

    ewoiisnr = zcat.oii_3727_2_ew[0]/zcat.oii_3727_2_ew[1]
    oiisnr = zcat.oii_3727_2_amp[0]/(zcat.oii_3727_2_amp[1]+$
      (zcat.oii_3727_2_amp[1] eq 0))*(zcat.oii_3727_2_amp[1] ne 0)
    
; --------------------------------------------------
; gr vs rz - stellar contamination
    stars = mrdfits(targpath+'deep2egs-photstars.fits.gz',1)
    stars = stars[where(stars.ugriz[2] lt magcut1,nstar)]
    
    loz = where(allphot.zhelio lt 0.6 and allphot.ugriz[2] lt magcut1,nloz)
    hiz = where(allphot.zhelio gt 0.6 and allphot.zhelio lt 1.2 and $
      allphot.ugriz[2] lt magcut1,nhiz)
    vhiz = where(allphot.zhelio gt 1.2 and allphot.zhelio lt 1.6 and $
      allphot.ugriz[2] lt magcut1,nvhiz)
    
    w1 = where(allphot.zhelio gt 0.6 and allphot.ugriz[2] lt magcut1,nw1)
    w2 = where(allphot.zhelio gt 1.2 and allphot.ugriz[2] lt magcut1,nw2)
    splog, 'Fraction of z>0.6 in grz box: ', n_elements(desi_get_hizelg($
      allphot[w1].ugriz,magcut=magcut1))/float(nw1)
    splog, 'Fraction of z>1.2 in grz box: ', n_elements(desi_get_hizelg($
      allphot[w2].ugriz,magcut=magcut1))/float(nw2)
    
    psfile = cdrpath+'deep2-elg-stars-grz.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
      xmargin=[1.5,0.4], charsize=1.7
    
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='r - z', ytitle='g - r', xrange=[-0.5,2.2], yrange=[-0.3,2]
    djs_oplot, allphot[loz].ugriz[2]-allphot[loz].ugriz[4], $
      allphot[loz].ugriz[1]-allphot[loz].ugriz[2], psym=symcat(16), symsize=0.3
    djs_oplot, allphot[hiz].ugriz[2]-allphot[hiz].ugriz[4], $
      allphot[hiz].ugriz[1]-allphot[hiz].ugriz[2], psym=symcat(9,thick=2), $
      color=cgcolor('firebrick'), symsize=0.3
    djs_oplot, allphot[vhiz].ugriz[2]-allphot[vhiz].ugriz[4], $
      allphot[vhiz].ugriz[1]-allphot[vhiz].ugriz[2], psym=symcat(6), $
      color=cgcolor('forest green'), symsize=0.5
    djs_oplot, stars.ugriz[2]-stars.ugriz[4], symsize=0.15, $
      stars.ugriz[1]-stars.ugriz[2], psym=symcat(7), color='blue'
       
;; Mostek       
;    rzaxis = range(0.2,1.2,500)
;    djs_oplot, rzaxis, poly(rzaxis,[-0.08,0.68]), line=5, thick=6
;    djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[-0.08,0.68])], line=5, thick=6
;    djs_oplot, 1.3*[1,1], [!y.crange[0],poly(1.3,[-0.08,0.68])], line=5, thick=6
;    
;; proposed
;    rzaxis = range(0.1,1.1,500)
;    int = -0.08 & slope = 1.0
;    djs_oplot, rzaxis, poly(rzaxis,[int,slope]), line=0, thick=6
;    djs_oplot, 0.1*[1,1], [!y.crange[0],poly(0.1,[int,slope])], line=0, thick=6
;    djs_oplot, [1.1,1.8], poly(1.1,[int,slope])*[1,1], line=0, thick=6
;    djs_oplot, 1.8*[1,1], [!y.crange[0],poly(1.1,[int,slope])], line=0, thick=6
    
    im_legend, ['Stars','z<0.6','0.6<z<1.2','1.2<z<1.6'], $
      /left, /top, box=0, psym=[7,16,9,6], position=[0.22,0.86], /norm, $
      color=['blue','','firebrick','forest green'], charsize=1.5
    im_legend, ['18.5<r<23'], spacing=2.0, /left, /top, box=0, margin=0
    
;   im_legend, ['Mostek','Proposed'], /right, /bottom, box=0, $
;     line=[0,5], pspacing=1.7, thick=6
    im_plotconfig, psfile=psfile, /psclose, /pdf


    
stop       

; --------------------------------------------------
; Figure 3.4 - example ELG spectrum
    version = desi_deep2_template_version()
    templatepath = getenv('IM_PROJECTS_DIR')+'/desi/templates/'+version+'/'

    templatefile = templatepath+'desi_deep2elg_templates_'+version+'.fits.gz'
    templates = mrdfits(templatefile,0,hdr)
    wave = 10D^make_wave(hdr)
    info = mrdfits(templatefile,1)

;   ww = where(info.z gt 0.99 and info.z lt 1.01 and info.oii_3727 gt 1D-16 and $
;     info.sigma_kms gt 50 and info.sigma_kms lt 85,nww)
;   for ii = 0, nww-1 do begin
;      djs_plot, wave, templates[*,ww[ii]], xsty=3, ysty=3, $ ; /ylog, $
;        psym=-8, xrange=[3707,3747]      ; , xrange=[2000,10000]
;      splog, info[ww[ii]].z, info[ww[ii]].oii_3727, info[ww[ii]].oii_3727_ew, info[ww[ii]].logmstar, $
;        info[ww[ii]].objno, info[ww[ii]].sigma_kms
;      cc = get_kbrd(1)
;   endfor

    this = where(info.objno eq 13057544L)
    flux = templates[*,this]
    flux = flux/interpol(flux,wave,5500)
;   desi_quicksim, wave, flux, model='elg', simdat=simdat

;   im_lineid_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;     xrange=[1500,1E4], yrange=[0.3,200], xtitle='Rest-Frame Wavelength (\AA)', $
;     ytitle='Relative Flux', /ylog
    psfile = cdrpath+'ELG-deep2-example.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=4.5, width=6.7, $
      xmargin=[1.3,0.5], thick=4
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=[1500,1E4], yrange=[0.35,200], xtitle='Rest-Frame Wavelength (\AA)', $
      ytitle='Relative Flux', /ylog
    djs_oplot, wave, flux, psym=10, color=cgcolor('firebrick'), thick=3

;   zrange = where(wave*(1+info[this].z) gt 6550 and wave*(1+info[this].z) lt 9800)
;   djs_oplot, wave[zrange], flux[zrange], psym=10, color=cgcolor('dodger blue'), thick=3
    
    xyouts, 3710, 115, textoidl('[OII]'), /data, align=0.5, charsize=1.3
    xyouts, 4890, 45, textoidl('H\beta+[OIII]'), /data, align=0.5, charsize=1.3
    xyouts, 6560, 90, textoidl('H\alpha+[NII]'), /data, align=0.5, charsize=1.3
    
    xr = 3727+10*[-1,1] & get_element, wave, xr, ww
    djs_plot, [0], [0], /nodata, /noerase, position=[0.21,0.67,0.35,0.9], $ ; /noerase, 
      xsty=9, ysty=9, xrange=xr, yrange=[0,1.1], /norm, xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), xticks=0, yticks=0
    djs_oplot, wave, flux/max(flux[ww[0]:ww[1]]), psym=10, color=cgcolor('dodger blue')
    
    im_plotconfig, psfile=psfile, /psclose, /pdf
    
stop

; --------------------------------------------------
; Figure 3.6 - grz color-color plot coded by redshift

    gr = allphot.ugriz[1]-allphot.ugriz[2]
    rz = allphot.ugriz[2]-allphot.ugriz[4]

    grrange = [-0.5,1.8]
    rzrange = [-0.5,2.2]
    
    grzhist = hist_nd(transpose([[rz],[gr]]),0.1,min=[rzrange[0],grrange[0]],$
      max=[rzrange[1],grrange[1]],reverse_ind=ri)
    sgrzhist = smooth(grzhist,3)
    zmap = grzhist*0.0;-1
    for ii = 0, n_elements(grzhist)-1 do begin
       	if ri[ii] ne ri[ii+1] then begin
           if sgrzhist[ii] gt 3 then $
             zmap[ii] = djs_median(allphot[ri[ri[ii]:ri[ii+1]-1]].zhelio)
        endif
    endfor

    psfile = cdrpath+'grz-zweight.ps'
    cgPS_Open, psfile
     
    cgLoadCT, 16
    TVLCT, cgColor('white', /Triple), 0
    TVLCT, r, g, b, /Get
    palette = [ [r], [g], [b] ]

    zmax = 1.5
    zcrap = where(zmap lt 0,comp=zgood)
    zimage = bytscl(zmap,min=0.0,max=zmax)
;   zimage[zcrap] = cgcolor('white')
    
    cgImage, zimage, XRange=rzrange, YRange=grrange, $
      /Axes, Palette=palette, XTitle='r - z', YTitle='g -r', $
      Position=[0.125, 0.125, 0.9, 0.8]
;   cgContour, grzhist, LEVELS=max(sgrzhist)*[0.25,0.5,0.75,0.9,0.95], /OnImage, $
;     C_Colors=['Tan','Tan', 'Brown'], $
;     C_Annotation=['0.25', '0.5', '0.75'], $
;     C_Thick=thick, C_CharThick=thick
    cgColorbar, Position=[0.125, 0.875, 0.9, 0.925], Title='Mean Redshift', $
      Range=[0,zmax], NColors=254, Bottom=1, $ ; OOB_Low='white', $
      TLocation='Top'
    cgloadct, 0

; Mostek       
    rzaxis = range(0.2,1.3,500)
    djs_oplot, rzaxis, poly(rzaxis,[-0.08,0.68]), line=0, thick=6
    djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[-0.08,0.68])], line=0, thick=6
    djs_oplot, 1.3*[1,1], [!y.crange[0],poly(1.3,[-0.08,0.68])], line=0, thick=6

    cgPS_Close
    
    
stop    
    
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='r - z', ytitle='g - r', xrange=[-0.5,2.2], yrange=[-0.3,1.8]
    djs_oplot, zcat[hiz].ugriz[2]-zcat[hiz].ugriz[4], $
      zcat[hiz].ugriz[1]-zcat[hiz].ugriz[2], psym=symcat(5), $
      color=cgcolor('firebrick'), symsize=0.2


    
; --------------------------------------------------
; Figure 3.10 - [OII] flux vs redshift from DEEP2


return
end
    

    


    
;;; --------------------------------------------------
;;; gr vs rz - stellar contamination
;;    stars = mrdfits(targpath+'deep2egs-photstars.fits.gz',1)
;;    stars = stars[where(stars.ugriz[2] lt magcut1,nstar)]
;;    
;;    loz = where(zcat.zbest lt 0.6 and zcat.ugriz[2] lt magcut1,nloz)
;;    hiz = where(zcat.zbest gt 0.6 and zcat.zbest lt 1.2 and $
;;      zcat.ugriz[2] lt magcut1,nhiz)
;;    vhiz = where(zcat.zbest gt 1.2 and zcat.zbest lt 1.6 and $
;;      zcat.ugriz[2] lt magcut1,nvhiz)
;;    
;;    w1 = where(zcat.zbest gt 0.6 and zcat.ugriz[2] lt magcut1,nw1)
;;    w2 = where(zcat.zbest gt 1.2 and zcat.ugriz[2] lt magcut1,nw2)
;;    splog, 'Fraction of z>0.6 in grz box: ', n_elements(desi_get_hizelg($
;;      zcat[w1].ugriz,magcut=magcut1))/float(nw1)
;;    splog, 'Fraction of z>1.2 in grz box: ', n_elements(desi_get_hizelg($
;;      zcat[w2].ugriz,magcut=magcut1))/float(nw2)
;;    
;;    psfile = cdrpath+'deep2-elg-stars-grz.ps'
;;    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
;;      xmargin=[1.5,0.4], charsize=1.7
;;    
;;    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;;      xtitle='r - z', ytitle='g - r', xrange=[-0.5,2.2], yrange=[-0.3,1.8]
;;    djs_oplot, zcat[loz].ugriz[2]-zcat[loz].ugriz[4], $
;;      zcat[loz].ugriz[1]-zcat[loz].ugriz[2], psym=symcat(16), symsize=0.3
;;    djs_oplot, zcat[hiz].ugriz[2]-zcat[hiz].ugriz[4], $
;;      zcat[hiz].ugriz[1]-zcat[hiz].ugriz[2], psym=symcat(5), $
;;      color=cgcolor('firebrick'), symsize=0.2
;;    djs_oplot, zcat[vhiz].ugriz[2]-zcat[vhiz].ugriz[4], $
;;      zcat[vhiz].ugriz[1]-zcat[vhiz].ugriz[2], psym=symcat(6), $
;;      color=cgcolor('forest green'), symsize=0.5
;;    djs_oplot, stars.ugriz[2]-stars.ugriz[4], symsize=0.1, $
;;      stars.ugriz[1]-stars.ugriz[2], psym=symcat(7), color='blue'
;;       
;;;; Mostek       
;;;    rzaxis = range(0.2,1.2,500)
;;;    djs_oplot, rzaxis, poly(rzaxis,[-0.08,0.68]), line=5, thick=6
;;;    djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[-0.08,0.68])], line=5, thick=6
;;;    djs_oplot, 1.3*[1,1], [!y.crange[0],poly(1.3,[-0.08,0.68])], line=5, thick=6
;;;    
;;;; proposed
;;;    rzaxis = range(0.1,1.1,500)
;;;    int = -0.08 & slope = 1.0
;;;    djs_oplot, rzaxis, poly(rzaxis,[int,slope]), line=0, thick=6
;;;    djs_oplot, 0.1*[1,1], [!y.crange[0],poly(0.1,[int,slope])], line=0, thick=6
;;;    djs_oplot, [1.1,1.8], poly(1.1,[int,slope])*[1,1], line=0, thick=6
;;;    djs_oplot, 1.8*[1,1], [!y.crange[0],poly(1.1,[int,slope])], line=0, thick=6
;;    
;;    im_legend, ['Stars','z<0.6','0.6<z<1.2','1.2<z<1.6'], $
;;      /left, /top, box=0, psym=[7,16,5,6], position=[0.22,0.86], /norm, $
;;      color=['blue','','firebrick','forest green'], charsize=1.5
;;    im_legend, ['18.5<r<23'], spacing=2.0, /left, /top, box=0, margin=0
;;    
;;;   im_legend, ['Mostek','Proposed'], /right, /bottom, box=0, $
;;;     line=[0,5], pspacing=1.7, thick=6
;;    im_plotconfig, psfile=psfile, /psclose, /pdf
