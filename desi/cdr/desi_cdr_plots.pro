function get_dndm, rmag, weight=weight, faintcut=faintcut, $
  brightcut=brightcut, magaxis=magaxis
; get the number counts    
    area = 0.4342 ; deg^2
    binsize = 0.2
    nbins = ceil((faintcut-brightcut)/binsize)
    magaxis = lindgen(nbins)*binsize+brightcut+binsize/2.0
    dndm = hogg_histogram(rmag,[brightcut,faintcut],$ ; #/0.2 mag/deg^2
      nbins,weight=weight)/area
;   dndm = im_hist1d(rmag,weight,min=brightcut,$ 
;     max=faintcut,binsize=binsize)/area/binsize
return, dndm
end

pro desi_cdr_plots
; jm14jun16siena - build plots for the DESI/CDR section on ELGs

    catpath = deep2_path(/cat)
    cdrpath = getenv('IM_SVNREPOS')+'/desi/cdr/4_Targets/plots/'
    targpath = getenv('IM_PROJECTS_DIR')+'/desi/targeting/'

    brightcut = 18.5
    faintcut = 24.0

;   magcut2 = 23.0
    oiisnrcut = 1.0
    oiicut1 = 8D-17 ; [erg/s/cm2]
    area = 0.4342 ; deg^2

    phot = mrdfits(targpath+'deep2egs-photparent.fits.gz',1)
    zcat = mrdfits(targpath+'deep2egs-zcatparent.Q34.fits.gz',1)

;   oiisnr = zcat.oii_3727_2_amp[0]/(zcat.oii_3727_2_amp[1]+$
;     (zcat.oii_3727_2_amp[1] eq 0))*(zcat.oii_3727_2_amp[1] ne 0)

; -------------------------
; assign [OII] fluxes to galaxies with no [OII] measured, but only
; between z=0.8-1.45    
    zmin = 0.6                  ; 0.8
    zmax = 1.6

    refindx = where(zcat.zbest gt zmin and zcat.zbest lt zmax and $
      zcat.oii_3727[1] ne -2 and $
      zcat.oii_3727_2_ew[0]/zcat.oii_3727_2_ew[1] ge 1.0,nrefindx)
    rref = zcat[refindx].ugriz[2]
    rzref = zcat[refindx].ugriz[2]-zcat[refindx].ugriz[4]
    grref = zcat[refindx].ugriz[1]-zcat[refindx].ugriz[2]
    zref = zcat[refindx].zbest
    oiiref = zcat[refindx].oii_3727

    noneindx = where(zcat.zbest gt zmin and zcat.zbest lt zmax and $
      (zcat.oii_3727[1] eq -2 or $
      zcat.oii_3727_2_ew[0]/zcat.oii_3727_2_ew[1] lt 1.0),nnoneindx)
    rnone = zcat[noneindx].ugriz[2]
    rznone = zcat[noneindx].ugriz[2]-zcat[noneindx].ugriz[4]
    grnone = zcat[noneindx].ugriz[1]-zcat[noneindx].ugriz[2]
    znone = zcat[noneindx].zbest

    for ii = 0, nnoneindx-1 do begin
       dist = sqrt((rref-rnone[ii])^2+(rzref-rznone[ii])^2+$
         (grref-grnone[ii])^2+(zref-znone[ii])^2)
       mindist = min(dist,thisref)
;      print, rref[thisref], rnone[ii], rzref[thisref], rznone[ii], $
;        grref[thisref], grnone[ii], zref[thisref], znone[ii]

       zcat[noneindx[ii]].oii_3727 = zcat[refindx[thisref]].oii_3727           ; flux
       zcat[noneindx[ii]].oii_3727_2_ew = zcat[refindx[thisref]].oii_3727_2_ew ; EW
    endfor

; --------------------------------------------------
; gr vs rz coded by [OII] strength
    zmin = 0.6 ; 0.8
    zmax = 1.6
    magcut1 = 23.2
    
    all = where(zcat.ugriz[2] lt magcut1,nall)
    loz = where(zcat.zbest lt zmin and zcat.ugriz[2] lt magcut1,nloz)
    oiibright = where(zcat.zbest ge zmin and $ ; zcat.zbest lt zmax and $
      zcat.ugriz[2] lt magcut1 and zcat.oii_3727[1] ne -2.0 and $
      zcat.oii_3727_2_ew[0]/zcat.oii_3727_2_ew[1] gt 1.0 and zcat.oii_3727[0] gt oiicut1,noiibright)
    oiifaint = where(zcat.zbest ge zmin and $ ; zcat.zbest lt zmax and $
      zcat.ugriz[2] lt magcut1 and zcat.oii_3727[1] ne -2.0 and $
      zcat.oii_3727_2_ew[0]/zcat.oii_3727_2_ew[1] gt 1.0 and zcat.oii_3727[0] lt oiicut1,noiifaint)
    oiinone = where(zcat.zbest ge zmin and $ ; zcat.zbest lt zmax and $
      zcat.ugriz[2] lt magcut1 and (zcat.oii_3727[1] eq -2 or $
      zcat.oii_3727_2_ew[0]/zcat.oii_3727_2_ew[1] le 1.0),noiinone)

; for testing
    oiibright_vhiz = where(zcat.zbest ge 1.2 and $ ; zcat.zbest lt zmax and $
      zcat.ugriz[2] lt magcut1 and zcat.oii_3727[1] ne -2.0 and $
      zcat.oii_3727_2_ew[0]/zcat.oii_3727_2_ew[1] gt 1.0 and zcat.oii_3727[0] gt oiicut1)
    
    splog, nall, nloz, noiibright, noiifaint, noiinone, $
      nloz+noiibright+noiifaint+noiinone
    
    psfile = cdrpath+'deep2-elg-grz-oii.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
      xmargin=[1.5,0.4], charsize=1.7
    
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='r - z', ytitle='g - r', xrange=[-0.5,2.2], yrange=[-0.3,2]
    djs_oplot, zcat[loz].ugriz[2]-zcat[loz].ugriz[4], $
      zcat[loz].ugriz[1]-zcat[loz].ugriz[2], psym=symcat(16), symsize=0.3
    djs_oplot, zcat[oiifaint].ugriz[2]-zcat[oiifaint].ugriz[4], $
      zcat[oiifaint].ugriz[1]-zcat[oiifaint].ugriz[2], psym=symcat(6), $
      color=cgcolor('forest green'), symsize=0.5
    djs_oplot, zcat[oiibright].ugriz[2]-zcat[oiibright].ugriz[4], $
      zcat[oiibright].ugriz[1]-zcat[oiibright].ugriz[2], psym=symcat(5), $
      color=cgcolor('firebrick'), symsize=0.3
    djs_oplot, zcat[oiibright_vhiz].ugriz[2]-zcat[oiibright_vhiz].ugriz[4], $
      zcat[oiibright_vhiz].ugriz[1]-zcat[oiibright_vhiz].ugriz[2], psym=symcat(7), $
      color=cgcolor('blue'), symsize=0.3
;   djs_oplot, zcat[oiinone].ugriz[2]-zcat[oiinone].ugriz[4], $
;     zcat[oiinone].ugriz[1]-zcat[oiinone].ugriz[2], psym=symcat(7), $
;     color=cgcolor('blue'), symsize=0.3

    hiz_all = desi_get_hizelg(zcat.ugriz,magcut=magcut1)
    hiz_oiibright = desi_get_hizelg(zcat[oiibright].ugriz);,magcut=magcut1)
    hiz_oiifaint = desi_get_hizelg(zcat[oiifaint].ugriz);,magcut=magcut1)
;   hiz_oiinone = desi_get_hizelg(zcat[oiinone].ugriz);,magcut=magcut1)
    hiz_loz = desi_get_hizelg(zcat[loz].ugriz);,magcut=magcut1)
    splog, 'All ', n_elements(hiz_all), $
      ' Bright ', n_elements(hiz_oiibright), 1.0*n_elements(hiz_oiibright)/n_elements(hiz_all), $
      ' Faint ', n_elements(hiz_oiifaint), 1.0*n_elements(hiz_oiifaint)/n_elements(hiz_all), $
;     ' None ', 1.0*n_elements(hiz_oiinone)/n_elements(hiz_all), $
      ' Low-z ', n_elements(hiz_loz), 1.0*n_elements(hiz_loz)/n_elements(hiz_all)
;   help, hiz_all, hiz_oiibright, hiz_oiifaint, hiz_oiinone, loz
;   djs_oplot, zcat[hiz_all].ugriz[2]-zcat[hiz_all].ugriz[4], $
;     zcat[hiz_all].ugriz[1]-zcat[hiz_all].ugriz[2], psym=7, symsize=0.2, $
;     color='orange'
;    djs_oplot, zcat[loz[hiz_loz]].ugriz[2]-zcat[loz[hiz_loz]].ugriz[4], $
;      zcat[loz[hiz_loz]].ugriz[1]-zcat[loz[hiz_loz]].ugriz[2], psym=7, symsize=1.1, $
;      color='blue'
    
;; Mostek       
;    rzaxis = range(0.2,1.2,500)
;    djs_oplot, rzaxis, poly(rzaxis,[-0.08,0.68]), line=5, thick=6
;    djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[-0.08,0.68])], line=5, thick=6
;    djs_oplot, 1.3*[1,1], [!y.crange[0],poly(1.3,[-0.08,0.68])], line=5, thick=6
    
;; proposed
;    rzaxis = range(0.1,1.1,500)
;    int = 0.0 & slope = 0.9
;    djs_oplot, [!x.crange[0],0.3], 0.1*[1,1], line=0, thick=6
;    djs_oplot, [0.3,1.3], poly([0.3,1.3],[int,slope]), line=0, thick=6
;    djs_oplot, 1.3*[1,1], [!y.crange[0],poly(1.3,[int,slope])], line=0, thick=6

; proposed
    rzaxis1 = range(0.2,0.9,100)
;   rzaxis1 = range(0.2,1.4,100)
    rzaxis2 = range(0.9,1.4,100)
    int1 = -0.1 & slope1 = 1.0
    int2 = 1.7 & slope2 = -1.0
    djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[int1,slope1])], thick=6
;   djs_oplot, [!x.crange[0],0.2], poly(0.2,[int1,slope1])*[1,1], line=0, thick=6
;   djs_oplot, [-0.2,0.3], poly(0.3,[int1,slope1])*[1,1], line=0, thick=6
;   djs_oplot, 0.1*[1,1], [!y.crange[0],poly(0.1,[int1,slope1])], line=0, thick=6
    djs_oplot, rzaxis1, poly(rzaxis1,[int1,slope1]), line=0, thick=6
;   djs_oplot, rzaxis1, poly(rzaxis1,[int1,slope1])<poly(0.9,[int1,slope1]), line=0, thick=6
    djs_oplot, rzaxis2, poly(rzaxis2,[int2,slope2]), line=0, thick=6
;   djs_oplot, [0.9,!x.crange[1]], poly(0.9,[int1,slope1])*[1,1], line=0, thick=6
;   djs_oplot, [0.9,1.4], poly(0.9,[int1,slope1])*[1,1], line=0, thick=6
    djs_oplot, 1.4*[1,1], [!y.crange[0],poly(1.4,[int2,slope2])], line=0, thick=6

;   djs_oplot, zcat[oiibright[hiz]].ugriz[2]-zcat[oiibright[hiz]].ugriz[4], $
;     zcat[oiibright[hiz]].ugriz[1]-zcat[oiibright[hiz]].ugriz[2], psym=symcat(16), $
;     symsize=2, color='green'
    
;; Johan's proposed cuts    
;    rzaxis = range(0.2,1.25,500)
;    int = 0.1 & slope = 0.55/1.25
;    djs_oplot, rzaxis, slope*rzaxis+int, line=5, thick=6
;    djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[int,slope])], line=5, thick=6
;    djs_oplot, 1.25*[1,1], [!y.crange[0],poly(1.25,[int,slope])], line=5, thick=6
    
;   rmag-zmag>0.2
;   rmag-zmag<1.25
;   gmag-rmag<0.55*(rmag-zmag)/1.25+0.1
;   gmag-rmag>0

    magstr = string(magcut1,format='(F4.1)')
    im_legend, ['18.5<r<'+magstr], spacing=2.0, /left, /top, box=0, $
      position=[0.17,0.92], /norm, charsize=1.7

;   im_legend, ['0.75<z<1.45'], spacing=2.0, /left, /top, box=0, $
;     position=[0.18,0.855], /norm, charsize=1.4
;   im_legend, ['[OII]>8\times10^{-17}',$
;     '[OII]<8\times10^{-17}','[OII] unmeasured'], /left, /top, $
;     box=0, psym=[5,6,7], position=[0.21,0.8], /norm, $
;     color=['firebrick','forest green','blue'], charsize=1.2
;   im_oplot_box, 0.9, 0.5, 0.0, xoff=0.0, yoff=1.48

    zstr = string(zmin,format='(F3.1)')
    im_legend, ['z<'+zstr+'','z>'+zstr+'; [OII]>8\times10^{-17}',$
      'z>'+zstr+'; [OII]<8\times10^{-17}'], /left, /top, $
      box=0, psym=[16,5,6], position=[0.2,0.86], /norm, $
      color=['','firebrick','forest green'], charsize=1.5
;   im_legend, ['z<0.8','z>0.8; [OII]>8\times10^{-17}',$
;     'z>0.8; [OII]<8\times10^{-17}',$
;     'z>0.8; [OII] unmeasured'], $
;     /left, /top, $
;     box=0, psym=[16,5,6,7], position=[0.2,0.86], /norm, $
;     color=['','firebrick','forest green','blue'], charsize=1.3
;   im_legend, ['z<0.75','0.75<z<1.45; [OII]>8\times10^{-17}',$
;     '0.75<z<1.45; [OII]<8\times10^{-17}',$
;     '0.75<z<1.45; [OII] unmeasured'], /left, /top, $
;     box=0, psym=[16,5,6,7], position=[0.2,0.86], /norm, $
;     color=['','firebrick','forest green','blue'], charsize=1.4
    
    xyouts, 1.4, 0.35, 'In color box:', align=0.0, charsize=1.5, /data
    im_legend, [$
      'N='+strtrim(n_elements(hiz_loz),2)+' ('+$
      string(round(100.0*n_elements(hiz_loz)/n_elements(hiz_all)),format='(I0)')+'%)',$
      'N='+strtrim(n_elements(hiz_oiibright),2)+' ('+$
      string(round(100.0*n_elements(hiz_oiibright)/n_elements(hiz_all)),format='(I0)')+'%)',$
      'N='+strtrim(n_elements(hiz_oiifaint),2)+' ('+$
      string(round(100.0*n_elements(hiz_oiifaint)/n_elements(hiz_all)),format='(I0)')+'%)'], $
      /left, /bottom, box=0, psym=[16,5,6], position=[0.72,0.25], /norm, $
      color=['','firebrick','forest green'], charsize=1.4

    im_plotconfig, psfile=psfile, /psclose, /pdf

; --------------------------------------------------
; cumulative number counts

; assume that the objects with formal flux limits above our [OII] cut
; are detections (this is a small number of objects)    
    hiz = desi_get_hizelg(zcat.ugriz)
    oiibright = where(zcat[hiz].oii_3727[1] ne -2.0 and $ ; oiisnr gt oiisnrcut and $
      zcat[hiz].oii_3727_2_ew[0]/zcat[hiz].oii_3727_2_ew[1] gt 1.0 and $
      zcat[hiz].oii_3727[0] gt oiicut1,noiibright)
    oiifaint = where(zcat[hiz].oii_3727[1] ne -2.0 and $ ; oiisnr gt oiisnrcut and $
      zcat[hiz].oii_3727_2_ew[0]/zcat[hiz].oii_3727_2_ew[1] gt 1.0 and $
      zcat[hiz].oii_3727[0] lt oiicut1,noiifaint)
    oiinone = where(zcat[hiz].oii_3727[1] eq -2 or $
      zcat[hiz].oii_3727_2_ew[0]/zcat[hiz].oii_3727_2_ew[1] le 1.0,noiinone)

    dndm_phot = get_dndm(phot.ugriz[2],faintcut=faintcut,$
      brightcut=brightcut,magaxis=magaxis)
    dndm = get_dndm(zcat.ugriz[2],weight=zcat.final_weight,$
      faintcut=faintcut,brightcut=brightcut,magaxis=magaxis)
    dndm_hiz = get_dndm(zcat[hiz].ugriz[2],weight=zcat[hiz].final_weight,$
      faintcut=faintcut,brightcut=brightcut)
    dndm_oiibright = get_dndm(zcat[hiz[oiibright]].ugriz[2],$
      weight=zcat[hiz[oiibright]].final_weight,$
      faintcut=faintcut,brightcut=brightcut)
    
; get the *ratio* of the number of bright-to-faint [OII] sources so
; that we can correct for the objects with missing [OII] measurements
; (see plot below)     
    dndm_oiifaint = get_dndm(zcat[hiz[oiifaint]].ugriz[2],$
      weight=zcat[hiz[oiifaint]].final_weight,$
      faintcut=faintcut,brightcut=brightcut)
    dndm_oiinone = get_dndm(zcat[hiz[oiinone]].ugriz[2],$
      weight=zcat[hiz[oiinone]].final_weight,$
      faintcut=faintcut,brightcut=brightcut)

;   denom = dndm_oiibright+dndm_oiifaint
;   dndm_oiicor = dndm_oiinone*dndm_oiibright/(denom+(denom eq 0))*(denom ne 0)
;   dndm_oiibright_cor = dndm_oiibright + dndm_oiicor    
    
    psfile = cdrpath+'deep2-elg-dndm.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
      xmargin=[1.5,0.4], charsize=1.7
    
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='r (AB mag)', ytitle='log N(<r) (gal / deg^{2})', $
      xrange=[brightcut-0.5,faintcut+0.5], yrange=alog10([5,2E5]) ;, /ylog
    djs_oplot, magaxis, alog10(total(dndm_phot,/cumu)), line=0, thick=6, color='grey'
;     line=2, thick=6, color=cgcolor('dark grey')
;   djs_oplot, magaxis, alog10(total(dndm,/cumu)), line=0, thick=6
    
    ww = where(total(dndm_hiz,/cumu) gt 0)
    djs_oplot, magaxis[ww], alog10((total(dndm_hiz,/cumu))[ww]), line=0, thick=8
    
    ww = where(total(dndm_oiibright,/cumu) gt 0)
    djs_oplot, magaxis[ww], alog10((total(dndm_oiibright,/cumu))[ww]), color=cgcolor('firebrick'), $
      line=5, thick=8           ; psym=symcat(6,thick=4), symsize=1.3

    ww = where(total(dndm_oiifaint,/cumu) gt 0)
    djs_oplot, magaxis[ww], alog10((total(dndm_oiifaint,/cumu))[ww]), color=cgcolor('forest green'), $
      line=4, thick=8           ; psym=symcat(6,thick=4), symsize=1.3

;   ww = where(total(dndm_oiinone,/cumu) gt 0)
;   djs_oplot, magaxis[ww], alog10((total(dndm_oiinone,/cumu))[ww]), color=cgcolor('dark grey'), $
;     line=1, thick=6           ; psym=symcat(6,thick=4), symsize=1.3

;    ww = where(total(dndm_oiibright_cor,/cumu) gt 0)
;    djs_oplot, magaxis[ww], alog10((total(dndm_oiibright_cor,/cumu))[ww]), $
;      color=cgcolor('red'), line=5, thick=8 ; psym=symcat(6,thick=4), symsize=1.3
    
    numcut = 3000
    rcut = interpol(magaxis,total(dndm_hiz,/cumu),numcut)
;   numcut = 2400
;   rcut = interpol(magaxis,total(dndm_oiibright_cor,/cumu),numcut)
;   rcut = 22.7
    oiinumcut = long(interpol(total(dndm_oiibright,/cumu),magaxis,rcut))
    splog, rcut, numcut, oiinumcut, 1.0*oiinumcut/numcut
;   oiinumcut = long(interpol(total(dndm_oiibright_cor,/cumu),magaxis,rcut))
    djs_oplot, [!x.crange[0],rcut], alog10(numcut)*[1,1], line=1, thick=6
    djs_oplot, [!x.crange[0],rcut], alog10(oiinumcut)*[1,1], line=1, thick=6
    djs_oplot, rcut*[1,1], [!y.crange[0],alog10(numcut)], line=1, thick=6
    
    xyouts, 18.3, alog10(numcut*1.2), textoidl(string(numcut,$
      format='(I0)')+' gal/deg^{2}'), align=0.0, /data, charsize=1.5
    xyouts, 18.3, alog10(oiinumcut*0.6), textoidl(string(oiinumcut,$
      format='(I0)')+' gal/deg^{2}'), align=0.0, /data, charsize=1.5
    xyouts, rcut+0.1, 2.3, textoidl('r='+string(rcut,$
      format='(F4.1)')), align=0.0, orientation=270, /data, charsize=1.5
    
    im_legend, ['All Galaxies'], /right, /top, box=0, thick=6, $
      charsize=1.3, spacing=2.0, margin=0, line=0, pspacing=1.5, color='grey'
    xyouts, 0.21, 0.88, 'In grz color box:', /norm, align=0.0, charsize=1.3
    im_legend, ['All galaxies','[OII]>8\times10^{-17}','[OII]<8\times10^{-17}'], $
      color=['black','firebrick','forest green'], /left, /top, box=0, thick=6, $
      charsize=1.3, spacing=2.0, margin=0, line=[0,5,4], pspacing=1.5, $
      position=[0.22,0.85], /norm
;   im_legend, ['PhotParent','zcatParent','grz Color Box',$
;     'grz & F([OII])>8\times10^{-17}'], $
;     color=['dark grey','','orange','blue'], /left, /top, box=0, thick=6, $
;     charsize=1.3, spacing=2.0, margin=0, line=[2,0,3,5], pspacing=1.5
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

; --------------------------------------------------
; redshift histogram of sources selected using my grz color-cuts;
; remove the 
    magcut1 = 22.84
    
    hiz = desi_get_hizelg(zcat.ugriz,magcut=magcut1,sigma_kms=zcat.sigma_kms)
    oiibright = where(zcat[hiz].oii_3727[1] ne -2.0 and $ ; oiisnr gt oiisnrcut and $
      zcat[hiz].oii_3727_2_ew[0]/zcat[hiz].oii_3727_2_ew[1] gt 1.0 and $
      zcat[hiz].oii_3727[0] gt oiicut1,noiibright)
    oiifaint = where(zcat[hiz].oii_3727[1] ne -2.0 and $ ; oiisnr gt oiisnrcut and $
      zcat[hiz].oii_3727_2_ew[0]/zcat[hiz].oii_3727_2_ew[1] gt 1.0 and $
      zcat[hiz].oii_3727[0] lt oiicut1,noiifaint)
    oiinone = where(zcat[hiz].oii_3727[1] eq -2 or $
      zcat[hiz].oii_3727_2_ew[0]/zcat[hiz].oii_3727_2_ew[1] le 1.0,noiinone)
    
    zmin = 0.0
    zmax = 2.0
    zbin = 0.1
    nzbin = ceil((zmax-zmin)/zbin)
    zhist = lindgen(nzbin)*zbin+zmin+zbin/2.0
    nz = hogg_histogram(zcat.zbest,[zmin,zmax],nzbin,$
      weight=zcat.final_weight)/area
    nz_hiz = hogg_histogram(zcat[hiz].zbest,[zmin,zmax],$
      nzbin,weight=zcat[hiz].final_weight)/area
    nz_oiibright = hogg_histogram(zcat[hiz[oiibright]].zbest,[zmin,zmax],$
      nzbin,weight=zcat[hiz[oiibright]].final_weight)/area
    nz_oiifaint = hogg_histogram(zcat[hiz[oiifaint]].zbest,[zmin,zmax],$
      nzbin,weight=zcat[hiz[oiifaint]].final_weight)/area

; extrapolate linearly
    anchor = where(zhist gt 1.1 and zhist lt 1.4)
    extrap = where(zhist gt 1.4 and zhist lt 1.8)

;   plot, zhist, alog10(nz_hiz), psym=8, xr=[1.0,1.6], symsize=3
    nz_hiz_extrap = nz_hiz
;   nz_hiz_extrap[extrap] = interpol(nz_hiz[anchor],$
;     zhist[anchor],zhist[extrap])>nz_hiz[extrap]
    nz_hiz_extrap[extrap] = 10^poly(zhist[extrap],$
      linfit(zhist[anchor],alog10(nz_hiz[anchor])))
    niceprint, zhist, nz_hiz_extrap, nz_hiz & print

;   plot, zhist, alog10(nz_oiibright), psym=8, xr=[1.0,1.6], symsize=3
    nz_oiibright_extrap = nz_oiibright
;   nz_oiibright_extrap[extrap] = interpol(nz_oiibright[anchor],$
;     zhist[anchor],zhist[extrap])>nz_oiibright[extrap]
    nz_oiibright_extrap[extrap] = 10^poly(zhist[extrap],$
      linfit(zhist[anchor],alog10(nz_oiibright[anchor])))
    niceprint, zhist, nz_oiibright_extrap, nz_oiibright & print
    
    nz_oiifaint_extrap = nz_oiifaint
    nz_oiifaint_extrap[extrap] = interpol(nz_oiifaint[anchor],$
      zhist[anchor],zhist[extrap])>nz_oiifaint[extrap]
;   niceprint, zhist, nz_oiifaint_extrap, nz_oiifaint

    splog, total(nz_hiz_extrap), total(nz_oiibright_extrap)
    
;; correct for the missing [OII] sources       
;   nz_oiinone = hogg_histogram(zcat[hiz[oiinone]].zbest,[zmin,zmax],$
;     nzbin,weight=zcat[hiz[oiinone]].final_weight)/area
    
;   denom = nz_oiibright+nz_oiifaint
;   nz_oiicor = nz_oiinone*nz_oiibright/(denom+(denom eq 0))*(denom ne 0)
;   nz_oiibright_cor = nz_oiibright + nz_oiicor    
    
    psfile = cdrpath+'deep2-elg-dndz.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
      xmargin=[1.5,0.4], charsize=1.7
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='Redshift', ytitle='dN/dz (gal / 0.1dz / deg^{2})', $
      xrange=[zmin,zmax], yrange=[0.0,1100]
;     xrange=[zmin,zmax], yrange=[0.0,6500]
;   djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
;     [0,nz,0], psym=10, thick=6, line=0
    djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
      [0,nz_hiz_extrap,0], psym=10, thick=8, line=0, color='black'
    djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
      [0,nz_oiibright_extrap,0], psym=10, thick=6, line=5, color=cgcolor('firebrick')
;   djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
;     [0,nz_oiibright_cor,0], psym=10, thick=8, line=5, color='blue'
    djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
      [0,nz_oiifaint_extrap,0], psym=10, thick=6, line=4, color=cgcolor('forest green')

    magstr = string(magcut1,format='(F4.1)')
    im_legend, ['18.5<r<'+magstr], spacing=2.0, /left, /top, box=0, $
      position=[0.17,0.92], /norm, charsize=1.7

    xyouts, 0.63, 0.88, 'In grz color box:', /norm, align=0.0, charsize=1.3
    im_legend, ['All galaxies','[OII]>8\times10^{-17}','[OII]<8\times10^{-17}'], $
      color=['black','firebrick','forest green'], /left, /top, box=0, thick=6, $
      charsize=1.3, spacing=2.0, margin=0, line=[0,5,4], pspacing=1.5, $
      position=[0.65,0.85], /norm

;   im_legend, ['zcatParent','grz Sample','grz & F([OII])>8\times10^{-17}'], $
;     line=[0,3,5], pspacing=1.7, color=['','orange','blue'], /right, /top, box=0, $
;     charsize=1.4, spacing=2.0, margin=0, thick=6
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop       
       
; --------------------------------------------------
; gr vs rz - stellar contamination
    allphot = mrdfits(catpath+'deep2.pcat_ext.fits.gz',1)
    keep = where(allphot.zquality ge 3 and allphot.g gt 0 and $
      allphot.r gt 0 and allphot.z gt 0 and allphot.gerr lt 1 and $
      allphot.rerr lt 1 and allphot.zerr lt 1 and $
      allphot.badflag eq 0 and allphot.pgal ge 1)
    allphot = deep2_get_ugriz(allphot[keep])

    stars = mrdfits(targpath+'deep2egs-photstars.fits.gz',1)
    stars = stars[where(stars.ugriz[2] lt magcut1,nstar)]
    
    loz = where(allphot.zhelio lt 0.6 and allphot.ugriz[2] lt magcut1,nloz)
    hiz = where(allphot.zhelio gt 0.6 and allphot.zhelio lt 1.2 and $
      allphot.ugriz[2] lt magcut1,nhiz)
    vhiz = where(allphot.zhelio gt 1.2 and allphot.zhelio lt 1.5 and $
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
    
    im_legend, ['Stars','z<0.6','0.6<z<1.2','1.2<z<1.5'], $
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

;;    rmin = 19.5
;;    rmax = 24.0
;;    rbin = 0.2
;;    nrbins = long((rmax-rmin)/rbin+1)
;;    raxis = lindgen(nrbins)*rbin+rmin+rbin/2.0
;;    
;;    rzmin = -0.5
;;    rzmax = 2.2
;;    rzbin = 0.15
;;    nrzbins = long((rzmax-rzmin)/rzbin+1)
;;    rzaxis = lindgen(nrzbins)*rzbin+rzmin+rzbin/2.0
;;
;;    grmin = -0.5
;;    grmax = 2.2
;;    grbin = 0.15
;;    ngrbins = long((grmax-grmin)/grbin+1)
;;    graxis = lindgen(ngrbins)*grbin+grmin+grbin/2.0
;;
;;    zmin = 0.6                  ; 0.8
;;    zmax = 1.6
;;    zbin = 0.05
;;    nzbins = long((zmax-zmin)/zbin+1)
;;    zaxis = lindgen(nzbins)*zbin+zmin+zbin/2.0
;;    
;;    ewoiisnr = zcat.oii_3727_2_ew[0]/zcat.oii_3727_2_ew[1]
;;    refindx = where(zcat.zbest gt zmin and zcat.zbest lt zmax and $
;;      zcat.oii_3727[1] ne -2 and ewoiisnr ge 1.0,nrefindx)
;;    rref = zcat[refindx].ugriz[2]
;;    rzref = zcat[refindx].ugriz[2]-zcat[refindx].ugriz[4]
;;    grref = zcat[refindx].ugriz[1]-zcat[refindx].ugriz[2]
;;    zref = zcat[refindx].zbest
;;    oiiref = zcat[refindx].oii_3727
;;
;;    noneindx = where(zcat.zbest gt zmin and zcat.zbest lt zmax and $
;;      (zcat.oii_3727[1] eq -2 or ewoiisnr lt 1.0),nnoneindx)
;;    rnone = zcat[noneindx].ugriz[2]
;;    rznone = zcat[noneindx].ugriz[2]-zcat[noneindx].ugriz[4]
;;    grnone = zcat[noneindx].ugriz[1]-zcat[noneindx].ugriz[2]
;;    znone = zcat[noneindx].zbest
;;    
;;    histref = hist_nd(transpose([$
;;      [rref],$ ; r-band first
;;      [(rzref>rzmin)<rzmax],$ ; r-z
;;      [(grref>grmin)<grmax],$ ; g-r
;;      [zref]]),$              ; redshift
;;      [rbin,rzbin,grbin,zbin],min=[rmin,rzmin,grmin,zmin],$
;;      max=[rmax,rzmax,grmax,zmax],reverse_ind=riref)
;;    histnone = hist_nd(transpose([$
;;      [rnone],$ ; r-band first
;;      [(rznone>rzmin)<rzmax],$ ; r-z
;;      [(grnone>grmin)<grmax],$ ; g-r
;;      [znone]]),$              ; redshift
;;      [rbin,rzbin,grbin,zbin],min=[rmin,rzmin,grmin,zmin],$
;;      max=[rmax,rzmax,grmax,zmax],reverse_ind=rinone)
;;
;;;    for ii = 0, n_elements(histref)-1 do begin
;;;       if riref[ii] ne riref[ii+1] then begin
;;;print, ii & stop
;;;       endif
;;;    endfor
;;
;;;    for ii = 0, n_elements(histnone)-1 do begin
;;;       indx = array_indices(histnone,ii)
;;;       print, raxis[indx[0]], rzaxis[indx[1]], $
;;;         graxis[indx[2]], zaxis[indx[3]]
;;;    endfor
;;
;;    refgood = where(histref gt 0)
;;    for ii = 0, n_elements(histnone)-1 do begin
;;       if rinone[ii] ne rinone[ii+1] then begin
;;;         print, rinone[rinone[ii]:rinone[ii+1]-1]
;;;         print, rnone[rinone[rinone[ii]:rinone[ii+1]-1]], $
;;;           znone[rinone[rinone[ii]:rinone[ii+1]-1]], $
;;;           grnone[rinone[rinone[ii]:rinone[ii+1]-1]], $
;;;           rznone[rinone[rinone[ii]:rinone[ii+1]-1]]
;;;         keepgoing = 1
;;;         jj = ii
;;          if riref[ii] ne riref[ii+1] then begin
;;;            print, rref[riref[riref[ii]:riref[ii+1]-1]], $
;;;              zref[riref[riref[ii]:riref[ii+1]-1]], $
;;;              grref[riref[riref[ii]:riref[ii+1]-1]]
;;
;;; randomly pick one of the "reference" galaxies
;;             theseref = riref[riref[ii]:riref[ii+1]-1]
;;             for kk = 0, histnone[ii]-1 do begin
;;                thisnone = (rinone[rinone[ii]:rinone[ii+1]-1])[kk]
;;                thisref = theseref[fix(randomu(seed,1)*n_elements(theseref))]
;;                zcat[noneindx[thisnone]].oii_3727 = zcat[refindx[thisref]].oii_3727           ; flux
;;                zcat[noneindx[thisnone]].oii_3727_2_ew = zcat[refindx[thisref]].oii_3727_2_ew ; EW
;;;               if thisnone eq 16 then stop
;;             endfor
;;          endif else begin
;;             indx = array_indices(histnone,ii)
;;;            print, raxis[indx[0]], rzaxis[indx[1]], $
;;;              graxis[indx[2]], zaxis[indx[3]]
;;             
;;; find the "nearest" galaxy with a good [OII] measurement in the ND
;;; space                 
;;             for kk = 0, histnone[ii]-1 do begin
;;                thisnone = (rinone[rinone[ii]:rinone[ii+1]-1])[kk]
;;                dist = sqrt((rref-rnone[thisnone])^2+(rzref-rznone[thisnone])^2+$
;;                  (grref-grnone[thisnone])^2+(zref-znone[thisnone])^2)
;;
;;                mindist = min(dist,thisref)
;;
;;;                print, rref[thisref], raxis[indx[0]], rnone[thisnone]
;;;                print, rzref[thisref], rzaxis[indx[1]], rznone[thisnone]
;;;                print, grref[thisref], graxis[indx[2]], grnone[thisnone]
;;;                print, zref[thisref], zaxis[indx[3]], znone[thisnone]
;;;                print
;;
;;                zcat[noneindx[thisnone]].oii_3727 = zcat[refindx[thisref]].oii_3727           ; flux
;;                zcat[noneindx[thisnone]].oii_3727_2_ew = zcat[refindx[thisref]].oii_3727_2_ew ; EW
;;;               print, thisnone
;;;               if thisnone eq 16 then stop
;;;               if (rinone[rinone[ii]:rinone[ii+1]-1])[kk] eq 16 then stop
;;             endfor
;;
;;;             while keepgoing do begin
;;;                jj++
;;;                if riref[jj] ne riref[jj+1] then begin
;;;                   keepgoing = 0
;;;                endif
;;;             endwhile
;;          endelse
;;;         print, rnone[rinone[rinone[ii]:rinone[ii+1]-1]], rref[riref[riref[jj]:riref[jj+1]-1]], $
;;;           znone[rinone[rinone[ii]:rinone[ii+1]-1]], zref[riref[riref[jj]:riref[jj+1]-1]], $
;;;           grnone[rinone[rinone[ii]:rinone[ii+1]-1]], grref[riref[riref[jj]:riref[jj+1]-1]], $
;;;           rznone[rinone[rinone[ii]:rinone[ii+1]-1]], rzref[riref[riref[jj]:riref[jj+1]-1]]
;;;         print & print, jj, oiiref[*,riref[riref[jj]:riref[jj+1]-1]] & print
;;
;;;; randomly pick one of the "reference" galaxies
;;;          theseref = riref[riref[jj]:riref[jj+1]-1]
;;;          for kk = 0, histnone[ii]-1 do begin
;;;             thisnone = (rinone[rinone[ii]:rinone[ii+1]-1])[kk]
;;;             thisref = theseref[fix(randomu(seed,1)*n_elements(theseref))]
;;;             zcat[noneindx[thisnone]].oii_3727 = zcat[refindx[thisref]].oii_3727 ; flux
;;;             zcat[noneindx[thisnone]].oii_3727_2_ew = zcat[refindx[thisref]].oii_3727_2_ew ; EW
;;;;            if zcat[refindx[thisref]].oii_3727[1] eq -2 then stop
;;;
;;;; print this as a diagnostic             
;;;             print, zcat[noneindx[thisnone]].zbest, zcat[refindx[thisref]].zbest, $
;;;               zcat[noneindx[thisnone]].ugriz[2], zcat[refindx[thisref]].ugriz[2]
;;;
;;;;            if thisnone eq 100 then stop
;;;;            if (rinone[rinone[ii]:rinone[ii+1]-1])[kk] eq 100 then stop
;;;          endfor
;;;         if histref[jj] gt 1 then stop
;;;         if total(rinone[rinone[ii]:rinone[ii+1]-1] eq 100) gt 0 then stop
;;       endif
;;    endfor
