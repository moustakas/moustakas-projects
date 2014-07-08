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

pro qaplot_desi_deep2egs
; jm14mar01siena - make some QAplots from the output of
;   BUILD_DESI_DEEP2EGS for the TechNote

    templatepath = getenv('IM_PROJECTS_DIR')+'/desi/templates/v1.1/'
    isedfit_paramfile = templatepath+'desi_deep2_paramfile.par'
    targpath = getenv('IM_PROJECTS_DIR')+'/desi/targeting/'
    technotepath = getenv('IM_SVNREPOS')+'/desi/technotes/targeting/elg-deep2/trunk/'
    winpath = deep2_path(/window)
    catpath = deep2_path(/cat)
    
    area = 0.4342 ; deg^2
;   area = 0.6 ; [deg^2]

    brightcut = 18.5
    faintcut = 24.0
    oiicut1 = 8D-17 ; [erg/s/cm2]
    oiisnrcut = 3.0




    
; -------------------------
; plot the half-light radius       
       phot = mrdfits(catpath+'deep2.pcat_ext.fits.gz',1)
       acs = mrdfits(getenv('IM_DATA_DIR')+'/acs-gc/egs_v_i_public_catalog_V1.0.fits.gz',1)
;      spherematch, phot.ra_deep, phot.dec_deep, acs.ra, acs.dec, 1D/3600.0, m1, m2
       match, phot.objno, acs.objno, m1, m2
       phot = deep2_get_ugriz(phot[m1])
       acs = acs[m2]

       hiz_mostek = desi_get_hizelg(phot.ugriz,/mostek)
       these = where(acs[hiz_mostek].flag_galfit_hi eq 0 and $
         phot[hiz_mostek].ugriz[2] lt 23.0,ngal)

       requant = weighted_quantile(acs[hiz_mostek[these]].re_galfit_hi*0.03,quant=[0.5,0.25,0.75])
       nquant = weighted_quantile(acs[hiz_mostek[these]].n_galfit_hi,quant=[0.5,0.25,0.75])
;      zmed = weighted_quantile(zcat[hiz_mostek[these]].z,$
;        zcat[hiz_mostek[these]].final_weight,quant=0.5)
       
       psfile = technotepath+'egs_halflight.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.7, $
         xmargin=[1.4,0.4], charsize=1.7

       djs_plot, [0], [0], /nodata, position=pos, /xlog, /ylog, xsty=1, ysty=1, $
         xrange=[0.01,10], yrange=[0.1,20.0], xtitle='Half-light radius (arcsec)', $
         ytitle='Sersic n'
       im_legend, ['DEEP2/EGS','18.5<r<23','grz-selected'], /left, /top, box=0, margin=0
;      im_legend, ['DEEP2/EGS','18.5<r<24','1.0<z<1.4'], /left, /top, box=0, margin=0
       im_legend, ['<r_{e}>='+strtrim(string(requant[0],format='(F12.3)'),2)+'"',$
         '<n>='+strtrim(string(nquant[0],format='(F12.3)'),2)], $
;        '<z>='+strtrim(string(zmed,format='(F12.3)'),2)], 
         /left, /bottom, box=0, margin=0, charsize=1.8
       djs_oplot, acs[hiz_mostek[these]].re_galfit_hi*0.03, acs[hiz_mostek[these]].n_galfit_hi, $
         psym=symcat(16), color=cgcolor('dark grey'), symsize=0.8
       oploterror, requant[0], nquant[0], requant[1], nquant[1], psym=symcat(6,thick=5), $
         color=cgcolor('dodger blue'), errcolor=cgcolor('dodger blue'), /lobar, symsize=3, $
         errthick=8
       oploterror, requant[0], nquant[0], requant[2], nquant[2], psym=symcat(6,thick=1), $
         color=cgcolor('dodger blue'), errcolor=cgcolor('dodger blue'), /hibar, symsize=3, $
         errthick=8
       
       im_plotconfig, psfile=psfile, /psclose, /pdf

stop       

       phot = mrdfits(targpath+'deep2egs-photparent.fits.gz',1)
;      zcat = mrdfits(targpath+'deep2egs_parent_zcat.fits.gz',1)
       zcat = mrdfits(targpath+'deep2egs-zcatparent.Q34.fits.gz',1)

       ewoiisnr = zcat.oii_3727_2_ew[0]/zcat.oii_3727_2_ew[1]
       oiisnr = zcat.oii_3727_2_amp[0]/(zcat.oii_3727_2_amp[1]+$
         (zcat.oii_3727_2_amp[1] eq 0))*(zcat.oii_3727_2_amp[1] ne 0)

       magcut1 = 23.0
       magcut2 = 23.5

; Sersic index vs half-light radius
       hiz_mostek = desi_get_hizelg(zcat.ugriz,/mostek)
       these = where(zcat[hiz_mostek].flag_galfit_hi eq 0 and zcat[hiz_mostek].ugriz[2] lt 23.0,ngal)
;      these = where(zcat.n_galfit_hi gt 0 and zcat.z gt 1.0 and zcat.z lt 1.4,ngal)
       splog, ngal, weighted_quantile(zcat[hiz_mostek[these]].z,zcat[hiz_mostek[these]].final_weight,quant=0.5)

       requant = weighted_quantile(zcat[hiz_mostek[these]].re_galfit_hi*0.03,$
         zcat[hiz_mostek[these]].final_weight,quant=[0.5,0.25,0.75])
       nquant = weighted_quantile(zcat[hiz_mostek[these]].n_galfit_hi,$
         zcat[hiz_mostek[these]].final_weight,quant=[0.5,0.25,0.75])
       zmed = weighted_quantile(zcat[hiz_mostek[these]].z,$
         zcat[hiz_mostek[these]].final_weight,quant=0.5)
       
       psfile = technotepath+'egs_halflight.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.7, $
         xmargin=[1.4,0.4], charsize=1.7

       djs_plot, [0], [0], /nodata, position=pos, /xlog, /ylog, xsty=1, ysty=1, $
         xrange=[0.01,10], yrange=[0.1,20.0], xtitle='Half-light radius (arcsec)', $
         ytitle='Sersic n'
       im_legend, ['DEEP2/EGS','18.5<r<23','grz-selected'], /left, /top, box=0, margin=0
;      im_legend, ['DEEP2/EGS','18.5<r<24','1.0<z<1.4'], /left, /top, box=0, margin=0
       im_legend, ['<r_{e}>='+strtrim(string(requant[0],format='(F12.3)'),2)+'"',$
         '<n>='+strtrim(string(nquant[0],format='(F12.3)'),2),$
         '<z>='+strtrim(string(zmed,format='(F12.3)'),2)], /left, /bottom, box=0, $
         margin=0, charsize=1.8
       djs_oplot, zcat[hiz_mostek[these]].re_galfit_hi*0.03, zcat[hiz_mostek[these]].n_galfit_hi, $
         psym=symcat(16), color=cgcolor('dark grey'), symsize=0.8
       oploterror, requant[0], nquant[0], requant[1], nquant[1], psym=symcat(6,thick=5), $
         color=cgcolor('dodger blue'), errcolor=cgcolor('dodger blue'), /lobar, symsize=3, $
         errthick=8
       oploterror, requant[0], nquant[0], requant[2], nquant[2], psym=symcat(6,thick=1), $
         color=cgcolor('dodger blue'), errcolor=cgcolor('dodger blue'), /hibar, symsize=3, $
         errthick=8
       
       im_plotconfig, psfile=psfile, /psclose, /pdf
       
stop
       
; gr vs rz - stellar contamination
       stars = mrdfits(targpath+'deep2egs-photstars.fits.gz',1)
       stars = stars[where(stars.ugriz[2] lt magcut1,nstar)]

       loz = where(zcat.zbest lt 0.6 and zcat.ugriz[2] lt magcut1,nloz)
       hiz = where(zcat.zbest gt 0.6 and zcat.zbest lt 1.2 and $
         zcat.ugriz[2] lt magcut1,nhiz)
       vhiz = where(zcat.zbest gt 1.2 and zcat.zbest lt 1.6 and $
         zcat.ugriz[2] lt magcut1,nvhiz)
       
       w1 = where(zcat.zbest gt 0.6 and zcat.ugriz[2] lt magcut1,nw1)
       w2 = where(zcat.zbest gt 1.2 and zcat.ugriz[2] lt magcut1,nw2)
       splog, 'Fraction of z>0.6 in grz box: ', n_elements(desi_get_hizelg(zcat[w1].ugriz,$
         magcut=magcut1))/float(nw1)
       splog, 'Fraction of z>1.2 in grz box: ', n_elements(desi_get_hizelg(zcat[w2].ugriz,$
         magcut=magcut1))/float(nw2)
       
       psfile = technotepath+'egs_grz_stars.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
         xmargin=[1.5,0.4], charsize=1.7

       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='r - z', ytitle='g - r', xrange=[-0.5,2.2], yrange=[-0.3,1.8]
       djs_oplot, zcat[loz].ugriz[2]-zcat[loz].ugriz[4], $
         zcat[loz].ugriz[1]-zcat[loz].ugriz[2], psym=symcat(16), symsize=0.3
       djs_oplot, zcat[hiz].ugriz[2]-zcat[hiz].ugriz[4], $
         zcat[hiz].ugriz[1]-zcat[hiz].ugriz[2], psym=symcat(5), $
         color=cgcolor('firebrick'), symsize=0.2
       djs_oplot, zcat[vhiz].ugriz[2]-zcat[vhiz].ugriz[4], $
         zcat[vhiz].ugriz[1]-zcat[vhiz].ugriz[2], psym=symcat(6), $
         color=cgcolor('forest green'), symsize=0.5
       djs_oplot, stars.ugriz[2]-stars.ugriz[4], symsize=0.1, $
         stars.ugriz[1]-stars.ugriz[2], psym=symcat(7), color='blue'
       
; Mostek       
       rzaxis = range(0.2,1.2,500)
       djs_oplot, rzaxis, poly(rzaxis,[-0.08,0.68]), line=5, thick=6
       djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[-0.08,0.68])], line=5, thick=6
       djs_oplot, 1.3*[1,1], [!y.crange[0],poly(1.3,[-0.08,0.68])], line=5, thick=6

; proposed
       rzaxis = range(0.1,1.1,500)
       int = -0.08 & slope = 1.0
       djs_oplot, rzaxis, poly(rzaxis,[int,slope]), line=0, thick=6
       djs_oplot, 0.1*[1,1], [!y.crange[0],poly(0.1,[int,slope])], line=0, thick=6
       djs_oplot, [1.1,1.8], poly(1.1,[int,slope])*[1,1], line=0, thick=6
       djs_oplot, 1.8*[1,1], [!y.crange[0],poly(1.1,[int,slope])], line=0, thick=6
       
       im_legend, ['Stars (Pgal<0.2)','z<0.6','0.6<z<1.2','1.2<z<1.6'], $
         /left, /top, box=0, psym=[7,16,5,6], position=[0.22,0.86], /norm, $
         color=['blue','','firebrick','forest green'], charsize=1.5
       im_legend, ['18.5<r<23'], spacing=2.0, /left, /top, box=0, margin=0

;      im_legend, ['Mostek','Proposed'], /right, /bottom, box=0, $
;        line=[0,5], pspacing=1.7, thick=6
       im_plotconfig, psfile=psfile, /psclose, /pdf

; gr vs rz - [OII] strength
       hiz = desi_get_hizelg(zcat.ugriz)
       loz = where(zcat.zbest lt 0.6 and zcat.ugriz[2] lt magcut1,nloz)
       oiibright = where(zcat.zbest gt 0.6 and zcat.ugriz[2] lt magcut1 and $
         zcat.oii_3727_1[1] ne -2.0 and oiisnr gt oiisnrcut and $
         ewoiisnr gt 1.0 and zcat.oii_3727[0] gt oiicut1,noiibright)
       oiifaint = where(zcat.zbest gt 0.6 and zcat.ugriz[2] lt magcut1 and $
         zcat.oii_3727_1[1] ne -2.0 and oiisnr gt oiisnrcut and $
         ewoiisnr gt 1.0 and zcat.oii_3727[0] lt oiicut1,noiifaint)
       oiinone = where(zcat.zbest gt 0.6 and zcat.ugriz[2] lt magcut1 and $
         zcat.oii_3727[1] eq -2,noiinone)

;      noiibright_cor = noiibright+1.0*noiibright/(noiifaint+noiibright)*noiinone
;      frac_oiinone = noiinone/total(zcat.zbest gt 0.6 and zcat.ugriz[2] lt magcut1)
;      splog, 'Fraction of z>0.6 with strong [OII] in grz box: ', $
;        n_elements(desi_get_hizelg(zcat[oiibright].ugriz,magcut=magcut1))/float(n_elements(hiz))

       psfile = technotepath+'egs_grz_oii.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
         xmargin=[1.5,0.4], charsize=1.7

       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='r - z', ytitle='g - r', xrange=[-0.5,2.2], yrange=[-0.3,1.8]
       djs_oplot, zcat[loz].ugriz[2]-zcat[loz].ugriz[4], $
         zcat[loz].ugriz[1]-zcat[loz].ugriz[2], psym=symcat(16), symsize=0.3
       djs_oplot, zcat[oiibright].ugriz[2]-zcat[oiibright].ugriz[4], $
         zcat[oiibright].ugriz[1]-zcat[oiibright].ugriz[2], psym=symcat(5), $
         color=cgcolor('firebrick'), symsize=0.3
       djs_oplot, zcat[oiifaint].ugriz[2]-zcat[oiifaint].ugriz[4], $
         zcat[oiifaint].ugriz[1]-zcat[oiifaint].ugriz[2], psym=symcat(6), $
         color=cgcolor('forest green'), symsize=0.5
       djs_oplot, zcat[oiinone].ugriz[2]-zcat[oiinone].ugriz[4], $
         zcat[oiinone].ugriz[1]-zcat[oiinone].ugriz[2], psym=symcat(7), $
         color=cgcolor('blue'), symsize=0.3
       
; Mostek       
       rzaxis = range(0.2,1.2,500)
       djs_oplot, rzaxis, poly(rzaxis,[-0.08,0.68]), line=5, thick=6
       djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[-0.08,0.68])], line=5, thick=6
       djs_oplot, 1.3*[1,1], [!y.crange[0],poly(1.3,[-0.08,0.68])], line=5, thick=6

; proposed
       rzaxis = range(0.1,1.1,500)
       int = -0.08 & slope = 1.0
       djs_oplot, rzaxis, poly(rzaxis,[int,slope]), line=0, thick=6
       djs_oplot, 0.1*[1,1], [!y.crange[0],poly(0.1,[int,slope])], line=0, thick=6
       djs_oplot, [1.1,1.8], poly(1.1,[int,slope])*[1,1], line=0, thick=6
       djs_oplot, 1.8*[1,1], [!y.crange[0],poly(1.1,[int,slope])], line=0, thick=6
       
       im_legend, ['z<0.6','z>0.6; [OII]>8\times10^{-17}','z>0.6; [OII]<8\times10^{-17}',$
         'z>0.6; No [OII]'], /left, /top, box=0, psym=[16,5,6,7], $
         position=[0.22,0.86], /norm, $
         color=['','firebrick','forest green','blue'], charsize=1.4
       im_legend, ['18.5<r<23'], spacing=2.0, /left, /top, box=0, margin=0

       im_plotconfig, psfile=psfile, /psclose, /pdf
             
; [OII] flux vs redshift
       hiz = desi_get_hizelg(zcat.ugriz)
       hiz_mostek = desi_get_hizelg(zcat.ugriz,/mostek)

       psfile = technotepath+'egs_oii_redshift.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
         xmargin=[1.5,0.4], charsize=1.7

       factor = 1D17
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='Redshift', ytitle='F([O II] \lambda\lambda3726,29) '+$
         '(10^{-17} erg s^{-1} cm^{-2})', $
         xrange=[0.65,1.6], yrange=[0.1,300], /ylog;, title='grz Sample'

       good = where(zcat.oii_3727_1[1] ne -2.0 and oiisnr gt oiisnrcut and $
         ewoiisnr gt 1.0 and zcat.ugriz[2] lt magcut2,ngood)
       lim = where(zcat.oii_3727_1[1] ne -2.0 and zcat.ugriz[2] lt magcut2 and $
         ((oiisnr lt oiisnrcut) or (ewoiisnr lt 1.0)),nlim)
       none = where(zcat.oii_3727[1] eq -2 and zcat.ugriz[2] lt magcut2,nnone)
       splog, ngood, nlim, nnone
       
       djs_oplot, zcat[good].zbest, factor*zcat[good].oii_3727[0], $
         psym=symcat(16), symsize=0.6, color=cgcolor('dark grey')
       djs_oplot, zcat[lim].zbest, factor*zcat[lim].oii_3727_limit[0], $
         psym=symcat(11,thick=4), symsize=0.8, color=cgcolor('firebrick')
       djs_oplot, zcat[none].zbest, replicate(1.0,nnone), $
         psym=symcat(6), symsize=0.4, color=cgcolor('dodger blue')
       djs_oplot, !x.crange, factor*oiicut1*[1,1], line=0, thick=6
;      xyouts, 1.5, factor*4.5D-17, textoidl('8\times10^{-17}!cerg s^{-1} cm^{-2}'), $
       xyouts, 1.47, factor*4.8D-17, textoidl('8\times10^{-17}'), $
         align=0.0, charsize=1.3, /data
       im_legend, ['18.5<r<23.5'], spacing=2.0, /right, /top, box=0, margin=0
       
       im_legend, ['Well-measured','Upper limit (1\sigma)','No measurement'], $
         /right, /bottom, box=0, psym=[16,11,6], $
         color=['dark grey','firebrick','dodger blue'], $
         charsize=1.6
       im_plotconfig, psfile=psfile, /psclose, /pdf

; cumulative number counts

; assume that the objects with formal flux limits above our [OII] cut
; are detections (this is a small number of objects)    
       hiz = desi_get_hizelg(zcat.ugriz)
       oiibright = where(zcat[hiz].oii_3727_1[1] ne -2.0 and oiisnr gt oiisnrcut and $
         ewoiisnr gt 1.0 and zcat[hiz].oii_3727[0] gt oiicut1,noiibright)
       oiifaint = where(zcat[hiz].oii_3727_1[1] ne -2.0 and oiisnr gt oiisnrcut and $
         ewoiisnr gt 1.0 and zcat[hiz].oii_3727[0] lt oiicut1,noiifaint)
       oiinone = where(zcat[hiz].oii_3727[1] eq -2,noiinone)

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
       denom = dndm_oiibright+dndm_oiifaint
       dndm_oiicor = dndm_oiinone*dndm_oiibright/(denom+(denom eq 0))*(denom ne 0)
       dndm_oiibright_cor = dndm_oiibright + dndm_oiicor    
       
       psfile = technotepath+'egs_oii_dndm.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
         xmargin=[1.5,0.4], charsize=1.7

       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='r (AB mag)', ytitle='log N(<r) (gal / deg^{2})', $
         xrange=[brightcut-0.5,faintcut+0.5], yrange=alog10([10,1E5]) ;, /ylog
       djs_oplot, magaxis, alog10(total(dndm_phot,/cumu)), $
         line=2, thick=6, color=cgcolor('dark grey')
       djs_oplot, magaxis, alog10(total(dndm,/cumu)), line=0, thick=6
       
       ww = where(total(dndm_hiz,/cumu) gt 0)
       djs_oplot, magaxis[ww], alog10((total(dndm_hiz,/cumu))[ww]), $
         color='orange', line=3, thick=6

;       ww = where(total(dndm_oiibright,/cumu) gt 0)
;       djs_oplot, magaxis[ww], alog10((total(dndm_oiibright,/cumu))[ww]), color='red', $
;         line=1, thick=6 ; psym=symcat(6,thick=4), symsize=1.3

       ww = where(total(dndm_oiibright_cor,/cumu) gt 0)
       djs_oplot, magaxis[ww], alog10((total(dndm_oiibright_cor,/cumu))[ww]), $
         color=cgcolor('blue'), line=5, thick=8 ; psym=symcat(6,thick=4), symsize=1.3
       
       numcut = 2400
       rcut = interpol(magaxis,total(dndm_oiibright_cor,/cumu),numcut)
;      rcut = 22.7
;      numcut = long(interpol(total(dndm_oiibright_cor,/cumu),magaxis,rcut))
       djs_oplot, [!x.crange[0],rcut], alog10(numcut)*[1,1], line=1, thick=6
       djs_oplot, rcut*[1,1], [!y.crange[0],alog10(numcut)], line=1, thick=6

       xyouts, 18.3, alog10(numcut*1.2), textoidl(string(numcut,$
         format='(I0)')+' gal/deg^{2}'), align=0.0, /data, charsize=1.5
       xyouts, rcut+0.1, 2.5, textoidl('r='+string(rcut,$
         format='(F4.1)')), align=0.0, orientation=270, /data, charsize=1.5
       
       im_legend, ['PhotParent','zcatParent','grz Sample',$
         'grz & F([OII])>8\times10^{-17}'], $
         color=['dark grey','','orange','blue'], /left, /top, box=0, thick=6, $
         charsize=1.3, spacing=2.0, margin=0, line=[2,0,3,5], pspacing=1.5
       
       im_plotconfig, psfile=psfile, /psclose, /pdf

; redshift histogram of sources selected using my grz color-cuts 
       hiz = desi_get_hizelg(zcat.ugriz)
       oiibright = where(zcat[hiz].oii_3727_1[1] ne -2.0 and oiisnr gt oiisnrcut and $
         ewoiisnr gt 1.0 and zcat[hiz].oii_3727[0] gt oiicut1,noiibright)
       oiifaint = where(zcat[hiz].oii_3727_1[1] ne -2.0 and oiisnr gt oiisnrcut and $
         ewoiisnr gt 1.0 and zcat[hiz].oii_3727[0] lt oiicut1,noiifaint)
       oiinone = where(zcat[hiz].oii_3727[1] eq -2,noiinone)

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

; correct for the missing [OII] sources       
       nz_oiinone = hogg_histogram(zcat[hiz[oiinone]].zbest,[zmin,zmax],$
         nzbin,weight=zcat[hiz[oiinone]].final_weight)/area
       nz_oiifaint = hogg_histogram(zcat[hiz[oiifaint]].zbest,[zmin,zmax],$
         nzbin,weight=zcat[hiz[oiifaint]].final_weight)/area

       denom = nz_oiibright+nz_oiifaint
       nz_oiicor = nz_oiinone*nz_oiibright/(denom+(denom eq 0))*(denom ne 0)
       nz_oiibright_cor = nz_oiibright + nz_oiicor    

       psfile = technotepath+'egs_oii_zhist.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
         xmargin=[1.5,0.4], charsize=1.7
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='Redshift', ytitle='dN/dz (gal / 0.1dz / deg^{2})', $
         xrange=[zmin,zmax], yrange=[0.0,6500]
       djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
         [0,nz,0], psym=10, thick=6, line=0
       djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
         [0,nz_hiz,0], psym=10, thick=8, line=3, color='orange'
;      djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
;        [0,nz_oiibright,0], psym=10, thick=6, line=5, color='orange'
       djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
         [0,nz_oiibright_cor,0], psym=10, thick=8, line=5, color='blue'

;       im_plotconfig, 6, pos2, height=[6.0,2.5], width=6.6, $
;         xmargin=[1.5,0.4], charsize=1.7
;       djs_plot, [0], [0], /nodata, position=pos2[*,0], xsty=1, ysty=1, $
;         xtitle='', ytitle='dN/dz (gal / 0.1dz / deg^{2})', $
;         xrange=[zmin,zmax], yrange=[0.0,5000], xtickname=replicate(' ',10)
;       djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
;         [0,nz,0], psym=10, thick=6, line=5
;       djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
;         [0,nz_hiz,0], psym=10, thick=6, line=0, color='red'
;
;       djs_plot, [0], [0], /nodata, /noerase, position=pos2[*,1], xsty=1, ysty=1, $
;         xtitle='Redshift', ytitle='Fraction', $
;         xrange=[zmin,zmax], yrange=[0.0,1.1], ytickinterval=0.5
;       djs_oplot, !x.crange, [1,1], line=1, thick=6
;       djs_oplot, zhist, nz_hiz/nz, psym=10, thick=6

       im_legend, ['zcatParent','grz Sample','grz & F([OII])>8\times10^{-17}'], $
         line=[0,3,5], pspacing=1.7, color=['','orange','blue'], /right, /top, box=0, $
         charsize=1.4, spacing=2.0, margin=0, thick=6
;      im_legend, ['Parent Sample (N_{eff}='+strtrim(long(total(nz)*area),2)+')',$
;        'grz Sample (N_{eff}='+strtrim(long(total(nz_hiz)*area),2)+')'], $
;        line=[5,0], pspacing=1.7, color=['','red'], /right, /top, box=0, $
;        charsize=1.4, spacing=2.0, margin=0, thick=6

       im_plotconfig, psfile=psfile, /psclose, /pdf
       
;; gr vs rz
;       psfile = technotepath+'egs_grz.ps'
;       im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
;         xmargin=[1.5,0.4], charsize=1.7
;
;       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;         xtitle='r - z', ytitle='g - r', xrange=[-0.5,2.5], yrange=[-0.5,2.0]
;       djs_oplot, zcat[loz].ugriz[2]-zcat[loz].ugriz[4], $
;         zcat[loz].ugriz[1]-zcat[loz].ugriz[2], psym=symcat(16), symsize=0.3
;       djs_oplot, zcat[hiz].ugriz[2]-zcat[hiz].ugriz[4], $
;         zcat[hiz].ugriz[1]-zcat[hiz].ugriz[2], psym=symcat(5), $
;         color=cgcolor('firebrick'), symsize=0.2
;       djs_oplot, zcat[vhiz].ugriz[2]-zcat[vhiz].ugriz[4], $
;         zcat[vhiz].ugriz[1]-zcat[vhiz].ugriz[2], psym=symcat(6), $
;         color=cgcolor('dodger blue'), symsize=0.4
;       im_legend, ['18.5<r<24'], spacing=2.0, /right, /top, box=0
;       
;       rzaxis = range(0.2,1.2,500)
;       djs_oplot, rzaxis, poly(rzaxis,[-0.08,0.68]), line=0, thick=6
;       djs_oplot, 0.2*[1,1], [!y.crange[0],poly(0.2,[-0.08,0.68])], line=0, thick=6
;       djs_oplot, [1.2,1.8], poly(1.2,[-0.08,0.68])*[1,1], line=0, thick=6
;       djs_oplot, 1.8*[1,1], [!y.crange[0],poly(1.2,[-0.08,0.68])], line=0, thick=6
;;      djs_oplot, 1.3*[1,1], [!y.crange[0],poly(1.3,[-0.08,0.68])], line=0, thick=6
;       
;       rzcut = 1.2
;       blue = where(rzaxis lt rzcut,comp=red)
;;      djs_oplot, rzaxis[blue], 0.7*rzaxis[blue]+0.05, line=5, thick=6
;       djs_oplot, rzaxis[red], 1.4*rzaxis[red]-0.79, line=5, thick=6
;       im_legend, ['z<0.7','0.7<z<1.2','z>1.2'], /left, /top, box=0, $
;         psym=[16,5,6], color=['','firebrick','dodger blue']
;       
;;      im_legend, ['Mostek','Proposed'], /right, /bottom, box=0, $
;;        line=[0,5], pspacing=1.7, thick=6
;       im_plotconfig, psfile=psfile, /psclose, /pdf
;
;stop       
;       
;; ug vs gr    
;       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;         xtitle='g - r', ytitle='u - g', xrange=[-0.5,2.0], yrange=[-0.5,3]
;       djs_oplot, zcat[loz].ugriz[1,*]-zcat[loz].ugriz[2,*], $
;         zcat[loz].ugriz[0,*]-zcat[loz].ugriz[1,*], psym=symcat(15), symsize=0.3
;       djs_oplot, zcat[hiz].ugriz[1,*]-zcat[hiz].ugriz[2,*], $
;         zcat[hiz].ugriz[0,*]-zcat[hiz].ugriz[1,*], psym=symcat(16), $
;         color='orange', symsize=0.3
;       im_legend, ['z<0.7','z>0.7'], /left, /top, box=0, $
;         psym=[15,16], color=['','orange']
       
    endif

; ###########################################################################
; build some QAplots of the parent samples for the technote 
    if keyword_set(qaplots_parent) then begin
       allphot = mrdfits(deep2_path(/cat)+'deep2.pcat_ext.fits.gz',1)    
       egs = where(strmid(strtrim(allphot.objno,2),0,1) eq '1',negs)
       allphot = allphot[egs]

       allzcat = read_deep2_zcat(/all)
       egs = where(strmid(strtrim(allzcat.objno,2),0,1) eq '1',negs)
       allzcat = allzcat[egs]
       
       phot = mrdfits(targpath+'deep2egs-photparent.fits.gz',1)
       zcat = mrdfits(targpath+'deep2egs-zcatparent.fits.gz',1)
       zcat_q34 = mrdfits(targpath+'deep2egs-zcatparent.Q34.fits.gz',1)

       weight = zcat.targ_weight
       weight_q34 = zcat_q34.final_weight

       nphot = n_elements(phot)
       nzcat = n_elements(zcat)
       nzcat_weighted = long(total(weight))

       nzcat_q34 = n_elements(zcat_q34)
       nzcat_q34_weighted = long(total(weight_q34))
       splog, nzcat, nzcat_weighted, nzcat_q34, nzcat_q34_weighted

; ra, dec       
       psfile = technotepath+'egs_radec.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
         xmargin=[1.5,0.4], charsize=1.7;, ypage=6.5, xpage=8.5

       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='R. A. (J2000, deg)', ytitle='Dec (J2000, deg)', $
         xrange=[216.5,213.0], yrange=[51.8,53.8], xtickinterval=1.0;, $
;        title='DEEP2/EGS (Area='+string(area,format='(F5.3)')+' deg^{2})'
       djs_oplot, allphot.ra_deep, allphot.dec_deep, psym=3, color='grey'
       djs_oplot, phot.ra, phot.dec, psym=6, symsize=0.15, color='red'
       djs_oplot, allzcat.ra, allzcat.dec, psym=7, symsize=0.1, color='cyan'
       djs_oplot, zcat.ra, zcat.dec, psym=symcat(16), symsize=0.1, color='blue'

       im_legend, ['PhotParent (N='+strtrim(nphot,2)+')',$
         'zcatParent (N='+strtrim(nzcat_q34,2)+')'],$
         psym=[6,16], color=['red','blue'], /left, /bottom, box=0, $
         charsize=1.4, spacing=2.0, margin=0

       im_legend, ['Matthews+13 (all photometry)','Newman+13 (all spectroscopy)'], $
         psym=[16,6], color=['grey','cyan'], symsize=0.5, /right, /top, box=0, $
         charsize=1.2, spacing=2.0, margin=0

;       im_legend, ['Matthews+13 (N='+strtrim(nphot,2)+')'], $
;         /left, /top, box=0, psym=6, color='red', $
;         position=[0.58,0.91], /norm, charsize=1.4

;      im_legend, ['18.5<r<24','w>0','P_{gal}>0.5'], spacing=2.0, $
;        im_legend, ['18.5<r<24','P_{gal}>0.5'], spacing=2.0, $
;          /left, /top, box=0, position=[0.62,0.87], /norm, charsize=1.4

;        im_legend, ['DEEP2 (N='+strtrim(nspec,2)+')'], $
;          /left, /bottom, box=0, psym=16, color='blue', $
;          position=[0.18,0.35], /norm, charsize=1.4
;        im_legend, ['18.5<r<24'], spacing=1.9, $
;          /left, /bottom, box=0, position=[0.22,0.28], /norm, charsize=1.4
;      im_legend, ['18.5<r<24','Q>=3'], spacing=1.9, $
;        /left, /bottom, box=0, position=[0.22,0.26], /norm, charsize=1.4
;      im_legend, ['18.5<r<24','w>0','Q>=3'], spacing=1.9, $
;        /left, /bottom, box=0, position=[0.22,0.22], /norm, charsize=1.4
       im_plotconfig, psfile=psfile, /psclose, /pdf

; cumulative number counts 
       dndm = get_dndm(phot.ugriz[2,*],magaxis=magaxis,faintcut=faintcut,$
         brightcut=brightcut)
       zcat_dndm_noweight = get_dndm(zcat_q34.ugriz[2,*],$
         faintcut=faintcut,brightcut=brightcut)
;      zcat_dndm = get_dndm(zcat.ugriz[2,*],weight=weight,$
;        faintcut=faintcut,brightcut=brightcut)
       zcat_dndm = get_dndm(zcat_q34.ugriz[2,*],weight=weight_q34,$
         faintcut=faintcut,brightcut=brightcut)
       
       psfile = technotepath+'egs_dndm.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
         xmargin=[1.5,0.4], charsize=1.7

       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='r (AB mag)', ytitle='log N(<r) (gal / deg^{2})', $
         xrange=[brightcut,faintcut], yrange=alog10([10,1E5]);, /ylog
       djs_oplot, magaxis, alog10(total(dndm,/cumu)), psym=symcat(16)
       djs_oplot, magaxis, alog10(total(zcat_dndm_noweight,/cumu)), $
         color='orange', psym=symcat(15)
       djs_oplot, magaxis, alog10(total(zcat_dndm,/cumu)), color='red', $
         psym=symcat(6,thick=4), symsize=1.3

       im_legend, ['PhotParent (N='+strtrim(nphot,2)+')',$
         'zcatParent (observed, N='+strtrim(nzcat_q34,2)+')',$
         'zcatParent (weighted, N_{eff}='+strtrim(nzcat_q34_weighted,2)+')'], $
         psym=[16,15,6], color=['','orange','red'], /left, /top, box=0, $
         charsize=1.4, spacing=2.0, margin=0
;      im_legend, ['Matthews+13 (N='+strtrim(nphot,2)+')',$
;        'DEEP2 (unweighted, N='+strtrim(nspec,2)+')',$
;        'DEEP2 (weighted, N_{eff}='+strtrim(nspec_weighted,2)+')'], $
;        psym=[16,15,6], color=['','orange','red'], /left, /top, box=0, $
;        charsize=1.4, spacing=2.0, margin=0
       im_plotconfig, psfile=psfile, /psclose, /pdf

; redshift histogram
       psfile = technotepath+'egs_zhist.ps'
       im_plotconfig, 0, pos, psfile=psfile, height=5.0, width=6.6, $
         xmargin=[1.5,0.4], charsize=1.7

       zmin = 0.0
       zmax = 1.8
       zbin = 0.1
       nzbin = ceil((zmax-zmin)/zbin)
       zhist = lindgen(nzbin)*zbin+zmin+zbin/2.0
       nn = hogg_histogram(zcat_q34.zbest,[zmin,zmax],nzbin)/area
       nn_cor = hogg_histogram(zcat_q34.zbest,[zmin,zmax],nzbin,$
         weight=weight_q34)/area
  
       djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
         xtitle='Redshift', ytitle='dN/dz (gal / 0.1dz / deg^{2})', $
         xrange=[0.0,1.8], yrange=[0.0,7000]
       djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
         [0,nn_cor,0], psym=10, thick=6, line=0, color='red'
       djs_oplot, [(zhist[0]-zbin/2.0),zhist,(zhist[nzbin-1]+zbin/2.0)], $
         [0,nn,0], psym=10, thick=6, line=5
;      xyouts, 0.7, 2400.0, 'x15', charsize=1.5, align=0.5, /data
       
       im_legend, ['zcatParent (unweighted, N='+strtrim(nzcat_q34,2)+')',$
         'zcatParent (weighted, N_{eff}='+strtrim(nzcat_q34_weighted,2)+')'], $
         line=[5,0], pspacing=1.7, color=['','red'], /right, /top, box=0, $
         charsize=1.3, spacing=2.0, margin=0       
       im_plotconfig, psfile=psfile, /psclose, /pdf
    endif

    
return
end
    
