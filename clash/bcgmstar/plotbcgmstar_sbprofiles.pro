pro plotbcgmstar_sbprofiles, pdf=pdf
; jm13oct19siena - plot the SB profiles

    paperpath = bcgmstar_path(/paper)
    ellpath = bcgmstar_path(/ellipse)
    sersicpath = bcgmstar_path(/sersic)

    sample = read_bcgmstar_sample(/zsort)
    ncl = n_elements(sample)

    pixscale = 0.065
    errfloor = 0.0D ; 0.02      ; magnitude error floor on my SB measurements 
    ncol = 3 ; number of columns
    
    nicecluster = repstr(repstr(strupcase(sample.dirname),'ABELL','Abell'),'_',' ')

; ---------------------------------------------------------------------------
; make 15 figures for the appendix which shows all the clusters and
; all the bands (see bcgmstar_sersicfit, /qaplot_sbprofiles
    modelrr = [0,range(0.01,200,500,/log)] ; equivalent radius [kpc]

    allsersic = 1
    dosersic2 = 0

    xrange = [0.3,120]
    
;   for ic = 1, 1 do begin
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)

       psfile = paperpath+'bcg_all_sbprofiles_'+cluster+'.ps'
       im_plotconfig, 0, pos, psfile=psfile, charsize=1.3

       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]

       if allsersic then begin
          sersic = mrdfits(sersicpath+cluster+'-allsersic.fits.gz',1,/silent)
          sersic_results = mrdfits(sersicpath+cluster+'-allsersic-results.fits.gz',1,/silent)
       endif else begin
          sersic = mrdfits(sersicpath+cluster+'-sersic.fits.gz',1,/silent)
       endelse
       modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
       nfilt = n_elements(modphot)

       nrow = ceil(nfilt/float(ncol))
;      if allsersic then nrow = ceil((nfilt+1)/float(ncol)) else $
;        nrow = ceil(nfilt/float(ncol))
       wid = 2.4
       hei = 3
       ym = [0.2,1.1]
       xm = [0.9,0.4]
       xpage = total(xm)+total(replicate(wid,ncol))
       ypage = total(ym)+total(replicate(hei,nrow))
       pos = im_getposition(nx=ncol,ny=nrow,yspace=0.0,xspace=0.0,$
         xmargin=xm,ymargin=ym,width=wid,height=hei,xpage=xpage,$
         ypage=ypage)
       count = 0

;      for ib = nfilt-1, 0, -1 do begin ; reverse order
       for ib = 0, nfilt-1 do begin
          band = strtrim(strupcase(modphot[ib].band),2)

          modgood = where(modphot[ib].majora*pixscale*arcsec2kpc le sersic[ib].amax_kpc and $
            modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0)
          modbad = where(modphot[ib].majora*pixscale*arcsec2kpc gt sersic[ib].amax_kpc and $
            modphot[ib].sb0fit gt 0 and modphot[ib].sb0fit_ivar gt 0)
          
          rr = modphot[ib].radius_kpc[modgood] ; [kpc]
          sb = -2.5*alog10(modphot[ib].sb0fit[modgood])
          sberr = 2.5/(modphot[ib].sb0fit[modgood]*sqrt(modphot[ib].sb0fit_ivar[modgood])*alog(10))

; just for the plot make the uncertainty have a minimum floor so that
; the SB profile shows up!          
          sberr = sberr>0.1
          
;         if count eq 1 then title = strupcase(cluster) else delvarx, title
          if count ge nfilt-3 then begin
             delvarx, xtickname
          endif else begin
             xtickname = replicate(' ',10)
          endelse
          if (count mod 3) eq 0 then begin
             delvarx, ytickname
          endif else begin
             ytickname = replicate(' ',10)
          endelse
          
          djs_plot, [0], [0], /nodata, /xlog, noerase=count gt 0, $
            xrange=xrange, xsty=1, yrange=[28,16], position=pos[*,count], $
            xtickname=xtickname, ytickname=ytickname, $ ; title=title, $
            symsize=0.5, ytickinterval=3, ysty=1
;         djs_oplot, rr, sb, color=cgcolor('dark grey'), thick=10
;         djs_oplot, rr, sb, psym=symcat(15), symsize=0.5
          polyfill, [rr,reverse(rr)], [sb-sberr,reverse(sb+sberr)], /fill, $
            color=cgcolor('medium grey'), noclip=0
          djs_oplot, rr, sb+sberr, line=0, color=cgcolor('dark grey')
          djs_oplot, rr, sb-sberr, line=0, color=cgcolor('dark grey')

          if count eq 0 then im_legend, nicecluster[ic], /right, $
            /top, box=0, charsize=1.5, margin=0, position=[pos[2,count]+0.005,pos[3,count]-0.02], /norm
;         if count eq 0 then im_legend, nicecluster[ic], /left, /bottom, box=0, charsize=1.5, margin=0
          
          if allsersic then begin
             label = [$
               '\mu_{e1}='+strtrim(string(-2.5*alog10(sersic[ib].sersic2_all_sbe1),format='(F12.1)'),2)+','+$
               'n_{1}='+strtrim(string(sersic[ib].sersic2_all_n1,format='(F12.2)'),2)+','+$
               'r_{e1}='+strtrim(string(sersic[ib].sersic2_all_re1,format='(F12.2)'),2)+' kpc',$
               '\mu_{e2}='+strtrim(string(-2.5*alog10(sersic[ib].sersic2_all_sbe2),format='(F12.1)'),2)+','+$
               'n_{2}='+strtrim(string(sersic[ib].sersic2_all_n2,format='(F12.2)'),2)+','+$
               'r_{e2}='+strtrim(string(sersic[ib].sersic2_all_re2,format='(F12.2)'),2)+' kpc']
          endif else begin
             label = [$
               '\chi^{2}_{\nu, single}='+$
               strtrim(string(sersic[ib].sersic_chi2/sersic[ib].sersic_dof,format='(F12.2)'),2),$
               '\mu_{e}='+strtrim(string(sersic[ib].sersic_sbe,format='(F12.1)'),2)+','+$
               'n='+strtrim(string(sersic[ib].sersic_n,format='(F12.2)'),2)+','+$
               'r_{e}='+strtrim(string(sersic[ib].sersic_re,format='(F12.1)'),2)+' kpc']
             if dosersic2 then begin
                if sersic[ib].sersic2_sbe1 eq 0.0 or sersic[ib].sersic2_sbe2 eq 0.0 then begin
                   label = [label,'Sersic-2 dropped']
                endif else begin
                   label = [label,$
                     '\chi^{2}_{\nu, double}='+strtrim(string(sersic[ib].sersic2_chi2/$
                     sersic[ib].sersic2_dof,format='(F12.2)'),2),$
                     '\mu_{e1}='+strtrim(string(-2.5*alog10(sersic[ib].sersic2_sbe1),format='(F12.1)'),2)+','+$
                     'n_{1}='+strtrim(string(sersic[ib].sersic2_n1,format='(F12.2)'),2)+','+$
                     'r_{e1}='+strtrim(string(sersic[ib].sersic2_re1,format='(F12.2)'),2)+' kpc',$
                     '\mu_{e2}='+strtrim(string(-2.5*alog10(sersic[ib].sersic2_sbe2),format='(F12.1)'),2)+','+$
                     'n_{2}='+strtrim(string(sersic[ib].sersic2_n2,format='(F12.2)'),2)+','+$
                     'r_{e2}='+strtrim(string(sersic[ib].sersic2_re2,format='(F12.2)'),2)+' kpc']
                endelse 
             endif
          endelse 
;         im_legend, label, /left, /bottom, box=0, margin=0, charsize=0.7, charthick=1.8

;          if modbad[0] ne -1 then begin
;             djs_oplot, modphot[ib].radius_kpc[modbad], -2.5*alog10(modphot[ib].sb0fit[modbad]), $
;               psym=symcat(9), color=cgcolor('medium grey'), symsize=0.5
;          endif
          
          im_legend, band, /left, /bottom, box=0, margin=0, charsize=1.2;, charthick=1.8
;         im_legend, band, /right, /top, box=0, margin=0, charsize=1.2;, charthick=1.8
          
          if allsersic then begin
;            djs_oplot, modelrr, -2.5*alog10(bcgmstar_sersic2_func(modelrr,params=sersic[0],/allbands)), $
;              color=cgcolor('forest green')
             djs_oplot, modelrr, -2.5*alog10(bcgmstar_sersic2_func(modelrr,params=sersic[ib],/allbands)), $
               color=cgcolor('tomato'), thick=6
             djs_oplot, modelrr, bcgmstar_sersic_func(modelrr,[-2.5*alog10(sersic[ib].sersic2_all_sbe1),$
               sersic[ib].sersic2_all_re1,sersic[ib].sersic2_all_n1]), color=cgcolor('forest green'), line=2
             djs_oplot, modelrr, bcgmstar_sersic_func(modelrr,[-2.5*alog10(sersic[ib].sersic2_all_sbe2),$
               sersic[ib].sersic2_all_re2,sersic[ib].sersic2_all_n2]), color=cgcolor('forest green'), line=2
          endif else begin
             djs_oplot, modelrr, bcgmstar_sersic_func(modelrr,params=sersic[ib]), $
               color=cgcolor('firebrick')
             if ib gt 0 then djs_oplot, modelrr, bcgmstar_sersic_func(modelrr,$
               params=sersic[0]), color=cgcolor('forest green')
             
             if dosersic2 then begin
                djs_oplot, modelrr, -2.5*alog10(bcgmstar_sersic2_func(modelrr,params=sersic[ib])), $
                  color=cgcolor('dodger blue')
                if sersic[ib].sersic2_sbe1 eq 0.0 or sersic[ib].sersic2_sbe2 eq 0.0 then begin
                   splog, '  '+band+': second Sersic dropped!'
                endif else begin
                   djs_oplot, modelrr, bcgmstar_sersic_func(modelrr,[-2.5*alog10(sersic[ib].sersic2_sbe1),$
                     sersic[ib].sersic2_re1,sersic[ib].sersic2_n1]), color=cgcolor('orange'), line=2
                   djs_oplot, modelrr, bcgmstar_sersic_func(modelrr,[-2.5*alog10(sersic[ib].sersic2_sbe2),$
                     sersic[ib].sersic2_re2,sersic[ib].sersic2_n2]), color=cgcolor('orange'), line=2
                endelse
             endif
          endelse
          
          djs_oplot, 10^!x.crange, (modphot[ib].sblimit-2.5*alog10(3.0))*[1,1], line=1 ; 3-sigma
;         djs_oplot, [70.0,10^!x.crange[1]], modphot[ib].sblimit*[1,1], line=0
          count++
       endfor 

;      if allsersic then begin
;         label = ['\chi^{2}_{\nu}='+strtrim(string(sersic_results.chi2/$
;           sersic_results.dof,format='(F12.2)'),2),$
;           '\nu='+strtrim(string(sersic_results.dof,format='(I0)'),2),$
;           '\alpha_{1}='+strtrim(string(sersic_results.alpha1,format='(F12.3)'),2)+'\pm'+$
;           strtrim(string(sersic_results.alpha1_err,format='(F12.3)'),2),$
;           '\alpha_{2}='+strtrim(string(sersic_results.alpha2,format='(F12.3)'),2)+'\pm'+$
;           strtrim(string(sersic_results.alpha2_err,format='(F12.3)'),2),$
;           '\beta_{1}='+strtrim(string(sersic_results.beta1,format='(F12.3)'),2)+'\pm'+$
;           strtrim(string(sersic_results.beta1_err,format='(F12.3)'),2),$
;           '\beta_{2}='+strtrim(string(sersic_results.beta2,format='(F12.3)'),2)+'\pm'+$
;           strtrim(string(sersic_results.beta2_err,format='(F12.3)'),2)]
;         djs_plot, [0], [0], /nodata, /noerase, position=pos[*,count], $
;           xsty=5, ysty=5
;         im_legend, label, /left, /bottom, box=0, margin=0, charsize=1, charthick=1.8
;      endif
          
       xyouts, min(pos[0,*])-0.06, (max(pos[3,*])-min(pos[1,*]))/2.0+min(pos[1,*]), $
         textoidl('\mu (mag arcsec^{-2})'), orientation=90, align=0.5, charsize=1.6, /norm
       xyouts, (max(pos[2,*])-min(pos[0,*]))/2.0+min(pos[0,*]), min(pos[1,*])-0.07, $
         textoidl('Equivalent Radius (kpc)'), align=0.5, charsize=1.6, /norm
       im_plotconfig, psfile=psfile, /psclose, /pdf
    endfor 

stop    

; ---------------------------------------------------------------------------
; make a plot for the main body of the paper which shows the surface
; brightness profiles and Sersic fits for three representative
; bandpasses     
    
; plot these bands    
    thesefilt = ['f475w','f814w','f160w']
    color1 = ['medium grey','powder blue','tomato']
    color2 = ['black','navy','firebrick']
    line = [2,3,0]

;   thesefilt = ['f160w','f814w','f475w']
;   color1 = ['tomato','powder blue','medium grey']
;   color2 = ['firebrick','navy','black']

;   color = ['orange','firebrick','navy']
;   color2 = ['orange red','orange','dodger blue']
    nthese = n_elements(thesefilt)

; overplot the galaxy photometry?  overplot Marc's SB profiles?
    plotgal = 0
    plotmarc = 0
    
; make the plot
;   xtitle = 'Semi-Major Axis (kpc)'    
;   ytitle = '\mu (mag arcsec^{-2})'
    xrange = [0.5,100]
;   xrange = [3*pixscale*arcsec2kpc,150]
    yrange = [28,16]

    psfile = paperpath+'bcg_example_sbprofiles.eps'
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.3
    pos = im_getposition(nx=5,ny=3,/landscape,yspace=0.1,xspace=0.1,$
      xmargin=[1.0,0.4])

;   for ic = 0, 0 do begin
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       print & splog, cluster, sample[ic].z
;      datapath = bcgmstar_path(/bcg)+cluster+'/'

       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
       
       sersic = mrdfits(sersicpath+cluster+'-allsersic.fits.gz',1,/silent)
;      sersic = mrdfits(sersicpath+cluster+'-sersic.fits.gz',1,/silent)
;      galphot = mrdfits(ellpath+cluster+'-ellipse-image.fits.gz',1,/silent)
       modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
       nfilt = n_elements(modphot)

;      xrange = [3*pixscale*arcsec2kpc,max(galphot[reffilt].radius_kpc)]
;      xrange = [pixscale*arcsec2kpc,max(galphot.radius_kpc)*1.5]
;      yrange = [max(modphot.sb0fit),min(modphot.sb0fit)]

       if ic gt 9 then begin
          delvarx, xtickname
       endif else begin
          xtickname = replicate(' ',10)
       endelse

       if (ic mod 5) eq 0 then begin
          delvarx, ytickname
       endif else begin
          ytickname = replicate(' ',10)
       endelse
       
       djs_plot, [0], [0], /nodata, position=pos[*,ic], noerase=ic gt 0, $
         xsty=1, ysty=1, yrange=yrange, xrange=xrange, /xlog, $
         xtitle=xtitle, ytitle=ytitle, ytickname=ytickname, xtickname=xtickname, $
         ytickinterval=4.0
       im_legend, strupcase(cluster), /right, /top, box=0, margin=0, charsize=1.0

       if ic eq 0 then begin
          im_legend, ['F475W','F814W','F160W'], /left, /bottom, box=0, $
            line=line, color=color2, pspacing=1.8, margin=-0.2, charsize=0.8, $
            charthick=3
       endif
       
;      for ii = 0, nfilt-1 do begin
;         modgood = where(modphot[ii].sb0fit gt 10^(-0.4*modphot[ii].sblimit))
;         djs_oplot, modphot[ii].radius_kpc[modgood], $
;           -2.5*alog10(modphot[ii].sb0fit[modgood])
;      endfor

       for ii = 0, nthese-1 do begin
          this = where(thesefilt[ii] eq strtrim(modphot.band,2))

;         if plotgal then begin
;            modgood = where(modphot[this].majora*pixscale*arcsec2kpc le sersic[this].amax_kpc and $
;              modphot[this].sb0fit gt 0 and modphot[this].sb0fit_ivar gt 0)
;            djs_oplot, modphot[this].radius_kpc[modgood], $
;              -2.5*alog10(modphot[this].sb0fit[modgood]), $
;              line=0, color=cgcolor('medium gray') ; color=cgcolor(color2[ii])
;            oploterror, modphot[this].radius_kpc[notzero], modphot[this].sb0fit[notzero], $
;              1.0/sqrt(modphot[this].sb0fit_ivar[notzero]), color=cgcolor('dodger blue'), $
;              psym=symcat(16), symsize=0.4
;         endif

          modgood = where(modphot[this].majora*pixscale*arcsec2kpc le sersic[this].amax_kpc and $
            modphot[this].sb0fit gt 0 and modphot[this].sb0fit_ivar gt 0)
;         modgood = where(modphot[this].sb0fit gt 10^(-0.4*modphot[this].sblimit) and $
;           modphot[this].sb0fit_ivar gt 0 and modphot[this].radius_kpc lt 100.0) ; r<100 kpc
          splog, thesefilt[ii], max(modphot[this].radius_kpc[modgood])

          rr = modphot[this].radius_kpc[modgood]
          sb = -2.5*alog10(modphot[this].sb0fit[modgood])
          sberr = 2.5/(modphot[this].sb0fit[modgood]*sqrt(modphot[this].sb0fit_ivar[modgood])*alog(10))

; just for the plot make the uncertainty have a minimum floor so that
; the SB profile shows up!          
          sberr = sberr>0.1

          polyfill, [rr,reverse(rr)], [sb-sberr,reverse(sb+sberr)], /fill, $
            color=cgcolor(color1[ii]), noclip=0
;         oploterror, rr, sb, sberr, color=cgcolor(color[ii]), line=line[ii]
;         djs_oplot, rr, sb, color=cgcolor(color[ii]), line=line[ii]
;         djs_oplot, 10^!x.crange, modphot[this].sblimit*[1,1], line=5, $
;           color=cgcolor(color[ii])

;; overplot the single Sersic function
;          rrser = [0,range(1E-3,200.0,500,/log)]
;;         rrser = [0,range(1E-3,max(rr),500,/log)]
;          djs_oplot, rrser, -2.5*alog10(bcgmstar_sersic_func(rrser,params=sersic[this])), $
;            color=cgcolor(color2[ii]), line=line[ii], thick=3

; overplot the Sersic+deVac function
          rrser = [0,range(1E-3,200.0,500,/log)]
;         rrser = [0,range(1E-3,max(rr),500,/log)]

          djs_oplot, rrser, -2.5*alog10(bcgmstar_sersic2_func(rrser,params=sersic[this],/allbands)), $
            color=cgcolor(color2[ii]), line=line[ii], thick=3
;         djs_oplot, rrser, -2.5*alog10(bcgmstar_sersic2_func(rrser,params=sersic[this])), $
;           color=cgcolor(color2[ii]), line=line[ii], thick=3
          
; overplot Marc's SB profile, as a consistency check
          if plotmarc then begin
             pp = read_bcg_profiles(cluster,these_filters=thesefilt[ii])
             if size(pp,/type) eq 8 then begin
                ww = where(pp.sma gt -90)
                djs_oplot, pp.sma[ww], pp.mu[ww], psym=symcat(15,thick=2), symsize=0.4
             endif
          endif
       endfor
;      cc = get_kbrd(1)
    endfor

    xyouts, min(pos[0,*])-0.05, (max(pos[3,*])-min(pos[1,*]))/2.0+min(pos[1,*]), $
      textoidl('\mu (mag arcsec^{-2})'), orientation=90, align=0.5, charsize=1.4, /norm
    xyouts, (max(pos[2,*])-min(pos[0,*]))/2.0+min(pos[0,*]), min(pos[1,*])-0.09, $
      textoidl('Equivalent Radius (kpc)'), $
;     textoidl('Equivalent Radius r=a\sqrt{1-\epsilon} (kpc)'), $
      align=0.5, charsize=1.4, /norm

    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
return
end
    
