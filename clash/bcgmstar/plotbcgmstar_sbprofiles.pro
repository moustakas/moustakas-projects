pro plotbcgmstar_sbprofiles, pdf=pdf
; jm13oct19siena - plot the SB profiles

    paperpath = bcgmstar_path(/paper)
    ellpath = bcgmstar_path(/ellipse)
    sersicpath = bcgmstar_path(/sersic)

    sample = read_bcgmstar_sample(/zsort)
    ncl = n_elements(sample)

    pixscale = 0.065
    errfloor = 0.0D ; 0.02      ; magnitude error floor on my SB measurements 
    ncol = 6 ; number of columns
    
    nicecluster = repstr(repstr(strupcase(sample.dirname),'ABELL','Abell'),'_',' ')

; ---------------------------------------------------------------------------
; make 15 figures for the appendix which shows all the clusters and
; all the bands (see bcgmstar_sersicfit, /qaplot_sbprofiles
    modelrr = [0,range(0.01,200,500,/log)] ; equivalent radius [kpc]
    xrange = [0.3,120]
    
;   for ic = 1, 1 do begin
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]

       sersic = read_bcgmstar_sersic(cluster,results=sersic_results,$
         radius=modelrr,model=modelsb,firstsersic=firstsersic,$
         secondsersic=secondsersic)
       modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
       nfilt = n_elements(modphot)

       psfile = paperpath+'bcg_all_sbprofiles_'+cluster+'.eps'
       im_plotconfig, 1, pos, psfile=psfile, charsize=1.3

       nrow = 2*ceil(nfilt/float(ncol))
;      if allsersic then nrow = ceil((nfilt+1)/float(ncol)) else $
;        nrow = ceil(nfilt/float(ncol))
       ym = [0.4,1.1]
       xm = [1.4,0.4]
       wid = replicate(2.6,ncol)
       hei = reform(cmreplicate([1.2,0.8],nrow/2),nrow)
       xsp = replicate(0.05,ncol-1)
       ysp = [reform(cmreplicate([0.0,0.05],(nrow-2)/2),nrow-2),0.0]
       
       xpage = total(xm)+total(wid)+total(xsp)
       ypage = total(ym)+total(hei)+total(ysp)
;      ypage = total(ym)+total(replicate(hei,nrow))

       pos = im_getposition(nx=ncol,ny=nrow,yspace=ysp,xspace=xsp,$
         xmargin=xm,ymargin=ym,width=wid,height=hei,xpage=xpage,$
         ypage=ypage,/landscape)
       pos = reform(pos,4,ncol,nrow)
       plotpos = reform(pos[*,*,2*lindgen(nrow/2)],4,ncol*nrow/2)
       residpos = reform(pos[*,*,2*lindgen(nrow/2)+1],4,ncol*nrow/2)

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
          if count ge nfilt-ncol then begin
             delvarx, xtickname
          endif else begin
             xtickname = replicate(' ',10)
          endelse
          if (count mod ncol) eq 0 then begin
             delvarx, ytickname
             plottitle = '\mu' ; 
             residtitle = '\Delta'+'\mu'
          endif else begin
             ytickname = replicate(' ',10)
             delvarx, plottitle, residtitle
          endelse
          
          djs_plot, [0], [0], /nodata, /xlog, noerase=count gt 0, $
            xrange=xrange, xsty=1, yrange=[26,16], position=plotpos[*,count], $
            xtickname=replicate(' ',10), ytickname=ytickname, $ ; title=title, $
            symsize=0.5, ytickinterval=3, ysty=1, ytitle=plottitle

;         djs_oplot, rr, sb, color=cgcolor('dark grey'), thick=10
;         djs_oplot, rr, sb, psym=symcat(15), symsize=0.5
          polyfill, [rr,reverse(rr)], [sb-sberr,reverse(sb+sberr)], /fill, $
            color=cgcolor('tomato'), noclip=0
;         djs_oplot, rr, sb+sberr, line=0, color=cgcolor('dark grey')
;         djs_oplot, rr, sb-sberr, line=0, color=cgcolor('dark grey')

;         if count eq 0 then im_legend, nicecluster[ic], /left, $
;           /bottom, box=0, charsize=1.5, margin=0;, $
;           position=[pos[2,count]+0.005,pos[3,count]-0.02], /norm
;         if count eq 0 then im_legend, nicecluster[ic], /left, /bottom, box=0, charsize=1.5, margin=0
          
;          if modbad[0] ne -1 then begin
;             djs_oplot, modphot[ib].radius_kpc[modbad], -2.5*alog10(modphot[ib].sb0fit[modbad]), $
;               psym=symcat(9), color=cgcolor('medium grey'), symsize=0.5
;          endif

; overplot the model          
          djs_oplot, modelrr, modelsb[*,ib], line=0
          if cluster eq 'a209' or cluster eq 'a2261' or cluster eq 'rxj2248' then begin
             djs_oplot, modelrr, firstsersic[*,ib], line=5
             djs_oplot, modelrr, secondsersic[*,ib], line=5
          endif

          im_legend, band, /right, /top, box=0, margin=0, charsize=1.1;, charthick=1.8
;         im_legend, band, /right, /top, box=0, margin=0, charsize=1.2;, charthick=1.8
                    
          djs_oplot, 10^!x.crange, (modphot[ib].sblimit-2.5*alog10(2.0))*[1,1], line=1 ; 3-sigma
;         djs_oplot, [70.0,10^!x.crange[1]], modphot[ib].sblimit*[1,1], line=0

; residuals
          resid = sb-interpol(modelsb[*,ib],modelrr,rr)
          djs_plot, [0], [0], /nodata, /xlog, /noerase, $
            xrange=xrange, xsty=1, yrange=0.8*[-1,1], position=residpos[*,count], $
            xtickname=xtickname, ytickname=ytickname, $ ; title=title, $
            ytickinterval=0.5, ysty=1, ytitle=residtitle

          polyfill, [rr,reverse(rr)], [resid-sberr,reverse(resid+sberr)], /fill, $
            color=cgcolor('tomato'), noclip=0
;         djs_oplot, rr, resid+sberr, line=0, color=cgcolor('dark grey')
;         djs_oplot, rr, resid-sberr, line=0, color=cgcolor('dark grey')
;         djs_oplot, rr, resid, psym=symcat(16), symsize=0.3

          djs_oplot, 10^!x.crange, [0,0], line=0, thick=2
          
          count++
       endfor 

;       xyouts, min(pos[0,*])-0.06, (max(pos[3,*])-min(pos[1,*]))/2.0+min(pos[1,*]), $
;         textoidl('\mu (mag arcsec^{-2})'), orientation=90, align=0.5, charsize=1.6, /norm
       xyouts, (max(residpos[2,*])-min(residpos[0,*]))/2.0+min(residpos[0,*]), $
         min(residpos[1,*])-0.1, textoidl('Galactocentric Radius (kpc)'), $
         align=0.5, charsize=1.6, /norm

       xyouts, (max(residpos[2,*])-min(residpos[0,*]))/2.0+min(residpos[0,*]), $
         max(plotpos[3,*])+0.02, nicecluster[ic], $
         align=0.5, charsize=1.6, /norm

       im_plotconfig, psfile=psfile, /psclose, /pdf
;stop
    endfor 

stop    
    
; ---------------------------------------------------------------------------
; make a plot for the main body of the paper which shows the surface
; brightness profiles and Sersic fits for three representative
; bandpasses     
    
; plot these bands    
    thesefilt = ['f475w','f814w','f160w']
    thesefilt2 = ['f555w','f814w','f160w']
    color1 = ['dark grey','dodger blue','tomato']
;   color2 = ['black','black','black']
    color2 = ['black','navy','firebrick']
    line = [2,1,0]

;   thesefilt = ['f160w','f814w','f475w']
;   color1 = ['tomato','powder blue','medium grey']
;   color2 = ['firebrick','navy','black']

;   color = ['orange','firebrick','navy']
;   color2 = ['orange red','orange','dodger blue']
    nthese = n_elements(thesefilt)

; overplot the galaxy photometry?  overplot Marc's SB profiles?
    plotgal = 0
    plotmarc = 0
    
    modelrr = [0,range(1E-3,200.0,500,/log)]
;   modelrr = [0,range(1E-3,max(rr),500,/log)]

; make the plot
;   xtitle = 'Semi-Major Axis (kpc)'    
;   ytitle = '\mu (mag arcsec^{-2})'
    xrange = [0.2,120]
;   xrange = [3*pixscale*arcsec2kpc,150]
    yrange = [26,16]

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
       sersic = read_bcgmstar_sersic(cluster,results=sersic_results,radius=modelrr,model=modelsb)

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
          im_legend, ['F160W','F814W','F475W or!c F555W'], /left, /bottom, box=0, $
            line=reverse(line), color=reverse(color1), pspacing=1.7, $ ; margin=-0.2, $
            charsize=0.7, charthick=2, position=[pos[0,ic]+0.003,pos[1,ic]+0.025], $
            /norm, thick=6
       endif
       
;      for ii = 0, nfilt-1 do begin
;         modgood = where(modphot[ii].sb0fit gt 10^(-0.4*modphot[ii].sblimit))
;         djs_oplot, modphot[ii].radius_kpc[modgood], $
;           -2.5*alog10(modphot[ii].sb0fit[modgood])
;      endfor

       for ii = 0, nthese-1 do begin
          case cluster of
             'macs0744': this = where(thesefilt2[ii] eq strtrim(modphot.band,2))
             'macs1149': this = where(thesefilt2[ii] eq strtrim(modphot.band,2))
             'macs2129': this = where(thesefilt2[ii] eq strtrim(modphot.band,2))
             else: this = where(thesefilt[ii] eq strtrim(modphot.band,2))
          endcase

          modgood = where(modphot[this].majora*pixscale*arcsec2kpc le sersic[this].amax_kpc and $
            modphot[this].sb0fit gt 0 and modphot[this].sb0fit_ivar gt 0)
;         modgood = where(modphot[this].sb0fit gt 10^(-0.4*modphot[this].sblimit) and $
;           modphot[this].sb0fit_ivar gt 0 and modphot[this].radius_kpc lt 100.0) ; r<100 kpc
          splog, thesefilt[ii], max(modphot[this].radius_kpc[modgood])

          rr = modphot[this].radius_kpc[modgood]
          sb = -2.5*alog10(modphot[this].sb0fit[modgood])
          sberr = 2.5/(modphot[this].sb0fit[modgood]*sqrt(modphot[this].sb0fit_ivar[modgood])*alog(10))

; overplot the model fit
          djs_oplot, modelrr, modelsb[*,this], color=cgcolor(color1[ii]), line=line[ii], thick=3
;         djs_oplot, modelrr, -2.5*alog10(bcgmstar_sersic2_func(modelrr,params=sersic[this])), $
;           color=cgcolor(color2[ii]), line=line[ii], thick=3
          
; just for the plot make the uncertainty have a minimum floor so that
; the SB profile shows up!          
          sberr = sberr>0.1

          polyfill, [rr,reverse(rr)], [sb-sberr,reverse(sb+sberr)], /fill, $
            color=cgcolor(color1[ii]), noclip=0
;         oploterror, rr, sb, sberr, color=cgcolor(color[ii]), line=line[ii]
;         djs_oplot, rr, sb, color=cgcolor(color[ii]), line=line[ii]

;          djs_oplot, 10^!x.crange, (modphot[this].sblimit-2.5*alog10(2.0))*[1,1], $
;            line=2, thick=2, color=cgcolor(color1[ii])     ; 3-sigma

;; overplot the single Sersic function
;          modelrr = [0,range(1E-3,200.0,500,/log)]
;;         modelrr = [0,range(1E-3,max(rr),500,/log)]
;          djs_oplot, modelrr, -2.5*alog10(bcgmstar_sersic_func(modelrr,params=sersic[this])), $
;            color=cgcolor(color2[ii]), line=line[ii], thick=3

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
      textoidl('Galactocentric Radius (kpc)'), $
;     textoidl('Equivalent Radius r=a\sqrt{1-\epsilon} (kpc)'), $
      align=0.5, charsize=1.4, /norm

    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
return
end
    
