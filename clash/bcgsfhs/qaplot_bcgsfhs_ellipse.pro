pro qaplot_bcgsfhs_ellipse
; jm13oct22siena - generate a QAplot for the output of BCGSFHS_ELLIPSE

    ellpath = bcgsfhs_path(/ellipse)
    qapath = ellpath+'qaplots/'

    sample = read_bcgsfhs_sample()
    ncl = n_elements(sample)

    pixscale = 0.065
    
; make the plot
;   xtitle = 'Semi-Major Axis (kpc)'    
;   ytitle = '\mu (mag arcsec^{-2})'

;   yrange = [30,16]
    xrange = [0.3,300]
;   xrange = [3*pixscale*arcsec2kpc,150]

    galcolor = 'dodger blue' ; 'dark gray'
    galline = 0
    modline = 0
    modcolor = 'navy'
    refmodcolor = 'firebrick'
    refmodline = 5

    showa3 = 0 ; optionally plot the a3, a4 coefficients

; overplot Marc's SB profiles?
    plotmarc = 1
    
;   for ic = 0, 6 do begin
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       print & splog, cluster, sample[ic].z

       psfile = qapath+'qa_ellipse_'+cluster+'.ps'
       im_plotconfig, 0, pos, psfile=psfile, charsize=1.3

       if showa3 then begin
          pos = im_getposition(nx=1,ny=5,yspace=0.0,xmargin=[1.3,2.2],width=5.0)
       endif else begin
          pos = im_getposition(nx=1,ny=3,yspace=0.0,$
            xmargin=[1.3,2.2],width=5.0)
       endelse
       
       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
       
       galphot = mrdfits(ellpath+cluster+'-ellipse-image.fits.gz',1,/silent)
       modphot = mrdfits(ellpath+cluster+'-ellipse-model.fits.gz',1,/silent)
       nfilt = n_elements(modphot)

; read Marc's stuff
       pp = read_bcg_profiles(cluster,these_filters=strtrim(modphot.band,2))
       if n_elements(pp) ne nfilt then message, 'Missing SB profiles!'
          
       for ii = 0, nfilt-1 do begin
          band = strtrim(strupcase(galphot[ii].band),2)
          
          galgood = where(galphot[ii].sb0fit gt 0)
          modgood = where(modphot[ii].sb0fit gt 0)
          refmodgood = where(modphot[0].sb0fit gt 0)

          amax = max(pp[ii].sma)
          splog, '   ', band, amax

;         galgood = where(galphot[ii].sb0fit gt 0)
;         modgood = where(modphot[ii].sb0fit gt 10^(-0.4*modphot[ii].sblimit))
;         refmodgood = where(modphot[0].sb0fit gt 10^(-0.4*modphot[0].sblimit))

          sbmin = 10^(-0.4*modphot[ii].sblimit)
;         rmax = interpolate(modphot[ii].radius_kpc[modgood],findex($
;           modphot[ii].sb0fit[modgood],sbmin))

; SB profile (data + model + Sersic fit)
          yrange = [min(modphot[ii].sb0fit[modgood])<min(modphot[0].sb0fit[refmodgood])<$
            min(galphot[ii].sb0fit[galgood]),$
            max(modphot[ii].sb0fit[modgood])>max(modphot[0].sb0fit[refmodgood])>$
            max(galphot[ii].sb0fit[galgood])]
          yrange = -2.5*alog10(yrange)

          djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=3, $
            yrange=yrange, xrange=xrange, /xlog, xtickname=replicate(' ',10), $
            title=strupcase(cluster)+'/'+band, yminor=3, ytickinterval=3, $
            ytitle='\mu (mag arcsec^{-2})'
          djs_oplot, amax*[1,1], !y.crange, line=0, color=cgcolor('grey')
          djs_oplot, [3.0,10^!x.crange[1]], -2.5*alog10(sbmin)*[1,1], $
            line=0, color=cgcolor('grey')

          if ii gt 0 then begin
             djs_oplot, modphot[0].majora[refmodgood]*pixscale*arcsec2kpc, $
               -2.5*alog10(modphot[0].sb0fit[refmodgood]), $
               color=cgcolor(refmodcolor), line=refmodline
             im_legend, ['F160W/Model',band+'/Model',band+'/Data'], /left, $
               /bottom, box=0, margin=0, line=[refmodline,modline,galline], pspacing=1.7, $
               color=cgcolor([refmodcolor,modcolor,galcolor]), thick=8, charsize=1.0
          endif else begin
             im_legend, ['F160W/Model','F160W/Data'], /left, /bottom, box=0, $
               margin=0, line=[modline,galline], pspacing=1.7, $
               color=cgcolor([modcolor,galcolor]), thick=8, charsize=1.0
          endelse

          djs_oplot, galphot[ii].majora[galgood]*pixscale*arcsec2kpc, $
            -2.5*alog10(galphot[ii].sb0fit[galgood]), $
            color=cgcolor(galcolor), line=galline, psym=-symcat(16,thick=2), symsize=0.5
          djs_oplot, modphot[ii].majora[modgood]*pixscale*arcsec2kpc, $
            -2.5*alog10(modphot[ii].sb0fit[modgood]), $
            color=cgcolor(modcolor), line=modline

; overplot Marc's measurements, as a consistency check
          if plotmarc then begin
             ww = where(pp[ii].sma gt -90)
             djs_oplot, pp[ii].sma[ww], pp[ii].mu[ww], color=cgcolor('orange'), $
               psym=symcat(15,thick=2), symsize=0.4
          endif
          
; ellipticity (model)
          djs_plot, [0], [0], /nodata, position=pos[*,1], /noerase, xsty=1, ysty=1, $
            yrange=[0.0,0.7], xrange=xrange, /xlog, xtickname=replicate(' ',10), $
            ytickinterval=0.2, ytitle='\epsilon', yminor=2
          djs_oplot, amax*[1,1], !y.crange, line=0, color=cgcolor('grey')

          if ii gt 0 then begin
             djs_oplot, modphot[0].majora[refmodgood]*pixscale*arcsec2kpc, $
               modphot[0].ellipticityfit[refmodgood], $
               color=cgcolor(refmodcolor), line=refmodline
          endif

          djs_oplot, galphot[ii].majora[galgood]*pixscale*arcsec2kpc, galphot[ii].ellipticityfit[galgood], $
            color=cgcolor(galcolor), line=galline, psym=-symcat(16,thick=2), symsize=0.5
          djs_oplot, modphot[ii].majora[modgood]*pixscale*arcsec2kpc, modphot[ii].ellipticityfit[modgood], $
            color=cgcolor(modcolor), line=modline

; overplot Marc's measurements, as a consistency check
          if plotmarc then begin
             ww = where(pp[ii].sma gt -90)
             djs_oplot, pp[ii].sma[ww], pp[ii].ell[ww], color=cgcolor('orange'), $
               psym=symcat(15,thick=2), symsize=0.4
          endif
          
; position angle (model)
          if showa3 then begin
             xtitle = ''
             xtickname = replicate(' ',10)
          endif else begin
             xtitle = 'Semi-Major Axis (kpc)'
             delvarx, xtickname
          endelse

          djs_plot, [0], [0], /nodata, position=pos[*,2], /noerase, xsty=1, ysty=1, $
            yrange=[0.0,180], xrange=xrange, /xlog, xtickname=xtickname, xtitle=xtitle, $
            ytitle='PA (degree)', yminor=2
          djs_oplot, amax*[1,1], !y.crange, line=0, color=cgcolor('grey')

          if ii gt 0 then begin
             djs_oplot, modphot[0].majora[refmodgood]*pixscale*arcsec2kpc, $
               modphot[0].pafit[refmodgood]*!radeg+90, $
               color=cgcolor(refmodcolor), line=refmodline
          endif

          djs_oplot, galphot[ii].majora[galgood]*pixscale*arcsec2kpc, galphot[ii].pafit[galgood]*!radeg+90, $
            color=cgcolor(galcolor), line=galline, psym=-symcat(16,thick=2), symsize=0.5
          djs_oplot, modphot[ii].majora[modgood]*pixscale*arcsec2kpc, modphot[ii].pafit[modgood]*!radeg+90, $
            color=cgcolor(modcolor), line=modline

; overplot Marc's measurements, as a consistency check
          if plotmarc then begin
             ww = where(pp[ii].sma gt -90)
             djs_oplot, pp[ii].sma[ww], 180.0-pp[ii].pa[ww], color=cgcolor('orange'), $
               psym=symcat(15,thick=2), symsize=0.4
          endif
          
          if showa3 then begin
; a3*100 (model)
             djs_plot, [0], [0], /nodata, position=pos[*,3], /noerase, xsty=1, ysty=1, $
               yrange=[-0.15,0.15], xrange=xrange, /xlog, xtickname=replicate(' ',10), $
               ytitle='100 a_{3}/a', ytickinterval=0.1, yminor=2
             djs_oplot, amax*[1,1], !y.crange, line=0, color=cgcolor('grey')
             
             if ii gt 0 then begin
                djs_oplot, modphot[0].majora[refmodgood]*pixscale*arcsec2kpc, $
                  100*modphot[0].a3fit[refmodgood]/modphot[0].majora[refmodgood], $
                  color=cgcolor(refmodcolor), line=refmodline
             endif
             
             djs_oplot, galphot[ii].majora[galgood]*pixscale*arcsec2kpc, $
               100*galphot[ii].a3fit[galgood]/galphot[ii].majora[galgood], $
               color=cgcolor(galcolor), line=galline, psym=-symcat(16,thick=2), symsize=0.5
             djs_oplot, modphot[ii].majora[modgood]*pixscale*arcsec2kpc, $
               100*modphot[ii].a3fit[modgood]/modphot[ii].majora[modgood], $
               color=cgcolor(modcolor), line=modline
             
; a4*100 (model)
             djs_plot, [0], [0], /nodata, position=pos[*,4], /noerase, xsty=1, ysty=1, $
               yrange=[-0.15,0.15], xrange=xrange, /xlog, xtitle='Semi-Major Axis (kpc)', $
               ytitle='100 a_{4}/a', ytickinterval=0.1, yminor=2
             djs_oplot, amax*[1,1], !y.crange, line=0, color=cgcolor('grey')

             if ii gt 0 then begin
                djs_oplot, modphot[0].majora[refmodgood]*pixscale*arcsec2kpc, $
                  100*modphot[0].a4fit[refmodgood]/modphot[0].majora[refmodgood], $
                  color=cgcolor(refmodcolor), line=refmodline
             endif
             
             djs_oplot, galphot[ii].majora[galgood]*pixscale*arcsec2kpc, $
               100*galphot[ii].a4fit[galgood]/galphot[ii].majora[galgood], $
               color=cgcolor(galcolor), line=galline, psym=-symcat(16,thick=2), symsize=0.5
             djs_oplot, modphot[ii].majora[modgood]*pixscale*arcsec2kpc, $
               100*modphot[ii].a4fit[modgood]/modphot[ii].majora[modgood], $
               color=cgcolor(modcolor), line=modline
          endif

       endfor
       im_plotconfig, psfile=psfile, /psclose, /pdf

    endfor 

stop    
    
return
end
    
