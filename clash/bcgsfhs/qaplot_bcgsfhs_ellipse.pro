pro qaplot_bcgsfhs_ellipse
; jm13oct22siena - generate a QAplot for the output of BCGSFHS_ELLIPSE

    paperpath = bcgsfhs_path()+'qaplots/'
    suffix = '.ps'

    sample = read_bcgsfhs_sample(/noa2261)
    ncl = n_elements(sample)

    pixscale = 0.065
    
; make the plot
;   xtitle = 'Semi-Major Axis (kpc)'    
;   ytitle = '\mu (mag arcsec^{-2})'

;   yrange = [30,16]
    xrange = [0.5,150]
;   xrange = [3*pixscale*arcsec2kpc,150]

    for ic = 4, 4 do begin
;   for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       print & splog, cluster, sample[ic].z
       datapath = bcgsfhs_path(/bcg)+cluster+'/'

       psfile = paperpath+'qa_ellipse_'+cluster+suffix
       im_plotconfig, 0, pos, psfile=psfile, charsize=1.3
       pos = im_getposition(nx=1,ny=5,yspace=0.0,xmargin=[1.3,0.4],width=6.0,xpage=7.7)
       
       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
       
       galphot = mrdfits(datapath+cluster+'-ellipse-image.fits.gz',1,/silent)
       modphot = mrdfits(datapath+cluster+'-ellipse-model.fits.gz',1,/silent)
       nfilt = n_elements(modphot)

       reffilt = where(strtrim(galphot.band,2) eq 'f160w') ; reference filter

       for ii = 0, nfilt-1 do begin
          band = strtrim(strupcase(galphot[ii].band),2)
          
          galgood = where(galphot[ii].sb0fit gt 0)
          modgood = where(modphot[ii].sb0fit gt 10^(-0.4*modphot[ii].sblimit))
          modgoodref = where(modphot[0].sb0fit gt 10^(-0.4*modphot[0].sblimit))

;         if ii eq 11 then stop

; SB profile (data + model + Sersic fit)
          yrange = [min(modphot[ii].sb0fit[modgood])<min(modphot[0].sb0fit[modgoodref])<$
            min(galphot[ii].sb0fit[galgood]),$
            max(modphot[ii].sb0fit[modgood])>max(modphot[0].sb0fit[modgoodref])>$
            max(galphot[ii].sb0fit[galgood])]
          yrange = -2.5*alog10(yrange)
          
          djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=3, $
            yrange=yrange, xrange=xrange, /xlog, xtickname=replicate(' ',10), $
            title=strupcase(cluster)+'/'+band, yminor=3, ytickinterval=3, $
            ytitle='\mu (mag arcsec^{-2})'

          djs_oplot, galphot[ii].radius_kpc[galgood], $
            -2.5*alog10(galphot[ii].sb0fit[galgood]), $
            psym=symcat(16,thick=2), symsize=0.5, color=cgcolor('grey');, thick=6

          if ii eq 0 then begin
             line = 5
             color = 'firebrick'
             im_legend, 'F160W', /right, /top, box=0, margin=0, line=5, $
               pspacing=1.7, color=cgcolor(color)
          endif else begin
             djs_oplot, modphot[0].radius_kpc[modgood], $
               -2.5*alog10(modphot[0].sb0fit[modgood]), $
               color=cgcolor('firebrick'), line=5
             line = 0
             color = 'navy'
             im_legend, ['F160W',band], /right, $
               /top, box=0, margin=0, line=[5,0], pspacing=1.7, $
               color=cgcolor(['firebrick','navy'])
          endelse

          djs_oplot, modphot[ii].radius_kpc[modgood], $
            -2.5*alog10(modphot[ii].sb0fit[modgood]), $
            color=cgcolor(color), line=line
          
; ellipticity (model)
          djs_plot, [0], [0], /nodata, position=pos[*,1], /noerase, xsty=1, ysty=1, $
            yrange=[0.0,0.7], xrange=xrange, /xlog, xtickname=replicate(' ',10), $
            ytickinterval=0.2, ytitle='\epsilon', yminor=2

          if ii eq 0 then begin
             line = 5
             color = 'firebrick'
          endif else begin
             djs_oplot, modphot[0].radius_kpc[modgood], modphot[0].ellipticityfit[modgood], $
               color=cgcolor('firebrick'), line=5
             line = 0
             color = 'navy'
          endelse

          djs_oplot, modphot[ii].radius_kpc[modgood], modphot[ii].ellipticityfit[modgood], $
            color=cgcolor(color), line=line

; position angle (model)
          djs_plot, [0], [0], /nodata, position=pos[*,2], /noerase, xsty=1, ysty=1, $
            yrange=[0.0,180], xrange=xrange, /xlog, xtickname=replicate(' ',10), $
            ytitle='PA (degree)', yminor=2

          if ii eq 0 then begin
             line = 5
             color = 'firebrick'
          endif else begin
             djs_oplot, modphot[0].radius_kpc[modgood], modphot[0].pafit[modgood]*!radeg+90, $
               color=cgcolor('firebrick'), line=5
             line = 0
             color = 'navy'
          endelse

          djs_oplot, modphot[ii].radius_kpc[modgood], modphot[ii].pafit[modgood]*!radeg+90, $
            color=cgcolor(color), line=line

; a3*100 (model)
          djs_plot, [0], [0], /nodata, position=pos[*,3], /noerase, xsty=1, ysty=1, $
            yrange=[-0.15,0.15], xrange=xrange, /xlog, xtickname=replicate(' ',10), $
            ytitle='100 a_{3}/a', ytickinterval=0.1, yminor=2

          if ii eq 0 then begin
             line = 5
             color = 'firebrick'
          endif else begin
             djs_oplot, modphot[0].radius_kpc[modgood], 100*modphot[0].a3fit[modgood]/modphot[0].majora, $
               color=cgcolor('firebrick'), line=5
             line = 0
             color = 'navy'
          endelse

          djs_oplot, modphot[ii].radius_kpc[modgood], 100*modphot[ii].a3fit[modgood]/modphot[ii].majora, $
            color=cgcolor(color), line=line

; a4*100 (model)
          djs_plot, [0], [0], /nodata, position=pos[*,4], /noerase, xsty=1, ysty=1, $
            yrange=[-0.15,0.15], xrange=xrange, /xlog, xtitle='Semi-Major Axis (kpc)', $
            ytitle='100 a_{4}/a', ytickinterval=0.1, yminor=2

          if ii eq 0 then begin
             line = 5
             color = 'firebrick'
          endif else begin
             djs_oplot, modphot[0].radius_kpc[modgood], 100*modphot[0].a4fit[modgood]/modphot[0].majora, $
               color=cgcolor('firebrick'), line=5
             line = 0
             color = 'navy'
          endelse

          djs_oplot, modphot[ii].radius_kpc[modgood], 100*modphot[ii].a4fit[modgood]/modphot[ii].majora, $
            color=cgcolor(color), line=line

       endfor
       im_plotconfig, psfile=psfile, /psclose, /pdf

    endfor 

stop    
    
return
end
    
