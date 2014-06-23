pro plotbcgmstar_massprofiles
; jm14jun18siena - plot the stellar mass profiles

    ellpath = bcgmstar_path(/ellipse)
    sersicpath = bcgmstar_path(/sersic)
    paperpath = bcgmstar_path(/paper)
    isedfit_dir = bcgmstar_path(/isedfit)

    prefix = 'bcgmstar'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

    sample = read_bcgmstar_sample()
    ncl = n_elements(sample)

    pixscale = 0.065
    
; make the plot
    xrange = [0.6,100]
    yrange = [9.8,11.8]
;   yrange = [10,12.2]

    psfile = paperpath+'bcg_massprofiles.eps'
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.3
    pos = im_getposition(nx=5,ny=3,/landscape,yspace=0.1,xspace=0.1,$
      xmargin=[1.0,0.4])

    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       splog, cluster

       phot = mrdfits(sersicpath+cluster+'-phot.fits.gz',1,/silent)
       rad = phot[0].photradius_kpc
       intindx = 0
       resindx = lindgen(n_elements(rad))+1 ; resolved masses (no integrated)
       ised1 = read_isedfit(isedfit_paramfile,outprefix=prefix+'_'+cluster,$
         isedfit_dir=isedfit_dir,thissfhgrid=1,/silent)
       ised2 = read_isedfit(isedfit_paramfile,outprefix=prefix+'_'+cluster,$
         isedfit_dir=isedfit_dir,thissfhgrid=2,/silent)
       ised3 = read_isedfit(isedfit_paramfile,outprefix=prefix+'_'+cluster,$
         isedfit_dir=isedfit_dir,thissfhgrid=3,/silent)
       ised4 = read_isedfit(isedfit_paramfile,outprefix=prefix+'_'+cluster,$
         isedfit_dir=isedfit_dir,thissfhgrid=4,/silent)

       good1 = where(ised1[resindx].chi2 lt 1E6 and total(ised1[resindx].ivarmaggies gt 0,1) ge 3)
       good2 = where(ised2[resindx].chi2 lt 1E6 and total(ised2[resindx].ivarmaggies gt 0,1) ge 3)
       good3 = where(ised3[resindx].chi2 lt 1E6 and total(ised3[resindx].ivarmaggies gt 0,1) ge 3)
       good4 = where(ised4[resindx].chi2 lt 1E6 and total(ised4[resindx].ivarmaggies gt 0,1) ge 3)

;      good1 = where(ised1.chi2 lt 1E6 and total(ised1.ivarmaggies gt 0,1) ge 5) ; minimum of 5 bands 
;      good2 = where(ised2.chi2 lt 1E6 and total(ised2.ivarmaggies gt 0,1) ge 5) ; minimum of 5 bands 
;      good3 = where(ised3.chi2 lt 1E6 and total(ised3.ivarmaggies gt 0,1) ge 5) ; minimum of 5 bands 
;      good4 = where(ised4.chi2 lt 1E6 and total(ised4.ivarmaggies gt 0,1) ge 5) ; minimum of 5 bands 

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
         ytickinterval=1
       im_legend, strupcase(cluster), /left, /top, box=0, margin=0, charsize=1.0

       oploterror, rad[good4], ised4[resindx[good4]].mstar_50, ised4[resindx[good4]].mstar_err, line=4, $
         psym=-symcat(16), color=cgcolor('firebrick'), errcolor=cgcolor('firebrick'), $
         symsize=0.8
       djs_oplot, rad[good1], ised1[resindx[good1]].mstar_50, line=0, color=cgcolor('black')
       djs_oplot, rad[good2], ised2[resindx[good2]].mstar_50, line=3, color=cgcolor('forest green')
       djs_oplot, rad[good3], ised3[resindx[good3]].mstar_50, line=5, color=cgcolor('dodger blue')
;      djs_oplot, rad[good4], ised4[resindx[good4]].mstar_50, line=4, color=cgcolor('firebrick')

       splog, ised1[intindx].mstar_50, ised2[intindx].mstar_50, ised3[intindx].mstar_50, $
         ised4[intindx].mstar_50
;      djs_oplot, 10^!x.crange, ised1[intindx].mstar_50*[1,1], line=0, color='grey'
       
;      cumumass = alog10(total(10.0^ised4[good4].mstar_50,/cumu))
;      polyfill, [rad[good],reverse(rad[good])], [cumumass+0.06,reverse(cumumass-0.06)], $
;        /fill, color=cgcolor('dodger blue')
    endfor

    xyouts, min(pos[0,*])-0.05, (max(pos[3,*])-min(pos[1,*]))/2.0+min(pos[1,*]), $
      textoidl('log(Stellar Mass) (M_{\odot})'), orientation=90, align=0.5, charsize=1.4, /norm
;   xyouts, min(pos[0,*])-0.05, (max(pos[3,*])-min(pos[1,*]))/2.0+min(pos[1,*]), $
;     textoidl('log(Cumulative Stellar Mass) (M_{\odot})'), orientation=90, align=0.5, charsize=1.4, /norm
    xyouts, (max(pos[2,*])-min(pos[0,*]))/2.0+min(pos[0,*]), min(pos[1,*])-0.09, $
      textoidl('Equivalent Radius (kpc)'), $
;     textoidl('Equivalent Radius r=a\sqrt{1-\epsilon} (kpc)'), $
      align=0.5, charsize=1.4, /norm

    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
return
end
    
