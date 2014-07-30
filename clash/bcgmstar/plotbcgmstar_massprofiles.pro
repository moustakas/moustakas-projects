pro plotbcgmstar_massprofiles
; jm14jun18siena - plot the stellar mass profiles

    ellpath = bcgmstar_path(/ellipse)
    massprofpath = bcgmstar_path(/massprofiles)
    paperpath = bcgmstar_path(/paper)
    isedfit_dir = bcgmstar_path(/isedfit)

    prefix = 'bcgmstar'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

    sample = read_bcgmstar_sample()
    ncl = n_elements(sample)

; --------------------------------------------------    
; plot the total stellar mass vs the virial mass

; first gather the BCG stellar masses    
    mbcg = fltarr(ncl)
    mbcg_err = fltarr(ncl)
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       prof = mrdfits(massprofpath+cluster+'-massprofile.fits.gz',1,/silent)
       mbcg[ic] = prof.totmstar
       mbcg_err[ic] = prof.totmstar_err
    endfor

; make the plot    
    xrange = alog10([3D13,3D15])
    yrange = [11.3,13.2]
    
    psfile = paperpath+'mstar_vs_mvir.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, xmargin=[1.3,0.4], $
      width=6.8

    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, $
      xtitle='log_{10} (M_{500} / M_{\odot})', $
      ytitle='log_{10} (M_{*,BCG} / M_{\odot})'
      
; overplot the power-law fit from Kravtsov (Table 2)    
    scat = 0.17
    xx = range(xrange[0],xrange[1],50)
    yy = poly(xx-14.5,[12.24,0.33])
    polyfill, [xx,reverse(xx)], [yy+scat,reverse(yy-scat)], $
      /fill, color=cgcolor('light grey'), noclip=0
;   djs_oplot, xx, yy, thick=4

; overlay Kravtsov
    mbcg_krav = [3.12,4.14,3.06,1.47,0.79,1.26,1.09,0.91,1.38]*1D12
    mbcg_krav_err = [0.36,0.3,0.3,0.13,0.05,0.11,0.06,0.05,0.14]*1D12
    m500_krav = [15.6,10.3,7,5.34,2.35,1.86,1.34,0.46,0.47]*1D14

    oploterror, alog10(m500_krav), alog10(mbcg_krav), mbcg_krav_err/mbcg_krav/alog(10), $
      psym=symcat(15), color=cgcolor('tomato'), errcolor=cgcolor('tomato'), $
      symsize=1.5

; overplot Gonzalez+13
    mbcg_gonz = [0.68,0.82,0.35,0.74,0.56,0.35,0.46,0.88,0.64,0.7,0.65,0.35]*1D13
    mbcg_gonz_err = [0.04,0.06,0.03,0.06,0.05,0.02,0.02,0.06,0.05,0.06,0.04,0.02]*1D13
    m500_gonz = [2.26,5.15,0.95,3.46,3.59,0.99,0.95,3.23,2.26,2.41,2.37,1.45]*1D14
    m500_gonz_err = [0.19,0.42,0.1,0.32,0.28,0.11,0.1,0.19,0.23,0.18,0.24,0.21]*1D14
    
;    oploterror, alog10(m500_gonz), alog10(mbcg_gonz), mbcg_gonz_err/mbcg_gonz/alog(10), $
;      psym=symcat(14), color=cgcolor('forest green'), errcolor=cgcolor('forest green'), $
;      symsize=1.8

; overplot CLASH    
    oploterror, sample.m500, mbcg, sample.m500_err, mbcg_err, psym=symcat(16), symsize=1.6

; make a legend
    im_legend, ['CLASH','Kravtsov+14'], /right, /bottom, box=0, $
      psym=[16,15], color=['black','tomato'], spacing=2.5
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
; --------------------------------------------------    
; plot the mass profiles
    xrange = [0.6,100]
;   yrange = [9.8,11.8]
    yrange = [10,13]

    psfile = paperpath+'bcg_massprofiles.eps'
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.3
    pos = im_getposition(nx=5,ny=3,/landscape,yspace=0.1,xspace=0.1,$
      xmargin=[1.0,0.4])

    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)

       prof = mrdfits(massprofpath+cluster+'-massprofile.fits.gz',1,/silent)
       rad = prof.photradius_kpc
       splog, cluster, prof.totmstar, prof.totmstar_err

       good = where(prof.mstar gt 0)
       good1 = where(prof.mstar_grid01 gt 0)
       good2 = where(prof.mstar_grid02 gt 0)
       good3 = where(prof.mstar_grid03 gt 0)
       good4 = where(prof.mstar_grid04 gt 0)
       
;       intindx = 0
;       resindx = lindgen(n_elements(rad))+1 ; resolved masses (no integrated)
;       ised1 = read_isedfit(isedfit_paramfile,outprefix=prefix+'_'+cluster,$
;         isedfit_dir=isedfit_dir,thissfhgrid=1,/silent)
;       ised2 = read_isedfit(isedfit_paramfile,outprefix=prefix+'_'+cluster,$
;         isedfit_dir=isedfit_dir,thissfhgrid=2,/silent)
;       ised3 = read_isedfit(isedfit_paramfile,outprefix=prefix+'_'+cluster,$
;         isedfit_dir=isedfit_dir,thissfhgrid=3,/silent)
;       ised4 = read_isedfit(isedfit_paramfile,outprefix=prefix+'_'+cluster,$
;         isedfit_dir=isedfit_dir,thissfhgrid=4,/silent)
;
;       good1 = where(ised1[resindx].chi2 lt 1E6 and total(ised1[resindx].ivarmaggies gt 0,1) ge 3)
;       good2 = where(ised2[resindx].chi2 lt 1E6 and total(ised2[resindx].ivarmaggies gt 0,1) ge 3)
;       good3 = where(ised3[resindx].chi2 lt 1E6 and total(ised3[resindx].ivarmaggies gt 0,1) ge 3)
;       good4 = where(ised4[resindx].chi2 lt 1E6 and total(ised4[resindx].ivarmaggies gt 0,1) ge 3)

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

       oploterror, rad[good], prof.mstar[good], prof.mstar_err[good], $
        psym=-symcat(16), color=cgcolor('orchid'), errcolor=cgcolor('orchid'), $
        symsize=0.6
       djs_oplot, rad[good1], prof.mstar_grid01[good1], line=0, color=cgcolor('black')
       djs_oplot, rad[good2], prof.mstar_grid02[good2], line=3, color=cgcolor('forest green')
       djs_oplot, rad[good3], prof.mstar_grid03[good3], line=5, color=cgcolor('dodger blue')
       djs_oplot, rad[good4], prof.mstar_grid04[good4], line=4, color=cgcolor('firebrick')

;      splog, ised1[intindx].mstar_50, ised2[intindx].mstar_50, ised3[intindx].mstar_50, $
;        ised4[intindx].mstar_50
       djs_oplot, 10^!x.crange, prof.totmstar*[1,1], line=0, color='grey'
       
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
    
