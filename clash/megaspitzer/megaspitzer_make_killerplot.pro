pro megaspitzer_make_killerplot, clobber=clobber
; jm12jul29siena - make the killer plot
    
    path = clash_path(/megaspitzer)
    prefix = 'megaspitzer'

; restore the model fits and posteriors
    isedpath = clash_path(/megaspitzer)+'isedfit/'
    isedfit_sfhgrid_dir = clash_path(/megaspitzer)+'montegrids/'
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/megaspitzer/megaspitzer_sfhgrid.par'
    paramfile = isedpath+prefix+'_supergrid01_isedfit.par'

    cat = read_megaspitzer()
    megaspitzer_to_maggies, cat, refmaggies, refivar
    refmag = maggies2mag(refmaggies,ivar=refivar,magerr=refmagerr)
    ngal = n_elements(cat)

    filt = megaspitzer_filterlist()
    isf160 = (where(strmatch(filt,'*f160w*',/fold)))[0]
    
; with irac    
    outprefix = 'killerplot'
    junk = isedfit_restore(paramfile,irac,iopath=isedpath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,/nomodels,outprefix=outprefix)
;   mstar = isedfit_reconstruct_posterior(paramfile,post=pirac,$
;     isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath,$
;     age=age,sfr0=sfr0,sfrage=sfrage,outprefix=outprefix)

; without irac
    outprefix = 'killerplot_noirac'
    junk = isedfit_restore(paramfile,noirac,iopath=isedpath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,/nomodels,$
      outprefix=outprefix)
;   nomstar = isedfit_reconstruct_posterior(paramfile,post=pnoirac,$
;     isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath,$
;     age=noage,sfr0=nosfr0,sfrage=nosfrage,outprefix=outprefix)
    
; make the fractional uncertainty plot a 4-panel plot based on dropout
; class - just show the 1D histogram of stellar mass
    psfile = path+'irac_is_good_histo.ps'
    im_plotconfig, 5, pos, psfile=psfile, xmargin=[1.2,0.3], $
      width=[3.5,3.5], height=[2.8,2.8], xspace=0.1, yspace=0.1, $
      charsize=1.8, charthick=5

    xrange = [0.15,1.7]
    yrange = [0,45] ; [0.0,1.2]
    csize1 = 1.7
    showmed = 0
    line1 = 0
    line2 = 1

    col1 = 'black'
    col2 = 'navy'
    fcol1 = 'grey40'
    fcol2 = 'sky blue'
    syms = 0.6
    bin = 0.04
    
; B Dropouts    
    match = where(strmatch(cat.galaxy,'*B-*'),nmatch)
    djs_plot, [0], [0], /nodata, xsty=5, ysty=5, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xtickname=replicate(' ',10), xtickinterval=0.5
    im_legend, 'z \sim 4', /right, /top, box=0, charsize=csize1, margin=-0.1

    plothist, irac[match].mass_err*alog(10), color=im_color(col1), line=line1, bin=bin, /over, peak=0, $
      /fill, fcolor=im_color(fcol1), thick=10
    plothist, noirac[match].mass_err*alog(10), color=im_color(col2), line=line2, bin=bin, /over, peak=0, $
      /fill, fcolor=im_color(fcol2), thick=10

    plothist, irac[match].mass_err*alog(10), color=im_color(col1), line=line1, bin=bin, /over, peak=0, thick=6
    plothist, noirac[match].mass_err*alog(10), color=im_color(col2), line=line2, bin=bin, /over, peak=0, thick=4

    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, /noerase, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xtickname=replicate(' ',10), xtickinterval=0.5

; V Dropouts
    match = where(strmatch(cat.galaxy,'*V-*'),nmatch)
    djs_plot, [0], [0], /nodata, xsty=5, ysty=5, /noerase, position=pos[*,1], $
      xrange=xrange, yrange=yrange, xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    im_legend, 'z \sim 5', /right, /top, box=0, charsize=csize1, margin=-0.1

    im_legend, ['Existing Shallow IRAC','Proposed Deep IRAC'], line=[line2,line1], $
      /left, /bottom, box=1, thick=8, color=[col2,col1], $
      charsize=1.3, symsize=1.8, position=[0.23,27.0], /data, pspacing=1.8

    plothist, irac[match].mass_err*alog(10), color=im_color(col1), line=line1, bin=bin, /over, peak=0, $
      /fill, fcolor=im_color(fcol1), thick=10
    plothist, noirac[match].mass_err*alog(10), color=im_color(col2), line=line2, bin=bin, /over, peak=0, $
      /fill, fcolor=im_color(fcol2), thick=10

    plothist, irac[match].mass_err*alog(10), color=im_color(col1), line=line1, bin=bin, /over, peak=0, thick=6
    plothist, noirac[match].mass_err*alog(10), color=im_color(col2), line=line2, bin=bin, /over, peak=0, thick=4

    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, /noerase, position=pos[*,1], $
      xrange=xrange, yrange=yrange, xtickname=replicate(' ',10), ytickname=replicate(' ',10)

; i Dropouts
    match = where(strmatch(cat.galaxy,'*I-*'),nmatch)
;   match = where(strmatch(cat.galaxy,'*I-*') or strmatch(cat.galaxy,'*Z-*'),nmatch)
    djs_plot, [0], [0], /nodata, xsty=5, ysty=5, /noerase, position=pos[*,2], $
      xrange=xrange, yrange=yrange, xtickinterval=0.5
    im_legend, 'z \sim 6', /right, /top, box=0, charsize=csize1, margin=-0.1

    plothist, irac[match].mass_err*alog(10), color=im_color(col1), line=line1, bin=bin, /over, peak=0, $
      /fill, fcolor=im_color(fcol1), thick=10
    plothist, noirac[match].mass_err*alog(10), color=im_color(col2), line=line2, bin=bin, /over, peak=0, $
      /fill, fcolor=im_color(fcol2), thick=10

    plothist, irac[match].mass_err*alog(10), color=im_color(col1), line=line1, bin=bin, /over, peak=0, thick=6
    plothist, noirac[match].mass_err*alog(10), color=im_color(col2), line=line2, bin=bin, /over, peak=0, thick=4

    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, /noerase, position=pos[*,2], $
      xrange=xrange, yrange=yrange, xtickinterval=0.5    
    
; zY Dropouts
    match = where(strmatch(cat.galaxy,'*Z-*') or strmatch(cat.galaxy,'*Y-*'),nmatch)

    djs_plot, [0], [0], /nodata, xsty=5, ysty=5, /noerase, position=pos[*,3], $
      xrange=xrange, yrange=yrange, ytickname=replicate(' ',10), xtickinterval=0.5
    im_legend, 'z \sim 7-10', /right, /top, box=0, charsize=csize1, margin=-0.1

    plothist, irac[match].mass_err*alog(10), color=im_color(col1), line=line1, bin=bin, /over, peak=0, $
      /fill, fcolor=im_color(fcol1), thick=10
    plothist, noirac[match].mass_err*alog(10), color=im_color(col2), line=line2, bin=bin, /over, peak=0, $
      /fill, fcolor=im_color(fcol2), thick=10

    plothist, irac[match].mass_err*alog(10), color=im_color(col1), line=line1, bin=bin, /over, peak=0, thick=6
    plothist, noirac[match].mass_err*alog(10), color=im_color(col2), line=line2, bin=bin, /over, peak=0, thick=4

    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, /noerase, position=pos[*,3], $
      xrange=xrange, yrange=yrange, ytickname=replicate(' ',10), xtickinterval=0.5

; add WeiOne
    ww = where(cat.z gt 9.0)
    yy = 18
    arrow, noirac[ww].mass_err*alog(10), yy, irac[ww].mass_err*alog(10)+0.035, yy, $
      /data, thick=8.0, hsize=-0.2, hthick=8.0, color=im_color('dark red')
    plots, noirac[ww].mass_err*alog(10), yy, psym=symcat(6,thick=8), symsize=1.5, color=im_color('dark red')
    plots, irac[ww].mass_err*alog(10), yy, psym=symcat(15), symsize=1.5, color=im_color('dark red')
    xyouts, noirac[ww].mass_err*alog(10)-0.1, yy+5, textoidl('MACS1149-JD (z\sim9.6)'), $
      align=0.5, charsize=1.3, color=im_color('dark red')
        
    xyouts, pos[0,0]-0.09, djs_mean([pos[1,0],pos[3,2]]), 'Number of Galaxies', $
      align=0.5, orientation=90, /normal
    xyouts, djs_mean([pos[2,2],pos[0,3]]), pos[1,3]-0.1, 'Fractional Uncertainty in Stellar Mass', $
      align=0.5, /normal
    
    im_plotconfig, psfile=psfile, /psclose, /pdf
    
stop    

; make the fractional uncertainty plot a 4-panel plot based on dropout
; class 
    psfile = path+'irac_is_good.ps'
    im_plotconfig, 5, pos, psfile=psfile, xmargin=[1.2,0.3], $
      width=[3.5,3.5], height=[3.0,3.0], xspace=0.1, yspace=0.1, $
      charsize=1.8, charthick=5

    xrange = [0.0,0.9]
    yrange = [0,1.7]
    csize1 = 1.7
    showmed = 0

;   col1 = 'grey40'
    col1 = 'black'
    col2 = 'dodger blue'
    syms = 0.6
    
; B Dropouts    
    match = where(strmatch(cat.galaxy,'*B-*'),nmatch)
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xtickname=replicate(' ',10), ytickinterval=0.5
    im_legend, '<z> = 3.8', /left, /top, box=0, charsize=csize1, margin=-0.1
;   im_legend, 'B Dropouts (N='+strtrim(nmatch,2)+')', /right, /bottom, box=0, charsize=csize1, margin=-0.1
;   im_legend, 'B Dropouts Dropouts (N='+strtrim(nmatch,2)+')', /right, /bottom, box=0, charsize=csize1, margin=0

    djs_oplot, irac[match].sfrage_err/irac[match].sfrage_50, irac[match].mass_err*alog(10), $
      psym=symcat(15), symsize=syms, color=im_color(col1)
    djs_oplot, noirac[match].sfrage_err/noirac[match].sfrage_50, noirac[match].mass_err*alog(10), $
      psym=symcat(16), symsize=syms, color=im_color(col2)

    if showmed then begin
       oploterror, djs_median(irac[match].sfrage_err/irac[match].sfrage_50), $
         djs_median(irac[match].mass_err*alog(10)), 0*djsig(irac[match].sfrage_err/irac[match].sfrage_50), $
         0*djsig(irac[match].mass_err*alog(10)), psym=symcat(16), symsize=3.5, $
         color='black', errcolor='black', errthick=10, line=0
       oploterror, djs_median(noirac[match].sfrage_err/noirac[match].sfrage_50), $
         djs_median(noirac[match].mass_err*alog(10)), 0*djsig(noirac[match].sfrage_err/noirac[match].sfrage_50), $
         0*djsig(noirac[match].mass_err*alog(10)), psym=symcat(15), symsize=3.5, $
         color=im_color(col2), errcolor=im_color(col2), errthick=10, line=5
    endif

; V Dropouts
    match = where(strmatch(cat.galaxy,'*V-*'),nmatch)
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, /noerase, position=pos[*,1], $
      xrange=xrange, yrange=yrange, xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    im_legend, '<z> = 4.9', /left, /top, box=0, charsize=csize1, margin=-0.1
;   im_legend, 'V Dropouts (N='+strtrim(nmatch,2)+')', /right, /bottom, box=0, charsize=csize1, margin=-0.1
    djs_oplot, irac[match].sfrage_err/irac[match].sfrage_50, irac[match].mass_err*alog(10), $
      psym=symcat(15), symsize=syms, color=im_color(col1)
    djs_oplot, noirac[match].sfrage_err/noirac[match].sfrage_50, noirac[match].mass_err*alog(10), $
      psym=symcat(16), symsize=syms, color=im_color(col2)

    if showmed then begin
       oploterror, djs_median(irac[match].sfrage_err/irac[match].sfrage_50), $
         djs_median(irac[match].mass_err*alog(10)), 0*djsig(irac[match].sfrage_err/irac[match].sfrage_50), $
         0*djsig(irac[match].mass_err*alog(10)), psym=symcat(16), symsize=3.5, $
         color='black', errcolor='black', errthick=10, line=0
       oploterror, djs_median(noirac[match].sfrage_err/noirac[match].sfrage_50), $
         djs_median(noirac[match].mass_err*alog(10)), 0*djsig(noirac[match].sfrage_err/noirac[match].sfrage_50), $
         0*djsig(noirac[match].mass_err*alog(10)), psym=symcat(15), symsize=3.5, $
         color=im_color(col2), errcolor=im_color(col2), errthick=10, line=5
    endif

; i Dropouts
    match = where(strmatch(cat.galaxy,'*I-*'),nmatch)
;   match = where(strmatch(cat.galaxy,'*I-*') or strmatch(cat.galaxy,'*Z-*'),nmatch)
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, /noerase, position=pos[*,2], $
      xrange=xrange, yrange=yrange, ytickinterval=0.5, xtickinterval=0.2
    im_legend, ['Shallow IRAC','Deep IRAC (This Proposal)'], psym=[16,15], $
      /left, /bottom, box=0, line=[5,0], thick=8, color=[col2,col1], $
      charsize=1.4, symsize=1.8, margin=-0.1

    im_legend, '<z> = 5.9', /left, /top, box=0, charsize=csize1, margin=-0.1
;   im_legend, 'i Dropouts (N='+strtrim(nmatch,2)+')', /right, /bottom, box=0, charsize=csize1, margin=-0.1
    djs_oplot, irac[match].sfrage_err/irac[match].sfrage_50, irac[match].mass_err*alog(10), $
      psym=symcat(15), symsize=syms, color=im_color(col1)
    djs_oplot, noirac[match].sfrage_err/noirac[match].sfrage_50, noirac[match].mass_err*alog(10), $
      psym=symcat(16), symsize=syms, color=im_color(col2)

    if showmed then begin
       oploterror, djs_median(irac[match].sfrage_err/irac[match].sfrage_50), $
         djs_median(irac[match].mass_err*alog(10)), 0*djsig(irac[match].sfrage_err/irac[match].sfrage_50), $
         0*djsig(irac[match].mass_err*alog(10)), psym=symcat(16), symsize=3.5, $
         color='black', errcolor='black', errthick=10, line=0
       oploterror, djs_median(noirac[match].sfrage_err/noirac[match].sfrage_50), $
         djs_median(noirac[match].mass_err*alog(10)), 0*djsig(noirac[match].sfrage_err/noirac[match].sfrage_50), $
         0*djsig(noirac[match].mass_err*alog(10)), psym=symcat(15), symsize=3.5, $
         color=im_color(col2), errcolor=im_color(col2), errthick=10, line=5
    endif

; zY Dropouts
    match = where(strmatch(cat.galaxy,'*Z-*') or strmatch(cat.galaxy,'*Y-*'),nmatch)
;   match = where(strmatch(cat.galaxy,'*Y-*'),nmatch)
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, /noerase, position=pos[*,3], $
      xrange=xrange, yrange=yrange, ytickname=replicate(' ',10), xtickinterval=0.2

    im_legend, '<z> = 6.8-7.8', /left, /top, box=0, charsize=csize1, margin=-0.1
;   im_legend, 'z,Y Dropouts (N='+strtrim(nmatch,2)+')', /right, /bottom, box=0, charsize=csize1, margin=-0.1
    djs_oplot, irac[match].sfrage_err/irac[match].sfrage_50, irac[match].mass_err*alog(10), $
      psym=symcat(15), symsize=syms, color=im_color(col1)
    djs_oplot, noirac[match].sfrage_err/noirac[match].sfrage_50, noirac[match].mass_err*alog(10), $
      psym=symcat(16), symsize=syms, color=im_color(col2)

; add WeiOne
    ww = where(cat.z gt 9.0)
;   djs_oplot, [irac[ww].sfrage_err/irac[ww].sfrage_50,noirac[ww].sfrage_err/noirac[ww].sfrage_50], $
;     [irac[ww].mass_err*alog(10),noirac[ww].mass_err*alog(10)], line=0, thick=5
;   im_symbols, 115, thick=8, color=im_color(col2)
    plots, noirac[ww].sfrage_err/noirac[ww].sfrage_50, noirac[ww].mass_err*alog(10), $
;     psym=8, symsize=4.0, color=im_color(col2)
      psym=symcat(45,thick=10), symsize=3.5, color=im_color('dark red')

;   im_symbols, 115, /fill
    plots, irac[ww].sfrage_err/irac[ww].sfrage_50, irac[ww].mass_err*alog(10), $
;     psym=8, symsize=4.0, color=im_color(col1)
      psym=symcat(46), symsize=4.0, color=im_color('dark red')
;   arrow, 0.4, 0.3, 0.45, 0.4, /data, thick=4, hthick=10
    xyouts, 0.45, 0.25, 'MACS1149-JD (z=9.6)', align=0.5, charsize=1.3, color=im_color('dark red')

    if showmed then begin
       oploterror, djs_median(irac[match].sfrage_err/irac[match].sfrage_50), $
         djs_median(irac[match].mass_err*alog(10)), 0*djsig(irac[match].sfrage_err/irac[match].sfrage_50), $
         0*djsig(irac[match].mass_err*alog(10)), psym=symcat(16), symsize=3.5, $
         color='black', errcolor='black', errthick=10, line=0
       oploterror, djs_median(noirac[match].sfrage_err/noirac[match].sfrage_50), $
         djs_median(noirac[match].mass_err*alog(10)), 0*djsig(noirac[match].sfrage_err/noirac[match].sfrage_50), $
         0*djsig(noirac[match].mass_err*alog(10)), psym=symcat(15), symsize=3.5, $
         color=im_color(col2), errcolor=im_color(col2), errthick=10, line=5
    endif

    xyouts, pos[0,0]-0.09, djs_mean([pos[1,0],pos[3,2]]), 'Fractional Uncertainty in Stellar Mass', $
      align=0.5, orientation=90, /normal
    xyouts, djs_mean([pos[2,2],pos[0,3]]), pos[1,3]-0.1, 'Fractional Uncertainty in SFR-Weighted Age', $
      align=0.5, /normal
    
    im_plotconfig, psfile=psfile, /psclose, /pdf
    
stop    

; now some contour plots    
    psfile = path+'killerplot_v2.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.4,0.4], $
      width=6.7, height=5.0, charsize=1.7

    for ii = 0, ngal-1 do begin
       title = cat[ii].galaxy+', F160W = '+string(refmag[isf160,ii],format='(F4.1)')
;      title = cat[ii].galaxy+'(z='+string(irac[ii].zobj,format='(F3.1)')+$
;        '), F160W = '+string(refmag[isf160,ii],format='(F4.1)')
;      title = 'Obj #'+strtrim(ii+1,2)+' (z='+string(irac[ii].zobj,format='(F3.1)')+$
;        '), F160W='+string(refmag[isf160,ii],format='(F4.1)')

       hogg_scatterplot, noage[*,ii]*1E3, nomstar[*,ii], position=pos, $
         xsty=1, ysty=1, /nogrey, ccolor=im_color('dodger blue'), cthick=4.0, $
         xrange=[0,1500], yrange=[8.0,11.5], /internal, cline=5, $
         xtitle=textoidl('<t>_{SFR} (Myr)'), ytitle=textoidl('log (M_{*}/M'+sunsymbol()+')'), $
         levels=errorf([1.0,2.0]/sqrt(2)), title=title

       im_legend, ['Shallow IRAC','Deep IRAC'], $
         /right, /bottom, box=0, line=[5,0], thick=8, color=['dodger blue','black']

       hogg_scatterplot, age[*,ii]*1E3, mstar[*,ii], position=pos, $
         xsty=5, ysty=5, /noerase, /nogrey, ccolor=im_color('black'), cthick=8, $
         xrange=[50,1500], yrange=[8.0,11.5], /internal, $
         levels=errorf([1.0,2.0]/sqrt(2))
    endfor
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
return
end


;;; make the fractional uncertainty plot a 4-panel plot based on dropout
;;; class 
;;    psfile = path+'killerplot_uncert.ps'
;;    im_plotconfig, 5, pos, psfile=psfile, xmargin=[1.1,0.4], $
;;      width=[3.5,3.5], height=[3.0,3.0], xspace=0.1, yspace=0.1, $
;;      charsize=1.8, charthick=5
;;
;;    xrange=[0,0.9]
;;    yrange=[0,1.6]
;;    csize1 = 1.7
;;    
;;; B-band    
;;    match = where(strmatch(cat.galaxy,'*B-*'),nmatch)
;;    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, position=pos[*,0], $
;;      xrange=xrange, yrange=yrange, xtickname=replicate(' ',10), ytickinterval=0.5
;;    im_legend, 'B-band (N='+strtrim(nmatch,2)+')', /right, /bottom, box=0, charsize=csize1, margin=0
;;;   im_legend, 'B-band Dropouts (N='+strtrim(nmatch,2)+')', /right, /bottom, box=0, charsize=csize1, margin=0
;;    im_legend, ['Deep IRAC','Shallow IRAC'], psym=[16,15], $
;;      /left, /top, box=0, line=[0,5], thick=8, color=['black','dodger blue'], $
;;      charsize=1.5, symsize=1.5, margin=0
;;
;;    djs_oplot, irac[match].sfrage_err/irac[match].sfrage_50, irac[match].mass_err*alog(10), $
;;      psym=symcat(16), symsize=0.5, color='black'
;;    djs_oplot, noirac[match].sfrage_err/noirac[match].sfrage_50, noirac[match].mass_err*alog(10), $
;;      psym=symcat(16), symsize=0.5, color=im_color('dodger blue')
;;    
;;    oploterror, djs_median(irac[match].sfrage_err/irac[match].sfrage_50), $
;;      djs_median(irac[match].mass_err*alog(10)), djsig(irac[match].sfrage_err/irac[match].sfrage_50), $
;;      djsig(irac[match].mass_err*alog(10)), psym=symcat(16), symsize=3.5, $
;;      color='black', errcolor='black', errthick=10, line=0
;;    oploterror, djs_median(noirac[match].sfrage_err/noirac[match].sfrage_50), $
;;      djs_median(noirac[match].mass_err*alog(10)), djsig(noirac[match].sfrage_err/noirac[match].sfrage_50), $
;;      djsig(noirac[match].mass_err*alog(10)), psym=symcat(15), symsize=3.5, $
;;      color=im_color('dodger blue'), errcolor=im_color('dodger blue'), errthick=10, line=5
;;
;;; V-band
;;    match = where(strmatch(cat.galaxy,'*V-*'),nmatch)
;;    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, /noerase, position=pos[*,1], $
;;      xrange=xrange, yrange=yrange, xtickname=replicate(' ',10), ytickname=replicate(' ',10)
;;    im_legend, 'V-band (N='+strtrim(nmatch,2)+')', /right, /bottom, box=0, charsize=csize1, margin=0
;;    djs_oplot, irac[match].sfrage_err/irac[match].sfrage_50, irac[match].mass_err*alog(10), $
;;      psym=symcat(16), symsize=0.5, color='black'
;;    djs_oplot, noirac[match].sfrage_err/noirac[match].sfrage_50, noirac[match].mass_err*alog(10), $
;;      psym=symcat(16), symsize=0.5, color=im_color('dodger blue')
;;
;;    oploterror, djs_median(irac[match].sfrage_err/irac[match].sfrage_50), $
;;      djs_median(irac[match].mass_err*alog(10)), djsig(irac[match].sfrage_err/irac[match].sfrage_50), $
;;      djsig(irac[match].mass_err*alog(10)), psym=symcat(16), symsize=3.5, $
;;      color='black', errcolor='black', errthick=10, line=0
;;    oploterror, djs_median(noirac[match].sfrage_err/noirac[match].sfrage_50), $
;;      djs_median(noirac[match].mass_err*alog(10)), djsig(noirac[match].sfrage_err/noirac[match].sfrage_50), $
;;      djsig(noirac[match].mass_err*alog(10)), psym=symcat(15), symsize=3.5, $
;;      color=im_color('dodger blue'), errcolor=im_color('dodger blue'), errthick=10, line=5
;;
;;; Iz-band
;;    match = where(strmatch(cat.galaxy,'*I-*') or strmatch(cat.galaxy,'*Z-*'),nmatch)
;;    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, /noerase, position=pos[*,2], $
;;      xrange=xrange, yrange=yrange, ytickinterval=0.5
;;    im_legend, 'Iz-band (N='+strtrim(nmatch,2)+')', /right, /bottom, box=0, charsize=csize1, margin=0
;;    djs_oplot, irac[match].sfrage_err/irac[match].sfrage_50, irac[match].mass_err*alog(10), $
;;      psym=symcat(16), symsize=0.5, color='black'
;;    djs_oplot, noirac[match].sfrage_err/noirac[match].sfrage_50, noirac[match].mass_err*alog(10), $
;;      psym=symcat(16), symsize=0.5, color=im_color('dodger blue')
;;
;;    oploterror, djs_median(irac[match].sfrage_err/irac[match].sfrage_50), $
;;      djs_median(irac[match].mass_err*alog(10)), djsig(irac[match].sfrage_err/irac[match].sfrage_50), $
;;      djsig(irac[match].mass_err*alog(10)), psym=symcat(16), symsize=3.5, $
;;      color='black', errcolor='black', errthick=10, line=0
;;    oploterror, djs_median(noirac[match].sfrage_err/noirac[match].sfrage_50), $
;;      djs_median(noirac[match].mass_err*alog(10)), djsig(noirac[match].sfrage_err/noirac[match].sfrage_50), $
;;      djsig(noirac[match].mass_err*alog(10)), psym=symcat(15), symsize=3.5, $
;;      color=im_color('dodger blue'), errcolor=im_color('dodger blue'), errthick=10, line=5
;;
;;; Y-band
;;    match = where(strmatch(cat.galaxy,'*Y-*'),nmatch)
;;    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, /noerase, position=pos[*,3], $
;;      xrange=xrange, yrange=yrange, ytickname=replicate(' ',10)
;;    im_legend, 'Y-band (N='+strtrim(nmatch,2)+')', /right, /bottom, box=0, charsize=csize1, margin=0
;;    djs_oplot, irac[match].sfrage_err/irac[match].sfrage_50, irac[match].mass_err*alog(10), $
;;      psym=symcat(16), symsize=0.5, color='black'
;;    djs_oplot, noirac[match].sfrage_err/noirac[match].sfrage_50, noirac[match].mass_err*alog(10), $
;;      psym=symcat(16), symsize=0.5, color=im_color('dodger blue')
;;
;;    oploterror, djs_median(irac[match].sfrage_err/irac[match].sfrage_50), $
;;      djs_median(irac[match].mass_err*alog(10)), djsig(irac[match].sfrage_err/irac[match].sfrage_50), $
;;      djsig(irac[match].mass_err*alog(10)), psym=symcat(16), symsize=3.5, $
;;      color='black', errcolor='black', errthick=10, line=0
;;    oploterror, djs_median(noirac[match].sfrage_err/noirac[match].sfrage_50), $
;;      djs_median(noirac[match].mass_err*alog(10)), djsig(noirac[match].sfrage_err/noirac[match].sfrage_50), $
;;      djsig(noirac[match].mass_err*alog(10)), psym=symcat(15), symsize=3.5, $
;;      color=im_color('dodger blue'), errcolor=im_color('dodger blue'), errthick=10, line=5
;;
;;    xyouts, pos[0,0]-0.09, djs_mean([pos[1,0],pos[3,2]]), 'Fractional Uncertainty in Stellar Mass', $
;;      align=0.5, orientation=90, /normal
;;    xyouts, djs_mean([pos[2,2],pos[0,3]]), pos[1,3]-0.1, 'Fractional Uncertainty in SFR-Weighted Age', $
;;      align=0.5, /normal
;;    
;;    im_plotconfig, psfile=psfile, /psclose, /pdf
;;    
