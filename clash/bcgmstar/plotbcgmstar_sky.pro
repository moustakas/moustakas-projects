pro plotbcgmstar_sky
; jm14aug28siena - make plots of the sky-subtraction stuff

    sample = read_bcgmstar_sample()
    ncl = n_elements(sample)
    
    paperpath = bcgmstar_path(/paper)
    skyinfopath = bcgmstar_path()+'skyinfo/'

;   filt = bcgmstar_filterlist(weff=weff)

    psfile = paperpath+'clash_skymode.eps'
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.3
    pos = im_getposition(nx=5,ny=3,/landscape,yspace=0.1,xspace=0.1,$
      xmargin=[1.1,0.4])

    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       
       skyinfo = mrdfits(skyinfopath+'skyinfo-'+cluster+'.fits.gz',1,/silent)
       parallel = mrdfits(skyinfopath+'parallel-skyinfo-'+cluster+'.fits.gz',1,/silent)
       par1 = where(strmatch(parallel.file,'*par1*'))
       par2 = where(strmatch(parallel.file,'*par2*'))

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
         xsty=1, ysty=1, xrange=[0.3,1.8], yrange=[-0.02,0.01], $
         ytickname=ytickname, xtickname=xtickname, ytickinterval=0.01, $
         xtickinterval=0.4
       djs_oplot, !x.crange, [0,0], line=1
       djs_oplot, skyinfo.weff/1D4, skyinfo.mode/skyinfo.factor, $
         psym=symcat(16), symsize=1.0, color=cgcolor('grey')
       djs_oplot, skyinfo.weff/1D4, skyinfo.mode/skyinfo.factor, $
         psym=symcat(9), symsize=1.1, color=cgcolor('black')
       
       djs_oplot, 1.03*parallel[par1].weff/1D4, parallel[par1].mode/parallel[par1].factor, $
         psym=symcat(15), symsize=1.0, color=cgcolor('tomato')
       djs_oplot, 1.03*parallel[par1].weff/1D4, parallel[par1].mode/parallel[par1].factor, $
         psym=symcat(6), symsize=1.1, color=cgcolor('firebrick')
       
       djs_oplot, 0.97*parallel[par2].weff/1D4, parallel[par2].mode/parallel[par2].factor, $
         psym=symcat(14), symsize=1.3, color=cgcolor('powder blue')
       djs_oplot, 0.97*parallel[par2].weff/1D4, parallel[par2].mode/parallel[par2].factor, $
         psym=symcat(4), symsize=1.5, color=cgcolor('blue')

       im_legend, strupcase(cluster), /right, /top, box=0, margin=0, charsize=1.0

;       oploterror, skyinfo.weff/1D4, skyinfo.mode/skyinfo.factor, $
;         skyinfo.sigma/skyinfo.factor/sqrt(skyinfo.nsky), $
;         psym=symcat(16), symsize=3.0
;       oploterror, parallel[par1].weff/1D4, parallel[par1].mode/parallel[par1].factor, $
;         parallel[par1].sigma/parallel[par1].factor/sqrt(parallel[par1].nsky), psym=symcat(15), $
;         symsize=2.0, color=cgcolor('orange')
;       oploterror, parallel[par2].weff/1D4, parallel[par2].mode/parallel[par2].factor, $
;         parallel[par2].sigma/parallel[par2].factor/sqrt(parallel[par2].nsky), psym=symcat(15), $
;         symsize=2.0, color=cgcolor('forest green')

    endfor
    
    xyouts, min(pos[0,*])-0.07, (max(pos[3,*])-min(pos[1,*]))/2.0+min(pos[1,*]), $
      textoidl('Sky Mode (electrons sec^{-1})'), orientation=90, align=0.5, charsize=1.4, /norm
    xyouts, (max(pos[2,*])-min(pos[0,*]))/2.0+min(pos[0,*]), min(pos[1,*])-0.09, $
      textoidl('Bandpass Wavelength (\mu'+'m)'), $
;     textoidl('Equivalent Radius r=a\sqrt{1-\epsilon} (kpc)'), $
      align=0.5, charsize=1.4, /norm
       
    im_plotconfig, psfile=psfile, /psclose, /pdf

return
end
    
    
