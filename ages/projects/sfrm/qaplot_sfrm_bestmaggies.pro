pro qaplot_sfrm_bestmaggies
; jm10feb10ucsd - compare the observed and best-fitting photometry

    sfrmpath = ages_path(/projects)+'sfrm/'
    paperpath = ages_path(/papers)+'sfrm/'

    parent = read_sfrm_sample()

    filt = ages_filterlist()
;   filt = strtrim(parent[0].filterlist,2)
    nfilt = n_elements(filt)
    band = strarr(nfilt)
    for ii = 0, nfilt-1 do band[ii] = repstr(strmid(filt[ii],$
      strpos(filt[ii],'_',/reverse_search)+1),'.par','')

    psfile = sfrmpath+'qaplots/sfrm_bestmaggies.ps'
    im_plotconfig, 6, pos, psfile=psfile, xmargin=[1.4,0.3], $
      height=[4.5,3.0], width=6.8

    mrange = [$
      [18.0,26.0],$ ; FUV
      [18.0,27.0],$ ; NUV
      [17.5,25.0],$ ; Bw
      [16.5,22.5],$ ; R
      [16.0,21.5],$ ; I
;     [15.0,20.5],$ ; z
      [15.5,21.0],$ ; J
      [15.0,21.0],$ ; H
      [14.5,21.0],$ ; Ks
      [15.5,22.0],$ ; ch1
      [15.5,22.0],$ ; ch2
      [16.0,22.5],$ ; ch3
      [15.0,22.0]]  ; ch4
    rrange1 = 3.9*[-1,1]
    rrange2 = 1.1*[-1,1]

    ngal = n_elements(parent)
    for ii = 0, nfilt-1 do begin
       good = where((parent.maggies[ii] gt 0),ngood)
       xx = -2.5*alog10(parent[good].maggies[ii])
       yy = -2.5*alog10(parent[good].bestmaggies[ii])

       hogg_scatterplot, xx, yy, position=pos[*,0], xsty=1, ysty=1, $
         xrange=mrange[*,ii], yrange=mrange[*,ii], $
         psym=symcat(6,thick=2), symsize=0.5, xtitle='', xtickname=replicate(' ',10), $
         ytitle=band[ii]+' (synthesized, AB mag)', /outlier, outcolor=djs_icolor('grey'), $
         levels=[0.5,0.75,0.9], /internal
;      djs_plot, xx, yy, position=pos[*,0], xsty=1, ysty=1, $
;        xrange=mrange[*,ii], yrange=mrange[*,ii], $
;        psym=symcat(6,thick=2), symsize=0.5, xtitle='', xtickname=replicate(' ',10), $
;        ytitle=band[ii]+' (synthesized, AB mag)'
       djs_oplot, !x.crange, !y.crange, line=0, thick=3, color='red'
       im_legend, 'N='+string(ngood,format='(I0)')+'/'+$
         string(ngal,format='(I0)'), /right, /bottom, $
         box=0, charsize=1.5, margin=0

       im_legend, '<\Delta'+band[ii]+'> = '+$
         im_string_stats(yy-xx,sigrej=3.0), /left, /top, $ ;/right, /bottom, $
         box=0, charsize=1.6, margin=0
       if (ii le 1) then rrange = rrange1 else rrange = rrange2
       hogg_scatterplot, xx, yy-xx, position=pos[*,1], /noerase, $
         xsty=1, ysty=1, xrange=mrange[*,ii], yrange=rrange, $
         psym=symcat(6,thick=2), symsize=0.5, xtitle=band[ii]+' (observed, AB mag)', $
         ytitle='Residuals (AB mag)', /outlier, outcolor=djs_icolor('grey'), $
         levels=[0.5,0.75,0.9], /internal
;      djs_plot, xx, yy-xx, position=pos[*,1], /noerase, $
;        xsty=1, ysty=1, xrange=mrange[*,ii], yrange=rrange, $
;        psym=symcat(6,thick=2), symsize=0.5, xtitle=band[ii]+' (observed, AB mag)', $
;        ytitle='Residuals (AB mag)'
       djs_oplot, !x.crange, [0,0], line=0, thick=3, color='red'
       !p.multi = 0
    endfor
    im_plotconfig, /psclose, psfile=psfile, /gzip
;   spawn, 'rsync -auv '+psfile+'.gz ~/', /sh

stop    
return
end
    
