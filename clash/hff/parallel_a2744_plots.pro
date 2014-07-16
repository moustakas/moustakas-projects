pro parallel_a2744_plots
; jm14jul16siena
    
    isedfit_dir = getenv('CLASH_PROJECTS')+'/hff/parallel_a2744/'
    montegrids_dir = isedfit_dir+'montegrids/'

    filt = hff_filterlist(pivotwave=weff,width=hwhm,/useirac)
    weff = weff/1D4
    hwhm = hwhm/1D4

; --------------------------------------------------
; QAplot: P(z)'s
    isedfit_paramfile = isedfit_dir+'parallel_a2744_photoz_paramfile.par'

    rr = read_isedfit(isedfit_paramfile,params=pp,isedfit_post=post)
    cat = read_parallel_a2744()
    ngal = n_elements(cat)

;   for ii = 0, ngal-1 do post[ii].pofz = alog(post[ii].pofz/total(post[ii].pofz))
    for ii = 0, ngal-1 do post[ii].pofz = post[ii].pofz/total(post[ii].pofz)
;   for ii = 0, ngal-1 do post[ii].pofz = post[ii].pofz/$
;     im_integral(pp.redshift,post[ii].pofz)

;   zmin = isedfit_find_zmin(pp.redshift,-2.0*alog(pofz_final>1D-20),nmin=nmin)
;   struct_print, zmin
    
    xrange = [0,12]
;   yrange = [-15,0]
    yrange = [0,max(post.pofz)*1.02]
    
    psfile = isedfit_dir+'parallel_a2744_pofz.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, height=4.5, $
      xmargin=[1.4,0.4], width=6.7
    for ii = 0, ngal-1 do begin
;      yrange = [0,max(post[ii].pofz)*1.05]
       djs_plot, [0], [0], /nodata, position=pos, xsty=5, ysty=5, $
         xrange=xrange, yrange=yrange ;, /ylog
       polyfill, [pp.redshift,reverse(pp.redshift)], $
;        [post[ii].pofz,post[ii].pofz*0+!y.crange[0]], $
         [post[ii].pofz,post[ii].pofz*0], $
         /fill, color=cgcolor('powder blue'), noclip=0
       djs_oplot, pp.redshift, post[ii].pofz, $
         thick=6, line=0, color=cgcolor('navy') ;, psym=10
       djs_plot, [0], [0], /nodata, position=pos, /noerase, $
         xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
         ytitle='Posterior Probability', $
;        ytitle='log (Likelihood)', $
         xtitle='Redshift'
       im_legend, cat[ii].galaxy, /left, /top, box=0, charsize=1.5, $
         margin=0, color=color, thick=8, pspacing=1.7
    endfor
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
return
end

