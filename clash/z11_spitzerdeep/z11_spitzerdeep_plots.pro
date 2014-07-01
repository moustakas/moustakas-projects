pro z11_spitzerdeep_plots

    isedfit_dir = getenv('CLASH_PROJECTS')+'/z11_spitzerdeep/'
    montegrids_dir = isedfit_dir+'montegrids/'

    filt = hff_filterlist(pivotwave=weff,width=hwhm,/useirac)
    weff = weff/1D4
    hwhm = hwhm/1D4

; --------------------------------------------------
; paper plot: P(z)'s
    isedfit_paramfile = isedfit_dir+'z11_spitzerdeep_photoz_paramfile.par'

    rr = read_isedfit(isedfit_paramfile,params=pp,isedfit_post=post,index=index)
    cat = read_z11_spitzerdeep()
    ngal = n_elements(cat)
    color = ['dodger blue','tan','orange']
    lcolor = ['navy','brown','red']
    
    for ii = 0, ngal-1 do post[ii].pofz = post[ii].pofz/total(post[ii].pofz)
;   for ii = 0, ngal-1 do post[ii].pofz = post[ii].pofz/$
;     im_integral(pp.redshift,post[ii].pofz)
    
    xrange = [0,12]
    yrange = [0,max(post.pofz)*1.05]
    
    psfile = isedfit_dir+'z11_spitzerdeep_pofz.eps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, height=4.5, $
      xmargin=[1.1,0.4], width=7.0
    djs_plot, [0], [0], /nodata, position=pos, xsty=5, ysty=5, $
      xrange=xrange, yrange=yrange
    
    for ii = 0, ngal-1 do begin
       polyfill, [pp.redshift,reverse(pp.redshift)], $
         [post[ii].pofz,post[ii].pofz*0], $
         /fill, color=cgcolor(color[ii]), noclip=0
    endfor

    for ii = 0, ngal-1 do djs_oplot, pp.redshift, post[ii].pofz, $
      thick=6, line=0, color=cgcolor(lcolor[ii]);, psym=10
;   djs_oplot, pp.redshift, pofz_final, line=0, thick=8, $
;     color=cgcolor('black')
    
    djs_plot, [0], [0], /nodata, position=pos, /noerase, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      ytitle='Redshift Probability', xtitle='Redshift', $
      ytickname=replicate(' ',10)

    im_legend, cat.galaxy, /left, /top, box=0, charsize=1.5, $
      margin=0, color=color, line=0, thick=8, pspacing=1.7
;   im_legend, ['JD1A','JD1B','JD1A x JD1B'], /left, /top, box=0, $
;     margin=0, color=[color,'black'], line=0, thick=8, pspacing=1.7
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

return
end

