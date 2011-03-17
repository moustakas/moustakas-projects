pro visualize_ages_archetypes, binary=binary
; jm08sep18nyu - visualize the distribution of Hogg's AGES
;   archetypes in physical parameter space
    
    aa = read_ages(/ispec)

    chistr = '1.0'
    if keyword_set(binary) then $
      filename = 'ages_binary_program.'+chistr+'.fits' else $
      filename = 'agesId.'+chistr+'.fits'
    bb = mrdfits(filename,1)

    ww = lonarr(n_elements(bb))
    for jj = 0, n_elements(bb)-1 do ww[jj] = where(bb[jj].ages_id eq aa.ages_id)
;   ages_display_spectrum, aa[ww], /postscript, /pdf, $
;     plottype=2, labeltype=3, specfit=ss

    symsize = 0.5*sqrt(bb.responsibility*1.0)

    im_plotfaves, /post
    dfpsplot, repstr(filename,'.fits','.ps'), /color, /square
    djs_plot, aa[ww].d4000_narrow_model[0], aa[ww].h_alpha_ew[0], $
      xtitle='D_{n}(4000)', ytitle='EW(H\alpha) (\AA)', /ylog, $
      ps=symcat(9,thick=!p.thick), symsize=symsize, xsty=3, ysty=3, $
      xrange=[1.0,2.30], yrange=[0.3,300], title=repstr(filename,'_','/')
    dfpsclose
    im_plotfaves

return
end
