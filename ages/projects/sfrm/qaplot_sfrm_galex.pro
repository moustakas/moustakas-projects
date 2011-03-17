pro qaplot_sfrm_galex
; jm10feb10ucsd - QAplots showing the quality of our GALEX
; K-corrections as a function of various properties

    sfrmpath = ages_path(/projects)+'sfrm/'
    paperpath = ages_path(/papers)+'sfrm/'

; require good photometry    
    parent = read_sfrm_sample()
    isolated = where(parent.i_segflags_aper_10 eq 0,comp=crowded)
;   isolated = where(parent.nfriends-1 le 0,comp=crowded)

    psfile = sfrmpath+'qaplots/sfrm_galex.ps'
    im_plotconfig, 5, pos, psfile=psfile, xspace=0.0, yspace=1.2

    band = ['FUV','NUV']
    for ii = 0, 1 do begin ; fuv, nuv
       good = where((parent.k_maggies[ii] gt 0),ngood)
       obsmag = -2.5*alog10(parent[good].k_maggies[ii])
       modelmag = -2.5*alog10(parent[good].k_bestmaggies[ii])
       resid = modelmag-obsmag
       
       imag = -2.5*alog10(parent[good].k_maggies[4])
       redshift = parent[good].z
       gr = parent[good].ugriz_absmag[1]-parent[good].ugriz_absmag[2]

       ytitle = band[ii]+' Residuals (Model minus Observed, AB mag)'

; UV magnitude       
       hogg_scatterplot, obsmag, resid, position=pos[*,0], xsty=1, ysty=1, $
         xrange=[18,25.7], yrange=[-4,4], xtitle=band[ii]+' (AB mag)', $
         ytitle='', /outlier, outcolor=djs_icolor('grey'), $
         levels=[0.5,0.75,0.9], /internal
       djs_oplot, obsmag[isolated], resid[isolated], psym=symcat(6,thick=3), color='red', symsize=0.12
       
; optical magnitude       
       hogg_scatterplot, imag, resid, position=pos[*,1], /noerase, xsty=1, ysty=1, $
         xrange=[15.5,21], yrange=[-4,4], xtitle='I (AB mag)', $
         ytitle='', /outlier, outcolor=djs_icolor('grey'), $
         levels=[0.5,0.75,0.9], /internal, ytickname=replicate(' ',10)
       djs_oplot, imag[isolated], resid[isolated], psym=symcat(6,thick=3), color='red', symsize=0.12
       
; redshift
       hogg_scatterplot, redshift, resid, position=pos[*,2], /noerase, xsty=1, ysty=1, $
         xrange=[0.0,0.59], yrange=[-4,4], xtitle='Redshift', $
         ytitle='', /outlier, outcolor=djs_icolor('grey'), $
         levels=[0.5,0.75,0.9], /internal
       djs_oplot, redshift[isolated], resid[isolated], psym=symcat(6,thick=3), color='red', symsize=0.12

; optical color
       hogg_scatterplot, gr, resid, position=pos[*,3], /noerase, xsty=1, ysty=1, $
         xrange=[-0.0,1.1], yrange=[-4,4], xtitle=textoidl('^{0.1}(g-r)'), $
         ytitle='', /outlier, outcolor=djs_icolor('grey'), $
         levels=[0.5,0.75,0.9], /internal, ytickname=replicate(' ',10)
       djs_oplot, gr[isolated], resid[isolated], psym=symcat(6,thick=3), color='red', symsize=0.12
       !p.multi = 0

       xyouts, 0.06, 0.54, ytitle, align=0.5, orientation=90, /norm
    endfor    
    im_plotconfig, /psclose, psfile=psfile, /gzip
    

return
end
    
