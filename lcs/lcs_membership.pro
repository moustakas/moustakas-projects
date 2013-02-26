pro lcs_membership, clobber=clobber
; jm13feb11siena - build a sample of galaxies from the NSA belonging
; to the LCS cluster sample

    common lcs_nsa, nsa
    if n_elements(nsa) eq 0L then nsa = read_nsa()

    lcspath = getenv('LCS_DATA')

    light = im_light(/km)
    H0 = 70.0
    scale = 1D3
    radius_deg = 3.0

    mproj_mpc = 3.0 ; perpendicular to LoS [Mpc]
    mz_mpc = 1.5 ; along LoS [Mpc]
    mz_kms = mz_mpc*H0
    
    mproj = H0*mproj_mpc/light
    mz = H0*mz_mpc/light
    splog, mproj_mpc, mz_mpc, mz_kms, mproj, mz

    cl = rsex(getenv('LCS_DIR')+'/lcs_sample.cat')
    ncl = n_elements(cl)
    struct_print, cl

    plot, nsa.ra, nsa.dec, ps=3, xr=[150,270.0], yr=[0,50]
    djs_oplot, cl.ra, cl.dec, psym=6, thick=8, color='orange'

;   im_window, 0, yr=0.8
    psfile = 'qa_lcs_membership.ps'
    im_plotconfig, 6, pos, yspace=0.9, psfile=psfile
    
;   for ii = 0, 1 do begin
    for ii = 0, ncl-1 do begin
       these = where(nsa.ra gt cl[ii].ra-radius_deg and nsa.ra lt cl[ii].ra+radius_deg and $
         nsa.dec gt cl[ii].dec-radius_deg and nsa.dec lt cl[ii].dec+radius_deg,ngal)
       ngal = n_elements(these)

       ing = zgroup(nsa[these].zdist,nsa[these].ra,nsa[these].dec,mproj,$
         mz,mult=mult,first=first,next=next)
       hh = histogram(ing,bin=1)
       peak = max(hh,indx)
       mem = where(ing eq indx,nmem,comp=nonmem)
       sig = djsig(nsa[these[mem]].zdist*light)
       czbin = sig*0.2
       splog, sig
       
       dist = djs_diff_angle(nsa[these].ra,nsa[these].dec,cl[ii].ra,cl[ii].dec)
       djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=3, ysty=3, $
         xtitle='Projected Clustercentric Distance (degrees)', $
         ytitle='Redshift (1000 km s^{-1})', yrange=minmax(nsa[these].zdist)*light/scale, $
         xrange=minmax(dist), title=strtrim(cl[ii].cluster,2)+': '+strtrim(nmem,2)+' members'
       djs_oplot, dist, nsa[these].zdist*light/scale, psym=8
       djs_oplot, dist[mem], nsa[these[mem]].zdist*light/scale, psym=8, $
         color=im_color('tan')
       djs_oplot, dist[nonmem], nsa[these[nonmem]].zdist*light/scale, psym=8
       djs_oplot, !x.crange, cl[ii].vel*[1,1]/scale, line=0, color=im_color('firebrick')

       im_plothist, nsa[these].zdist*light/scale, bin=czbin/scale, /noplot, xb, yb, $
         xsty=1, ysty=1
       im_plothist, nsa[these].zdist*light/scale, bin=czbin/scale, position=pos[*,1], $
         /noerase, xtitle='Redshift (1000 km s^{-1})', ytitle='Number of Galaxies', $
         xrange=(cl[ii].vel+10*sig*[-1,1])/scale, yrange=[0,max(yb)*1.05], xsty=1, ysty=1
       im_plothist, nsa[these[mem]].zdist*light/scale, bin=czbin/scale, /overplot, /fill, $
         fcolor=im_color('tan'), $ ; color=im_color('firebrick'), $
         position=pos[*,1], xsty=1, ysty=1, $
         xrange=(cl[ii].vel+10*sig*[-1,1])/scale, yrange=[0,max(yb)*1.05]
;      im_legend, strtrim(cl[ii].cluster,2)+': '+strtrim(nmem,2)+' members', $
;        /left, /top, box=0
       djs_oplot, cl[ii].vel*[1,1]/scale, !y.crange, color=im_color('firebrick')
;      cc = get_kbrd(1)
       
; write out
       nsafile = lcspath+strlowcase(cl[ii].cluster)+'_nsa.fits'
       im_mwrfits, nsa[these], nsafile, clobber=clobber
       
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf

return
end
    

    
