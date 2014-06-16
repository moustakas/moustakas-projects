pro lcs_membership, clobber=clobber
; jm13feb11siena - build a sample of galaxies from the NSA belonging
; to the LCS cluster sample

    common lcs_nsa, nsa
    if n_elements(nsa) eq 0L then nsa = read_nsa()

    lcspath = getenv('LCS_DATA')+'/membership/'

    light = im_light(/km)
    H0 = 70.0
    scale = 1D3
    radius_deg = 3.0

    mproj_mpc = 1.5 ; perpendicular to LoS [Mpc]
    mz_mpc = 3.0    ; along LoS [Mpc]
    mz_kms = mz_mpc*H0
    
    mproj = H0*mproj_mpc/light
    mz = H0*mz_mpc/light
    splog, mproj_mpc, mz_mpc, mz_kms, mproj, mz

    cl = rsex(getenv('LCS_DIR')+'/lcs_sample.cat')
    ncl = n_elements(cl)
    struct_print, cl

    niceprint, dangular(cl.vel/light,/kpc)/206265.0D, cl.cluster
    
;   plot, nsa.ra, nsa.dec, ps=3, xr=[150,270.0], yr=[0,50]
;   djs_oplot, cl.ra, cl.dec, psym=6, thick=8, color='orange'

;   im_window, 0, yr=0.8
    psfile = 'qa_lcs_membership.ps'
    im_plotconfig, 4, pos, yspace=[0.8,0.8], height=2.6*[1,1,1], $
      psfile=psfile, charsize=1.6, ymargin=[0.4,1.2]
    
;   for ii = 0, 1 do begin
    for ii = 0, ncl-1 do begin
       these = where(nsa.ra gt cl[ii].ra-radius_deg and nsa.ra lt cl[ii].ra+radius_deg and $
         nsa.dec gt cl[ii].dec-radius_deg and nsa.dec lt cl[ii].dec+radius_deg,ngal)
       ngal = n_elements(these)

       ing = zgroup(nsa[these].z,nsa[these].ra,nsa[these].dec,mproj,$
         mz,mult=mult,first=first,next=next)
       hh = histogram(ing,bin=1,locations=xx)
       srt = reverse(sort(hh))
       mem = where(ing eq xx[srt[0]],nmem)
       mem2 = where(ing eq xx[srt[1]],nmem2)
;      nonmem = where(ing ne xx[srt[0]])
       nonmem = where(ing ne xx[srt[0]] and ing ne xx[srt[1]])
       
       sig = djsig(nsa[these[mem]].z*light)
       czbin = sig*0.4
;      czbin = sig*0.2
       splog, sig, cl[ii].cluster
       
       dist = djs_diff_angle(nsa[these].ra,nsa[these].dec,cl[ii].ra,cl[ii].dec)
       djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=3, ysty=3, $
         xtitle='Projected Clustercentric Distance (degrees)', $
         ytitle='Redshift (1000 km s^{-1})', yrange=minmax(nsa[these].z)*light/scale, $
         xrange=minmax(dist);, title=strtrim(cl[ii].cluster,2)+': '+strtrim(nmem,2)+' members'
       djs_oplot, dist, nsa[these].z*light/scale, psym=8
       djs_oplot, dist[mem], nsa[these[mem]].z*light/scale, psym=8, $
         color=im_color('tan')
       djs_oplot, dist[mem2], nsa[these[mem2]].z*light/scale, psym=8, $
         color=im_color('dodger blue')
       djs_oplot, dist[nonmem], nsa[these[nonmem]].z*light/scale, psym=8
       djs_oplot, !x.crange, cl[ii].vel*[1,1]/scale, line=0, color=im_color('firebrick')

       im_plothist, nsa[these].z*light/scale, bin=czbin/scale, /noplot, xb, yb, $
         xsty=1, ysty=1
       yrange = [0,max(yb)*1.05]
;      xrange = (cl[ii].vel+10*sig*[-1,1])/scale
       xrange = minmax(nsa[these].z*light)/scale
       im_plothist, nsa[these].z*light/scale, bin=czbin/scale, position=pos[*,1], $
         /noerase, xtitle='Redshift (1000 km s^{-1})', ytitle='Number of Galaxies', $
         xrange=xrange, yrange=yrange, xsty=1, ysty=1
       im_plothist, nsa[these[mem]].z*light/scale, bin=czbin/scale, /overplot, /fill, $
         fcolor=im_color('tan'), $
         xrange=xrange, yrange=yrange, xsty=1, ysty=1
       im_plothist, nsa[these[mem2]].z*light/scale, bin=czbin/scale, /overplot, /fill, $
         fcolor=im_color('dodger blue'), $
         xrange=xrange, yrange=yrange, xsty=1, ysty=1
       im_legend, strtrim(cl[ii].cluster,2)+': '+strtrim(nmem,2)+' members', $
         /left, /top, box=0, margin=0, charsize=1.6
       djs_oplot, cl[ii].vel*[1,1]/scale, !y.crange, color=im_color('firebrick')
;      cc = get_kbrd(1)

; ra,dec plot       
       djs_plot, nsa[these].ra, nsa[these].dec, psym=8, xsty=3, ysty=3, $
         xtitle='\alpha_{J2000} (deg)', ytitle='\delta_{J2000} (deg)', $
         position=pos[*,2], /noerase, yminor=3
       djs_oplot, nsa[these[mem]].ra, nsa[these[mem]].dec, $
         psym=8, color=im_color('tan')
       djs_oplot, nsa[these[mem2]].ra, nsa[these[mem2]].dec, $
         psym=8, color=im_color('dodger blue')
       
; write out
       nsafile = lcspath+strlowcase(cl[ii].cluster)+'_nsa.fits'
       out = struct_addtags(replicate({object_position: 0L},nmem),nsa[these[mem]])
       out.object_position = these[mem]
       im_mwrfits, out, nsafile, clobber=clobber
       
    endfor
    im_plotconfig, psfile=psfile, /psclose, /pdf

return
end
    

    
