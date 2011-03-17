pro oii_nondetections, postscript=postscript
; jm05jun09uofa - 

    if keyword_set(postscript) then begin
       postthick = 8.0 
    endif else begin
       postthick = 2.0
       im_window, 0, xratio=0.6, /square
    endelse

;   pspath = atlas_path(/papers)+'sfrs/FIG_SFRS/'
    pspath = atlas_path(/projects)+'sfrs/'
    
    psname = 'oii_nondetections'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=7.1, height=8.5, $
      xmargin=[1.1,0.3], ymargin=[0.3,1.1], xpage=8.5, ypage=8.5, $
      position=pos, /normal
   
    atlas = read_integrated(linefitnodust=atlasnodust,niikeep=niikeep,hakeep=hakeep,$
      hbkeep=hbkeep,oiikeep=oiikeep,oiiikeep=oiiikeep)
    nfgs = read_nfgs(linefitnodust=nfgsnodust,niikeep=nniikeep,hakeep=nhakeep,$
      hbkeep=nhbkeep,oiikeep=noiikeep,oiiikeep=noiiikeep)

    indx1 = cmset_op(cmset_op(cmset_op(hakeep,'AND',hbkeep),'AND',where(atlasnodust.ebv_hahb_err gt 0.0)),$
      'AND',oiikeep)
    indx2 = cmset_op(cmset_op(cmset_op(hakeep,'AND',hbkeep),'AND',$
      where(atlasnodust.ebv_hahb_err gt 0.0)),'AND',/not2,oiikeep)

    nindx1 = cmset_op(cmset_op(cmset_op(nhakeep,'AND',nhbkeep),'AND',where(nfgsnodust.ebv_hahb_err gt 0.0)),$
      'AND',oiikeep)
    nindx2 = cmset_op(cmset_op(cmset_op(nhakeep,'AND',nhbkeep),'AND',$
      where(nfgsnodust.ebv_hahb_err gt 0.0)),'AND',/not2,noiikeep)

    x = [atlas[indx1].M_B,nfgs[nindx1].M_B]
    y = [atlasnodust[indx1].ebv_hahb,nfgsnodust[nindx1].ebv_hahb]
    
    xnon = [atlas[indx2].M_B,nfgs[nindx2].M_B]
    ynon = [atlasnodust[indx2].ebv_hahb,nfgsnodust[nindx2].ebv_hahb]

    xtitle = 'M_{B}'
    ytitle = 'E(B-V)'

    xrange = [-11,-24]
    yrange = [-0.05,1.4]

    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
      charsize=1.8, charthick=postthick, xthick=postthick, ythick=postthick, $
      xsty=1, ysty=1, xtitle=xtitle, ytitle=ytitle, position=pos[*,0]
    plotsym, 0, 0.7, /fill & djs_oplot, x, y, ps=8, color='red'
    plotsym, 8, 1.2, /fill & djs_oplot, xnon, ynon, ps=8, color='blue'
    
    im_openclose, postscript=postscript, /close    

return
end
