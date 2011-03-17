pro mosaic_sexcat2xy

    rootpath = '/global/bias1/ioannis/mosaic/'
    datapath = rootpath+'data/'

;   catlist = file_search(datapath+'fobj125.cat')
    catlist = file_search(datapath+'fobj???.cat')
    fitslist = repstr(catlist,'.cat','.fits')
    reglist = repstr(catlist,'.cat','.reg')
;   for icat = 0L, 1L do begin
    for icat = 0L, n_elements(catlist)-1L do begin
       for iext = 1L, 8L do begin
          hdr = headfits(fitslist[icat],ext=iext)
          extast, hdr, astr
          cat = mrdfits(catlist[icat],2*iext,/silent)
;         stars = lindgen(n_elements(cat))
          stars = where(cat.flux_auto/cat.fluxerr_auto gt 10.0,nstars)
          xy2ad, cat[stars].xwin_image, cat[stars].ywin_image, astr, aa, dd
          xylist1 = replicate({x_world: 0.0, y_world: 0.0},nstars) & xylist1.x_world = aa & xylist1.y_world = dd
          if (iext eq 1L) then xylist = xylist1 else xylist = [xylist,xylist1]
;         reglist = repstr(catlist[icat],'.cat','')+'.ext'+string(iext,format='(I0)')+'.reg'
;         splog, 'Writing '+reglist
;         struct_print, struct_trimtags(cat[stars],select=['XWIN_IMAGE','YWIN_IMAGE']), file=reglist, /no_head
       endfor
       splog, 'Writing '+reglist[icat]
       struct_print, xylist, file=reglist[icat], /no_head
    endfor
    
return
end
    
