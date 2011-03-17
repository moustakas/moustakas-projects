pro lbt_sexcat2xy

    rootpath = '/global/bias1/ioannis/lbt/'
    datapath = rootpath+'data/'

    catlist = [file_search(datapath+'p94_*_?.cat'),file_search(datapath+'p94_*_10.cat')]
;   for icat = 0L, 1L do begin
    for icat = 0L, n_elements(catlist)-1L do begin
       for iext = 1L, 4L do begin
          cat = mrdfits(catlist[icat],2*iext,/silent)
;         stars = lindgen(n_elements(cat))
          stars = where(cat.flux_auto/cat.fluxerr_auto gt 10.0)
;         stars = where(cat.class_star gt 0.8)
          reglist = repstr(catlist[icat],'.cat','')+'.ext'+string(iext,format='(I0)')+'.reg'
          splog, 'Writing '+reglist
          struct_print, struct_trimtags(cat[stars],select=['XWIN_IMAGE','YWIN_IMAGE']), file=reglist, /no_head
       endfor
    endfor
    
return
end
    
