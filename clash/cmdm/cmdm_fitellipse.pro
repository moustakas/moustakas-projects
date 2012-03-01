pro cmdm_fitellipse
; jm12feb27ucsd - fit ellipse to everything
    
    path = clash_path()+'projects/cmdm/'

    cl = 'macs1206'
    cat = mrdfits(path+cl+'_cat.fits.gz',1)
    ngal = n_elements(cat)
    
    imfile = file_search(path+'data/cutouts/macs1206_'+string(cat.id,format='(I5.5)')+'.fits')
    weightfile = repstr(imfile,'.fits','.weight.fits')

;   for ii = 1, ngal-1 do begin
    for ii = 0L, ngal-1 do begin
       splog, file_basename(imfile[ii])
       im = mrdfits(imfile[ii],0,hdr)
       ivar = 1.0/im
;      ivar = mrdfits(weightfile[ii],0)

       namax = ceil(2*cat[ii].flux_radius)
;      namax = ceil(5*sqrt(cat[ii].awin_image*cat[ii].bwin_image))
       cmdm_ellipse, im, imivar=ivar, ellipse=ell, namax=namax
       bgt_ellipse_sersic, ell, outellipse=final
;      bgt_ellipse_show, im, final
;      cc = get_kbrd(1)
    endfor
    
stop    

return
end
    
