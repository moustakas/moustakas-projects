pro cmdm_cutouts
; jm12feb27ucsd - get cutouts of all the cluster members
    
    path = clash_path()+'projects/cmdm/'

    cl = 'macs1206'
    cat = mrdfits(path+cl+'_cat.fits.gz',1)
    ngal = n_elements(cat)
    
; get the cutouts    
    pngfile = path+'data/MACSJ1206_BVRIZ.png'
    imfile = path+'data/macsj1206_2009_IC_O1.fits'
    weightfile = path+'data/macsj1206_2009_IC_O1.weight.fits'
    
;   png = read_png(pngfile)
    im = mrdfits(imfile,0,hdr,/silent)
    weight = mrdfits(weightfile,0,weighthdr,/silent)
    sz = size(im,/dim)
    extast, hdr, astr
    pixsize = sqrt(abs(determ(astr.cd)))*3600

; output image files    
    outimfile = path+'data/cutouts/macs1206_'+string(cat.id,format='(I5.5)')+'.fits'
    outweightfile = repstr(outimfile,'.fits','.weight.fits')
    
    for ii = 0L, ngal-1 do begin
       diam = ceil(1.5*2.0*cat[ii].flux_radius2)>10L
;      diam = ceil(10*sqrt(cat[ii].awin_image*cat[ii].bwin_image))
;      diam = (ceil(2*cat[ii].flux_radius))
;      diam = (ceil(2*sqrt(cat[ii].area/!pi)))>20 ; [pixel]
       ad2xy, cat[ii].xwin_world, cat[ii].ywin_world, astr, xx, yy
;      ad2xy, cat[ii].ra, cat[ii].dec, astr, xx, yy
       x0 = long(xx-diam/2.0)>0L
       x1 = long(xx+diam/2.0)<(sz[0]-1)
       y0 = long(yy-diam/2.0)>0
       y1 = long(yy+diam/2.0)<(sz[1]-1)
       splog, diam, x0, x1, y0, y1

; image       
       hextract, im, hdr, newim, newhdr, x0, x1, y0, y1, /silent
       mwrfits, newim, outimfile[ii], newhdr, /silent, /create

; weight map       
       hextract, weight, weighthdr, newweight, newweighthdr, x0, x1, y0, y1, /silent
       mwrfits, newweight, outweightfile[ii], newweighthdr, /silent, /create
    endfor

return
end
    
