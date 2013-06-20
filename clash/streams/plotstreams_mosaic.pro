pro plotstreams_mosaic
; jm13jun10siena - show the Subaru + HST imaging for MACS1206

    path = streams_path()
    paperpath = streams_path(/paper)
    
    cluster = 'macs1206' ; start with one cluster
    info = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    info = info[where(cluster eq strtrim(info.shortname,2))]
    arcsec2kpc = dangular(info.z,/kpc)/206265D ; [kpc/arcsec]

; --------------------------------------------------
; build the HST color
    outfile = paperpath+'hst_mosaic.eps'

; read the HST color image and astrometric header    
    rgbfile = path+'macs1206.png'
    splog, 'Reading '+rgbfile
    png = read_png(rgbfile,/verbose)
    hdr = headfits(path+'macs1206.hdr.fits.gz')
    extast, hdr, astr

; read the WFC3/IR footprint region files
    wfc3regfile = path+'macs1206_f160w.reg'
    wfc3 = ds9polygon_indices(wfc3regfile,header=hdr,$
      xvert=wxvert,yvert=wyvert)
;   xywfc3 = array_indices(fltarr(astr.naxis),wfc3)

    png0 = bytarr(astr.naxis)+255
    png1 = bytarr(astr.naxis)+255
    png2 = bytarr(astr.naxis)+255
    png0[wfc3] = (png[0,*,*])[wfc3]
    png1[wfc3] = (png[1,*,*])[wfc3]
    png2[wfc3] = (png[2,*,*])[wfc3]

    newpng = transpose([[[png0]],[[png1]],[[png2]]],[2,0,1])
    
; zoom into the BCG    
    pixscale = sqrt(abs(determ(astr.cd)))*3600D ; [arcsec/pixel]
;   ad2xy, info.ra_bcg, info.dec_bcg, astr, xcen, ycen
;   xcen = long(xcen) & ycen = long(ycen)
;   
;   xdd = 1000D
;   ydd = 500D
;   xdiam = long(xdd/arcsec2kpc/pixscale) ; [DD Mpc in pixels]
;   ydiam = long(ydd/arcsec2kpc/pixscale) ; [DD Mpc in pixels]

    xdiam = long(max(wxvert)-min(wxvert))+10
    ydiam = long(max(wyvert)-min(wyvert))+10
    xcen = xdiam/2+long(min(wxvert))
    ycen = ydiam/2+long(min(wyvert))

    x0 = xcen-xdiam/2
    x1 = xcen+xdiam/2
    y0 = ycen-ydiam/2
    y1 = ycen+ydiam/2
;   image = png[*,x0:x1,y0:y1]
    image = newpng[*,x0:x1,y0:y1]
    sz = size(image,/dim)
    
    ps_start, outfile, /encapsulated
    cgimage, image, /save, /keep_aspect

    nv = n_elements(wxvert)
    wxv = [wxvert[0:nv-1],wxvert[0]]-(xcen-xdiam/2)
    wyv = [wyvert[0:nv-1],wyvert[0]]-(ycen-ydiam/2)
    for ii = 0, nv-1 do cgoplot, [wxv[ii],wxv[ii+1]], $
      [wyv[ii],wyv[ii+1]], line=0, color='red', thick=7
    ps_end, /png

stop
    
; --------------------------------------------------
; build the Subaru color mosaic
    outfile = paperpath+'subaru_mosaic.eps'
    
; read the Subaru color image and astrometric header    
    rgbfile = path+'MACSJ1206_BVRIZ.png'
    splog, 'Reading '+rgbfile
    png = read_png(rgbfile,/verbose)
    hdr = headfits(path+'MACSJ1206_BVRIZ.hdr.fits.gz')
    extast, hdr, astr

; read the WFC3/IR and ACS footprint region files
    wfc3regfile = path+'macs1206_f160w.reg'
    wfc3 = ds9polygon_indices(wfc3regfile,header=hdr,$
      xvert=wxvert,yvert=wyvert)
    
    acsregfile = path+'macs1206_f814w.reg'
    acs = ds9polygon_indices(acsregfile,header=hdr,$
      xvert=axvert,yvert=ayvert)
    
; extract the center of the mosaic centered on the BCG
    pixscale = sqrt(abs(determ(astr.cd)))*3600D ; [arcsec/pixel]
    ad2xy, info.ra_bcg, info.dec_bcg, astr, xcen, ycen
    xcen = long(xcen) & ycen = long(ycen)
    
    dd = 2000D ; [=2 Mpc]
    diam_deg = dd/arcsec2kpc/3600D ; [2 Mpc in degrees]
    diam = long(dd/arcsec2kpc/pixscale) ; [2 Mpc in pixels]

    x0 = xcen-diam/2
    x1 = xcen+diam/2
    y0 = ycen-diam/2
    y1 = ycen+diam/2
    image = png[*,x0:x1,y0:y1]
    sz = size(image,/dim)

; get the 2MASS stars in the field
;   star = sdss_sweep_circle(181.55082D,-8.8009182D,0.1D,type='star')
    star = im_read_tbl(path+'macs1206-2mass-stars.tbl')
    bright = where(star.k_m lt 14.0,nstar)
    ad2xy, star[bright].ra, star[bright].dec, astr, xstar, ystar
    xstar = xstar-(xcen-diam/2)
    ystar = ystar-(ycen-diam/2)

    keep = where(xstar ge 0 and ystar ge 0 and xstar le sz[1]-1 and ystar le sz[2]-1,nstar)
    xstar = xstar[keep]
    ystar = ystar[keep]
    kmag = star[bright[keep]].k_m
    jmag = star[bright[keep]].j_m
    splog, 'Kmag and J-K of the brightest stars'
    niceprint, kmag, jmag-kmag
    
; render the plot
    ps_start, outfile, /encapsulated
    cgimage, image, /save

    for ii = 0, nstar-1 do tvcircle, 50.0, xstar[ii], $
      ystar[ii], /data, color=cgcolor('white')

    nv = n_elements(wxvert)
    wxv = [wxvert[0:nv-1],wxvert[0]]-(xcen-diam/2)
    wyv = [wyvert[0:nv-1],wyvert[0]]-(ycen-diam/2)
    for ii = 0, nv-1 do cgoplot, [wxv[ii],wxv[ii+1]], $
      [wyv[ii],wyv[ii+1]], line=0, color='red', thick=7
    
    nv = n_elements(axvert)
    axv = [axvert[0:nv-1],axvert[0]]-(xcen-diam/2)
    ayv = [ayvert[0:nv-1],ayvert[0]]-(ycen-diam/2)
    for ii = 0, nv-1 do cgoplot, [axv[ii],axv[ii+1]], $
      [ayv[ii],ayv[ii+1]], line=0, color='yellow', thick=7
    
    xoff = 10
    yoff = 100
    cgarrow, sz[1]/2, yoff, sz[1]-xoff, yoff, /data, hthick=3, thick=5, /solid, color='white'
    cgarrow, sz[1]/2, yoff, xoff, yoff, /data, hthick=3, thick=5, /solid, color='white'
    cgtext, sz[1]/2, yoff*1.5, '2 Mpc', align=0.5, color='white'
    ps_end, /png
    

return
end
    

