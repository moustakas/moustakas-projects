function streams_mge, image, badpixels=badpixels, pixscale=pixscale, $
  model=model, twist=twist, plot=plot, minlevel=minlevel, $
  bulge_disk=bulge_disk, sector_width=sector_width
; fit the BCG using MGE; ToDo: fit the PSF using an MGE

    if n_elements(pixscale) eq 0 then pixscale = 1.0
    if n_elements(minlevel) eq 0 then minlevel = 0.0
    
    forward_function multi_gauss, multi_gauss_twist
    resolve_routine, 'mge_print_contours', /compile_full_file, /no_recompile
    resolve_routine, 'mge_print_contours_twist', /compile_full_file, /no_recompile

; find the galaxy              
    find_galaxy, image, size, ellipticity, posangle, xcen, ycen, $
      xcen_lum, ycen_lum, fraction=fraction, index=index, level=level, $
      nblob=nblob, plot=plot;, /quiet

stop    
    
; get the photometry
    if keyword_set(twist) then begin
       sectors_photometry_twist, image, ellipticity, xcen, ycen, radius, $
         phi, counts, n_sectors=n_sectors, sector_width=sector_width, $
         badpixels=badpixels, minlevel=minlevel
    endif else begin
       sectors_photometry, image, ellipticity, posangle, xcen_lum, ycen_lum, $
         radius, phi, counts, n_sectors=n_sectors, sector_width=sector_width, $
         badpixels=badpixels, minlevel=minlevel
    endelse

; do the MGE fitting
    if keyword_set(twist) then begin
       mge_fit_sectors_twist, radius, phi, counts, ellipticity, ngauss=ngauss, $
         scale=pixscale, negative=negative, normpsf=normpsf, print=print, $
         qbounds=qbounds, rbounds=rbounds, sigmapsf=sigmapsf, $
         outer_slope=outer_slope, sol=sol, absdev=absdev
;      mge_print_contours_twist, image, ang, xpeak, ypeak, sol, model=model
       model = multi_gauss_twist(sol,image,0.0,1.0,xcen,ycen,posangle)
       model1 = multi_gauss_twist(sol,image,0.0,1.0,xcen,ycen,90) ; aligned left-right
    endif else begin
       mge_fit_sectors, radius, phi, counts, ellipticity, ngauss=ngauss, $
         scale=pixscale, bulge_disk=bulge_disk, fastnorm=fastnorm, $
         linear=linear, negative=negative, normpsf=normpsf, print=print, $
         qbounds=qbounds, rbounds=rbounds, /quiet, sigmapsf=sigmapsf, $
         outer_slope=outer_slope, sol=sol, absdev=absdev
;      mge_print_contours, cutimage, ang, xpeak, ypeak, sol, model=model
;      model = multi_gauss(sol,image,sigmapsf,normpsf,xcen,ycen,posangle)
       model = float(multi_gauss(sol,image,0.0,1.0,xcen,ycen,posangle))
       model1 = multi_gauss(sol,image,0.0,1.0,xcen,ycen,90) ; aligned left-right
    endelse
    
;; get the surface brightness profile of the galaxy and the model along
;; the semi-major and semi-minor axes
;    sz = size(image,/dim)
;    raxis = [0.0,range(pixscale,1.5*size,200,/log)] ; position along the major,minor axes [pixels]
;    if keyword_set(twist) then begin
;       splog, 'Code me!'
;    endif else begin
;; data
;       mjr = where(phi eq 0.0,nmjr)
;       mnr = where(phi eq 90.0,nmnr)
;       sb_major = interpolate(counts[mjr],findex(radius[mjr],raxis),missing=0.0)
;       sb_minor = interpolate(counts[mnr],findex(radius[mnr],raxis),missing=0.0)
;; model
;       sbmodel_major = reform(interpolate(model1[xcen:sz[0]-1,ycen],$
;         findex(range(xcen,sz[0]-1,sz[0]-xcen),xcen+raxis)))
;       sbmodel_minor = reform(interpolate(model1[xcen,ycen:sz[1]-1],$
;         findex(range(ycen,sz[1]-1,sz[1]-ycen),ycen+raxis)))
;
;       djs_plot, raxis, sb_major, psym=8, /xlog, /ylog, $
;         xrange=[raxis[1]<min(radius[mjr]),max(raxis)], $
;         yrange=[minlevel*0.9,max(counts)], xsty=3, ysty=3
;       djs_oplot, raxis, sbmodel_major, color='orange'
;       djs_oplot, 10^!x.crange, minlevel*[1,1], line=0
;stop       
;    endelse

    mge = {$
      size:        size*pixscale,$ ; approximate galaxy "size" [arcsec]
      ellipticity: ellipticity,$   ; average ellipticity = 1-b/a
      posangle:    posangle,$      ; position angle measured from the image Y-axis
      xcen:        xcen,$          ; galaxy x,y center [integer pixels]
      ycen:        ycen,$   
      xcen_lum:    xcen_lum,$      ; luminosity-weighted x,y centroid 
      ycen_lum:    ycen_lum,$
      absdev:      absdev,$        ; mean absolute deviation of the data from the model
      sol:         sol};,$

;     raxis:       raxis*pixscale,$ ; [arcsec]
;     sb_minor:      float(sb_minor),$
;     sb_major:      float(sb_major),$
;     sbmodel_minor: float(sbmodel_minor),$
;     sbmodel_major: float(sbmodel_major)}
;   mge.sbprofile = sol[0,*]/(2.0*!dpi*sol[1,*]^2*sol[2,*])

return, mge
end
    
