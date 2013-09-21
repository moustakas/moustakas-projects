pro streams_find_bcg, debug=debug
; jm13sep04siena - find the BCGs in the CLASH clusters; get the
; central coordinates and the basic parameters (ellipticity, position
; angle, etc.) that we'll need to do MGE modeling later 

; note! images in units of [10^-12 erg/s/cm^2/Hz] (pico-maggies)
    
; for no particular reason, find the BCGs the full sample of CLASH 
; clusters, not just the STREAMS project ones    
;   sample = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    sample = rsex(streams_path(/propath)+'streams_sample.sex')
    struct_print, sample
    ncl = n_elements(sample) 

    pixscale = 0.065D ; [arcsec/pixel]
    rmaxkpc = 40D     ; [kpc]    

    outfile = streams_path()+'bcg_info.fits'

; wrap on each cluster    
    for ic = 0, ncl-1 do begin
;   for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       splog, 'Working on cluster '+cluster
       skypath = streams_path(/skysub)+cluster+'/'
       
; get the preliminary coordinates of the BCG(s) and the physical
; scale; CLJ1226 is affected by a bright neighbor; also note that the
; coordinates for the merging clusters are somewhat arbitrary
       rmaxkpc = 40D
       if cluster eq 'clj1226' then rmaxkpc = 10D
       if cluster eq 'macs2129' then rmaxkpc = 8
       bcg_ra1 = 15D*hms2dec(sample[ic].ra)
       bcg_dec1 = hms2dec(sample[ic].dec)
       
       arcsec2kpc = dangular(sample[ic].z,/kpc)/206265D ; [kpc/arcsec]
       rmax = ceil(rmaxkpc/arcsec2kpc/pixscale)         ; [pixels]
       rmax = rmax+odd(rmax)                            ; force even

; read the skyinfo structure to get the filters
       skyinfo = mrdfits(streams_path(/skysub)+'skyinfo-'+$
         cluster+'.fits.gz',1,/silent)
       short = strtrim(skyinfo.band,2)
       reffilt = where(short eq 'f160w') ; reference filter
       nfilt = n_elements(skyinfo)
    
; read the reference image
       imfile = skypath+cluster+'-'+short[reffilt]+'.fits.gz'
       splog, 'Reading '+imfile
       rimage = mrdfits(imfile,0,rhdr,/silent)
       rinvvar = mrdfits(imfile,1,ivarrhdr,/silent)

; get a RMAXKPC by RMAXKPC cutout centered on the BCG so that we can
; find it; iterate twice
       for iter = 0, 1 do begin
          adxy, rhdr, bcg_ra1, bcg_dec1, xcen1, ycen1
          xcen1 = fix(xcen1) & ycen1 = fix(ycen1)
          hextract, rimage, rhdr, cutrimage, cutrhdr, xcen1-rmax, $
            xcen1+rmax, ycen1-rmax, ycen1+rmax, /silent

          find_galaxy, cutrimage, size, ellipticity, posangle, xcen, ycen, $
            xcen_lum, ycen_lum, fraction=fraction, index=index, level=level, $
            nblob=nblob, plot=plot, /quiet

;         dpeaks, cutrimage, xcen=xx, ycen=yy, maxnpeaks=10
;         xx = xcen1-rmax+xx & yy = ycen1-rmax+yy
;         dist = min(sqrt((xx-xcen1)^2+(yy-ycen1)^2),this)

          xyad, rhdr, xcen_lum+xcen1-rmax, $
            ycen_lum+ycen1-rmax, bcg_ra1, bcg_dec1
       endfor

       drefine, cutrimage, xcen_lum, ycen_lum, $
         xr=xcen_final, yr=ycen_final, box=11L

;      adxy, rhdr, bcg_ra1, bcg_dec1, xcen1, ycen1
;      xcen1 = fix(xcen1) & ycen1 = fix(ycen1)
;      hextract, rimage, rhdr, cutrimage, cutrhdr, xcen1-rmax, $
;        xcen1+rmax, ycen1-rmax, ycen1+rmax, /silent
;      hextract, rinvvar, ivarrhdr, cutrinvvar, cutivarrhdr, $
;        xcen1-rmax, xcen1+rmax, ycen1-rmax, ycen1+rmax, /silent
;
;      adxy, rhdr, bcg_ra1, bcg_dec1, xcen1, ycen1
;      xcen1 = fix(xcen1) & ycen1 = fix(ycen1)
;      hextract, rimage, rhdr, cutrimage, cutrhdr, xcen1-rmax, $
;        xcen1+rmax, ycen1-rmax, ycen1+rmax, /silent
;
;      find_galaxy, cutrimage, size, ellipticity, posangle, xcen, ycen, $
;        xcen_lum, ycen_lum, fraction=fraction, index=index, level=level, $
;        nblob=nblob, plot=plot, /quiet
;      drefine, cutrimage, xcen_lum, ycen_lum, xr=xcen_final, yr=ycen_final, box=3L

       if keyword_set(debug) then begin
          delvarx, pp
          cgimage, cutrimage, stretch=2, /save, /keep_aspect, $
            margin=0, position=pp, clip=3
          plots, xcen_lum, ycen_lum, psym=7, color=cgcolor('green')
          plots, xcen_final, ycen_final, psym=8, color=cgcolor('red')
          xyouts, pp[0]+0.03, pp[3]-0.1, strupcase(cluster), /norm

          adxy, rhdr, 15D*hms2dec(sample[ic].ra), hms2dec(sample[ic].dec), xx, yy
          plots, xx-xcen1+rmax, yy-ycen1+rmax, psym=6, color=cgcolor('blue')
       endif

       ycen += ycen1-rmax
       xcen_lum += xcen1-rmax
       ycen_lum += ycen1-rmax
       xcen_final += xcen1-rmax
       ycen_final += ycen1-rmax
       xyad, rhdr, xcen_final[0], ycen_final[0], bcg_ra, bcg_dec
       
       info1 = {$
         cluster:   cluster,$
         z:    sample[ic].z,$
         ra:         bcg_ra,$
         dec:       bcg_dec,$
         xcen:     float(xcen_final[0]),$   ; luminosity-weighted x,y centroid 
         ycen:     float(ycen_final[0]),$
         size:        float(size*pixscale),$ ; approximate galaxy "size" [arcsec]
         size_kpc:    float(size*pixscale*arcsec2kpc),$ ; approximate galaxy "size" [kpc]
         ellipticity: float(ellipticity),$   ; average ellipticity = 1-b/a
         posangle:    float(posangle)}       ; position angle measured from the image Y-axis

       if keyword_set(debug) then begin
          help, info1, /str
          cc = get_kbrd(1)
       endif

       if n_elements(info) eq 0 then info = info1 else info = [info,info1]
    endfor                      ; close cluster loop

    im_mwrfits, info, outfile, /clobber
    
return
end
    
