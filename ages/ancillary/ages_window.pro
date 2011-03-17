;+
; NAME:
;   AGES_WINDOW
; PURPOSE:
;   Build the window function for AGES.
; INPUTS: 
;   None.
; KEYWORD PARAMETERS: 
; OUTPUTS: 
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 02, UCSD
;
; Copyright (C) 2010, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro ages_window, noplot=noplot, clobber=clobber

    catpath = ages_path(/catalogs)
    windowpath = ages_path(/window)
    scheme = '-Ps0,10'

; some catalogs we need/want    
    weight = mrdfits(catpath+'catalog.spectweight.fits.gz',1)
    usno = mrdfits(catpath+'catalog.usno.fits.gz',1)
    codes = mrdfits(catpath+'catalog.codes.fits.gz',1)
    sdss = mrdfits(catpath+'catalog.sdss_photometry.fits.gz',1)
    kcorr = mrdfits(catpath+'catalog.kcorr.v3.fits.gz',1)
    ages = mrdfits(catpath+'catalog.cat.noguidestars.fits.gz',1)
    ages.ra = ages.ra*15D

;;    win = im_is_in_window(windowpath+'ages_fields.ply',ra=ages.ra,dec=ages.dec)
;;    istar = (sdss.type eq 6) and (sdss.psfmag_r le 19.0)
;;    imain = (codes.gshort gt 0) and (codes.gshort ne 2048) and $
;;      (codes.bgood eq 1) and (codes.rgood eq 1) and (win eq 1) and $
;;      (istar eq 0) and (kcorr.i_obs le 19.95)
;;    w1 = where((imain eq 1) and (weight.main_weight eq 0))
;;    w2 = where((imain eq 0) and (weight.main_weight eq 1))
;;
;;    djs_plot, ages[where(imain)].ra, ages[where(imain)].dec, ps=3, ysty=3 ;, xr=[217,218]
;;    djs_oplot, ages[w1].ra, ages[w1].dec, ps=6, sym=2, color='red'
;;    djs_oplot, ages[w2].ra, ages[w2].dec, ps=6, sym=2, color='orange'
;;    
;;    ww = where(weight[main].main_weight eq 0)
;;    ww = where(weight[main].main_weight eq 0)
    
; --------------------------------------------------
; build polygons for the 2004 field centers; the following file gives
; the field centers of radius 0.49 deg
    fieldfile = windowpath+'ages_fields.ply'
    if (file_test(fieldfile) eq 0) or (keyword_set(clobber) eq 1) then begin
       readcol, ages_path(/catalogs)+'bootes_standard_fields.dat', $
         field, ra, dec, format='I,D,D', /silent, comment='#'
       nfield = n_elements(field)
       radius = 0.49D           ; [deg]

       circles = construct_polygon(ncaps=1,nelem=nfield)
       for ii = 0, nfield-1 do begin 
          circles[ii].caps = ptr_new(circle_cap(ra=ra[ii],dec=dec[ii],radius))
          circles[ii].ncaps = 1 
          circles[ii].use_caps = 1L
          circles[ii].str = garea(circles[ii])
          circles[ii].weight = 1.0
       endfor
       junkfile = '/tmp/junk.ply'
       write_mangle_polygons, junkfile, circles
       destruct_polygon, circles
       
       im_call_mangle, junkfile, outfile=fieldfile, /pixelize, $
         /snap, /balkanize, scheme=scheme

;      if (keyword_set(noplot) eq 0) then begin
;         read_mangle_polygons, fieldfile, fieldpoly
;         plot_poly, fieldpoly
;      endif
    endif

; --------------------------------------------------
; now build and write out the bright star mask
    bsfile = windowpath+'ages_bsmask.ply'
    if (file_test(bsfile) eq 0) or (keyword_set(clobber) eq 1) then begin

;; --------------------------------------------------       
;; the following is good code, but it's not what was used for
;; the targeting! 
;       usnofile = catpath+'ages_usno.fits.gz' ; see AGES_GET_USNO() 
;       splog, 'Reading '+usnofile
;       usno = mrdfits(usnofile,1)
;       keep = where(usno.rmag lt 15.0,nkeep)
;       usno = usno[keep]
;       srt = sort(usno.rmag)
;       usno = usno[srt]
;       ing = spheregroup(usno.ra,usno.dec,30.0/3600.0,$
;         firstg=firstg,multg=multg,nextg=nextg) ; remove duplicates    
;       firstg = firstg[0L:max(ing)]
;
;; choose the radius based on the R-band magnitude (see the 2005 AGES
;; documentation)     
;       usno = usno[firstg]
;       radius = 10.0+2.5*(15.0-usno.rmag) ; [arcsec]
;       nstar = n_elements(radius)

       rad = 10.0+2.5*(15.0-usno.mstar) ; [arcsec]
       stars = where((usno.dstar lt rad) and (usno.bstar eq 1),nstar)

; ra,dec and radius-of-influence of each USNO star       
       starradius = rad[stars] ; [arcsec]
       starra = ages[stars].ra+usno[stars].dra/3600.0
       stardec = ages[stars].dec+usno[stars].ddec/3600.0

;      spherematch, ra, dec, ages.ra, ages.dec, 0.5/3600.0, m1, m2
;      djs_plot, ra, dec, ps=6, xsty=3, ysty=3, color='blue'
;      djs_oplot, ages[stars].ra, ages[stars].dec, psym=6, sym=0.5, color='yellow'
;      djs_oplot, ra[m1], dec[m1], ps=7, color='red', sym=2

; build the polygons with a weight of zero
       delvarx, keep
       circles = construct_polygon(ncaps=1,nelem=nstar)
       for ii = 0L, nstar-1L do begin 
          circles[ii].caps = ptr_new(circle_cap(ra=starra[ii],$
            dec=stardec[ii],starradius[ii]/3600.0))
          circles[ii].ncaps = 1 
          circles[ii].use_caps = 1L
          circles[ii].str = garea(circles[ii])
          circles[ii].weight = 1.0
;;; identify all ages objects around this purported bright star; if any
;;; of them were observed and found to have a non-zero redshift, then
;;; get rid of this bright star polygon
;;          flag = where(is_in_window(ra=ages.ra,$
;;            dec=ages.dec,circles[ii]) eq 1,nflag)
;;          check = where((usno[flag].match eq 1) and $
;;            (weight[flag].redshift ne -1.0),ncheck)
;;          if (ncheck eq 0) then if (n_elements(keep) eq 0) then $
;;            keep = ii else keep = [keep,ii]
;;;         if (ncheck ne 0) then begin
;;;            struct_print, codes[flag]
;;;            struct_print, weight[flag]
;;;            struct_print, usno[flag]
;;;            stop
;;;         endif
       endfor
;      circles = circles[keep]

; run mangle
       junkfile = '/tmp/junk.ply'
       write_mangle_polygons, junkfile, circles
       destruct_polygon, circles

       im_call_mangle, junkfile, outfile=bsfile, /pixelize, $
         /snap, /balkanize, scheme=scheme

;      if (keyword_set(noplot) eq 0) then begin
;         read_mangle_polygons, bsfile, bspoly
;         plot_poly, bspoly
;      endif

;      djs_plot, ages.ra, ages.dec, ps=3, xsty=3, ysty=3
;      for ii=0L,nstar-1 do tvcircle, 10*radius[ii]/3600.0, $
;        ra[ii], dec[ii], color=djs_icolor('yellow'), /data
       
    endif

; --------------------------------------------------
; build the final window function by combining the field and
; bright-star polygons; put the bright-star mask last to ensure that
; those polygons contribute a weight of zero to the output window 
    
    finalfile = windowpath+'ages_window_final.ply'
    if (file_test(finalfile) eq 0) or (keyword_set(clobber) eq 1) then begin

; first combine the field and bright star polygons
       tempfile = repstr(finalfile,'.ply','.tmp.ply')
       im_call_mangle, [fieldfile,bsfile], outfile=tempfile, $
         /pixelize, /snap, /balkanize, scheme=scheme

; next, read the field and bright-star polygons back in and retain
; polygons that are within the field boundaries but that are *outside*
; the bright-star polygons, and then run mangle one last time
       read_mangle_polygons, tempfile, temppoly
       wxyz = vmid(temppoly,ra=wra,dec=wdec)

       field = im_is_in_window(fieldfile,ra=wra,dec=wdec,/silent) ; good
       bstar = im_is_in_window(bsfile,ra=wra,dec=wdec,/silent)    ; reject
       final = where((field gt 0) and (bstar eq 0),nfinal)
       write_mangle_polygons, tempfile, temppoly[final]

       im_call_mangle, tempfile, outfile=finalfile, $
         /pixelize, /snap, /balkanize, scheme=scheme
       spawn, '/bin/rm -f '+tempfile, /sh
       
; some numbers...
;      read_mangle_polygons, finalfile, finalpoly
;      print, total(finalpoly.str)*(180.0/!dpi)^2
    endif
       
return
end
    
