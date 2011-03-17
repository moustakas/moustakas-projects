;+
; NAME:
;   VAGC_GALEX_WINDOW
;
; PURPOSE:
;   Build the union of the GALEX and the VAGC window function for
;   various LSS subsamples.
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;   sample, letter, poststr - VAGC/LSS sample (default: sample='dr72',
;     letter='bsafe', poststr='25')
;
; KEYWORD PARAMETERS: 
;   rebuild_galex_window - rebuild the window function for the full
;     set of GALEX GR4/GR5 tiles
;
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 May 03, UCSD
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

pro galex_window, polyfile
; internal support routine to build the GALEX GR4/GR5 window function
; *once* 
    
    common vagc_galex_win, tiles1, galex

    if (n_elements(galex) eq 0L) then $
      galex = mrdfits(vagcpath+'object_galex_gr45.fits.gz',1)
    if (n_elements(tiles1) eq 0L) then $
      tiles1 = galex_read_gr45_tileinfo()

    splog, 'Building GALEX GR4/5 window function'
    t0 = systime(1)

; select the tiles that overlap with the SDSS
    good = where(galex.tile ne -999)
    alltiles = galex[good].tile
    tiles = alltiles[uniq(alltiles,sort(alltiles))]

    match, tiles1.tile_num, tiles, m1, m2
    tileinfo = tiles1[m1]
    ntiles = n_elements(tileinfo)

; the GALEX FOV is a circle that is 0.61 deg in radius
    radius = 0.61D
    circles = construct_polygon(ncaps=1,nelem=ntiles)
    for ii = 0, ntiles-1 do begin 
       circles[ii].caps = ptr_new(circle_cap(ra=tileinfo[ii].tile_ra,$
         dec=tileinfo[ii].tile_dec,radius))
       circles[ii].ncaps = 1 
       circles[ii].use_caps = 1L
       circles[ii].str = garea(circles[ii])
       circles[ii].weight = 1.0
    endfor

    tempfile = repstr(polyfile,'.ply','.tmp.ply')
    write_mangle_polygons, tempfile, circles
    destruct_polygon, circles
    
    spawn, 'pixelize '+tempfile+' '+tempfile, /sh
    spawn, 'snap '+tempfile+' '+tempfile, /sh
    spawn, 'balkanize '+tempfile+' '+tempfile, /sh
    spawn, 'unify '+tempfile+' '+tempfile, /sh
    spawn, 'sed 2,4d '+tempfile+' > '+polyfile, /sh
    spawn, '/bin/rm -f '+tempfile, /sh
    splog, 'Total time (minutes) = ', (systime(1)-t0)/60.0

    splog, 'Total area (deg^2) = ', total(balkans.str)*(180/!dpi)^2
;   read_mangle_polygons, polyfile, balkans
;   plot_poly, balkans
;   destruct_polygon, balkans

return
end

pro vagc_galex_window, sample=sample, letter=letter, poststr=poststr, $
  rebuild_galex_window=rebuild_galex_window, doplot=doplot
; jm10may03ucsd - build the GALEX/VAGC window function

    vagcpath = getenv('VAGC_REDUX')+'/'
    galexpolyfile = vagcpath+'galex/galex_gr45.ply'
    if keyword_set(rebuild_galex_window) then $
      galex_window, galexpolyfile
    
    if (n_elements(sample) eq 0) then sample = 'dr72'
    if (n_elements(letter) eq 0) then letter = 'bsafe'
    if (n_elements(poststr) eq 0) then poststr = '25'

    lsspath = getenv('LSS_REDUX')+'/'+sample+'/'+letter+'/'+poststr+'/lss/' 
    windowfile = lsspath+'window.'+sample+letter+poststr+'.ply'
    maskfile = lsspath+'mask.'+sample+letter+poststr+'.ply'

; output file    
    polyfile = vagcpath+'galex/'+'galex_gr45.'+sample+letter+poststr+'.ply' 
    polyfile1 = polyfile+'1'

; first combine the GALEX and VAGC window polygons
    splog, 'Reading '+galexpolyfile
    read_mangle_polygons, galexpolyfile, galexpoly
    splog, 'Reading '+windowfile
    read_mangle_polygons, windowfile, windowpoly

    tempfile = repstr(polyfile1,'.ply1','.tmp.ply1')
    write_mangle_polygons, tempfile, [galexpoly,windowpoly]
    spawn, 'pixelize '+tempfile+' '+tempfile, /sh
    spawn, 'snap '+tempfile+' '+tempfile, /sh
    spawn, 'balkanize '+tempfile+' '+tempfile, /sh
    spawn, 'sed 2,4d '+tempfile+' > '+polyfile1, /sh
    spawn, '/bin/rm -f '+tempfile, /sh

; next, add (balkanize) the mask polygons; keep polygons that are
; within the field boundaries but that are *outside* the mask polygons
    splog, 'Reading '+maskfile
    read_mangle_polygons, maskfile, maskpoly
    splog, 'Reading '+polyfile1
    read_mangle_polygons, polyfile1, poly1

    tempfile = repstr(polyfile1,'.ply1','.tmp.ply1')
    write_mangle_polygons, tempfile, [poly1,maskpoly]
    spawn, 'pixelize '+tempfile+' '+tempfile, /sh
    spawn, 'snap '+tempfile+' '+tempfile, /sh
    spawn, 'balkanize '+tempfile+' '+tempfile, /sh
    spawn, 'sed 2,4d '+tempfile+' > '+polyfile1, /sh
    spawn, '/bin/rm -f '+tempfile, /sh

    read_mangle_polygons, polyfile1, balkans
    wxyz = vmid(balkans,ra=wra,dec=wdec)
    isgood = is_in_window(poly1,ra=wra,dec=wdec)    ; good
    masked = is_in_window(maskpoly,ra=wra,dec=wdec) ; masked
    final = where((isgood gt 0) and (masked eq 0),nfinal)
    write_mangle_polygons, tempfile, balkans[final]

; run through mangle one last time and write out the final polygon
; file 
    spawn, 'pixelize '+tempfile+' '+tempfile, /sh
    spawn, 'snap '+tempfile+' '+tempfile, /sh
    spawn, 'balkanize '+tempfile+' '+tempfile, /sh
    spawn, 'sed 2,4d '+tempfile+' > '+polyfile, /sh
    spawn, '/bin/rm -f '+tempfile, /sh
    spawn, '/bin/rm -f '+polyfile1, /sh
    
return
end
    
