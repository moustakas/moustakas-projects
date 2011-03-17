
pro plot_smooth, $
                 setup_file, data_file, $
                 nshell=nshell, vshell=vshell, $
                 ncontours=nc, $
                 title=map_title, $
                 no_zero=no_zero, $
                 no_slabel=no_slabel, $
                 no_zlabel=no_zlabel, $
                 add_noise=add_noise, $
                 no_grid=no_grid, $
                 glinestyle=glinestyle, $
                 charsize=charsize, $
                 help=help

;+
; PLOT_SMOOTH
;   Plots a contour map of the smoothed overdensity field 
;
; INPUTS
;   setup_file -- File containing cell coordinates
;   data_file -- File containing the smoothed data array
;
; KEYWORDS
;   Note one of either nshell or vshell must be set
;   nshell -- Plot this radial shell (integer)
;   vshell -- Plot a shell at this velocity, interpolating between
;             shells at two neighbouring velocities
;   nc -- Approximate number of contours
;   title -- Map title
;   no_zero -- Do not draw contour at zero
;   no_slabel -- No statistics label
;   no_zlabel -- No redshift label
;   no_grid -- No map grid
;   glinestyle -- linestyle for grid (default 0)
;   add_noise -- Work around for buggy contour plotting
;
; HISTORY
;   1996oct09  JEB
;   1998aug25  JEB  Fixed up contours for plots with galactic 
;                   anti-center at the origin!
;-

if n_params() lt 2 or keyword_set( help ) then begin
    doc_library, 'plot_smooth'
    retall
endif

if not keyword_set( nshell ) and not keyword_set( vshell ) then begin
    doc_library, 'plot_smooth'
    retall
endif

if not keyword_set( glinestyle ) then glinestyle = 0
if not keyword_set( charsize ) then charsize = 1

; --Read the cell coordinates and data array
read_ors_bins, setup_file, v, mu, wmu, lon, nv, nmu, nphi
read_smooth_data, data_file, data, nv, nmu, nphi
v = v( 0:nv-1 )

lat = asin( mu )

; --Select radial shell number nshell
if keyword_set( nshell ) then begin
    if nshell lt 1 or nshell gt nv then begin
        print, format='(" Shell number must be between 1 and ",i0)', nv
        retall
    endif
    data = data( nshell-1, *, * )
    vshell = v( nshell-1 )

; --Linearly interpolate between two shells nearest the chosen velocity
endif else begin
    if vshell gt v( nv-1 ) then vshell = v( nv-1 )
    if vshell lt v( 0 ) then vshell = 0 
; --Find the lower shell
    nlow = 0
    while v( nlow+1 ) lt vshell do begin
        nlow = nlow + 1
    endwhile
; --Find the higher shell
    nhi = nlow
    while v( nhi ) lt vshell do begin
        nhi = nhi + 1
    endwhile
; --Set linear weights
    if ( vshell eq v( nhi )) then wlow = 0. $
    else if ( vshell eq v( nlow )) then wlow = 1. $
    else wlow = ( v( nhi ) - vshell )/( v( nhi ) - v( nlow ))
; --Interpolate
    data = wlow*data( nlow, *, * ) + ( 1. - wlow )*data( nhi, *, * )
; --Write shell numbers, velocities and weights
    fmt = '(" Shell number ", i0, ": v=", i0, " km/s,  wt=", f5.3)'
    if wlow gt 0 then print, format=fmt, nlow, fix( v( nlow )+0.5 ), wlow
    if wlow lt 1 then print, format=fmt, nhi, fix( v( nhi )+0.5 ), 1.-wlow
endelse

data = reform( data )

; --Now find the mean over the sky by doing Gauss-Legendre weighted sum
sum_phi = total( data, 2 )
mean = 1./( 2.*nphi ) * total( wmu * sum_phi )

; --Note that data are written (r, mu, phi), but we want (lon, lat), so 
; --we need to take the transpose
data = transpose( data )

; --Assume data are read in with lon=[0,360], wrap to [-180,180]
;nlon = nphi
;nlon2 = nlon/2
;lon = [ lon[nlon2:nlon-1], lon[0:nlon2-1] ]
;tmp = data
;tmp[ nlon2:nlon-1, * ] = data[ 0:nlon2-1, * ]
;tmp[ 0:nlon2-1, * ] = data[ nlon2:nlon-1, * ]
;data = tmp

; --Add an extra phi point to wrap with the first
nlon = nphi + 1
nlat = nmu
lon = [ lon, lon[0]+360 ]
tmp = fltarr( nlon, nlat )
tmp[ 0:nlon-2, * ] = data
tmp[ nlon-1, * ] = data[ 0, * ]
data = tmp

dmin = min( data )
dmax = max( data )
print, format='( " Data: min=",f10.3,",  max=",f10.3)', dmin, dmax
print, format='(" Mean over sky=", f10.3)', mean

; --Heavy contour at 0, dot-dashed under 0.  No 0 contour if no_zero set
set_contours, levs, dmin, dmax, nc=nc
indxlt0 = where( levs lt 0., nlt0 )
indxeq0 = where( levs eq 0., neq0 )
if keyword_set( no_zero ) and neq0 gt 0 then begin
    levs = levs( where( levs ne 0. ))
    neq0 = 0
endif
contsty = intarr( n_elements( levs ))
if nlt0 gt 0 then contsty( indxlt0 ) = 3
contthick = intarr( n_elements( levs )) + 1
if neq0 gt 0 then contthick( indxeq0 ) = 2
if !p.thick ne 0 then contthick = 0.5 * !p.thick * contthick
if !p.multi( 2 ) gt 1 then charsize = charsize * ( 0.7 )^( !p.multi( 2 ) - 1 )

; --Add a little noise to work around IDL bug
if keyword_set( add_noise ) then $
  data = data + 0.001*randomn( seed, nv, nmu, nphi )

; --Draw the map
indx = where( !p.multi(1:2) gt 1, count )
if count eq 0 then erase
aitoff_position, pos=pos, xmargin=[1,1], ymargin=[1,2]
if keyword_set( no_grid ) then begin
    map_set, /advance, /aitoff, /noborder, glinesty=glinestyle, $
      glinethick=1, title=map_title, pos=pos, xmargin=[1,1], ymargin=[1,2]
endif else begin
    map_set, /advance, /aitoff, /grid, /noborder, glinesty=glinestyle, $
      glinethick=1, title=map_title, pos=pos, xmargin=[1,1], ymargin=[1,2]
    map_label, charsize=1.5*charsize
endelse
draw_lon_line, 179.99, linestyle=glinestyle, thick=1
;draw_lat_line, 20, linesty=5, thick=2
;draw_lat_line, -20, linesty=5, thick=2

contour, data, 180-lon, lat, $
  /overplot, $
  levels=levs, $
  c_linestyle=contsty, c_thick=contthick, $
  c_charsize=charsize*0.75

; --Labels
if not keyword_set( no_slabel ) then begin
    xn = 0.01 + fltarr( 4 )
    yn = 0.12 - 0.03 * findgen( 4 )
    reg2dev, xn, yn, x1, y1
    xn = 0.8 + fltarr( 4 )
    reg2dev, xn, yn, x2, y2
    x1 = x1( 0 )
    x2 = x2( 0 )
    xyouts, x1, y1( 0 ), /dev, 'galactic', charsize=charsize*0.7
    if dmax lt 20 then begin
        fmt = '( "(\delta_{min}, \delta_{max}) = (", f6.3, ", ", f6.3, ")" )'
        fmt2 = '( "<\delta> = ", f6.3 )'
        fmt3 = '("contour spacing = ", f5.3)'
    endif else begin
        fmt = '( "(v_{min}, v_{max}) = (", i0, ", ", i0, ")" )'
        fmt2 = '( "<v> = ", i0 )'
        fmt3 = '("contour spacing = ", i0)'
    endelse
    str = string( format=fmt3, levs(1)-levs(0) )
    xyouts, x1, y1( 1 ), /dev, str, charsize=charsize*0.7
    str = textoidl( string( format=fmt, dmin, dmax ))
    xyouts, x1, y1( 2 ), /dev, str, charsize=charsize*0.7
    str = textoidl( string( format=fmt2, mean ))
    xyouts, x1, y1( 3 ), /dev, str, charsize=charsize*0.7
endif
if not keyword_set( no_slabel ) then begin
    fmt = '( "cz = ", i0, " km s^{-1}" )'
    str = textoidl( string( format=fmt, fix( vshell + 0.5 ))) 
    xyouts, x2, y2( 2 ), /dev, str, charsize=charsize
endif
print, " Contour spacing: ", levs(1)-levs(0)

end
