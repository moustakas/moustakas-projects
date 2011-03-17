
pro plot_neighbors, file, symbol=sym, symsize=ssize, nolabels=nolabels, $
                    nocolor=nocolor, colors=clr, columns=col, help=help
;+
; NAME:
;    plot_neighbors
;
; PURPOSE:
;    Plots galaxy positions over a pre-existing map.  Different symbols
;    (and colors) are used depending on the value of a flag.  Color table
;    39 yields good results
;
; CATEGORY:
;    ORS -- TRGB
;
; CALLING SEQUENCE:
;    plot_neighbors(, file[, symbol, nolabels, nocolor, colors])
; 
; INPUTS:
;
; OPTIONAL INPUTS:
;    file    Filename containing galaxy (l, b) coordinates and a flag 
;            with value 0, 1, or 2.  Default file name is 
;            '/coma8/jbaker/TRGB/neighbor.lis'
;	
; KEYWORD PARAMETERS:
;    symbol -- Three-element vector containing plotting symbols to use
;              for flag=0, 1, 2
;    symsize -- Three-element vector giving symbol sizes
;    nolabels -- Set this keyword to suppress labels
;    nocolor -- Set this keyword to use one color only
;    colors -- A 3-element vector with colors for different flags
;    columns -- A 3-element vector giving index of l, b, flag in input file 
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;    1996sep10  JEB
;-

if keyword_set( help ) then begin
    retall
    doc_library, 'plot_neighbors'
endif

; --Try to find an input file
if n_elements( file ) eq 0 then file = '/coma8/jbaker/TRGB/neighbor.dat'
if n_elements( col ) eq 0 then col = [ 4, 5, 8 ]

; --Set plotting parameters
if not keyword_set( sym ) then sym = [ 5, 6, 2, 1, 4]
if not keyword_set( ssize ) then ssize = [ 0.5, 1, 1, 2, 2]
th = 1 
if !p.thick ne 0 then th = 0.5*!p.thick
if keyword_set( nocolor ) then begin
    clr = intarr( 5 )
    if !d.name ne 'PS' then clr = clr + !d.n_colors - 1
endif else if not keyword_set( clr ) then $
  clr = [ 1., 0.661, 0.767, 1., 1.] * (!d.n_colors - 1)

; --Read the data
;read_data, file, x, col=[ 1, 2, 5 ]
read_data, file, x, col = col

; --Plot galaxies not to be observed
indx = where( x(*, 2) eq 0, count )
if count gt 0 then begin
    gl = x( indx, 0 )
    gb = x( indx, 1 )
    oplot, 180-gl, gb, psym=sym(0), symsize=ssize(0), thick=th, color=clr(0)
endif

; --Galaxies to be observed by Keck
indx = where( x(*, 2) eq 1, count )
if count gt 0 then begin
    gl = x( indx, 0 )
    gb = x( indx, 1 )
    oplot, 180-gl, gb, psym=sym(1), symsize=ssize(1), thick=th, color=clr(1)
endif

; --Galaxies to be observed in Cycle 8
indx = where( x(*, 2) eq 2, count )
if count gt 0 then begin
    gl = x( indx, 0 )
    gb = x( indx, 1 )
    oplot, 180-gl, gb, psym=sym(2), symsize=ssize(2), thick=th, color=clr(2)
endif

; --Galaxies in HST Archive
indx = where( x(*, 2) eq 3, count )
if count gt 0 then begin
    gl = x( indx, 0 )
    gb = x( indx, 1 )
    oplot, 180-gl, gb, psym=sym(3), symsize=ssize(3), thick=th, color=clr(3)
endif

; --Galaxies in HST Archive
indx = where( x(*, 2) eq 4, count )
if count gt 0 then begin
    gl = x( indx, 0 )
    gb = x( indx, 1 )
    oplot, 180-gl, gb, psym=sym(4), symsize=ssize(4), thick=th, color=clr(4)
endif

; --Legend
if not keyword_set( nolabels ) then begin
    xl = 0.03
    yl = 0.11 - 0.025*indgen( 4 )
    plots, xl, yl( 0 ) + 0.005, /norm, psym=sym(1), color=clr(1)
    plots, xl, yl( 1 ) + 0.005, /norm, psym=sym(2), color=clr(2)
    plots, xl, yl( 2 ) + 0.005, /norm, psym=sym(3), color=clr(3)
    plots, xl, yl( 3 ) + 0.005, /norm, psym=sym(4), color=clr(4)
    xyouts, xl + 0.015, yl( 0 ), /norm, 'observed by Keck'
    xyouts, xl + 0.015, yl( 1 ), /norm, 'observed by HST in Cycle 8'
    xyouts, xl + 0.015, yl( 2 ), /norm, 'Archive object'
    xyouts, xl + 0.015, yl( 3 ), /norm, 'observed by Seitzer in Cycle 8'
endif

end

