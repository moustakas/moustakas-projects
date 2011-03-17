
pro map_label, lon_left=lon_left, lonlab=lonlab, latlab=latlab, $
               geographic=geo, charsize=charsize

;+
; MAP_LABEL
;   Draws labels on lon-lat grids
;
; KEYWORDS
;   lon_left -- leftmost longitude line (not labelled)
;   lonlab -- the latitutde of the longitude line labels
;   latlab -- the longitude of the latitude line labesl
;   geo -- set this keyword for longitude increasing left to right
;   char_size -- change relative size of labels (default=1)
;-

if not keyword_set( lon_left ) then lon_left = 360
if not keyword_set( lonlab ) then lonlab = 1
if not keyword_set( latlab ) then latlab = 180
if not keyword_set( charsize ) then charsize = 1

; --Label longitude lines
lsp = 90
for lmap = -180 + lsp, 180, lsp do begin
    lon = 180 - lmap - lon_left
    if keyword_set( geo ) then lon = -lon
    if lon lt 0 then lon = lon + 360
    xyouts, lmap, lonlab, string( format='(i0)', lon ), align=0.5, $
      size=charsize
endfor

; --Label latitude lines
lmap = 180 - latlab - lon_left
if keyword_set( geo ) then lmap = -lmap
xyouts, lmap, 31, '30', align=1, size=charsize
xyouts, lmap, -29, '-30', align=1, size=charsize

; --Mark NP
xyouts, 0, 90, 'N', align=0.5, size=0.7*charsize

end
