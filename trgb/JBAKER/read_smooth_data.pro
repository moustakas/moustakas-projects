
pro read_smooth_data, file, data, n1, n2, n3, complex=complex, help=help

;+
; READ_SMOOTH_DATA
;   Reads f77 unformatted data file.  First line contains dimensions of
;   the array
;
; INPUTS
;   file -- input file name
;
; OPTIONAL OUTPUTS
;   data -- the array of data
;   n1, n2, n3 -- dimensions of the array
;
; KEYWORDS
;   complex -- If set, array is complex; otherwise, float
;-

if keyword_set( help ) or n_params() lt 1 then begin
    doc_library, 'read_smooth_data'
    retall
endif

print, '...Reading file ', file, '...'
openr, fp, file, /f77_unformatted, /get_lun
n = lonarr( 3 )
readu, fp, n
n1 = n( 0 )
n2 = n( 1 )
n3 = n( 2 )
if keyword_set( complex ) then data = complexarr( n1, n2, n3 ) $
else data = fltarr( n1, n2, n3 )
readu, fp, data
free_lun, fp

end

