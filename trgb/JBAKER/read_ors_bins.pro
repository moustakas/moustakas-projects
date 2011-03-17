
pro read_ors_bins, file, v, mu, wmu, phi, nv, nmu, nphi, help=help

;+
; READ_ORS_BINS
;   Reads ORS bin coordinates from a file.  The first line contains the 
;   numbers of v, mu, phi bins.  This is followed by the v bins, one per
;   line.  Each mu bin has a corresponding weight on the same line. 
;   Last are the phi bins.
;
; INPUTS
;   file -- Name of input file
;
; OUTPUTS
;   v, mu, phi -- coordinates
;   wmu -- mu Gauss-Legendre weights
;
; OPTIONAL OUTPUTS
;   nv, nmu, nphi -- numbers of bins
;
; HISTORY
;   1996oct09  JEB
;-

if n_params() lt 5 or keyword_set( help ) then begin
    doc_library, 'read_ors_bins'
    retall
endif
print, '...Reading file ', file, '...'
read_data, file, x, nrows=1
nv = x( 0 )
nmu = x( 1 )
nphi = x( 2 )
read_data, file, v, row_start=2, nrows=nv
read_data, file, mu, wmu, row_start=nv+2, nrows=nmu
read_data, file, phi, row_start=2+nv+nmu

end

