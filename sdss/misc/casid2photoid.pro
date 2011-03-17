
;+
; NAME:
;  casid2photoid
;
; PURPOSE:
;  Convert a set of CAS ids (skyserver ids) to a photoids.
;
; CATEGORY:
;  SDSS
;
; CALLING SEQUENCE:
;  pid = casid2photoid(casid)
;
; INPUTS:
;  casid: A scalar or array of CAS ids.
;
; OUTPUTS:
;  photoid(s).  See the program photoid.pro
;
; OPTIONAL OUTPUTS:
;  sky_version, first_field. Extra info kept in CAS ids not 
;    kept in photoids.
;
; MODIFICATION HISTORY:
;  Created: 2007-04-03. Erin Sheldon, NYU
;
;-

function casid2photoid, casid, sky_version=sky_version, first_field=first_field

    if n_elements(casid) eq 0 then begin
        print,'-Syntax: casid2photoid(casid, sky_version=, first_field=)'
        on_error, 2
        message,'Halting'
    endif

    casid_extract, casid, run, rerun, camcol, field, id, sky_version=sky_version, first_field=first_field
    return, photoid(run, rerun, camcol, field, id)

end
