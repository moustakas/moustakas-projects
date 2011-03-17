
;+
; NAME:
;  photoid2casid
;
; PURPOSE:
;  Convert a set of photoids to CAS ids (skyserver ids)
;
; CATEGORY:
;  SDSS
;
; CALLING SEQUENCE:
;  cid = photoid2casid(pid)
;
; INPUTS:
;  pid: A scalar or array of photoids as created by the photoid.pro
;   program.
;
; OPTIONAL OUTPUTS:
;  sky_version, first_field. Extra info kept in CAS ids not 
;    kept in photoids.  See the casid() program for more info.
;
; OUTPUTS:
;  casid(s).  
;
; MODIFICATION HISTORY:
;  Created: 2007-04-03. Erin Sheldon, NYU
;
;-

function casid2photoid, pid, sky_version=sky_version, first_field=first_field

    if n_elements(casid) eq 0 then begin
        print,'-Syntax: photoid2casid(photoids, sky_version=, first_field=)'
        on_error, 2
        message,'Halting'
    endif

    photoid_extract, pid, run, rerun, camcol, field, id
    return, casid(run, rerun, camcol, field, id, sky_version=sky_version, first_field=first_field)

end
