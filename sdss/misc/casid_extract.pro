;+
; NAME:
;   CASID_EXTRACT
;
; PURPOSE:
;   Extract run,rerun,camcol,field,id (and sky, first_field) 
;   from a casid (long64) from the cas (skyserver) or as created
;   by the casid.pro program.
;
; CATEGORY:
;   SDSS specific routine
;
; CALLING SEQUENCE:
;   casid_extract, super, run, rerun, camcol, field, id, 
;      sky_version=, first_field=
;
; INPUTS:
;   super: A super id
;
; OUTPUTS:
;   run,rerun,camcol,field,id
;
; OPTIONAL OUTPUTS:
;   sky_version, first_field
;
; EXAMPLE:
;   superid = casid(run,rerun,camcol,field,id)
;   casid_extract, superid, run, rerun, camcol, field, id
;
; MODIFICATION HISTORY:
;   Created: ?/?/? Ryan Scranton
;
;-


PRO casid_extract, casid, run, rerun, camcol, field, id, sky_version=sky_version, first_field=first_field

    if n_elements(casid) eq 0 then begin
        on_error, 2
        print,'-Syntax: extract_casid, casid, run, rerun, camcol, field, id, sky_version=, first_field='
        message,'Halting'
    endif


    run =    ishft(casid,-32) AND '0000FFFF'XL
    rerun =  ishft(casid,-48) AND '000007FF'XL
    camcol = ishft(casid,-29) AND '00000007'XL
    field =  ishft(casid,-16) AND '000000FFF'XL
    id =     casid AND '0000FFFF'XL

    if arg_present(sky_version) then sky_version = ishft(casid,-59) AND '0000000F'XL 
    if arg_present(first_field) then first_field = ishft(casid,-28) AND '000000001'XL

END
