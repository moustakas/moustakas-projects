function casid, run, rerun, camcol, field, obj_id, $
    first_field=first_field, sky_version=sky_version, verbose=verbose

    on_error, 2
    if n_params() lt 5 then begin
        print,'-Syntax: cid = casid(run, rerun, camcol, field, id, sky_version=, first_field=, /verbose)'
        print,'  defaults: sky_version=0, first_field=0'
        message,'Halting'
    endif

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; 
    ;; Takes run, rerun, camcol, field, id [and sky_version and 
    ;; first_field] to create a 64-bit integer that should match the 
    ;; unique id given in the SDSS SQL data base.
    ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


    id = lon64arr(n_elements(run))

    if n_elements(sky_version) eq 0 then begin
        if keyword_set(verbose) then begin
            print,'Assuming that sky_version is RUNS (15)'
        endif
        sky_version = intarr(n_elements(run)) + 15
    endif

    if sky_version(0) eq -1 then sky_version = intarr(n_elements(run)) + 15

    if n_elements(first_field) eq 0 then begin
        if keyword_set(verbose) then begin
            print,'Assuming that first_field is 0'
        endif
        first_field = intarr(n_elements(run))
    endif

    w = where(sky_version lt 0 or sky_version gt 15,total)

    if total gt 0 then message,'sky_version out of range (0 <= sky_version <= 15).'
    id = long64(sky_version)

    w = where(rerun LT 0 OR rerun GT 2047,total)

    if total gt 0 then message,'rerun out of range (0 <= rerun <= 2047).'
    id = ishft(id,11)
    id = id + long64(rerun)

    w = where(run lt 0 or run gt 65535,total)

    if total gt 0 then message,'run out of range (0 <= run <= 65535).'
    id = ishft(id,16)
    id = id + long64(run)

    w = where(camcol lt 0 or camcol gt 6,total)

    if total gt 0 then message,'camcol out of range (0 <= camcol <= 6).'
    id = ishft(id,3)
    id = id + long64(camcol)

    w = where(first_field lt 0 or first_field gt 1,total)

    if total gt 0 then message,'first_field out of range (first_field = 0,1).'
    id = ishft(id,1)
    id = id + long64(first_field)

    w = where(field lt 0 or field gt 4095,total)

    if total gt 0 then message,'field out of range (0 <= field <= 4095).'
    id = ishft(id,12)
    id = id + long64(field)

    w = where(obj_id lt 0 or obj_id gt 65535,total)

    if total gt 0 then message,'obj_id out of range (0 < obj_id < 65535).'
    id = ishft(id,16)
    id = id + long64(obj_id)

    return, id

end


