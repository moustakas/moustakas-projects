function desi_bgs_templates_version, isedfit=isedfit
; isedfit version
    if keyword_set(isedfit) then begin
       version = 'v1.0'
    endif else begin
; template version
;   * original effort
       version = 'v1.0'
    endelse
return, version
end
