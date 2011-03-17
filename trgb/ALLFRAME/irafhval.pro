; jm00mar15ucb
function irafhval, parfile, pname
; function to pull out information from the setimpars output from IRAF
; passed as a string array.  returns the value corresponding to pname
; (parameter name)  

        good = where(strpos(parfile,pname) ne -1)
        gstring = parfile[good]
        if n_elements(gstring) eq 1 then $
          value = float(strmid(gstring[0],strpos(gstring[0],'=')+1,10))

return, value
end

