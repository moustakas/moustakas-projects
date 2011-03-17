function deep2_version, kcorr=kcorr, ispec=ispec
; jm08mar31nyu - created
;
; kcorr
;   v1.0 - earlier effort
;   v2.0 - mainly updated k-corrections
;   v2.1 - matches the sample corresponding to the v2.0 spectra 
; ispec
;   v1.0 - earlier effort
;   v2.0 - remeasured spectra with updated ISPECLINEFIT() code [jm08sep04nyu]

    version = -1.0
    
    if keyword_set(kcorr) then version = 'dr3_v2.1'
    if keyword_set(ispec) then version = 'dr3_v2.0'

return, version
end
