function ages_upweight_field, reject_pass, field=field
; jm09dec01ucsd - for several applications we have to reject certain
;   plates because of e.g., poor spectrophotometry or other issues;
;   each fiber is assigned to a field number, where fields 1-15
;   correspond to the "main" survey area; this routine computes the
;   weight that the observations corresponding to each field must be
;   upweighted in order to retain the statistical fidelity of the
;   sample
;
; for example, the fibers in pass=104 span the field numbers 0, 3, 4,
; and 8; so we call this routine:
;   weight = ages_upweight_field(104,field=field)
;   niceprint, field, weight
;       0       0.997951
;       1        1.00000
;       2        1.00000
;       3       0.983501
;       4       0.847826
;       5        1.00000
;       6        1.00000
;       7        1.00000
;       8       0.981300
;       9        1.00000
;      10        1.00000
;      11        1.00000
;      12        1.00000
;      13        1.00000
;      14        1.00000
;      15        1.00000
;
; in other words, fibers corresponding to field numbers 0, 3, 4, and 8
; must be *upweighted* (i.e., divided by WEIGHT) by the amounts listed
; above 
    phot = read_ages(/phot)
    main = where((phot.main_weight ge 1) and $
      (phot.spec_yesno ge 1) and (phot.z_yesno ge 1))
    phot = phot[main]

    allfield = phot.field
    allpass = phot.pass
    pass = allpass[uniq(allpass,sort(allpass))]
    npass = n_elements(pass)

    field = indgen(16) ; includes null field 0
    weight = field*0.0+1.0
    nfield = n_elements(field)

    for ii = 0, npass-1 do begin
       isrej = where(pass[ii] eq reject_pass,nisrej)
       if (nisrej ne 0) then begin
          for jj = 0, nfield-1 do begin
             nrej = total((allpass eq pass[ii]) and (allfield eq jj))
             npossible = total(allfield eq jj+1)
             weight[jj] = 1.0-nrej/(npossible+(npossible eq 0))*(npossible ne 0)
          endfor
       endif
    endfor

return, weight
end
