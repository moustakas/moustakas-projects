function deep2_grid_axis, params, indx, midbin=midbin
; jm10sep13ucsd - simple support routine for the completeness
; parameter structure
    nbins = params.nbins[indx]
    binsz = params.binsz[indx]
    binsz_midbin = params.binsz_midbin[indx]
    hmin = params.hmin[indx]
    if keyword_set(midbin) then $
      axis = findgen(nbins)*binsz_midbin+hmin+binsz_midbin/2.0 else $
        axis = findgen(nbins)*binsz+hmin
return, axis
end

