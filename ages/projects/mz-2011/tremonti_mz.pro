function tremonti_mz, massaxis=massaxis, pivotmass=pivotmass, $
  coeff=coeff, kk04=kk04
; jm09mar28nyu - Tremonti+04 MZ relation; setting /KK04 adds +0.05
;   dex, putting the relation on the KK04 abundance scale
    
    if (n_elements(massaxis) eq 0L) then massaxis = im_array(8.8,11.3,0.01)
    if (n_elements(pivotmass) eq 0L) then pivotmass = 10.5

;   coeff = [9.10284,0.16154,-0.08026] ; T04 with a 10.5 pivot mass
;   oh = poly(massaxis-pivotmass,coeff)

    coeff = [-1.492,1.847,-0.08026] ; T04 published
    if keyword_set(kk04) then coeff[0] = coeff[0] + 0.05
    oh = poly(massaxis,coeff)

return, oh
end
