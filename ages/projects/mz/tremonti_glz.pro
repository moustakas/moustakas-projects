function tremonti_glz, magaxis=magaxis, pivotmag=pivotmag, kk04=kk04
; jm09mar28nyu - Tremonti+04 g-band LZ relation; setting /KK04 adds
;   +0.05 dex, putting the relation on the KK04 abundance scale
    
    if (n_elements(magaxis) eq 0L) then magaxis = im_array(-24.0,-17.0,0.05)
    if (n_elements(pivotmag) eq 0L) then pivotmag = -20.5

;   glzcoeff = [5.195,-0.186] ; T04 published
    glzcoeff = [9.008,-0.186] ; T04 with a -20.5 pivot magnitude 
    if keyword_set(kk04) then glzcoeff[0] = glzcoeff[0]+0.05
    oh = poly(magaxis-pivotmag,glzcoeff)

return, oh
end
