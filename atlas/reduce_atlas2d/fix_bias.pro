pro fix_bias
; jm01oct24uofa
; repair a cosmic ray in the night 7 bias frame

    bias = readfits('a.2201.fits',h)

    xbad = 650L
    ybad = 29L

    box = 5L

    medval = median(bias[xbad-box:xbad+box,ybad-box:ybad+box])
    bias[xbad,ybad] = medval

    writefits, 'a.2201.fits', bias, h
    
return
end
