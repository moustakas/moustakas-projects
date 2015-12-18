pro rowfix_0415, wfits=wfits
; jm02mar2uofa
; pixels [0:171,93:96] were hit by a massive CR: repair individually.
; [these pixels should be added to the bad pixel mask]

    datapath = '/home/ioannis/kennicutt/data/98apr/'
    imname = 'a.0415.fits'

    image = float(readfits(imname,header,/silent))
    imsize = size(image,/dimension)
    
    mask = make_array(imsize[0],imsize[1],/byte)
    mask[0:171,93:96] = 1B
       
    cleanimage = djs_maskinterp(image,mask,iaxis=1L) 
       
;   sxaddhist, 'Pixels [0:171,93:96] repaired '+im_today()+'.', header

    if keyword_set(wfits) then begin
       print, 'Writing '+imname+'.'
       writefits, imname, cleanimage, header
    endif
       
return
end
