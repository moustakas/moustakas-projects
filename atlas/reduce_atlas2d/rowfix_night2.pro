pro rowfix_night2, wfits=wfits
; jm02jan7uofa
; repair rows 49 through 51 in the night 2 (1999 April 21) data
; [these pixels should be added to the bad pixel mask]

    datapath = '/home/ioannis/kennicutt/data/99apr/'
    flist = findfile(datapath+'a.17*.fits',count=fcount)

    for k = 0L, fcount-1L do begin

       image = float(readfits(flist[k],header,/silent))
       imsize = size(image,/dimension)

       mask = make_array(imsize[0],imsize[1],/byte)
       mask[*,49:51] = 1B
       
       cleanimage = djs_maskinterp(image,mask,iaxis=1L) 

       sxdelpar, header, 'HISTORY'
       sxaddhist, 'Rows 49 through 51 repaired '+im_today()+'.', header

       if keyword_set(wfits) then begin
          print, 'Writing '+flist[k]+'.'
          writefits, flist[k], cleanimage, header
       endif
       
    endfor

return
end
