function mmt_almanac
; http://sculptor.as.arizona.edu/foltz/www/mmtsite.html

    lstmidstr = ['11:08:16'] ; ST at midnight on March 14/15
    lstmid = im_hms2dec_arr(lstmidstr)
    
return, lstmid
end

pro pisces_astrom
; jm01jul6uofa

; the data headers need: LST, LHA, RA, DEC, altitude, azimuth,
; airmass, epoch, parallactic angle, position angle (if possible),
; plate scale (optionally)

    root = 'qj35a'
    path = pisces_path()

    flist = findfile(path)
    fgood = where(strmatch(flist,root+'*.fit') eq 1B,fcount) ; good matches
    flist = flist[fgood]

    lstmid = (mmt_almanac())[0]
    rapoint = '14:26:00.8' & ra = im_hms2dec(rapoint)
    decpoint = '+35:31:31' & dec = im_hms2dec(decpoint)
    latitude = '+31:41:19.6' & lat = im_hms2dec(latitude)
    longitude = '110:53:04.4' & longit = im_hms2dec(longitude)
    
    for j = 0L, fcount-1L do begin

       head = headfits(path+flist[j])

       sxaddpar, head, 'RA', rapoint
       sxaddpar, head, 'DEC', decpoint
       sxaddpar, head, 'EPOCH', 'J2000'
       
       time = fxpar(head,'TIME')
       time = im_hms2dec(time)

       lst = lstmid - (float(12)-time) ; local apparent sidereal time
       aytriangle, lst, lat, ra, dec, zd, ha, alt, az, pangle

       niceprint, flist[j], im_dec2hms(time), ha, alt, az, pangle

    endfor

stop       


    readfast, 'source_list.txt', fixim, skip=1
    readfast, 'temp_list.txt', rotim, skip=1

    x_0 = fixim[0,*]-511.5
    y_0 = fixim[1,*]-511.5
    x_i = rotim[0,*]-511.5
    y_i = rotim[1,*]-511.5
    
    recpol, x_0, y_0, r_0, a_0, /degrees
    recpol, x_i, y_i, r_i, a_i, /degrees

    window, 0, xs=450, ys=450
    plot, a_0, r_0, ps=2
    oplot, a_i, r_i, ps=4

    jm = offset_from_pairs(r_0,a_0,r_i,a_i,dmax=5.0)    
    
;    xim = 360L
;    yim = ceil(max(r_0)>max(r_i))
;    im_0 = fltarr(xim,yim) & im_i = fltarr(xim,yim)
;
;    xy2image, x_0, y=y_0, im_0, /float, bin=1.0
;    xy2image, x_i, y=y_i, im_i, /float, bin=1.0
;
;    coalign, im_0, im_i, xshi=xshift, yshi=yshift

    

return
end
