; jm00june11ucb

; make pretty hst/wfpc2 images with cosmic rays rejected and bad
; pixels filled.  creates a montage of all four chips.

pro trgb_hst_makeim, objname, niceim

	spawn, ['mkdir -m 0755 COLOR']

        files = findfile(objname+'*.fits')
        fdata = files[where(strpos(files,'_dq.fits') eq -1L)]
        fcount = n_elements(fdata)

        mimage = {name: ' ', $
                  image: fltarr(1600,1600,fcount), $
                  header: strarr(500,fcount)}

        for i = 0L, fcount-1L do begin
            print & print, 'Reading '+fdata[i]+' . . . '
            wfpc2_read, fdata[i], imdata, dhead, /batwing
            mimage.image[*,*,i] = imdata
            mimage.header[0:n_elements(dhead)-1L] = dhead
        endfor
             
        print & print, 'Rejecting cosmic rays . . . '

	rdnoise = 1.00
        gain = 7.00

	cr_reject, mimage.image, rdnoise, 0, gain, 0.1, niceim, $
          combined_noise, combined_npix, $ ; exptime=exp, $
          nsig = [8,6,4], input_mask=mask, dilation=1, dfactor=0.5, $
          /noskyadjust, /noclearmask, mask_cube=cr_mask ;, /verbose

;        print & print, 'Interpolating over bad pixels . . . '
;        fill_missing, niceim, float(0), 2 ; interpolate over bad pixels            

        loadct, 3
        display, niceim < (median(niceim)+stdev(niceim)), min=0, /aspect
        writefits, 'COLOR/'+objname+'_montage.fits', niceim, dhead ; write out the new image

return
end
