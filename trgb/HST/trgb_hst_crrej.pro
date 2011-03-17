; jm00june9ucb

; generate a mean image with the cosmic rays rejected.  called by
; trgb_hst_redux.

; cr_mask is the mask of bad cosmic rays flagged by cr_reject (0/1 is bad/good)

pro trgb_hst_crrej, imcube, cr_mask, combined_image

	rdnoise = 1.00
        gain = 7.00
;       exp = sxpar(imcube.header[*,0],'UEXPODUR')	; seconds

        cube = imcube.image
        mask = imcube.mask

; reject cosmic rays

	cr_reject, cube, rdnoise, 0, gain, 0.1, combined_image, $
          combined_noise, combined_npix, $ ; exptime=exp, $
          nsig = [8,6,4], input_mask=mask, dilation=1, dfactor=0.5, $
          /noskyadjust, /noclearmask, mask_cube=cr_mask ;, /verbose

; make the median image

;       medarr, imcube.chip, medim

return
end
