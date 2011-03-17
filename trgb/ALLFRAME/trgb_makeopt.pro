function daophotopt, fwhm, psf, highbad, threshold, var

	daoopt = strarr(12)

        daoopt[0] = 'VAR = '+var
        daoopt[1] = 'READ = 6.30'
        daoopt[2] = 'LOW = 6.'
        daoopt[3] = 'FWHM = '+fwhm
        daoopt[4] = 'WATCH = 0'
        daoopt[5] = 'PSF = '+psf
        daoopt[6] = 'GAIN = 1.97'
        daoopt[7] = 'HIGH = '+highbad
        daoopt[8] = 'THRESHOLD = '+threshold
        daoopt[9] = 'FIT = '+fwhm
        daoopt[10] = 'EX = 3'
        daoopt[11] = 'AN = 1'
        
        daoopt = transpose(daoopt)
        return, daoopt

end

function photoopt, fwhm, is, os

	photopt = strarr(3)

        photopt[0] = 'A1 = '+fwhm
        photopt[1] = 'IS = '+is
        photopt[2] = 'OS = '+os
        
	return, photopt
end

function allstaropt, fwhm, is, os

	allsopt = strarr(10)

        allsopt[0] = 'FI = '+fwhm
        allsopt[1] = 'CE = 4.0'
        allsopt[2] = 'RE = 1'
        allsopt[3] = 'CR = 2.5'
        allsopt[4] = 'WA = 0'
        allsopt[5] = 'MAX = 50.'
        allsopt[6] = 'PER = 0.35'
        allsopt[7] = 'PRO = 3.'
        allsopt[8] = 'IS = '+is
        allsopt[9] = 'OS = '+os

	return, allsopt
end

function allframeopt, is, os

	allfopt = strarr(10)

        allfopt[0] = 'CE = 4.'
        allfopt[1] = 'CR = 2.5'
        allfopt[2] = 'GEO = 6'
        allfopt[3] = 'WA = 0'
        allfopt[4] = 'MIN = 5'
        allfopt[5] = 'MAX = 100.'
        allfopt[6] = 'IS = '+is
        allfopt[7] = 'OS = '+os
        allfopt[8] = 'PER = 0.35'
        allfopt[9] = 'PRO = 3.'

	return, allfopt
end

pro trgb_makeopt, parfile, psfhead
;+
; NAME:
;	TRGB_MAKEOPT
;
; PURPOSE:
;	Create all the .opt files for DAOPHOT for an image.
;
; INPUTS:
;	parfile : a scalar of the full pathname to the IRAF parameter
;	          file corresponding to the current image
;
; OUTPUTS:
;	Creates a daophot.opt, photo.opt, allstar.opt, and
;	allframe.opt in the current directory.
;
; PROCEDURE:
;	Reads in the specified parameter file and extracts the
;	necessary parameters to create the .opt files.
;
; PROCEDURES USED:
;	RDTXT(), IRAFHVAL(), 
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 March 15, UCB
;-

	spawn, ['pwd'], datapath	; current directory
        datapath = datapath[0]

	npar = n_params()
        if npar eq 2 then begin
            var = sxpar(psfhead,'VARORDER')
            var = strcompress(string(var,format='(F2.0)'),/remove_all)
        endif else var = strn(2)

        params = rdtxt(parfile)
        
        fwhm = irafhval(params,'datapars.fwhmpsf')
        fwhm = strcompress(string(fwhm,format='(F4.2)'),/remove_all)

        psf = irafhval(params,'daopars.psfrad')
        psf = strcompress(string(psf,format='(F5.2)'),/remove_all)

        highbad = irafhval(params,'datapars.datamax')
        highbad = strcompress(string(highbad,format='(F7.1)'),/remove_all)

        threshold = irafhval(params,'findpars.threshold')
        threshold = strcompress(string(threshold,format='(F2.0)'),/remove_all)

        annulus = irafhval(params,'fitskypars.annulus')
        dannulus = irafhval(params,'fitskypars.dannulus')
        osky = annulus+dannulus

        annulus = strcompress(string(annulus,format='(F5.2)'),/remove_all)
        dannulus = strcompress(string(dannulus,format='(F5.2)'),/remove_all)
        osky = strcompress(string(osky,format='(F5.2)'),/remove_all)

;       ap2 = 20.0 ; (pixels; second aperture for photo.opt)
;       ap2 = strcompress(string(ap2,format='(F4.1)'),/remove_all)
        
; create the daophot.opt file

        daoopt = daophotopt(fwhm, psf, highbad, threshold, var)
        openw, lun1, datapath+'/daophot.opt', /get_lun
        for i = 0L, n_elements(daoopt)-1L do $
          printf, lun1, daoopt[i]

; create the photo.opt file

        photopt = photoopt(fwhm, annulus, osky);, ap2)
        openw, lun2, datapath+'/photo.opt', /get_lun
        for i = 0L, n_elements(photopt)-1L do $
          printf, lun2, photopt[i]

; create the allstar.opt file

        allsopt = allstaropt(fwhm,annulus,osky)
        openw, lun3, datapath+'/allstar.opt', /get_lun
        for i = 0L, n_elements(allsopt)-1L do $
          printf, lun3, allsopt[i]

; create the allframe.opt file

        allfopt = allframeopt(annulus,osky)
        openw, lun4, datapath+'/allframe.opt', /get_lun
        for i = 0L, n_elements(allfopt)-1L do $
          printf, lun4, allfopt[i]

        free_lun, lun1, lun2, lun3, lun4

return
end


