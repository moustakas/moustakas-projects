pro trgb_hst_makeopt, pc=pc
;+
; NAME:
;	TRGB_HST_MAKEOPT
; PURPOSE:
;	Create all the .opt files for DAOPHOT for an image.
; INPUTS:
;
; OUTPUTS:
;
; KEYWORDS:
;	pc:  planetary camera keyword to distinguish between the
;	different hst chips
; PROCEDURE:
;
; PROGRAMS CALLED:
;
; MODIFICATION HISTORY:
;	Written, J. Moustakas, 2000 June 7
;-

	spawn, ['pwd'], datapath	; current directory
        datapath = datapath[0]

	daoopt = strarr(18)	; daophot.opt

        daoopt[0] = 'VAR = 2.00'
        daoopt[1] = 'READ = 1.00'
        if keyword_set(pc) then daoopt[2] = 'LOW = 40.00' else daoopt[2] = 'LOW = 100.00'
        if keyword_set(pc) then daoopt[3] = 'FWHM = 2.00' else daoopt[3] = 'FWHM = 1.50'
        daoopt[4] = 'WATCH = 0.00'
        if keyword_set(pc) then daoopt[5] = 'PSF = 14.00' else daoopt[5] = 'PSF = 10.00'
        daoopt[6] = 'GAIN = 7.00'
        daoopt[7] = 'HIGH = 8000.00'
        daoopt[8] = 'THRESHOLD = 5.00'
        daoopt[9] = 'FIT = 2.00'
        daoopt[10] = 'EX = 9.00'
        daoopt[11] = 'AN = -3.00'
        if keyword_set(pc) then daoopt[12] = 'LS = 0.10' else daoopt[12] = 'LS = 0.40'
        if keyword_set(pc) then daoopt[13] = 'HS = 1.00' else daoopt[13] = 'HS = 1.50'
        if keyword_set(pc) then daoopt[14] = 'LR = -0.80' else daoopt[14] = 'LR = -0.70'
        if keyword_set(pc) then daoopt[15] = 'HR = 0.80' else daoopt[15] = 'HR = 0.70'
        daoopt[16] = 'AN = -3.00'
        if keyword_set(pc) then daoopt[17] = 'USE = 1.00' else daoopt[17] = 'USE = 0.00'
        
        
        openw, lun1, datapath+'/daophot.opt', /get_lun
        for i = 0L, n_elements(daoopt)-1L do $
          printf, lun1, daoopt[i]

	photopt = strarr(14)	; photo.opt

        if keyword_set(pc) then photopt[0] = 'A1 = 3.261' else photopt[0] = 'A1 = 1.500'
        photopt[1] = 'A2 = 0.00'
        photopt[2] = 'A3 = 0.00'
        photopt[3] = 'A4 = 0.00'
        photopt[4] = 'A5 = 0.00'
        photopt[5] = 'A6 = 0.00'
        photopt[6] = 'A7 = 0.00'
        photopt[7] = 'A8 = 0.00'
        photopt[8] = 'A9 = 0.00'
        photopt[9] = 'AA = 0.00'
        photopt[10] = 'AB = 0.00'
        photopt[11] = 'AC = 0.00'
        if keyword_set(pc) then photopt[12] = 'IS = 20.00' else photopt[12] = 'IS = 14.00'
        if keyword_set(pc) then photopt[13] = 'OS = 40.00' else photopt[13] = 'OS = 35.00'
        
        openw, lun2, datapath+'/photo.opt', /get_lun
        for i = 0L, n_elements(photopt)-1L do $
          printf, lun2, photopt[i]

	allsopt = strarr(10)		; allstar.opt

        allsopt[0] = 'FI = 2.00'
        allsopt[1] = 'CE = 1.00'
        allsopt[2] = 'RE = 1.00'
        allsopt[3] = 'CR = 2.50'
        allsopt[4] = 'WA = 0.00'
        allsopt[5] = 'MAX = 50.00'
        allsopt[6] = 'PER = 0.00'
        if keyword_set(pc) then allsopt[7] = 'PRO = 2.00' else allsopt[7] = 'PRO = 1.01'
        if keyword_set(pc) then allsopt[8] = 'IS = 2.00' else allsopt[8] = 'IS = 0.70'
        if keyword_set(pc) then allsopt[9] = 'OS = 21.00' else allsopt[9] = 'OS = 15.00'

        openw, lun3, datapath+'/allstar.opt', /get_lun
        for i = 0L, n_elements(allsopt)-1L do $
          printf, lun3, allsopt[i]

	allfopt = strarr(10)		; allframe.opt

        allfopt[0] = 'CE = 8.00'
        allfopt[1] = 'CR = 6.00'
        allfopt[2] = 'GEO = 6.00'
        allfopt[3] = 'WA = 0.00'
        allfopt[4] = 'MIN = 5.00'
        allfopt[5] = 'MAX = 75.'
        if keyword_set(pc) then allfopt[6] = 'IS = 2.00' else allfopt[6] = 'IS = 1.00'
        if keyword_set(pc) then allfopt[7] = 'OS = 20.00' else allfopt[7] = 'OS = 14.00'
        allfopt[8] = 'PER = 0.00'
        allfopt[9] = 'PRO = 0.00'

        openw, lun4, datapath+'/allframe.opt', /get_lun
        for i = 0L, n_elements(allfopt)-1L do $
          printf, lun4, allfopt[i]

        free_lun, lun1, lun2, lun3, lun4

return
end


