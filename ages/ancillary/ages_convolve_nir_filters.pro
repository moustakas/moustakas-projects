pro ages_convolve_nir_filters, write=write
; jm06feb22uofa - convolve the AGES near-IR filter response functions
;                 (NDWFS/ONIS K and FLAMEX JKs) with the KPNO
;                 atmospheric transmission function

    filterpath = '/home/ioannis/ages/filters/'
    datapath = filterpath+'ORIGINAL/'

    im_window, 2, /square
    
; use the atmospheric transmission function provided by M. Brown

    readcol, '/home/ioannis/bin/observing/atmo_hybrid.dat', atmwave, atmresp, $
      format='F,X,F', /silent

; ONIS/DWFS K-band filter

    filtername = 'dwfs_onis_K.par'
    readcol, datapath+'dwfs_onis_K.dat', wave, resp, format='F,F', /silent
    resp = resp/100.0

    new_atmresp = interpol(atmresp,atmwave,wave)
    new_resp = resp*new_atmresp
    
    if keyword_set(write) then begin

       openw, lun, filterpath+filtername, /get_lun
       printf, lun, '# NDWFS/ONIS K filter curve, convolved with the atmospheric transmission function.'
       printf, lun, '# See AGES_CONVOLVE_NIR_FILTERS for details.'
       printf, lun, ''
       printf, lun, 'typedef struct {'
       printf, lun, '  double lambda;'
       printf, lun, '  double pass;'
       printf, lun, '} KFILTER;'
       printf, lun, ''
       for ii = 0L, n_elements(wave)-1L do printf, lun, 'KFILTER', wave[ii], $
         new_resp[ii], format='(A10,F12.2,F12.5)'
       free_lun, lun

    endif

    plot, atmwave, atmresp, xsty=3, ysty=3, yrange=[0,1], xtitle='Wavelength', $
      ytitle='Response', charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0, $
      title='ONIS K'
    djs_oplot, !x.crange, [0,0], color='blue'
    djs_oplot, wave, resp, color='red'
    djs_oplot, wave, new_resp, color='green'
    cc = get_kbrd(1)

; FLAMINGOS/FLAMEX J-band filter
    
    filtername = 'flamex_flamingos_J.par'
    readcol, datapath+'FLAMINGOS.BARR.J.MAN240.ColdWitness.txt', wave, resp, format='F,F', /silent
    resp = resp/100.0 & wave = wave*10.0

    new_atmresp = interpol(atmresp,atmwave,wave)
    new_resp = resp*new_atmresp
    
    if keyword_set(write) then begin

       openw, lun, filterpath+filtername, /get_lun
       printf, lun, '# FLAMINGOS/FLAMEX J filter curve, convolved with the atmospheric transmission function.'
       printf, lun, '# See AGES_CONVOLVE_NIR_FILTERS for details.'
       printf, lun, ''
       printf, lun, 'typedef struct {'
       printf, lun, '  double lambda;'
       printf, lun, '  double pass;'
       printf, lun, '} JFILTER;'
       printf, lun, ''
       for ii = 0L, n_elements(wave)-1L do printf, lun, 'JFILTER', wave[ii], $
         new_resp[ii], format='(A10,F12.2,F12.5)'
       free_lun, lun

    endif

    plot, atmwave, atmresp, xsty=3, ysty=3, yrange=[0,1], xtitle='Wavelength', $
      ytitle='Response', charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0, $
      title='FLAMINGOS/FLAMEX J'
    djs_oplot, !x.crange, [0,0], color='blue'
    djs_oplot, wave, resp, color='red'
    djs_oplot, wave, new_resp, color='green'
    cc = get_kbrd(1)

; FLAMINGOS/FLAMEX Ks-band filter
    
    filtername = 'flamex_flamingos_Ks.par'
    readcol, datapath+'FLAMINGOS.BARR.Ks.MAN306A.ColdWitness.txt', wave, resp, format='F,F', /silent
    resp = resp/100.0 & wave = wave*10.0

    new_atmresp = interpol(atmresp,atmwave,wave)
    new_resp = resp*new_atmresp
    
    if keyword_set(write) then begin

       openw, lun, filterpath+filtername, /get_lun
       printf, lun, '# FLAMINGOS/FLAMEX Ks filter curve, convolved with the atmospheric transmission function.'
       printf, lun, '# See AGES_CONVOLVE_NIR_FILTERS for details.'
       printf, lun, ''
       printf, lun, 'typedef struct {'
       printf, lun, '  double lambda;'
       printf, lun, '  double pass;'
       printf, lun, '} KsFILTER;'
       printf, lun, ''
       for ii = 0L, n_elements(wave)-1L do printf, lun, 'KsFILTER', wave[ii], $
         new_resp[ii], format='(A10,F12.2,F12.5)'
       free_lun, lun

    endif

    plot, atmwave, atmresp, xsty=3, ysty=3, yrange=[0,1], xtitle='Wavelength', $
      ytitle='Response', charsize=1.5, charthick=2.0, xthick=2.0, ythick=2.0, $
      title='FLAMINGOS/FLAMEX Ks'
    djs_oplot, !x.crange, [0,0], color='blue'
    djs_oplot, wave, resp, color='red'
    djs_oplot, wave, new_resp, color='green'
    cc = get_kbrd(1)

stop    

    info = im_filterspecs(filterlist=['dwfs_'+['Bw','I','R','onis_K'],'flamex_flamingos_'+['J','Ks']]+'.par',/verbose)


return
end
    
