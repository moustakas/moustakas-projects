pro write_04maier, outdata, debug=debug
; jm05jan01uofa

    snrcut = 3.0
    Bmsun = k_solar_magnitudes(filterlist='bessell_B.par')
    Bvega2ab = (k_vega2ab(filterlist='bessell_B.par',/kurucz))[0]
;   Bmsun = 5.42
    hahb = 2.86
    
    root = '04maier'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    raw = rsex(path+root+'_original.dat')
    ngalaxy = n_elements(raw)

; initialize the output data structure    

    data = init_cat_linefit(ngalaxy=ngalaxy)
    
; ---------------------------------------------------------------------------    
; store relevant galaxy properties
; ---------------------------------------------------------------------------    

    data.id = 'M04_'+string(lindgen(ngalaxy)+1L,format='(I2.2)')
    data.galaxy = strtrim(raw.galaxy,2)
    data.z_obj = raw.z
    data.distance = dluminosity(data.z_obj,/Mpc)

    data.M_B = raw.M_B ; Vega
    data.M_B_err = 0.1 ; reference?

    lum = 10.0D^(0.4*Bmsun) * 10.0^(-0.4*data.M_B)   ; [L_sun]
    lum_err = alog(10.0) * 0.4D * lum * data.M_B_err ; [L_sun]
       
    data.B_lum_err = lum_err/lum/alog(10.0) ; log L_sun
    data.B_lum = alog10(lum)                ; log L_sun
    
; ---------------------------------------------------------------------------    
; fluxes; convert from 10^-20 W/m2 by multiplying by 10-17
; ---------------------------------------------------------------------------    

    data.oii_3727[0] = raw.oii_3727*1D-17
    data.oii_3727[1] = raw.oii_3727_err*1D-17
    
    data.h_beta[0] = raw.h_beta*1D-17
    data.h_beta[1] = raw.h_beta_err*1D-17
    
    data.oiii_5007[0] = raw.oiii_5007*1D-17
    data.oiii_5007[1] = raw.oiii_5007_err*1D-17
    
;   data.h_alpha = HaHb*data.h_beta ; assume no reddening! [E(B-V)=0]

;   data.h_beta[0] = 1D-17               ; arbitrarily set!!
;   data.h_beta[1] = data.h_beta[0]*0.01 ; assume 1% error

;   data.oii_3727[0] = data.h_beta[0] * raw.oiii_hb / raw.oiii_oii
;   data.oii_3727[1] = data.h_beta[0] * im_compute_error(raw.oiii_hb,raw.oiii_hb_err,$
;     raw.oiii_oii,raw.oiii_oii_err,/quotient)

;   data.oiii_5007[0] = data.h_beta[0] * raw.oiii_hb
;   data.oiii_5007[1] = data.h_beta[0] * raw.oiii_hb_err

    data.oiii_4959 = data.oiii_5007 / 3.0D ; assume this is true

; ---------------------------------------------------------------------------    
; compute the reddening, metallicity, and other miscellaneous
; properties, then write out 
; ---------------------------------------------------------------------------    

    idata = iunred_linedust(data,snrcut=snrcut,/silent,/nopropagate)
    mz_log12oh, data, idata, data, ohdust, ohnodust, abund, abundnodust, $
      r23branch=raw.branch, branchmethod=5L, snrcut=snrcut

    outdata = struct_addtags(data,ohdust)

    lums = compute_linelums(outdata,select_lines=['oii_3727','oiii_5007','h_beta','h_alpha'])
    outdata = struct_addtags(outdata,lums)
    sfrs = compute_sfrs(outdata)

    class = iclassification(outdata,/kauffmann,/silent,doplot=debug)
    outdata = struct_addtags(struct_addtags(outdata,sfrs),class)

    splog, 'Writing '+path+root+'.fits.'
    mwrfits, outdata, path+root+'.fits', /create
    spawn, ['gzip -f ']+path+root+'.fits', /sh

; ---------------------------------------------------------------------------    
; write out statistics
; ---------------------------------------------------------------------------    

    indx = where((outdata.m_b gt -900.0) and (outdata.zstrong_12oh_kk04 gt -900.0),nindx)
    
    print, '##################################################'
    splog, 'Maier et al. (2004): '+string(nindx,format='(I0.0)')+' galaxies.'
    stats = im_stats(outdata[indx].z_obj)
    splog, strtrim(string(stats.min,format='(F12.2)'),2)+'<z<'+strtrim(string(stats.max,format='(F12.2)'),2)+$
      ', <z>='+strtrim(string(stats.mean,format='(F12.2)'),2)
    stats = im_stats(outdata[indx].m_b)
    splog, strtrim(string(stats.min,format='(F12.2)'),2)+'<M_B<'+strtrim(string(stats.max,format='(F12.2)'),2)+$
      ', <M_B>='+strtrim(string(stats.mean,format='(F12.2)'),2)
    print, '##################################################'
    struct_print, struct_trimtags(outdata,select=['GALAXY','M_B','ZSTRONG_12OH_KK04*','R23BRANCH_KK04'])
    print, '##################################################'
       
; ---------------------------------------------------------------------------    
; make some plots
; ---------------------------------------------------------------------------    

    if keyword_set(debug) then begin
       
       plot, outdata.m_b, outdata.zstrong_12oh_kk04, xrange=[-23,-15], yrange=[8.0,9.2], $
         ps=4, xsty=3, ysty=3, sym=2
       cc = get_kbrd(1)

    endif
       
return
end    
