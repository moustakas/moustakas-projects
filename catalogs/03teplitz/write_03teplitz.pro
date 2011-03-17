pro write_03teplitz, outdata
; jm04nov30uofa

    root = '03teplitz'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    table3 = im_read_fmr(path+root+'_table3.dat')
    table4 = im_read_fmr(path+root+'_table4.dat')

; only keep objects with a redshift flag greater than zero    

    keep = where(table3.f_z gt 0.0,ngalaxy)
    table3 = table3[keep]

; initialize the output data structure    

    data = init_cat_linefit(ngalaxy=ngalaxy)
    
; ---------------------------------------------------------------------------    
; store relevant galaxy properties
; ---------------------------------------------------------------------------    

    data.galaxy = table3.name
    data.z_obj = table3.z
    data.distance = dluminosity(data.z_obj,/Mpc)

; loop on each object

    for i = 0L, ngalaxy-1L do begin

       match = where(strtrim(table3[i].name,2) eq strtrim(table4.name,2),nmatch)
       if (nmatch eq 0L) then message, 'Houston, we have a problem.'

       for j = 0L, nmatch-1L do begin

          info = table4[match[j]]
          line = repstr(repstr(info.line,'[',''),']','') ; this is required to use STRMATCH
          
          if strmatch(line,'*OII*') then begin
             if (info.snr gt 0.0) then data[i].oii_3727     = info.flux*[1.0,1.0/info.snr]*1D-16
             if (info.csnr gt 0.0) then data[i].oii_3727_ew = info.eqwid*[1.0,1.0/info.csnr] / (1+data[i].z_obj) ; rest
          endif

          if strmatch(line,'*H{gamma}*') then begin
             if (info.snr gt 0.0) then data[i].h_gamma     = info.flux*[1.0,1.0/info.snr]*1D-16
             if (info.csnr gt 0.0) then data[i].h_gamma_ew = info.eqwid*[1.0,1.0/info.csnr] / (1+data[i].z_obj) ; rest
          endif

          if strmatch(line,'*H{beta}*') then begin
             if (info.snr gt 0.0) then data[i].h_beta     = info.flux*[1.0,1.0/info.snr]*1D-16
             if (info.csnr gt 0.0) then data[i].h_beta_ew = info.eqwid*[1.0,1.0/info.csnr] / (1+data[i].z_obj) ; rest
          endif

          if strmatch(line,'*OIII5007*') then begin
             if (info.snr gt 0.0) then data[i].oiii_5007     = info.flux*[1.0,1.0/info.snr]*1D-16
             if (info.csnr gt 0.0) then data[i].oiii_5007_ew = info.eqwid*[1.0,1.0/info.csnr] / (1+data[i].z_obj) ; rest
          endif

          if strmatch(line,'*H{alpha}*') then begin
             if (info.snr gt 0.0) then data[i].h_alpha     = info.flux*[1.0,1.0/info.snr]*1D-16
             if (info.csnr gt 0.0) then data[i].h_alpha_ew = info.eqwid*[1.0,1.0/info.csnr] / (1+data[i].z_obj) ; rest
          endif

       endfor

    endfor

    good = where(data.oiii_5007[1] gt 0.0,ngood)
    if (ngood ne 0L) then begin
       data[good].oiii_4959    = data[good].oiii_5007 / 3.0
       data[good].oiii_4959_ew = data[good].oiii_5007_ew / 3.0
    endif
    
; ---------------------------------------------------------------------------    
; continuum fluxes
; ---------------------------------------------------------------------------    

; ---------------------------------------------------------------------------    
; de-redden and compute abundances
; ---------------------------------------------------------------------------    

    adata = im_abundance(data,snrcut=1.0)
    outdata = struct_addtags(data,adata)

    splog, 'Writing '+path+root+'.fits.'
    mwrfits, outdata, path+root+'.fits', /create
    spawn, ['gzip -f ']+path+root+'.fits', /sh

; reddening-corrected fluxes    
    
    idata = iunred_linedust(data,snrcut=1.0,/silent)
    adata = im_abundance(idata,snrcut=1.0)
    outdata = struct_addtags(idata,adata)
    
    splog, 'Writing '+path+root+'_nodust.fits.'
    mwrfits, outdata, path+root+'_nodust.fits', /create
    spawn, ['gzip -f ']+path+root+'_nodust.fits', /sh
    
return
end    
