pro write_01kewley, out, out_nodust, _extra=extra
; jm07aug10nyu - classify the Kewley et al. galaxies

; SAMPLE AND DATA DON'T HAVE THE SAME NUMBER OF OBJECTS!!!    
    
    sample = im_read_fmr('table1.txt')
    ngalaxy = n_elements(sample)

    alldata = im_read_fmr('table3.txt')
    data = im_empty_structure(alldata,ncopies=ngalaxy)
    for ii = 0L, ngalaxy-1L do begin
       match = where(strtrim(sample[ii].name,2) eq strtrim(strmid(alldata.iras,0,10),2),nmatch)
       if (nmatch ne 0L) then data[ii] = alldata[match]
    endfor

STOP
    
    oratio = 2.984
    nratio = 3.054

; use ICLASSIFICATION to classify the galaxies in this sample; the
; ratios are *observed*, not reddening-corrected!

    line = {$
      galaxy:           '',        $
      alt_galaxy:       '',        $
      linename:         ['H_BETA','OIII_4959','OIII_5007','OI_6300',$
      'NII_6548','H_ALPHA','NII_6584','SII_6716','SII_6731'], $
      h_beta_wave:      4861.33, $
      oiii_4959_wave:   4958.91, $
      oiii_5007_wave:   5006.84, $
      oi_6300_wave:     6300.304,$
      nii_6548_wave:    6548.04, $
      h_alpha_wave:     6562.80, $
      nii_6584_wave:    6583.46, $
      sii_6716_wave:    6716.14, $
      sii_6731_wave:    6730.81, $
      h_beta:        [0.0,-2.0], $
      oiii_4959:     [0.0,-2.0], $
      oiii_5007:     [0.0,-2.0], $
      oi_6300:       [0.0,-2.0], $
      nii_6548:      [0.0,-2.0], $
      h_alpha:       [0.0,-2.0], $
      nii_6584:      [0.0,-2.0], $
      sii_6716:      [0.0,-2.0], $
      sii_6731:      [0.0,-2.0]}
    line = replicate(line,ngalaxy)

    line.alt_galaxy = 'IRAS'+strcompress(sample.name,/remove)
    line.galaxy = strcompress(sample.cname,/remove)
    empty = where(strtrim(line.galaxy,2) eq '')
    line[empty].galaxy = line[empty].alt_galaxy
    
    line.h_beta = [1.0,0.05] ; assume 5%
    
    g = where(data.ha_hb gt -900.0,ng)
    if (ng ne 0L) then begin
       line[g].h_alpha[0] = 10.0^data[g].ha_hb
       line[g].h_alpha[1] = line[g].h_alpha[0]*0.05
    endif

    g = where(data.oiii_hb gt -900.0,ng)
    if (ng ne 0L) then begin
       line[g].oiii_5007[0] = 10.0^data[g].oiii_hb
       line[g].oiii_5007[1] = line[g].oiii_5007[0]*0.05
       line[g].oiii_4959 = line[g].oiii_5007/oratio
    endif

    g = where((data.nii_ha gt -900.0) and (data.ha_hb gt -900.0),ng)
    if (ng ne 0L) then begin
       line[g].nii_6584[0] = 10.0^(data[g].nii_ha+data[g].ha_hb)
       line[g].nii_6584[1] = line[g].nii_6584[0]*0.05
       line[g].nii_6548 = line[g].nii_6584/nratio
    endif

    g = where((data.oi_ha gt -900.0) and (data.ha_hb gt -900.0),ng)
    if (ng ne 0L) then begin
       line[g].oi_6300[0] = 10.0^(data[g].oi_ha+data[g].ha_hb)
       line[g].oi_6300[1] = line[g].oi_6300[0]*0.05
    endif

    class = iclassification(line,ratios=r,snrcut=0.0)
    struct_print, class
    
;   line_nodust = iunred_linedust(line,snrcut=1.0)

stop    
return
end    
