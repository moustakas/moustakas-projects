function read_bpz_probs, file, redshift=zz
; jm13dec11siena - read the BPZ "probs" file
    
    txt = djs_readlines(file,nhead=1,head=head)
    ngal = n_elements(txt)

    pofz1 = strsplit(txt[0],' ',/extract)
    pofz1 = pofz1[1:n_elements(pofz1)-1]
    nzz = n_elements(pofz1)

    zinfo = float(strsplit(strmid(head,strpos(head,$
      'arange(')+7,21),',',/extract))
    zminmax = zinfo[0:1]
    dz = zinfo[2]
    zz = range(zminmax[0],zminmax[1]-dz ,nzz)

    data = replicate({id: 0L, pofz: fltarr(nzz)},ngal)
    for ii = 0, ngal-1 do begin
       pofz1 = strsplit(txt[ii],' ',/extract)
       data[ii].id = pofz1[0]
       data[ii].pofz = pofz1[1:n_elements(pofz1)-1]
    endfor
        
return, data    
end