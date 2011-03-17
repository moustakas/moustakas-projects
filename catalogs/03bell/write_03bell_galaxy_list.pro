pro write_ned_galaxy_list
; jm04nov29uofa - write the NED galaxy list

    path = getenv('CATALOGS_DIR')+'/bell_radiofir/'
    file = 'tableA1.dat'
    rdata = read_fmr(path+file)

; generate a list of galaxy names that are NED-compatible, then write
; out the galaxy list 

    data = struct_addtags(data,replicate({galaxy_ned: ''},ngalaxy))
    data.galaxy = strupcase(strtrim(data.galaxy,2))
    data.galaxy_ned = data.galaxy

    rc2 = where((strmatch(data.galaxy,'A?*') eq 1B) and (strmatch(data.galaxy,'*ARP*',/fold) ne 1B),nrc2)
    if (nrc2 ne 0L) then data[rc2].galaxy_ned = '[RC2]'+data[rc2].galaxy_ned

    ant = where(strmatch(data.galaxy,'*4038/9*') eq 1B,nant) ; The Antennae = Arp 244
    if (nant ne 0L) then data[ant].galaxy_ned = 'ARP244'
    
    vv = where(strmatch(data.galaxy,'*4567/8*') eq 1B,nvv) ; VV219
    if (nvv ne 0L) then data[vv].galaxy_ned = 'VV219'
    
    openw, lun, path+'galaxy_list.txt', /get_lun
    niceprintf, lun, data.galaxy_ned
    free_lun, lun
    
return
end    
