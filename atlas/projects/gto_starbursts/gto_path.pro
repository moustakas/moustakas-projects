function gto_path, ancillary=ancillary, dss=dss
; jm06dec07nyu - written

    home = getenv('HOME')
    datapath = getenv('RESEARCHPATH')+'/projects/gto_starbursts/'
    if keyword_set(ancillary) then datapath = datapath+'ancillary/'
    if keyword_set(dss) then datapath = datapath+'dss/'
    
return, datapath
end
