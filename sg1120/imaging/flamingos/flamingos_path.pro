function flamingos_path, original=original, mar06=mar06, $
  mosaics=mosaics, catalogs=catalogs
; jm07jan17nyu
; jm08jun14nyu - data removed from external drive

    datapath = getenv('RESEARCHPATH')+'/projects/sg1120/flamingos/'

    if keyword_set(dec03) then return, datapath+'03dec/'
    if keyword_set(feb06) then return, datapath+'06feb/'
    if keyword_set(mosaics) then return, datapath+'mosaics/'
    if keyword_set(catalogs) then return, datapath+'catalogs/'

    if keyword_set(original) then begin
       if file_test('/global/data/scr/ioannis/sg1120/flamingos/',/directory) then $
         datapath = '/global/data/scr/ioannis/sg1120/flamingos/' else begin
          message, 'Raw FLAMINGOS data not available on this machine.'
       endelse 
       return, datapath
    endif

return, datapath    
end
