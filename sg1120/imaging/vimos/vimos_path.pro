function vimos_path, original=original, dec03=dec03, feb06=feb06, $
  mosaics=mosaics, catalogs=catalogs, qc=qc, qaplots=qaplots
; jm07jan15nyu
; jm08jun09nyu - data removed from external drive

    datapath = getenv('RESEARCHPATH')+'/projects/sg1120/vimos/'

    if keyword_set(dec03) then return, datapath+'03dec/'
    if keyword_set(feb06) then return, datapath+'06feb/'
    if keyword_set(mosaics) then return, datapath+'mosaics/'
    if keyword_set(catalogs) then return, datapath+'catalogs/'
    if keyword_set(qc) then return, datapath+'qc/'
    if keyword_set(qaplots) then return, datapath+'qaplots/'

    if keyword_set(original) then begin
       if file_test('/mount/moon1/ioannis/sg1120/vimos/',/directory) then $
         datapath = '/mount/moon1/ioannis/sg1120/vimos/' else begin
          message, 'Raw VIMOS data not available on this machine.'
       endelse 
       return, datapath
    endif

return, datapath    
end
