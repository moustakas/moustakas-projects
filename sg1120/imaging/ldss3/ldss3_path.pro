function ldss3_path, original=original, feb06=feb06, sex=sex, $
  mosaics=mosaics, catalogs=catalogs
; jm06nov06nyu; ORIGINAL is the path name to the data straight off the
; DVD's; to access the same data, after running UNPACK_LDSS3, run
; ldss3_path(/feb06)+'raw/'
; jm08jun14nyu - data removed from external drive

    datapath = getenv('RESEARCHPATH')+'/projects/sg1120/ldss3/'

    if keyword_set(feb06) then return, datapath+'06feb/'
    if keyword_set(mosaics) then return, datapath+'mosaics/'
    if keyword_set(catalogs) then return, datapath+'catalogs/'

    if keyword_set(original) then begin
       if file_test('/mount/moon1/ioannis/sg1120/ldss3/ktran.27feb2006/',/directory) then $
         datapath = '/mount/moon1/ioannis/sg1120/ldss3/ktran.27feb2006/' else begin
          message, 'Raw LDSS3 data not available on this machine'
       endelse 
       return, datapath
    endif

return, datapath    
end
