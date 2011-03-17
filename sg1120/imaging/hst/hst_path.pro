function hst_path, mosaics=mosaics, catalogs=catalogs
; jm07jan17nyu
; jm09mar24nyu - rewritten

    datapath = getenv('RESEARCHPATH')+'/projects/sg1120/hst/'
    if keyword_set(mosaics) then return, datapath+'mosaics/'
    if keyword_set(catalogs) then return, datapath+'catalogs/gonzalez_deblend_0.015/'

return, datapath    
end
