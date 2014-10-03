pro build_hff_hiz1_photometry
; jm14sep19siena - build the photometric catalog for the HFF/hi-z1
; paper

    hffpath = getenv('CLASH_PROJECTS')+'/hff/hff_hiz1/'
    catpath = hffpath+'catalogs/'

    cat = {galaxy: '', z: 0.0, 
    
; Oesch+14 - luminous LBGs
    o14 = rsex(catpath+'14oesch_candels.dat')

    
    
; Zitrin+14 - triply-imaged z~10 candidate in A2744
    z14 = rsex(catpath+'a2744-z10.cat')

    
stop    
    
return
end
