pro atlas_specfit, test=test, integrated=integrated, nuclear=nuclear
; jm02jul14uofa - fit the integrated and nuclear galaxy spectra 
; jm05jul22uofa - updated
; jm07sep24nyu - updated to use ATLAS_SPECLINEFIT and to fit a broad
;                component to a small subset of the AGN in our sample 

    datapath = atlas_path(/atlas1d)
    specfitpath = atlas_path(/specfit)
    atlas = atlas_read_info()

    bc03version = atlas_version(/templates)
    eigendir = atlas_path(/specfit)
    eigenfile = 'BC03_Z02_salpeter_atlas_'+bc03version+'_templates.fits'

    specres = 7.9 ; 8.0

; fit the integrated spectra

    if keyword_set(integrated) then begin
       
       if keyword_set(test) then suffix = 'integrated_test' else $
         suffix = 'integrated_atlas'

;      int = where(atlas.drift and (atlas.drift_agnflag eq 1L))
       int = where(atlas.drift and (atlas.drift_agnflag eq 0L))
       speclist = strtrim(atlas[int].drift_file,2)

       specdata = atlas_ispeclinefit(speclist,specres=specres,snrcut=0.0,dustmodel=0,zsnrcross=5.0,$
         datapath=datapath,linepath=specfitpath,eigendir=eigendir,eigenfile=eigenfile,suffix=suffix,$
         Zmulti=Zmulti,/charlot,/zcrosscor,/postscript,/write,vmaxshift=500.0,/nologfile,starvdisp=100.0,$
         /integrated,/positive_emission)

       atlas_parse_specfit
       
    endif
       
; fit the nuclear spectra; VMAXSHIFT=600 km/s is needed for NGC1003;
; having trouble fitting NGC1275 (need better masking of the "other"
; emission lines); NGC1068 is OK, not great

    if keyword_set(nuclear) then begin

       if keyword_set(test) then suffix = 'nuclear_test' else $
         suffix = 'nuclear_atlas'

;      nuc = where(atlas.nuclear and (atlas.nuclear_agnflag eq 1L))
       nuc = where(atlas.nuclear and (atlas.nuclear_agnflag eq 0L))
       speclist = strtrim(atlas[nuc].nuclear_file,2)

       specdata = atlas_ispeclinefit(speclist,specres=specres,snrcut=0.0,dustmodel=0,zsnrcross=5.0,$
         datapath=datapath,linepath=specfitpath,eigendir=eigendir,eigenfile=eigenfile,suffix=suffix,$
         Zmulti=Zmulti,/charlot,/zcrosscor,/postscript,/write,vmaxshift=500.0,/nologfile,starvdisp=100.0,$
         /nuclear,/positive_emission)

       atlas_parse_specfit, /nuclear

    endif
       
stop    
    
return
end    
