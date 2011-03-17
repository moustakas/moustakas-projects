pro trgb_readata, objname, datapath, data, infobase, hst=hst, halo=halo, core=core, ccut=ccut
;+
; NAME:
;	TRGB_READATA
;
; PURPOSE:
;	Read in calibrated Keck or HST photometry data.
;
; INPUTS:
;	objname  : string name of the galaxy
;
; KEYWORD PARAMETERS:
;	hst  	: keyword specifying an HST object
;	halo 	: use the objname_starlist_halo.dat photometry file
;	          created in TRGB_REGIONS
;	core	: use the objname_starlist_core.dat photometry file
;	          created in TRGB_REGIONS
;	ccut	: use the objname_starlist_ccut.dat photometry file
;		  created in TRGB_COLORCUT
;
; OUTPUTS:
;	datapath : the full datapath to the data
;	data     : the data array
;	infobase : TRGB_OBJECT_DATA information structure on the
;		   selected object 
;
; COMMON BLOCKS:
;
; EXAMPLE:
;	trgb_readata, 'SextansB', datapath, data, infobase, /halo
;				OR
;	trgb_readata, 'ugc07577', datapath, data, infobase, /hst
;
; PROCEDURES USED:
;	READFAST, SREAD(), TRGB_DATAPATH()
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 26, UCB
;-

	paths = trgb_datapath()

        if keyword_set(hst) then begin ; HST data

            pushd, paths[0]+objname
            datapath = paths[0]+objname

            if (keyword_set(ccut) and keyword_set(halo)) then $
              filename = datapath+'/'+objname+'_starlist_halo_ccut.dat' else $
              if (keyword_set(ccut) and keyword_set(core)) then $
              filename = datapath+'/'+objname+'_starlist_core_ccut.dat' else $
              if keyword_set(ccut) then filename = datapath+'/'+objname+'_starlist_ccut.dat' else $
              if keyword_set(halo) then filename = datapath+'/'+objname+'_starlist_halo.dat' else $
              if keyword_set(core) then filename = datapath+'/'+objname+'_starlist_core.dat' else $
              filename = datapath+'/'+objname+'_starlist.dat'

            print & print, 'Reading '+filename+'.' & print

; retrieve the data structure            
            
            info = sread('/deep1/ioannis/trgb/results/hst_object.dat')
            indx = where(strupcase(info.object) eq strupcase(objname))
            infobase = info[indx]
            if infobase.color_info eq 'YES' then ncol = 8 else ncol = 5
        
            readfast, filename, data, ncols=ncol, skip=2	; read the data

        endif else begin ; KECK data

            objdata = ['holmbergii','holmbergix','i342','ngc2366', 'ngc1560', $
                       'ngc2903','ngc2976','ngc3109','sextansb']
            check = where(strupcase(objdata) eq strupcase(objname),count)
            if count gt 0L then date = '22dec97/' else date = '23dec97/'

            pushd, paths[2]+date+strlowcase(objname)
            datapath = paths[2]+date+strlowcase(objname)

            if (keyword_set(ccut) and keyword_set(halo)) then $
              filename = datapath+'/'+objname+'_IVmags_halo_ccut.dat' else $
              if (keyword_set(ccut) and keyword_set(core)) then $
              filename = datapath+'/'+objname+'_IVmags_core_ccut.dat' else $
              if keyword_set(ccut) then filename = datapath+'/'+objname+'_IVmags_ccut.dat' else $
              if keyword_set(halo) then filename = datapath+'/'+objname+'_IVmags_halo.dat' else $
              if keyword_set(core) then filename = datapath+'/'+objname+'_IVmags_core.dat' else $
              filename = datapath+'/'+objname+'_IVmags.dat'

            print & print, 'Reading '+filename+'.' & print

; retrieve the data structure            

            info = sread('/deep1/ioannis/trgb/results/keck_object.dat')
            indx = where(strupcase(info.object) eq strupcase(objname))
            infobase = info[indx]

            readfast, filename, data, skip=2, ncols=8

        endelse

        popd

return
end
