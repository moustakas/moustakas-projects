pro hcn_plots_shane, hcndust, hcnnodust, atlasdust1, atlasnodust1, $
  encapsulated=encapsulated, postscript=postscript, cleanpng=cleanpng, $
  nametag=nametag,lir_ahatot=lir_ahatot,lir_aha25=lir_aha25,$
  lir_ahaha=lir_ahaha,fhcn_fha24=fhcn_fha24,$
  lhcnco_lha24co=lhcnco_lha24co,lhcn_lha24=lhcn_lha24,$
  lir_lha24=lir_lha24,lhcn_lir=lhcn_lir,lhcn_lhacor=lhcn_lhacor
; jm06jul20uofa - written
; jm06dec11nyu - update
; sb07aug25uofa - update: divide into procedures

; rsync -auv hcn* howdy:"public_html/research/"

; read the samples 

path = hcn_path()

if (n_elements(hcndust) eq 0L) then begin
       hcndust = mrdfits(path+'atlas_hcn_speclinefit.fits.gz',1,/silent)
       hcnnodust = mrdfits(path+'atlas_hcn_speclinefit_nodust.fits.gz',$
         1,/silent)
;      hcndust = mrdfits('../analysis/atlas_hcn_speclinefit.fits.gz',1,/silent)
;      hcnnodust = mrdfits('../analysis/atlas_hcn_speclinefit_nodust.fits.gz',$
;        1,/silent)
endif
    
;if (n_elements(atlasdust1) eq 0L) then begin
;       atlasdust1 = mrdfits('../analysis/integrated_atlas_speclinefit.fits.gz',$
;         1,/silent)
;       atlasnodust1 = mrdfits('../analysis/integrated_atlas_speclinefit_'+$
;         'nodust.fits.gz',1,/silent)
;endif
;    
;indx = cmset_op(strtrim(atlasdust1.galaxy,2),'and',/not2,strtrim(hcndust.galaxy,2),/index)
;atlasdust = atlasdust1[indx] & atlasnodust = atlasnodust1[indx]

path = hcn_path()
;path = '../analysis/'
hcn = rsex(path+'04gao_apj.dat')

lsun = 3.826D33
snrcut = 2.0
lhcnerr = 0.1

lhcnaxis = findgen((10.0-6.0)/0.01+1)*0.01+6.0
;   lhcnaxis = findgen((!x.crange[1]-!x.crange[0])/0.05+1)*0.05+!x.crange[0]

; plotting and webpage variables    

htmlbase = 'html'

html_path = path
;pspath = path+htmlbase+'/'
pspath = path

if keyword_set(postscript) then begin
	postthick = 5.0 
	postthick2 = 8.0 
endif else begin
	postthick = 2
	postthick2 = 2
	window, 0, xsize=600, ysize=600, retain=2
endelse

@'xyrange_hcn'

hcnsfsym = 108 & hcnsfcolor = 'red' & hcnsfpsize = 1.5 & hcnsffill = 1
hcnsfsym_ul = 108 & hcnsfcolor_ul = 'red' & hcnsfpsize_ul = 4.0 
hcnsffill_ul = 1

hcnagnsym = 105 & hcnagncolor = 'blue' & hcnagnpsize = 1.8 & hcnagnfill = 0
hcnagnsym_ul = 105 & hcnagncolor_ul = 'blue' & hcnagnpsize_ul = 4.0 
hcnagnfill_ul = 1

hcnagnsfsym = 106 & hcnagnsfcolor = 'green' & hcnagnsfpsize = 1.8 
hcnagnsffill = 0 & hcnagnsfsym_ul = 106 & hcnagnsfcolor_ul = 'green' 
hcnagnsfpsize_ul = 4.0 & hcnagnsffill_ul = 1

hcnallsym = 108 & hcnallcolor = 'black' & hcnallpsize = 1.0 & hcnallfill = 1
hcnsym = 108 & hcncolor = 'red' & hcnpsize = 1.0 & hcnfill = 1
atlassym = 106 & atlascolor = 'grey' & atlaspsize = 1.0 & atlasfill = 1

; clean up old PS and PNG files

if keyword_set(cleanpng) then begin
	splog, 'Deleting all PNG and PS files in '+pspath
	spawn, ['rm -f '+pspath+'*.png*'], /sh
	spawn, ['rm -f '+pspath+'*ps'], /sh
endif

;indx1 = where((hcndust.gao_lhcn gt 0.0) and (hcndust.sfr[0] gt -900.0),nindx1)

; define B-V excess from hcnnodust structure
ebv = hcnnodust.ebv_hahb

; put L_IR, L_{Ha,obs}, L_{Ha,cor}, L_{IR,tot} in cgs units
lir = 10.^(hcnnodust.iras_25_lum[0])*lsun
haobs = 10.^(hcndust.h_alpha_lum[0])*lsun
hacor = 10.^(hcnnodust.h_alpha_lum[0])*lsun
lirtot = 10.^(hcnnodust.ir_lum[0])*lsun

; calculate A_Ha from Ha,obs/Ha,cor ratio
A_Ha = -2.5*alog10(haobs/hacor)

; set constant scaling factor to apply to L_IR
a = .031

; define A_IR, A_{IR,tot}
A_IR = 2.5*alog10(1. + a*lir/haobs)
A_IR_tot = 2.5*alog10(1. + a*lirtot/haobs)

; compute filter correction for IRAS25 --> Spitzer24
readcol,path+'25_24.dat',alpha,f25_24,f60_100
;readcol,'../analysis/25_24.dat',alpha,f25_24,f60_100
f25_24spline = interpol(f25_24,1000,/spline)
f60_100spline = interpol(f60_100,1000,/spline)
delta = 0.005

; filter out IRAS and Halpha non-detections 
indx = where((hcndust.iras_25_lum[0] gt -900.0) and $
  (hcndust.h_alpha_lum[0] gt -900.0),nindx)

; proceed if there is at least one common detection
if (nindx ne 0L) then begin

data25_24 = fltarr(nindx)

; speed of light in angstroms/sec
light = 3d18

; cm in one megaparsec
mpc2cm = 3.086D24

; iras 60 micron flux
f60 = hcndust[indx].iras_60*60D4*1D-23/(60.0D4)^1.0         

; iras 100 micron flux
f100 = hcndust[indx].iras_100*100D4*1D-23/(100.0D4)^1.0         

; surface area at distance of galaxy
area = 4.0*!dpi*(hcndust[indx].distance*mpc2cm)^2.0

; calculate IRAS 60,100 micron luminosity
iras_60_lnu = f60*area/lsun
iras_100_lnu = f100*area/lsun

; calculate 60/100 luminosity ratio
data60_100 = iras_60_lnu/iras_100_lnu

; cycle through all galaxies with detections
for i=0,nindx-1 do begin

	; determine index where 60/100 model ratio most closely matches data
	importance = where(f60_100spline lt data60_100[i]+delta and $
	  f60_100spline gt data60_100[i]-delta,nimp)
	
	; determine 25/24 ratio at that index
	f25_24importance = f25_24spline[importance]

	; store in data25_24 variable
	data25_24[i] = mean(f25_24importance)

endfor

print,minmax(data25_24)

; compute the combined SFR and associated error
hcndust[indx].sfr[0] = alog10(10.^(hcndust[indx].h_alpha_lum[0])*$
  lsun*5.3D-42 + 10.^(hcndust[indx].iras_25_lum[0])/data25_24*lsun*a*5.3D-42)
hcndust[indx].sfr[1] = sqrt(hcndust[indx].h_alpha_lum[1]^2 + $
  hcndust[indx].iras_25_lum[1]^2)

; old, incorrect equation to compute SFR
;hcndust[indx].sfr[0] = hcndust[indx].h_alpha_lum[0] + alog10(7.9D-42) + alog10(lsun) + $
;         hcndust[indx].iras_25_lum[0] + alog10(0.035*7.9D-42) + alog10(lsun) + alog10(0.5)
;       hcndust[indx].sfr[1] = sqrt(hcndust[indx].h_alpha_lum[1]^2 + hcndust[indx].iras_25_lum[1]^2)

; determine indx where hcn and Halpha dust-full data exist
indx = where((hcndust.gao_lhcn gt 0.0) and (hcndust.sfr[0] gt -900.0),nindx)

; determine indx containing upper limits
indx_ul = where((hcndust.gao_lhcn lt 0.0) and (hcndust.sfr[0] gt -900.0),$
  nindx_ul)

; put SFR from above in luminosity units: log(Lsun)
ha25 = alog10(10.^(hcndust[indx].sfr[0])/5.3d-42/lsun)

; same as above but using different method
ha25a = alog10(10.^(hcndust[indx].h_alpha_lum[0]) + $
 10.^(hcndust[indx].iras_25_lum[0])*a)

; compute same luminosity but with 24micron instead of 25micron
ha24a = alog10(10.^(hcndust[indx].h_alpha_lum[0]) + $
 10.^(hcndust[indx].iras_25_lum[0])/data25_24*a)

; compute same luminosity but with 24micron instead of 25micron
ha24a_ul = alog10(10.^(hcndust[indx_ul].h_alpha_lum[0]) + $
 10.^(hcndust[indx_ul].iras_25_lum[0])/data25_24*a)

endif

; determine indx where IR and Halpha dust-corrected exist
indx = where((hcnnodust.iras_25_lum[0] gt -900.0) and $
  (hcnnodust.h_alpha_lum[0] gt -900.0),nindx)

; if there's some data, compute the SFR from L_Ha,cor
if (nindx ne 0L) then begin

	hcnnodust[indx].sfr_h_alpha[0] = alog10(10.^($
	  hcnnodust[indx].h_alpha_lum[0])*lsun*5.3D-42)
	;hcnnodust[indx].sfr_h_alpha[1] = sqrt(hcnnodust[indx].h_alpha_lum[1]^2 + hcnnodust[indx].iras_25_lum[1]^2)

endif

;---------------------------------------------------
; Plots dealing with extinction correction A(Ha)

; PLOT: L_IR vs. A(Ha) [3-1000um]
;if keyword_set(lir_ahatot) then begin
;	lir_ahatot
;endif

; PLOT: L_IR vs. A(Ha) [25um]
;if keyword_set(lir_aha25) then begin
;	lir_aha25
;endif

; PLOT: L_IR vs. A(Ha) [Ha]
;if keyword_set(lir_ahaha) then begin
;	lir_ahaha
;endif

;----------------------------------------------------
; Distance corrected plots

; PLOT: Flux[HCN] vs Flux[(Ha)_obs+a*(24)]
;if keyword_set(fhcn_fha24) then begin
;	fhcn_fha24
;endif

; PLOT:  L(HCN)/L(CO) vs [L(Ha)_obs+L(24)]/L(CO)
if keyword_set(lhcnco_lha24co) then begin
	lhcnco_lha24co,pspath,hcndust,lhcnerr,lhcncorange,lharange4,ha24a,$
	  nametag=nametag,postscript=postscript,encapsulated=encapsulated,$
	  postthick,postthick2
endif

;-----------------------------------------------------
; Main plot!!

; PLOT: L(HCN) vs L(Ha)+a*L(24)
if keyword_set(lhcn_lha24) then begin
	lhcn_lha24,hcndust,pspath,postscript=postscript,$
	  encapsulated=encapsulated,lhcnrange,sfrhcnrange,sfrrange3,lharange3,$
	  ha24a,ha24a_ul,lhcnerr,nametag=nametag,postthick,postthick2
endif

; PLOT: L(HCN) vs SFR[L(Ha)_cor]
if keyword_set(lhcn_lhacor) then begin
	lhcn_lhacor,pspath,hcndust,lhcnerr,lhcnrange,lharange3,hcnnodust,$
	  nametag=nametag,sfrrange2,sfrhcnrange,$
	  postscript=postscript,encapsulated=encapsulated,$
	  postthick,postthick2
endif

;; PLOT: L(HCN) vs L(Ha)_obs
;if keyword_set(lhcn_lhaobs) then begin
;	lhcn_lhaobs
;endif

; PLOT: L(HCN) vs L(IR) only for subsample with Ha & HCN(1-0)
if keyword_set(lhcn_lir) then begin
	lhcn_lir,hcndust,hcn,pspath,postscript=postscript,$
	  encapsulated=encapsulated,lhcnerr,lhcnaxis,$
	  nametag=nametag,postthick,postthick2,lhcnrange2,lirrange2
endif

; PLOT: L(IR) vs L(Ha,obs)
;if keyword_set(lir_lhaobs) then begin
;	lir_lhaobs
;endif

; PLOT: L(IR) vs L(Ha,cor)
;if keyword_set(lir_lhacor) then begin
;	lir_lhacor
;endif

; PLOT L(IR) vs. L(Ha)+a*L(24)
if keyword_set(lir_lha24) then begin
	lir_lha24,hcndust,pspath,postscript=postscript,$
	  encapsulated=encapsulated,$
	  lirrange3,sfrhcnrange,sfrrange3,lharange3,ha24a,ha24a_ul,lhcnerr,$
	  nametag=nametag,postthick,postthick2

endif

; this could be important too
; PLOT: L(HCN) vs L(IR) GS04 reproduction
;if keyword_set(lhcn_lir_gao) then begin
;	lhcn_lir_gao
;endif

; PLOT: BPT Diagram
;if keyword_set(bpt) then begin
;	bpt
;endif

; --------------------------------------------------    
; GENERATE THE WEB PAGE
; --------------------------------------------------    

stop

if keyword_set(postscript) and keyword_set(encapsulated) and (not keyword_set(blackwhite)) then begin
im_ps2html, htmlbase, html_path=html_path, cleanpng=0, npscols=3, _extra=extra
endif

return
end

