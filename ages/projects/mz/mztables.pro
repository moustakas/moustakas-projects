function get_mzlzfit_info, fluxcor=fluxcor
    mzpath = mz_path()

    if keyword_set(fluxcor) then begin
       suffix = 'fluxcor' 
       calib = ['KK04','M91','T04','MPA-JHU']
    endif else begin
       suffix = 'ews'
       calib = ['KK04','M91','T04']
    endelse
    ncalib = n_elements(calib)

    info = replicate({calib: '', empty1: ' ', $
      ohstar: '', mstar: '', gamma: '', mzscatter: '', empty2: ' ', $
      c0: '', c1: '', lzscatter: ''},ncalib)
    info.calib = calib

    for ii = 0, ncalib-1 do begin
       calib1 = repstr(strlowcase(calib[ii]),'-','')

; MZ/closed box
       mzfit = mrdfits(mzpath+'mzlocal_sdss_'+suffix+'_'+calib1+'.fits.gz',1)
       info[ii].ohstar = '$'+im_sigfigs(mzfit.coeff[0],4)+'$'
       info[ii].mstar  = '$'+im_sigfigs(mzfit.coeff[1],4)+'$'
       info[ii].gamma = '$'+im_sigfigs(mzfit.coeff[2],3)+'$'
       info[ii].mzscatter = '$'+string(mzfit.scatter,format='(F4.2)')+'$'
;      info[ii].mzscatter = '$'+im_sigfigs(mzfit.scatter,3)+'$'

; LZ
       lzfit = mrdfits(mzpath+'lzlocal_sdss_'+suffix+'_'+calib1+'.fits.gz',1)
       info[ii].c0 = '$'+im_sigfigs(lzfit.coeff[0],4)+'$'
       info[ii].c1 = '$'+im_sigfigs(lzfit.coeff[1],3)+'$'
       info[ii].lzscatter = '$'+im_sigfigs(lzfit.scatter,3)+'$'
    endfor
    struct_print, info

return, info
end

function mzget_colhead, head, nobreak=nobreak
    if keyword_set(nobreak) then br = '' else br = '\\'
    nhead = n_elements(head)
    colhead = '\colhead{'+[head[0:nhead-2]+'} & ',head[nhead-1]+'} '+br]
    multi = where(strmatch(colhead,'*multicolumn*',/fold))
    if (multi[0] ne -1) then colhead[multi] = repstr(repstr(colhead[multi],'\colhead{',''),'}}','}')
return, colhead
end

pro mztables, preprint=preprint
; jm07apr11nyu - based on WRITE_ATLAS_TABLES
; jm10oct12ucsd - major rewrite

; by default build emulateapj format tables    
    if keyword_set(preprint) then filesuffix = '_preprint' else $
      filesuffix = '_apj'

    mzpath = ages_path(/projects)+'mz/'
    paperpath = ages_path(/papers)+'mz/'

; ---------------------------------------------------------------------------    
; Table: comparison with the literature

    table = {sample: '', ngal: '', area: '', zrange: '', $
      diagnostic: '', mlrange: '', dlogoh: '', dz: ''};, remark: ''}
    ntags = n_tags(table)

; caputi08
; rupke08
    
    sample = [$
; LZ relations
      '\citet{kobulnicky99b}',$
;     'CL01',$ ; Carollo & Lilly+01
      '\citet{lilly03a}\tablenotemark{h}',$
      '\citet{kobulnicky03b}',$
      '\citet{kobulnicky04a}',$
      '\citet{maier04a}',$ ; 'Maier+04',$
      '\citet{liang04a}',$ ; 'Liang+04',$
      '\citet{savaglio05a}\tablenotemark{i}',$ ; 'Savaglio+05',$
      '\citet{maier05a}',$ ; 'Maier+05',$
      '\citet{hoyos05a}',$ ; 'Hoyos+05',$
      '\citet{mouhcine06b}',$ ; 'Mouhcine+06a (cluster)',$
      '\citet{mouhcine06a}',$ ; 'Mouhcine+06b (field)',$
      '\citet{lama06b, lama06a}',$ ; 'Lamareille+06',B-band
      '\citet{lama09a}\tablenotemark{k}',$ ; VVDS-DEEP
      '\citet{lama09a}\tablenotemark{k}',$ ; VVDS-WIDE
      '\citet{zahid11a}',$ ; DEEP2
; MZ relations
      '\citet{savaglio05a}\tablenotemark{i}',$ ; 'Savaglio+05\tablenotemark{g}',$ ; MZ
;     'LHF06',$ ; Liang+06
      '\citet{rodrigues08a}\tablenotemark{m}',$
      '\citet{cowie08a}',$
      '\citet{lama09a}\tablenotemark{k}',$ ; VVDS-DEEP
      '\citet{lama09a}\tablenotemark{k}',$ ; VVDS-WIDE
      '\citet{zahid11a}'] ; DEEP2
      
    ntable = n_elements(sample)
    table = replicate(table,ntable)
    table.sample = sample

    table.ngal = [$
; LZ relations
      '$14$',$   ; KZ99
;     '$15$',$   ; CL01
      '$66$',$   ; Lilly+03
      '$64$',$   ; Kobulnicky+03
      '$204$',$  ; KK04
      '$16$',$   ; Maier+04
      '$46$',$   ; Liang+04 [Table 6, col 12 plus Table 7, col 8]
      '$56$',$   ; Savaglio+05 - LZ
      '$30$',$   ; Maier+05
      '$15$',$   ; Hoyos+05
      '$15$',$   ; Mouhcine+06a
      '$39$',$   ; Mouhcine+06b
      '$117$',$  ; Lama+06 ; B-band (see Sec 2.1.2); only 117 in '06lama.dat' with good M_B and O/H
      '$870$',$  ; Lama+09
      '$2191$',$ ; Lama+09
      '$1350$',$ ; Zahid+11
; MZ relations
      '$56$',$   ; Savaglio+05
;     '$34$',$ ; Liang+06
      '$88$',$ ; Rodrigues+08 (58 with lines; 7 are AGN, but they don't say which; 35 have good O/H; 30 from Liang06)
      '$\sim400$',$ ; Cowie+08 (rough estimate based on Fig 32)
      '$870$',$  ; Lama+09
      '$2191$',$ ; Lama+09
      '$1350$']  ; Zahid+11

; Area surveyed by Lama+09 (numbers are approximate): 
;   Deep = F02+CDFS = 0.61 deg^2 (Le Fevre+04)
;   Wide = F22+F10+F14+F02+CDFS = 3.0+0.6+0.9+0.5+0.11 (Le Fevre+04;
;     Garilli+08) = 
    
    table.area = [$
; LZ relations
      '?',$       ; KZ99
;     '?',$       ; CL01
      '$0.03$',$   ; Lilly+03 [Lilly+95, Sec 2.2]
      '$0.04$',$   ; Kobulnicky+03
      '$0.09$',$   ; KK04 [Wirth+04]
      '$0.11$',$   ; Maier+04 [105+98+107+103; see Maier+03]
      '?',$       ; Liang+04
      '$0.06$',$   ; Savaglio+05 (GDDS was 4x30 arcmin^2 + 112 from CFRS)
      '?',$       ; Maier+05
      '$0.09$',$   ; Hoyos+05 [TKRS]
      '?',$       ; Mouhcine+06a
      '?',$       ; Mouhcine+06b
      '?',$       ; Lama+06
      '$0.61$',$  ; Lama+09
      '$5.1$',$ ; Lama+09
      '$3.5$',$ ; Zahid+11
; MZ relations
      '0.06',$     ; Savaglio+05 (GDDS was 4x30 arcmin^2 + 112 from CFRS)
;     '?',$       ; Liang+06
      '$0.05$',$   ; Rodrigues+08 (four 6.8x6.8 arcmin^2 fields)
      '$0.04$',$   ; Cowie+08
      '$0.61$',$   ; Lama+09
      '$5.1$',$  ; Lama+09
      '$3.5$']   ; Zahid+11

;    table.area = [$
;; LZ relations
;      '?',$       ; KZ99
;;     '?',$       ; CL01
;      '$112$',$   ; Lilly+03 [Lilly+95, Sec 2.2]
;      '$127$',$   ; Kobulnicky+03
;      '$320$',$   ; KK04 [Wirth+04]
;      '$413$',$   ; Maier+04 [105+98+107+103; see Maier+03]
;      '?',$       ; Liang+04
;      '$232$',$   ; Savaglio+05 (GDDS was 4x30 arcmin^2 + 112 from CFRS)
;      '?',$       ; Maier+05
;      '$320$',$   ; Hoyos+05 [TKRS]
;      '?',$       ; Mouhcine+06a
;      '?',$       ; Mouhcine+06b
;      '?',$       ; Lama+06
;      '$2200$',$  ; Lama+09
;      '$18400$',$ ; Lama+09
;; MZ relations
;      '232',$     ; Savaglio+05 (GDDS was 4x30 arcmin^2 + 112 from CFRS)
;;     '?',$       ; Liang+06
;      '$185$',$   ; Rodrigues+08 (four 6.8x6.8 arcmin^2 fields)
;      '$145$',$   ; Cowie+08
;      '$2200$',$   ; Lama+09
;      '$18400$']  ; Lama+09

    table.zrange = [$
; LZ relations
      '$0.1-0.5$',$  ; KZ99
;     '$0.60-0.99$',$  ; CL01
      '$0.5-0.9$',$ ; Lilly+03
      '$0.3-0.8$',$ ; Kobulnicky+03
      '$0.3-0.9$',$ ; KK04
      '$0.4-0.6$',$ ; Maier+04
      '$0.1-0.9$',$ ; Liang+04
      '$0.5-1.0$',$ ; Savaglio+05 - LZ
      '$0.5-0.9$',$ ; Maier+05
      '$0.5-0.9$',$ ; Hoyos+05
      '$0.3-0.6$',$ ; Mouhcine+06a
      '$0.2-0.7$',$ ; Mouhcine+06b
      '$0.1-1.0$',$ ; Lama+06
      '$0.1-0.9$',$ ; Lama+09
      '$0.1-0.9$',$ ; Lama+09
      '$0.75-0.82$',$ ; Zahid+11
; MZ relations
      '$0.5-1.0$',$ ; Savaglio+05
;     '$-$',$   ; Liang+06
      '$0.4-1.0$',$   ; Rodrigues+08
      '$0.05-0.9$',$   ; Cowie+08
      '$0.1-0.9$',$   ; Lama+09
      '$0.1-0.9$',$   ; Lama+09
      '$0.75-0.82$'] ; Zahid+11

;    table.zrange = [$
;; LZ relations
;      '$0.11-0.50$',$  ; KZ99
;;     '$0.60-0.99$',$  ; CL01
;      '$0.47-0.92$',$ ; Lilly+03
;      '$0.26-0.82$',$ ; Kobulnicky+03
;      '$0.30-0.93$',$ ; KK04
;      '$0.40-0.64$',$ ; Maier+04
;      '$0.09-0.88$',$ ; Liang+04
;      '$0.47-0.96$',$ ; Savaglio+05 - LZ
;      '$0.47-0.90$',$ ; Maier+05
;      '$0.51-0.85$',$ ; Hoyos+05
;      '$0.31-0.59$',$ ; Mouhcine+06a
;      '$0.20-0.71$',$ ; Mouhcine+06b
;      '$0.07-1.00$',$ ; Lama+06
;      '$0.1-0.9$',$ ; Lama+09
;      '$0.1-0.9$',$ ; Lama+09
;; MZ relations
;      '$0.47-0.96$',$ ; Savaglio+05
;;     '$-$',$   ; Liang+06
;      '$0.40-0.96$',$   ; Rodrigues+08
;      '$0.05-0.9$',$   ; Cowie+08
;      '$0.1-0.9$',$   ; Lama+09
;      '$0.1-0.9$']    ; Lama+09

    table.diagnostic = [$
; LZ relations
      '$T_{e}$, \pagel',$ ; KZ99 (ZKH94)
;     '\pagel',$   ; CL01 (M91)
      '\pagel',$   ; Lilly+03 (M91)
      '\ewpagel',$ ; Kobulnicky+03 (M91)
      '\ewpagel',$ ; KK04 (M91+KD02)
      '\pagel',$   ; Maier+04 (M91)
      '\pagel',$   ; Liang+04 (M91?)
      '\pagel',$   ; Savaglio+05 (M91+KD02)
      '\oii, \hb, \oiii, \ha, \nii\tablenotemark{j}',$  ; Maier+05 (KD02)
      '$T_{e}$',$  ; Hoyos+05
      '\ewpagel',$ ; Mouhcine+06a
      '\ewpagel',$ ; Mouhcine+06b
      '\ewpagel',$ ; Lama+06
      'Various\tablenotemark{l}',$ ; Lama+09 (N2, O3N2, M91)
      'Various\tablenotemark{l}',$ ; Lama+09 (N2, O3N2, M91)
      '\ewpagel',$ ; Zahid+11
;     'EW(\nii/\ha); EW(\oiii/\nii); EW(\pagel)',$ ; Lama+09 (N2, O3N2, M91)
;     'EW(\nii/\ha); EW(\oiii/\nii); EW(\pagel)',$ ; Lama+09 (N2, O3N2, M91)
; MZ relations
      '\pagel',$   ; Savaglio+05
;     '\pagel',$   ; Liang+06
      '\pagel',$   ; Rodrigues+08 (T04)
      '\pagel',$   ; Cowie+08 ; T04 mostly, but they used various
      'Various\tablenotemark{l}',$ ; Lama+09 (N2, O3N2, M91)
      'Various\tablenotemark{l}',$ ; Lama+09 (N2, O3N2, M91)
      '\ewpagel'] ; Zahid+11
;     'EW(\nii/\ha); EW(\oiii/\nii); EW(\pagel)',$ ; Lama+09 (N2, O3N2, M91)
;     'EW(\nii/\ha); EW(\oiii/\nii); EW(\pagel)']  ; Lama+09 (N2, O3N2, M91)

    table.mlrange = [$
; LZ relations
      '$-22.2<\mb<-15.9$',$ ; KZ99
      '$-22.1<\mb<-19.4$',$ ; Lilly+03
      '$-22.7<\mb<-16.5$',$ ; Kobulnicky+03
      '$-22.3<\mb<-16.9$',$ ; KK04
      '$-19.3<\mb<-17.4$',$ ; Maier+04
      '$-22.4<\mb<-17.0$',$ ; Liang+04
      '$-22.1<\mb<-17.3$',$ ; Savaglio+05 - LZ
      '$-21.8<\mb<-19.6$',$ ; Maier+05
      '$-21.5<\mb<-18.6$',$ ; Hoyos+05
      '$-22.4<\mb<-20.1$',$ ; Mouhcine+06a
      '$-23.1<\mb<-19.2$',$ ; Mouhcine+06b
      '$-22.1<\mb<-15.4$',$ ; Lama+06
      '$-22\lesssim\mb\lesssim-16.5$',$ ; Lama+09 - DEEP
      '$-22\lesssim\mb\lesssim-16.5$',$ ; Lama+09 - WIDE
      '$-22.0<\mb<-19.5$',$ ; Zahid+11
; MZ relations
      '$8.3<\log(\mathcal{M}/\msun)<10.9$',$ ; Savaglio+05 (M_Kroupa = 1.2*M_Savaglio)
;     '\nodata',$ ; Liang+06
      '$9\lesssim\log(\mathcal{M}/\msun)\lesssim11$',$ ; Rodrigues+08
      '$9\lesssim\log(\mathcal{M}/\msun)\lesssim11$',$ ; Cowie+08
      '$8\lesssim\log(\mathcal{M}/\msun)\lesssim11$',$ ; Lama+09 - DEEP
      '$8\lesssim\log(\mathcal{M}/\msun)\lesssim11$',$ ; Lama+09 - WIDE
      '$9.2<\log(\mathcal{M}/\msun)<10.7$'] ; Zahid+11, Kroupa IMF

    table.dlogoh = [$
; LZ relations
      '$\sim0.1$',$     ; KZ99
      '$0.08\pm0.06$',$ ; Lilly+03
      '$\sim0.15$',$        ; Kobulnicky+03
      '$0.14\pm0.05$',$ ; KK04
      '$\sim0.3$',$ ; Maier+04
      '$\sim0.3$',$ ; Liang+04
      '$\sim0.4$',$ ; Savaglio+05 (estimated by eye from Fig. 11 at M_B=-20)
      '$\sim0.3$',$ ; Maier+05
      '$\sim0.4$',$ ; Hoyos+05 (estimated by eye from Fig 2a at M_B=-20)
      '\nodata',$ ; Mouhcine+06a
      '\nodata',$ ; Mouhcine+06a
      '$0.28-0.55$',$ ; Lama+06
      '$0.24-0.39$',$ ; Lama+09 - DEEP
      '$0.58-0.70$',$ ; Lama+09 - WIDE
      'numbers!',$ ; Zahid+11
; MZ relations
      '$\sim0.1$',$ ; Savaglio+05 (estimated by eye from Fig. 13 at 10^10 M_sun)
;     '\nodata',$   ; Liang+06
      '$0.39\pm0.09$',$   ; Rodrigues+08
      '$0.21\pm0.03$',$   ; Cowie+08
      '$0.12-0.23$',$   ; Lama+09 - DEEP
      '$0.23-0.36$',$    ; Lama+09 - WIDE
      'numbers!'] ; Zahid+11

    table.dz = [$
; LZ relations
      '$0.4$',$       ; KZ99
      '$0.75$',$      ; Lilly+03
      '$0.7$',$       ; Kobulnicky+03
      '$1.0$',$       ; KK04
      '$\sim0.5$',$  ; Maier+04 [actually, z~0.52]
      '$\sim0.65$',$  ; Liang+04
      '$\sim0.7$',$   ; Savaglio+05 - LZ
      '$\sim0.7$',$   ; Maier+05
      '$\sim0.7$',$   ; Hoyos+05
      '\nodata',$ ; Mouhcine+06a
      '\nodata',$ ; Mouhcine+06a
      '$0.9$',$     ; Lama+06
      '$0.8$',$     ; Lama+09 - DEEP
      '$0.8$',$     ; Lama+09 - WIDE
      '$0.8$',$     ; Zahid+11
; MZ relations
      '$\sim0.7$',$ ; Savaglio+05 (estimated by eye from Fig. 13 at 10^10 M_sun)
;     '\nodata',$ ; Liang+06
      '$0.85$',$ ; Rodrigues+08
      '$0.75$',$ ; Cowie+08
      '$0.8$',$ ; Lama+09 - DEEP
      '$0.8$',$  ; Lama+09 - WIDE
      '$0.8$']   ; Zahid+11

;    table.ref = [$
;; LZ relations
;      '1',$ ; KZ99
;      '2',$ ; Lilly+03
;      '3',$ ; Kobulnicky+03
;      '4',$ ; KK04
;      '5',$ ; Maier+04
;      '6',$ ; Liang+04
;      '7',$ ; Savaglio+05
;      '8',$ ; Maier+05
;      '9',$ ; Hoyos+05
;      '10',$ ; Mouhcine+06a
;      '11',$ ; Mouhcine+06b
;      '12',$ ; Lama06
;      '13',$ ; Lama09
;      '13',$ ; Lama09
;; MZ relations
;      '7',$ ; Savaglio+05
;;     '12',$ ; Liang+06
;      '14',$ ; Rodrigues+08
;      '15',$ ; Cowie+08
;      '13',$ ; Lama09
;      '13']  ; Lama09
    
;   table.remark = [$
;     'stuff',$
;     'Assume $\langle E(B-V)=0.3 \rangle$',$
;     'IR-Balmer correction',$
;     'Mean extinction correction',$
;     'stuff',$
;     'stuff',$
;     'stuff',$
;     'stuff',$
;     'stuff',$
;     'stuff']
    
;    refs = '(1) \citet{kobulnicky99b}; '+$
;      '(2) \citet{lilly03a}; '+$
;      '(3) \citet{kobulnicky03b}; '+$
;      '(4) \citet{kobulnicky04a}; '+$
;      '(5) \citet{maier04a}; '+$
;      '(6) \citet{liang04a}; '+$
;      '(7) \citet{savaglio05a}; '+$
;      '(8) \citet{maier05a}; '+$
;      '(9) \citet{hoyos05a}; '+$
;      '(10) \citet{mouhcine06b}; '+$
;      '(11) \citet{mouhcine06a}; '+$
;      '(12) \citet{lama06b, lama06a}; '+$
;      '(13) \citet{lama09a}; '+$
;      '(14) \citet{rodrigues08a}; '+$
;      '(15) \citet{cowie08a}.'

    colhead1 = $
      '\colhead{} & '+$
      '\colhead{} & '+$
      '\colhead{$\Omega$\tablenotemark{b}} & '+$
      '\colhead{Redshift} & '+$
      '\colhead{Abundance} & '+$
      '\colhead{Luminosity or} & '+$
      '\colhead{\dlogoh\tablenotemark{f}} & '+$
      '\colhead{} \\'
;     '\colhead{} & '+$
;     '\colhead{} \\'
    colhead2 = $
      '\colhead{Study} & '+$
      '\colhead{$N$\tablenotemark{a}} & '+$
      '\colhead{(deg$^2$)} & '+$
      '\colhead{Range\tablenotemark{c}} & '+$
      '\colhead{Diagnostics\tablenotemark{d}} & '+$
      '\colhead{Stellar Mass Range\tablenotemark{e}} & '+$
      '\colhead{(dex)} & '+$
      '\colhead{$\langle z \rangle$\tablenotemark{g}}'
;     '\colhead{Reference} & '+$
;     '\colhead{Remark}'

    texcenter = ['l',replicate('c',ntags-1)]
 
;   tablecomments = ['This table is not an exhaustive summary of results from the literature, but '+$
;     'we have attempted to include all the important studies.  We have not attempted to homogenize the '+$
;     'various methods used to derive rest-frame quantities such as luminosities and stellar masses oxygen abundances '+$
;     '(diagnostic used, stellar absorption, reddening correction), or other .']
    caption = 'Summary of Previously Published Chemical Evolution Measurements for '+$
      'Star-Forming Galaxies at $z<1$\label{table:litoh}'
    
    tablenotetext = [$
      '{a}{Number of galaxies in the sample with well-measured oxygen abundances, optical luminosities, and '+$
      'stellar masses.  In some studies the final number of objects used is not readily apparent from the '+$
      'information given in the paper, in which case we provide our best estimate of the correct number.}',$
      '{b}{Solid angle of the parent spectroscopic survey.  Heterogenously selected samples, or samples '+$
      'where we could not estimate the solid angle from the information given in the paper are indicated '+$
      'with a question mark.}',$
      '{c}{Approximate range of redshifts spanned by the sample of star-forming emission-line galaxies.}',$
      '{d}{See the corresponding paper and \S\ref{sec:oh} for the definitions of the abundance diagnostics '+$
      'indicated this column.  Electron temperature ($T_{e}$) abundance measurements are based on a '+$
      'detection of the temperature-sensitive \oiii~$\lambda4363$ line.}',$
      '{e}{We have converted the published luminosities and stellar masses '+$
      'to our adopted cosmology ($\Omega_{\rm m}=0.3$, $\Omega_{\Lambda}=0.7$, $H_{0}=70~\hubbleunits$), and to '+$
      'the AB magnitude system and the \citet{chabrier03a} initial mass function as necessary.}',$
      '{f}{Change in the mean oxygen abundance of star-forming galaxies since '+$
      '$\langle z \rangle$.  Note that a larger positive numbers indicates that galaxies were more metal-poor '+$
      'at higher redshift, which is opposite the convention used in \S\ref{sec:dlogoh} and Table~\ref{table:dlogohevol}.}',$
      '{g}{Mean redshift corresponding to \dlogoh.}',$
      '{h}{The absolute $B$-band magnitudes in \citet{lilly03a} were derived using a different cosmological model '+$
      'than what is stated in the paper; $0.45$~mag must be subtracted to transform them to a cosmological model '+$
      'with $\Omega_{M}=0.3$, $\Omega_{\Lambda}=0.7$, and $H_{0}=70$~\hubbleunits{} (C.~Maier 2006, private communication). '+$
      'In addition, the measurements reported for CFRS~22.0576 and CFRS~22.0599 must be switched in the '+$
      'published tables.}',$
      '{i}{Savaglio combined their sample from the Gemini Deep Deep Survey \citep[GDDS;][]{abraham04a} '+$
      'with previously published data from the Canada-France Redshift Survey \citep[CFRS;][]{lilly95a, lilly03a}.}',$
      '{j}{\citet{maier05a} modeled all five strong lines simultaneously using the \citet{kewley02a} photoionization models.}',$
      '{k}{The first and second row refers to the analysis by \citet{lama09a} of the VDDS-DEEP and '+$
      'VDDS-WIDE sample, respectively.}',$
      '{l}{\citet{lama09a} derived EW-based abundances using three different abundance diagnostics depending '+$
      'on redshift: \nii/\ha{} at $z<0.2$; (\oiii/\hb)/(\nii/\ha) \citep{pettini04a} at '+$
      '$z=0.2-0.5$; and \pagel{} \citep{mcgaugh91a} at $z=0.5-0.9$.  They then rescaled each '+$
      'abundance measurement to the \citet{tremonti04a} metallicity scale using the \citet{kewley08a} '+$
      'polynomial correction formulae.}',$
      '{m}{\citet{rodrigues08a} combines new spectroscopy with previously published measurements by '+$
      '\citet{liang06a}.}']

; write out
    texfile = paperpath+'mztable_litoh'+filesuffix+'.tex'
    splog, 'Writing '+texfile
    openw, lun, texfile, /get_lun
    if (keyword_set(preprint) eq 0) then begin
       printf, lun, '\LongTables'
       printf, lun, '\begin{landscape}'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;   printf, lun, '\tabletypesize{\tiny}'
    printf, lun, '\tablecaption{'+caption+'}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    printf, lun, '}'
    printf, lun, '\startdata'
    for ii = 0, ntable-1 do begin
       if (ii eq 0) then printf, lun, '\cutinhead{Luminosity-Metallicity Relations}'
       if (ii eq ntable-6) then printf, lun, '\cutinhead{Mass-Metallicity Relations}'
       line = strarr(ntags)
       for jj = 0L, ntags-1 do begin
          if (jj lt ntags-1) then suffix = ' & ' else if (ii lt ntable-1) then $
            suffix = ' \\ ' else suffix = ''
          line[jj] = table[ii].(jj)+suffix
       endfor
       printf, lun, line
    endfor
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
;   printf, lun, '\tablerefs{'+refs+'}'
    printf, lun, '\end{deluxetable}'
    printf, lun, '\clearpage'
    printf, lun, '\end{landscape}'
    free_lun, lun

stop    
    
; ---------------------------------------------------------------------------    
; mean metallicity in bins of mass and redshift, for each calibration
; (Fig 14)
    texfile = paperpath+'mztable_z_vs_oh_bymass'+filesuffix+'.tex'
    splog, 'Writing '+texfile
    
; read the fitting results
    zbins = mz_zbins(nz)
    sdss_zbins = mz_zbins(/sdss)
    massbins = mz_massbins(nmassbins)

    mass = '\log\,(\mass/\msun)'
    masslabel = ['$'+mass+'>11$','$10.5<'+mass+'<11$','$10<'+mass+'<10.5$','$9.5<'+mass+'<10$']
    
;   colhead1 = mzget_colhead(['Redshift','Median','\multicolumn{3}{c}{$12+\log\,(\textrm{O}/\textrm{H})$}'])
;   colhead2 = mzget_colhead(['Range','Redshift','KK04','T04','M91'],/nobreak)
;   colhead1 = mzget_colhead(['Redshift','','','','',''])
    colhead2 = mzget_colhead(['Redshift Range','$N$\tablenotemark{b}','$z_{\rm med}$\tablenotemark{c}',$
      '$\langle12+\log\,(\textrm{O}/\textrm{H})_{\rm KK04}\rangle$',$
      '$\langle12+\log\,(\textrm{O}/\textrm{H})_{\rm T04}\rangle$',$
      '$\langle12+\log\,(\textrm{O}/\textrm{H})_{\rm M91}\rangle$'],/nobreak)
    texcenter = ['c','c','c','c','c','c']

    tablenotetext = [$
      '{a}{Weighted mean metallicity of galaxies in multiple bins of stellar mass and redshift, '+$
      'based on the KK04, T04, and M91 calibrations.  Note that the metallicities in the first row of '+$
      'each stellar mass interval correspond to our SDSS sample, while the metallicities in the other rows are based on '+$
      'our AGES sample.}',$
      '{b}{Number of galaxies in this redshift interval and stellar mass bin.}',$
      '{c}{Median redshift of the subsample.}']
    
    openw, lun, texfile, /get_lun
    printf, lun, '\begin{deluxetable*}{'+strjoin(texcenter)+'}[!h]'
    printf, lun, '\tablecaption{Mean Oxygen Abundance in Bins of Stellar Mass \& '+$
      'Redshift\tablenotemark{a}\label{table:oh_bymass}}'   
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
;   niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    printf, lun, '}'
    printf, lun, '\startdata'

    calib = ['kk04','t04','m91']
    ncalib = n_elements(calib)
    
    for mm = 0, nmassbins-2 do begin
       printf, lun, '\multicolumn{6}{c}{'+masslabel[mm]+'} \\'
       printf, lun, '\cline{1-6}'

       info = replicate({zrange: '\nodata', ngal: '\nodata', medz: '\nodata', $
         oh_kk04: '\nodata', oh_t04: '\nodata', oh_m91: '\nodata'},nz+1) ; SDSS+AGES

; get the mean number of galaxies       
       allngal = lonarr(nz+1,ncalib)
       for ii = 0, ncalib-1 do begin
          mzevol = mrdfits(mzpath+'mzevol_'+calib[ii]+'.fits.gz',1,/silent)
          good = where(mzevol.ohmean_bymass[mm,*] gt -900.0)
;         print, [mzevol.sdss_ngal_bymass[mm],reform(mzevol.ngal_bymass[mm,*])]
          print, [mzevol.sdss_ngal_bymass[mm],reform(mzevol.ngal_bymass[mm,good])]
          allngal[[0,good+1],ii] = [mzevol.sdss_ngal_bymass[mm],reform(mzevol.ngal_bymass[mm,good])]
       endfor
       ngal = lonarr(nz+1)
       for zz = 0, nz do ngal[zz] = djs_median(allngal[zz,*])

       for ii = 0, ncalib-1 do begin
          mzevol = mrdfits(mzpath+'mzevol_'+calib[ii]+'.fits.gz',1,/silent)
          good = where(mzevol.ohmean_bymass[mm,*] gt -900.0)

          if (ii eq 0) then begin
             info.zrange = '$'+string([sdss_zbins.zlo,zbins.zlo],format='(F4.2)')+'-'+$
               string([sdss_zbins.zup,zbins.zup],format='(F4.2)')+'$'
             info[[0,good+1]].ngal = '$'+string(ngal[[0,good+1]],format='(I0)')+'$'
             info[[0,good+1]].medz = '$'+string([mzevol.sdss_medz_bymass[mm],$
               reform(mzevol.medz_bymass[mm,good])],format='(F4.2)')+'$'
          endif
          info[[0,good+1]].(ii+3) = '$'+string([mzevol.sdss_ohmean_bymass[mm],$ ; offset from the redshifts
            reform(mzevol.ohmean_bymass[mm,good])],format='(F5.3)')+$
            '\pm'+string([mzevol.sdss_ohmean_bymass_err[mm],$
            reform(mzevol.ohmean_bymass_err[mm,good])]>0.01,format='(F4.2)')+'$'
       endfor
       struct_print, info

       ntags = n_tags(info)
       for kk = 0, nz do begin
          for jj = 0, ntags-1 do begin
             if (jj eq ntags-1) then if (kk eq nz) then suffix = '' else $
               suffix = ' \\' else suffix = ' & '
             printf, lun, info[kk].(jj)+suffix
          endfor
       endfor
       printf, lun, '\\'
       if (mm lt nmassbins-2) then printf, lun, '\cline{1-5}'
    endfor

    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable*}'
    free_lun, lun

; ---------------------------------------------------------------------------    
; Table: samples and numbers of galaxies
    agesparent = read_mz_sample(/parent)
    agesemline = read_mz_sample(/mz_ispec)
    agesemline = agesemline[where(agesemline.oiii_5007[1] gt 0)]

    sdssparent = read_mz_sample(/parent,/sdss)
    sdssemline = read_mz_sample(/mz_ispec,/sdss)
    sdssemline = sdssemline[where(sdssemline.oiii_5007[1] gt 0)]
    
    agesmass = read_mz_sample(/mzhii_mass)
    agesancillary = read_mz_sample(/mzhii_ancillary)
    agesoh = read_mz_sample(/mzhii_log12oh)
    ageskk04 = mzlz_grab_info(agesoh,agesancillary,agesmass,/nolimit,/kk04)
    agest04 = mzlz_grab_info(agesoh,agesancillary,agesmass,/nolimit,/t04)
    agesm91 = mzlz_grab_info(agesoh,agesancillary,agesmass,/nolimit,/m91)
    
    sdssmass = read_mz_sample(/mzhii_mass,/sdss)
    sdssancillary = read_mz_sample(/mzhii_ancillary,/sdss)
    sdssoh = read_mz_sample(/mzhii_log12oh,/sdss)
    sdsskk04 = mzlz_grab_info(sdssoh,sdssancillary,sdssmass,/nolimit,/kk04)
    sdsst04 = mzlz_grab_info(sdssoh,sdssancillary,sdssmass,/nolimit,/t04)
    sdssm91 = mzlz_grab_info(sdssoh,sdssancillary,sdssmass,/nolimit,/m91)

; print some statistics for the paper    
; SDSS    
    nall = float(n_elements(sdssemline))
    nsf = total(sdssemline.agn eq 0)
    nagn = total(sdssemline.agn eq 1)
    splog, 'SDSS: Total star-forming, AGN'
    splog, nall, nsf, nsf/nall

; AGES
    nall = float(n_elements(agesemline))
    nsf = total(agesemline.agn eq 0)
    nagn = total(agesemline.agn eq 1)
    splog, 'AGESEMLINE: Total star-forming, AGN'
    splog, nall, nsf, nsf/nall

    splog, 'AGESEMLINE: Total classified as SF using the Yan diagram but *not* the BPT diagram'
    nyansf = total(agesemline.agn eq 0 and agesemline.bpt_agn eq -1 and agesemline.yan_agn eq 0)
    nyansf_hiz = total(agesemline.agn eq 0 and agesemline.bpt_agn eq -1 and $
      agesemline.yan_agn eq 0 and agesemline.z gt 0.4)    
    splog, nsf, nyansf, nyansf/nsf, nyansf_hiz, nyansf_hiz/nyansf

; now make the table    
    colhead1 = mzget_colhead(['Sample','$N$\tablenotemark{a}','Section\tablenotemark{b}'],/nobreak)
    colhead2 = '\cline{1-3} \\'
    texcenter = ['l','c','c']
    
    table = {sample: '', ngal: '', remark: ''}
    nages = 9 & nsdss = 9
    table = replicate(table,nages+nsdss)
    ntable = n_elements(table)
    ntags = n_tags(table)

    table.sample = [$
      'AGES',$
      '\hspace{0.2cm}Parent',$
      '\hspace{0.2cm}Emission-Line',$
      '\hspace{0.5cm}SF\tablenotemark{c}',$
      '\hspace{0.5cm}AGN',$
      '\hspace{0.2cm}Abundance\tablenotemark{d}',$
      '\hspace{0.5cm}KK04',$
      '\hspace{0.5cm}T04',$
      '\hspace{0.5cm}M91',$
      'SDSS',$
      '\hspace{0.2cm}Parent',$
      '\hspace{0.2cm}Emission-Line',$
      '\hspace{0.5cm}SF\tablenotemark{c}',$
      '\hspace{0.5cm}AGN',$
      '\hspace{0.2cm}Abundance\tablenotemark{d}',$
      '\hspace{0.5cm}KK04',$
      '\hspace{0.5cm}T04',$
      '\hspace{0.5cm}M91']

    table.ngal = [$
      '',$
      '$'+strtrim(n_elements(agesparent),2)+'$',$
      '$'+strtrim(n_elements(agesemline),2)+'$',$
      '$'+strtrim(long(total(agesemline.agn eq 0)),2)+'$',$
      '$'+strtrim(long(total(agesemline.agn eq 1)),2)+'$',$
      '',$
      '$'+strtrim(n_elements(ageskk04.oh),2)+'$',$
      '$'+strtrim(n_elements(agest04.oh),2)+'$',$
      '$'+strtrim(n_elements(agesm91.oh),2)+'$',$
      '',$
      '$'+strtrim(n_elements(sdssparent),2)+'$',$
      '$'+strtrim(n_elements(sdssemline),2)+'$',$
      '$'+strtrim(long(total(sdssemline.agn eq 0)),2)+'$',$
      '$'+strtrim(long(total(sdssemline.agn eq 1)),2)+'$',$
      '',$
      '$'+strtrim(n_elements(sdsskk04.oh),2)+'$',$
      '$'+strtrim(n_elements(sdsst04.oh),2)+'$',$
      '$'+strtrim(n_elements(sdssm91.oh),2)+'$']

    table.remark = [$
      '',$
      '\ref{sec:sample}',$
      '\ref{sec:selection}',$
      '\ref{sec:agn}',$
      '\ref{sec:agn}',$
      '',$
      '\ref{sec:oh}',$
      '\ref{sec:oh}',$
      '\ref{sec:oh}',$

      '',$
      '\ref{sec:sdss}',$
      '\ref{sec:selection}',$
      '\ref{sec:agn}',$
      '\ref{sec:agn}',$
      '',$
      '\ref{sec:oh}',$
      '\ref{sec:oh}',$
      '\ref{sec:oh}']

;   table.remark = [$
;     '',$
;     '$0.05<z<0.75$; $15.45<I<20.4$',$
;     '$\snr(\oii,\hb,\oiii)>4$; $\ewhb>2$~\AA',$
;     'BPT and \citet{yan10a} diagnostic diagrams; $\lesssim10\%$ may be AGN',$
;     'Identified using various optical, X-ray, mid-IR, and radio criteria',$
;     'Unambiguous; upper \pagel{} branch',$
;     '',$
;     '$0.033<z<0.25$; $14.5<r<17.6$',$
;     '$\snr(\oii,\hb,\oiii)>4$; $\ewhb>2$~\AA',$
;     'BPT diagram\tablenotemark{c}',$
;     'BPT diagram\tablenotemark{c}',$
;     'Unambiguous; upper \pagel{} branch']
    
;   tablecomments = ['AGES and SDSS sample definitions and sizes.']
    caption = 'AGES and SDSS Galaxy Samples\label{table:samples}'

    tablenotetext = [$
      '{a}{Number of galaxies in each sample.}',$
      '{b}{Section containing details regarding how we select each sample.}',$
      '{c}{Star-forming galaxies.}',$
      '{d}{Number of galaxies with well-measured oxygen abundances based on the '+$
      '\citet[M91]{mcgaugh91a}, \citet[T04]{tremonti04a}, and \citet[KK04]{kobulnicky04a} abundance calibrations.}'] 
      
;    tablenotetext = [$
;      '{a}{Number of galaxies in each sample.}',$
;      '{b}{Star-forming galaxies.}',$
;      '{c}{The BPT diagram \citep{baldwin81a} refers to the '+$
;      '\oiii/\hb{} vs.~\nii/\ha{} emission-line diagnostic diagram (see Fig.~\ref{fig:bpt}).}',$
;      '{d}{AGN in AGES were identified using various multiwavelength '+$
;      'criteria (see \S\ref{sec:agn}).}',$
;      '{e}{No classification possible.}',$
;      '{f}{Unclassifed SDSS galaxies are removed from subsequent analysis.}',$
;      '{g}{Number of galaxies with well-defined oxygen abundances (see \S\ref{sec:oh}).}']
      
; write out
    texfile = paperpath+'mztable_samples'+filesuffix+'.tex'
    splog, 'Writing '+texfile
    openw, lun, texfile, /get_lun
    if keyword_set(emulateapj) then begin
;      printf, lun, '\LongTables'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;   printf, lun, '\tabletypesize{\small}'
    printf, lun, '\tablecaption{'+caption+'}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    printf, lun, '}'
    printf, lun, '\startdata'
    for ii = 0, ntable-1 do begin
       line = strarr(ntags)
       for jj = 0L, ntags-1 do begin
          if (jj lt ntags-1) then suffix = ' & ' else if (ii lt ntable-1) then $
            suffix = ' \\ ' else suffix = ''
          line[jj] = table[ii].(jj)+suffix
       endfor
       printf, lun, line
       if (ii eq nages-1) then printf, lun, colhead2
    endfor
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    free_lun, lun

stop    
    
; ---------------------------------------------------------------------------    
; coefficients of the linear fits to the mean metallicity vs redshift
; plots (Fig 14) in bins of stellar mass and for each calibration 
    texfile = paperpath+'mztable_z_vs_oh_bymass_coeff'+filesuffix+'.tex'
    splog, 'Writing '+texfile
    
; read the fitting results
    mzavg = mrdfits(mzpath+'mzevol_avg.fits.gz',1,/silent)
    massbins = mz_massbins(nmassbins)

    calib = ['kk04','t04','m91']
    ncalib = n_elements(calib)

    mass = '\log\,(\mass/\msun)'
;   masslabel = ['\multicolumn{2}{c}{$'+mass+'>11$}','',$
;     '\multicolumn{2}{c}{$10.5<'+mass+'<11$}','',$
;     '\multicolumn{2}{c}{$10<'+mass+'<10.5$}','',$
;     '\multicolumn{2}{c}{$9.5<'+mass+'<10$}']
    masslabel = ['$'+mass+'>11$','$10.5<'+mass+'<11$','$10<'+mass+'<10.5$','$9.5<'+mass+'<10$']
    
    colhead1 = mzget_colhead(['','$\langle12+\log\,(\textrm{O}/\textrm{H})\rangle_{z=0.1}$',$
      '${\mathrm d}\log\,(\textrm{O}/\textrm{H})/{\mathrm d}z$','',$
      '$\langle12+\log\,(\textrm{O}/\textrm{H})\rangle_{z=0.1}$',$
      '${\mathrm d}\log\,(\textrm{O}/\textrm{H})/{\mathrm d}z$'])
    colhead2 = mzget_colhead(['Calibration','(dex)','(dex $z^{-1}$)','','(dex)','(dex $z^{-1}$)'],/nobreak)
    colhead3 = mzget_colhead(['','\multicolumn{2}{c}{'+masslabel[0]+'}','','\multicolumn{2}{c}{'+masslabel[1]+'}'])
    colhead4 = mzget_colhead(['','\multicolumn{2}{c}{'+masslabel[2]+'}','','\multicolumn{2}{c}{'+masslabel[3]+'}'])
    texcenter = replicate('c',6)
    
    tablenotetext = [$
      '{a}{The adopted linear model is given by: $\langle 12+\log\,(\textrm{O}/\textrm{H})\rangle = '+$
      '\langle 12+\log\,(\textrm{O}/\textrm{H})\rangle_{z=0.1} + {\mathrm d}\log\,(\textrm{O}/\textrm{H})/'+$
      '{\mathrm d}z\times(z-0.1)$, where $\langle12+\log\,(\textrm{O}/\textrm{H})\rangle_{z=0.1}$ is '+$
      'the mean metallicity at $z=0.1$ and ${\mathrm d}\log\,(\textrm{O}/\textrm{H})/{\mathrm d}z$ is the '+$
      'logarithmic rate of metallicity evolution.}']
    
    openw, lun, texfile, /get_lun
    printf, lun, '\begin{deluxetable*}{'+strjoin(texcenter)+'}[!h]'
    printf, lun, '\tablecaption{Mass-Dependent Evolution of the Mean Metallicity of '+$
      'Galaxies at $z=0.03-0.75$ \tablenotemark{a}\label{table:oh_bymass_coeff}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    printf, lun, '}'
    printf, lun, '\startdata'
    printf, lun, '\cline{1-6}'

; -------------------------
; M>11 and 10.5<M<11    
    printf, lun, colhead3
    printf, lun, '\cline{2-3}'
    printf, lun, '\cline{5-6}'

    for ii = 0, ncalib-1 do begin
       mzevol = mrdfits(mzpath+'mzevol_'+calib[ii]+'.fits.gz',1,/silent)
       printf, lun, strupcase(calib[ii])+' & '
; M>11
       printf, lun, $
         '$'+string(mzevol.coeffs_bymass[0,0],format='(F5.3)')+'\pm'+$
         string(mzevol.coeffs_bymass_err[0,0]>0.01,format='(F4.2)')+'$ & '+$
         '$'+string(mzevol.coeffs_bymass[1,0],format='(F7.3)')+'\pm'+$
         string(mzevol.coeffs_bymass_err[1,0],format='(F7.3)')+'$ & & '+$ ; space
; 10.5<M<11
         '$'+string(mzevol.coeffs_bymass[0,1],format='(F5.3)')+'\pm'+$
         string(mzevol.coeffs_bymass_err[0,1]>0.01,format='(F4.2)')+'$ & '+$
         '$'+string(mzevol.coeffs_bymass[1,1],format='(F7.3)')+'\pm'+$
         string(mzevol.coeffs_bymass_err[1,1],format='(F7.3)')+'$ \\ '
    endfor
    printf, lun, '\cline{1-6}'
    
; -------------------------
; 10<M<10.5 and 9.5<M<10
    printf, lun, colhead4
    printf, lun, '\cline{2-3}'
    printf, lun, '\cline{5-6}'

    for ii = 0, ncalib-1 do begin
       if (ii eq ncalib-1) then suffix = '' else suffix = '\\'
       mzevol = mrdfits(mzpath+'mzevol_'+calib[ii]+'.fits.gz',1,/silent)
       printf, lun, strupcase(calib[ii])+' & '
; 10<M<10.5
       printf, lun, $
         '$'+string(mzevol.coeffs_bymass[0,2],format='(F5.3)')+'\pm'+$
         string(mzevol.coeffs_bymass_err[0,2]>0.01,format='(F4.2)')+'$ & '+$
         '$'+string(mzevol.coeffs_bymass[1,2],format='(F7.3)')+'\pm'+$
         string(mzevol.coeffs_bymass_err[1,2],format='(F7.3)')+'$ & & '+$ ; space
; 9.5<M<10
         '$'+string(mzevol.coeffs_bymass[0,3],format='(F5.3)')+'\pm'+$
         string(mzevol.coeffs_bymass_err[0,3]>0.01,format='(F4.2)')+'$ & '+$
         '$'+string(mzevol.coeffs_bymass[1,3],format='(F7.3)')+'\pm'+$
         string(mzevol.coeffs_bymass_err[1,3],format='(F7.3)')+'$ '+suffix
    endfor

    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable*}'
    free_lun, lun

; ---------------------------------------------------------------------------    
; SDSS MZ/LZ relations
    texfile = paperpath+'mztable_mzlzlocal'+filesuffix+'.tex'
    splog, 'Writing '+texfile
    
; read the fitting results and gather the information we need based on
; both the reddening-corrected fluxes and the EWs
    fluxinfo = get_mzlzfit_info(/flux)
    ewinfo = get_mzlzfit_info()

    ncalib = n_elements(fluxinfo)
    ntags = n_tags(ewinfo)

    colhead1 = mzget_colhead(['','',$
      '\multicolumn{4}{c}{\mz\tablenotemark{b}}','',$
      '\multicolumn{3}{c}{\lz\tablenotemark{c}}'])
    colhead2 = '\cline{3-6}\cline{8-10}'
    colhead3 = mzget_colhead(['Calibration\tablenotemark{a}','',$
      '12+log\,\ohstar','$log\,(\mstar/10^{9}\, \msun$)','$\gamma$','$\sigma$\tablenotemark{d}','',$
      '$c_{0}$','$c_{1}$','$\sigma$\tablenotemark{d}'],/nobreak)
    texcenter = ['l','c',$
      replicate('c',4),'c',$
      replicate('c',3)]

    tablenotetext = [$
      '{a}{Results derived using the \citet[T04]{tremonti04a}, \citet[M91]{mcgaugh91a}, and '+$
      '\citet[KK04]{kobulnicky04a} calibrations of the \pagel{} parameter, using both '+$
      'reddening-corrected line-fluxes and equivalent widths (see \S\ref{sec:oh} for details). '+$
      'For comparison, we also derive the \mz{} and \lz{} relations using the oxygen abundances '+$
      'published by the MPA-JHU team (see \S\ref{sec:mzlzlocal}).}',$  
      '{b}{\mz{} relation given by '+$
      '$12+\log\,(\textrm{O}/\textrm{H}) = 12+\log\,(\textrm{O}/\textrm{H})^{\ast} -'+$
      '\log\, [1+(\mstar/10^{9}\mass)^{\gamma}]$.}',$
      '{c}{$B$-band \lz{} relation given by $\logoh=c_{0}+c_{1}(\mb+20.5)$.}',$
      '{d}{Residual scatter about the best-fitting relation.}']
    
    openw, lun, texfile, /get_lun
    printf, lun, '\begin{deluxetable*}{'+strjoin(texcenter)+'}[!h]'
    printf, lun, '\tablecaption{SDSS \mz{} and \lz{} Relations\label{table:mzlzlocal}}'   
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    niceprintf, lun, colhead3
    printf, lun, '}'
    printf, lun, '\startdata'

    printf, lun, '\multicolumn{2}{c}{} & \multicolumn{8}{c}{Reddening-Corrected Fluxes\tablenotemark{a}} \\'
    printf, lun, '\cline{1-10}'
    for ii = 0, ncalib-1 do begin
       for jj = 0, ntags-1 do begin
          if (jj eq ntags-1) then if (ii eq ncalib-1) then suffix = '' else $
            suffix = ' \\' else suffix = ' & '
          printf, lun, fluxinfo[ii].(jj)+suffix
       endfor
    endfor
    printf, lun, '\\'
    printf, lun, '\cline{1-10}'

    printf, lun, '\multicolumn{2}{c}{} & \multicolumn{8}{c}{Equivalent Widths\tablenotemark{a}} \\'
    printf, lun, '\cline{1-10}'
    for ii = 0, ncalib-2 do begin
       for jj = 0, ntags-1 do begin
          if (jj eq ntags-1) then if (ii eq ncalib-2) then suffix = '' else $
            suffix = ' \\' else suffix = ' & '
          printf, lun, ewinfo[ii].(jj)+suffix
       endfor
    endfor
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable*}'
    free_lun, lun

stop

stop
stop    
    
; ---------------------------------------------------------------------------    
; Table: delta-log(O/H) with redshift

; read the output from FIT_MZLZEVOL
    ohevol = mrdfits(mzpath+'ohevol.fits.gz',1)

    table = {z: '', dlogoh_q0: '', dlogoh_q0_cor: '', empty: '', $
      dlogoh_q15: '', dlogoh_q15_cor: ''}
    table = replicate(table,n_elements(ohevol))
;   table.z = '$'+string(ohevol.z,format='(F4.2)')+'$'
    table.z = '$'+['0.05-0.15','0.15-0.25','0.25-0.35','0.35-0.45','0.45-0.55','0.55-0.75']+'$'

; LZ    
    table1 = table
; Q=0
    pad = ['\phs','','','','','']
    table1.dlogoh_q0 = [pad+'$'+string(ohevol.ldlogoh_noevol,format='(F5.2)')+'\pm'+$
      string(ohevol.ldlogoh_noevol_err,format='(F5.2)')+'$']
    table1.dlogoh_q0_cor = [pad+'$'+string(ohevol.ldlogoh_noevol_cor,format='(F5.2)')+'\pm'+$
      string(ohevol.ldlogoh_noevol_cor_err,format='(F5.2)')+'$']
; Q=1.5
    table1.dlogoh_q15 = [pad+'$'+string(ohevol.ldlogoh_levol,format='(F5.2)')+'\pm'+$
      string(ohevol.ldlogoh_levol_err,format='(F5.2)')+'$']
    table1.dlogoh_q15_cor = [pad+'$'+string(ohevol.ldlogoh_levol_cor,format='(F5.2)')+'\pm'+$
      string(ohevol.ldlogoh_levol_cor_err,format='(F5.2)')+'$']

; MZ
    table2 = table
; Q=0
    pad = ['\phs','','','','','']
    table2.dlogoh_q0 = [pad+'$'+string(ohevol.mdlogoh,format='(F5.2)')+'\pm'+$
      string(ohevol.mdlogoh_err,format='(F5.2)')+'$']
    table2.dlogoh_q0_cor = [pad+'$'+string(ohevol.mdlogoh_noevol_cor,format='(F5.2)')+'\pm'+$
      string(ohevol.mdlogoh_noevol_cor_err,format='(F5.2)')+'$']
; Q=1.5
    table2.dlogoh_q15 = [pad+'$'+string(ohevol.mdlogoh,format='(F5.2)')+'\pm'+$
      string(ohevol.mdlogoh_err,format='(F5.2)')+'$']
    pad = ['\phs','','\phs','','','']
    table2.dlogoh_q15_cor = [pad+'$'+string(ohevol.mdlogoh_levol_cor,format='(F5.2)')+'\pm'+$
      string(ohevol.mdlogoh_levol_cor_err,format='(F5.2)')+'$']

    tags = tag_names(table1)
    ntable = n_elements(table1)
    ntags = n_tags(table1)
    
    colhead1 = $
      '\colhead{} & '+$
      '\multicolumn{2}{c}{\emph{L-Z}, $Q=0$} & '+$
      '\colhead{} & '+$
      '\multicolumn{2}{c}{\emph{L-Z}, $Q=1.5$} \\ '
    colhead2 = $
      '\colhead{} & '+$
      '\multicolumn{2}{c}{\emph{M-Z}, $Q=0$} & '+$
      '\colhead{} & '+$
      '\multicolumn{2}{c}{\emph{M-Z}, $Q=1.5$} \\ '
    colhead3 = '\cline{2-3}\cline{5-6}'
    colhead4 = $
      '\colhead{Redshift} & '+$
      '\colhead{$\langle \Delta\,\log\,({\rm O/H}) \rangle_{\rm obs}$} & '+$
      '\colhead{$\langle \Delta\,\log\,({\rm O/H}) \rangle_{\rm cor}$} & '+$
      '\colhead{} & '+$
      '\colhead{$\langle \Delta\,\log\,({\rm O/H}) \rangle_{\rm obs}$} & '+$
      '\colhead{$\langle \Delta\,\log\,({\rm O/H}) \rangle_{\rm cor}$} '
    colhead5 = '\cline{1-6}'
    texcenter = ['c','c','c','c','c','c']

    caption = 'Relative Change in Mean Oxygen Abundance '+$
      'with Redshift\label{table:dlogohevol}'
    
; write out
    texfile = paperpath+'mztable_dlogohevol_'+filesuffix+'.tex'
    splog, 'Writing '+texfile
    openw, lun, texfile, /get_lun
    if keyword_set(emulateapj) then begin
;      printf, lun, '\LongTables'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;   printf, lun, '\tabletypesize{\small}'
    printf, lun, '\tablecaption{'+caption+'}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead3
    niceprintf, lun, colhead4
    printf, lun, '}'
    printf, lun, '\startdata'
; LZ header and data
;   printf, lun, colhead1
;   printf, lun, colhead3
    for ii = 0, ntable-1 do begin
       line = strarr(ntags)
       for jj = 0L, ntags-1 do begin
          if (jj lt ntags-1) then suffix = ' & ' else suffix = ' \\ '
;         if (jj lt ntags-1) then suffix = ' & ' else if (ii lt ntable-1) then $
;           suffix = ' \\ ' else suffix = ''
          if strmatch(tags[jj],'*empty*',/fold) then line[jj] = suffix else $
            line[jj] = table1[ii].(jj)+suffix
       endfor
       printf, lun, line
    endfor
; MZ header and data
    printf, lun, colhead5
    printf, lun, colhead2
    printf, lun, colhead3
    printf, lun, colhead4+' \\'
    printf, lun, colhead5
    for ii = 0, ntable-1 do begin
       line = strarr(ntags)
       for jj = 0L, ntags-1 do begin
          if (jj lt ntags-1) then suffix = ' & ' else if (ii lt ntable-1) then $
            suffix = ' \\ ' else suffix = ''
          if strmatch(tags[jj],'*empty*',/fold) then line[jj] = suffix else $
            line[jj] = table2[ii].(jj)+suffix
       endfor
       printf, lun, line
    endfor
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
;   niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    free_lun, lun

stop
    
; ---------------------------------------------------------------------------    
; Table: local MZ relation

    ss = mrdfits(mzpath+'mzlocal_sdss_coeffs.fits.gz',1)
    aa = mrdfits(mzpath+'mzlocal_ages_coeffs.fits.gz',1)

    table = {sample: '', ngal: '', zrange: '', a1: '', b1: '', c1: '', scatter: ''}
    table = replicate(table,4)

    table.sample = ['AGES','SDSS','Tremonti+04\tablenotemark{d}','Adopted\tablenotemark{e}']
    table.ngal = ['$'+string([aa[1].ngal,ss[1].ngal,53400],format='(I0)')+'$','\nodata']

    table.zrange = ['$0.05<z<0.15$','$0.033<z<0.25$','$0.005<z<0.25$','$z\approx0.1$']
    table.a1 = [$
      '$'+string(aa[1].coeff[0],format='(F5.3)')+'\pm'+string(aa[1].coeff_err[0],format='(F5.3)')+'$',$
      '$'+string(ss[1].coeff[0],format='(F5.3)')+'\pm'+string(ss[1].coeff_err[0],format='(F5.3)')+'$',$
      '$'+string('9.10284',format='(F5.3)')+'$',$
      '$'+string(aa[2].coeff[0],format='(F5.3)')+'\pm'+string(aa[2].coeff_err[0],format='(F5.3)')+'$']
    table.b1 = [$
      '$'+string(aa[1].coeff[1],format='(F5.3)')+'\pm'+string(aa[1].coeff_err[1],format='(F5.3)')+'$',$
      '$'+string(ss[1].coeff[1],format='(F5.3)')+'\pm'+string(ss[1].coeff_err[1],format='(F5.3)')+'$',$
      '$'+string('0.16154',format='(F5.3)')+'$',$
      '$'+string(aa[2].coeff[1],format='(F5.3)')+'\pm'+string(aa[2].coeff_err[1],format='(F5.3)')+'$']
    table.c1 = [$
      '\nodata',$
      '$'+string(ss[1].coeff[2],format='(F7.4)')+'\pm'+string(ss[1].coeff_err[2],format='(F7.4)')+'$',$
      '$'+string('-0.08026',format='(F7.4)')+'$',$
      '$'+string(aa[2].coeff[2],format='(F7.4)')+'\pm'+string(aa[2].coeff_err[2],format='(F7.4)')+'$']

    table.scatter = ['$'+string([aa[1].scatter,ss[1].scatter,'0.10'],$
      format='(F4.2)')+'$','\nodata']
    struct_print, table

    ntable = n_elements(table)
    ntags = n_tags(table)
    
    colhead1 = $
      '\colhead{} & '+$
      '\colhead{} & '+$
      '\colhead{Redshift} & '+$
      '\multicolumn{3}{c}{Coefficients\tablenotemark{b}} & '+$
      '\colhead{} \\'
    colhead2 = '\cline{4-6}'
    colhead3 = $
      '\colhead{Sample} & '+$
      '\colhead{$N_{\rm gal}$\tablenotemark{a}} & '+$
      '\colhead{Range} & '+$
      '\colhead{$a$} & '+$
      '\colhead{$b$} & '+$
      '\colhead{$c$} & '+$
      '\colhead{Scatter\tablenotemark{c}}'
    texcenter = ['c','c','c','c','c','c','c']

    tablecomments = ['Each mass-metallicity relation was derived by computing '+$
      'the median oxygen abundance in $0.1$~dex wide bins of stellar mass and '+$
      'fitting a...']
    caption = 'Local Mass-Metallicity Relation\label{table:mzlocal}'

    tablenotetext = [$
      '{a}{Number of galaxies included in the fit.}',$
      '{b}{The data were fitted to a function of the form $12+\log\,({\rm O/H}) = a + b[\log\,(\mathcal{M})-10.5] + [\log\,(\mathcal{M})-10.5]^2$.}',$
      '{c}{$1\sigma$ dispersion in $\log\,({\rm O/H})$, relative to the polynomial fit.}',$
      '{d}{Mass-metallicity relation from \citet{tremonti04a}, shifted by $+0.05$~dex '+$
        'in $12+\log\,({\rm O/H})$ to match the zero-point of our adopted \citet{kobulnicky04a} '+$
        'abundance calibration.}',$
      '{e}{Adopted stellar mass-metallicity relation, taken to be the average '+$
      'of the SDSS and and \citet{tremonti04a} relations.}']
    
; write out
    texfile = paperpath+'mztable_mzlocal_'+filesuffix+'.tex'
    splog, 'Writing '+texfile
    openw, lun, texfile, /get_lun
    if keyword_set(emulateapj) then begin
;      printf, lun, '\LongTables'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;   printf, lun, '\tabletypesize{\small}'
    printf, lun, '\tablecaption{'+caption+'}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    niceprintf, lun, colhead3
    printf, lun, '}'
    printf, lun, '\startdata'
    for ii = 0, ntable-1 do begin
       line = strarr(ntags)
       for jj = 0L, ntags-1 do begin
          if (jj lt ntags-1) then suffix = ' & ' else if (ii lt ntable-1) then $
            suffix = ' \\ ' else suffix = ''
          line[jj] = table[ii].(jj)+suffix
       endfor
       printf, lun, line
    endfor
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    free_lun, lun

; ---------------------------------------------------------------------------    
; Table: local LZ relation

    ss = mrdfits(mzpath+'lzlocal_sdss_coeffs.fits.gz',1)
    aa = mrdfits(mzpath+'lzlocal_ages_coeffs.fits.gz',1)

    table = {sample: '', ngal: '', zrange: '', intercept: '', slope: '', scatter: ''}
    table = replicate(table,4)

    table.sample = ['AGES','SDSS','Tremonti+04\tablenotemark{d}','Adopted\tablenotemark{e}']
    table.ngal = ['$'+string([aa[1].ngal,ss[1].ngal,53400],format='(I0)')+'$','\nodata']

    table.zrange = ['$0.05<z<0.15$','$0.033<z<0.25$','$0.005<z<0.25$','$z\approx0.1$']
    table.slope = '$'+string([aa[1].coeff[1],ss[1].coeff[1],'-0.186',$
      aa[2].coeff[1]],format='(F6.3)')+'\pm'+string([aa[1].coeff_err[1],$
      ss[1].coeff_err[1],'0.001',aa[2].coeff_err[1]],format='(F5.3)')+'$'
    table.intercept = '$'+string([aa[1].coeff[0],ss[1].coeff[0],'9.058',$
      aa[2].coeff[0]],format='(F5.3)')+'\pm'+string([aa[1].coeff_err[0],$
      ss[1].coeff_err[0],'0.018',aa[2].coeff_err[0]],format='(F5.3)')+'$'
    table.scatter = ['$'+string([aa[1].scatter,ss[1].scatter,'0.16'],$
      format='(F4.2)')+'$','\nodata']
    struct_print, table

    ntable = n_elements(table)
    ntags = n_tags(table)
    
    colhead1 = $
      '\colhead{} & '+$
      '\colhead{} & '+$
      '\colhead{Redshift} & '+$
      '\multicolumn{2}{c}{Coefficients\tablenotemark{b}} & '+$
      '\colhead{} \\'
    colhead2 = '\cline{4-5}'
    colhead3 = $
      '\colhead{Sample} & '+$
      '\colhead{$N_{\rm gal}$\tablenotemark{a}} & '+$
      '\colhead{Range} & '+$
      '\colhead{$a$} & '+$
      '\colhead{$b$} & '+$
      '\colhead{Scatter\tablenotemark{c}}'
    texcenter = ['c','c','c','c','c','c']

    tablecomments = ['Each luminosity-metallicity relation was derived using an '+$
      'ordinary least-squares linear bisector fit to all objects with $-24.0<M_{0.1g}<-17.5$.']
    caption = 'Local $g$-band Luminosity-Metallicity Relation\label{table:glzlocal}'

    tablenotetext = [$
      '{a}{Number of galaxies included in the fit.}',$
      '{b}{The data were fitted to a function of the form $12+\log\,({\rm O/H}) = a + b(M_{0.1g}+20.5)$.}',$
      '{c}{$1\sigma$ dispersion in $\log\,({\rm O/H})$, relative to the best fit.}',$
      '{d}{$g$-band luminosity-metallicity relation from \citet{tremonti04a}, shifted by $+0.05$~dex '+$
        'in $12+\log\,({\rm O/H})$ to match the zero-point of our adopted \citet{kobulnicky04a} '+$
        'abundance calibration.}',$
      '{e}{Adopted $g$-band luminosity-metallicity relation, taken to be the average '+$
      'of the SDSS and AGES relations.}']
    
; write out
    texfile = paperpath+'mztable_glzlocal_'+filesuffix+'.tex'
    splog, 'Writing '+texfile
    openw, lun, texfile, /get_lun
    if keyword_set(emulateapj) then begin
;      printf, lun, '\LongTables'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;   printf, lun, '\tabletypesize{\small}'
    printf, lun, '\tablecaption{'+caption+'}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    niceprintf, lun, colhead3
    printf, lun, '}'
    printf, lun, '\startdata'
    for ii = 0, ntable-1 do begin
       line = strarr(ntags)
       for jj = 0L, ntags-1 do begin
          if (jj lt ntags-1) then suffix = ' & ' else if (ii lt ntable-1) then $
            suffix = ' \\ ' else suffix = ''
          line[jj] = table[ii].(jj)+suffix
       endfor
       printf, lun, line
    endfor
    printf, lun, '\enddata'
    printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    free_lun, lun

stop    
    
; ---------------------------------------------------------------------------    
; Table: measured properties

    ngalaxy = n_elements(agesdust)
    srt = sort(agesancillary.ages_id)

    ages = replicate({id: '', ra: '', dec: '', z: '', weight: '', $
      m_b: '', m_v: '', m_r: '', mass: '', $
      oii: '', hbeta: '', oiii: ''},ngalaxy)
    ntags = n_tags(ages)

    ages.id     = string(agesancillary[srt].ages_id,format='(I0)')
    ages.ra     = dec2hms(agesancillary[srt].ra/15D,/colon)
    ages.dec    = dec2hms(agesancillary[srt].dec,/colon)
    pos = where(agesancillary[srt].dec gt 0.0) & ages[pos].dec = '+'+ages[pos].dec

    ages.z      = string(agesancillary[srt].z,format='(F6.4)')
    ages.weight = string(agesancillary[srt].spec_weight,format='(F8.3)')

    ages.m_b  = '$'+string(agesancillary[srt].m_b,format='(F6.2)')+'\pm'+string(agesancillary[srt].m_b_err,format='(F4.2)')+'$'
    ages.m_v  = '$'+string(agesancillary[srt].m_v,format='(F6.2)')+'\pm'+string(agesancillary[srt].m_v_err,format='(F4.2)')+'$'
    ages.m_r  = '$'+string(agesancillary[srt].m_r,format='(F6.2)')+'\pm'+string(agesancillary[srt].m_r_err,format='(F4.2)')+'$'
    ages.mass = '$'+string(agesancillary[srt].kcorr_mass,format='(F5.2)')+'\pm'+string(agesancillary[srt].kcorr_mass_err>0.01,format='(F4.2)')+'$'
    
    format_flux_error, agesdust[srt].oii_3727_ew[0],  agesdust[srt].oii_3727_ew[1],  f_oii,  e_oii  & ages.oii   = '$'+f_oii +'\pm'+e_oii+'$'
    format_flux_error, agesdust[srt].oiii_5007_ew[0], agesdust[srt].oiii_5007_ew[1], f_oiii, e_oiii & ages.oiii  = '$'+f_oiii+'\pm'+e_oiii+'$'
    format_flux_error, agesdust[srt].h_beta_ew[0],    agesdust[srt].h_beta_ew[1],    f_hb,   e_hb   & ages.hbeta = '$'+f_hb  +'\pm'+e_hb+'$'

; replace negative [O III] values with upper limits

    neg = where(agesdust[srt].oiii_5007_ew[0] lt 0.0,nneg)
    if (nneg ne 0L) then ages[neg].oiii = '$>'+string(agesdust[srt[neg]].oiii_5007_ew_limit,format='(F4.1)')+'$' ; make sure the limits are positive!
;   format_flux_error, abs(agesdust[srt[neg]].oiii_5007_ew[0]), agesdust[srt[neg]].oiii_5007_ew[1], f, e
    
    colhead = strjoin('\colhead{'+['ID\tablenotemark{a}','$\alpha_{\rm J2000}$\tablenotemark{a}','$\delta_{\rm J2000}$\tablenotemark{a}',$
      'Redshift\tablenotemark{a}','Weight\tablenotemark{a}',$
      '$M_{B}$\tablenotemark{b}','$M_{V}$\tablenotemark{b}','$M_{R}$\tablenotemark{b}','$\log\,(M/M_{\sun})$\tablenotemark{b}',$
      'EW([O~{\sc ii}~$\lambda3727$])\tablenotemark{c}','EW(H$\beta$)\tablenotemark{c}','EW([O~{\sc iii}~$\lambda5007$])\tablenotemark{c}']+'}',' & ')+' \\'
    colunits = strjoin('\colhead{('+string(lindgen(ntags)+1L,format='(I0)')+')}',' & ') ;+' \\'
    colcenter = replicate('c',ntags)

    tablecomments = ['(1)','(2)']
    tablecomments = strjoin(tablecomments,' ')
    tablenotetext = [$
      '{a}{AGES identification numbers, celestial coordinates, spectroscopic redshifts, and statistical weights are '+$
      'taken from C.~S.~Kochanek et~al. (2007, in prep.) and D.~J. Eisenstein et~al. (2007, in prep.).}',$
      '{b}{Absolute magnitudes (Vega; $h=0.7$) and stellar masses (\citealt{salpeter55} IMF; $0.1-100~M_{\sun}$) have '+$
      'been derived using {\sc k-correct} \citep[ver.~4.1.3;][]{blanton06b}.}',$
      '{c}{Absorption-corrected, rest-frame emission-line equivalent widths in \AA.}']
    
    texfile = paperpath+'mztable_ages_mzdata_'+filesuffix+'.tex'
    splog, 'Writing '+texfile+'.'
    openw, lun, texfile, /get_lun
    printf, lun, '\begin{landscape}'
    printf, lun, '\begin{deluxetable}{'+strjoin(colcenter)+'}'
    printf, lun, '\tabletypesize{\tiny}'
;   printf, lun, '\rotate'
    printf, lun, '\tablecaption{AGES: Spectrophotometric Properties and Oxygen Abundances\label{table:ages_mzdata}}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead
    niceprintf, lun, colunits
    printf, lun, '}'
    printf, lun, '\startdata'

    for i = 0L, 24L do begin
;   for i = 0L, ngalaxy-1L do begin
       line = strarr(ntags)
       for j = 0L, ntags-1L do begin
          if (j lt ntags-1L) then suffix = ' & ' else begin
             if (i lt ngalaxy-1L) then suffix = ' \\ ' else suffix = ''
          endelse
          if (strcompress(ages[i].(j),/remove) eq '') then $
            line[j] = '\nodata'+suffix else line[j] = ages[i].(j)+suffix
       endfor
       printf, lun, line
    endfor
    
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
    niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    printf, lun, '\clearpage'
    printf, lun, '\end{landscape}'
    free_lun, lun

stop    

; ---------------------------------------------------------------------------    
; Table: delta-log(O/H) with redshift

; read the output from FIT_MZLZEVOL
    ohevol = mrdfits(mzpath+'ohevol.fits.gz',1)

    table = {z: '', dlogoh: '', dlogoh_cor: ''}
    table = replicate(table,n_elements(ohevol))
    table.z = '$'+string(ohevol.z,format='(F4.2)')+'$'

; L-Z, Q=0    
    table1 = table 
    table1.dlogoh = ['$'+string(ohevol.ldlogoh_noevol,format='(F6.2)')+'\pm'+$
      string(ohevol.ldlogoh_noevol_err,format='(F6.2)')+'$']
    table1.dlogoh_cor = ['$'+string(ohevol.ldlogoh_noevol_cor,format='(F6.2)')+'\pm'+$
      string(ohevol.ldlogoh_noevol_cor_err,format='(F6.2)')+'$']
; L-Z, Q=1.5
    table2 = table 
    table2.dlogoh = ['$'+string(ohevol.ldlogoh_levol,format='(F6.2)')+'\pm'+$
      string(ohevol.ldlogoh_levol_err,format='(F6.2)')+'$']
    table2.dlogoh_cor = ['$'+string(ohevol.ldlogoh_levol_cor,format='(F6.2)')+'\pm'+$
      string(ohevol.ldlogoh_levol_cor_err,format='(F6.2)')+'$']

    ntable = n_elements(table)
    ntags = n_tags(table)
    
    colhead1 = $
      '\colhead{} & '+$
      '\colhead{Observed} & '+$
      '\colhead{Corrected} \\'
    colhead2 = $
      '\colhead{Redshift} & '+$
      '\colhead{$\langle \Delta\,\log\,({\rm O/H}) \rangle$} & '+$
      '\colhead{$\langle \Delta\,\log\,({\rm O/H}) \rangle$} \\'
    texcenter = ['c','c','c']

    caption = 'Change in the Mean Oxygen Abundance\label{table:dlogohevol}'
    
; write out
    texfile = paperpath+'mztable_dlogohevol_'+filesuffix+'.tex'
    splog, 'Writing '+texfile
    openw, lun, texfile, /get_lun
    if keyword_set(emulateapj) then begin
;      printf, lun, '\LongTables'
    endif
    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;   printf, lun, '\tabletypesize{\small}'
    printf, lun, '\tablecaption{'+caption+'}'
    printf, lun, '\tablewidth{0pt}'
    printf, lun, '\tablehead{'
    niceprintf, lun, colhead1
    niceprintf, lun, colhead2
    printf, lun, '}'
    printf, lun, '\startdata'
; L-Z, Q=0    
    printf, lun, '\cutinhead{\emph{L-Z}, $Q=0$}'
    for ii = 0, ntable-1 do begin
       line = strarr(ntags)
       for jj = 0L, ntags-1 do begin
          if (jj lt ntags-1) then suffix = ' & ' else suffix = ' \\ '
;         if (jj lt ntags-1) then suffix = ' & ' else if (ii lt ntable-1) then $
;           suffix = ' \\ ' else suffix = ''
          line[jj] = table1[ii].(jj)+suffix
       endfor
       printf, lun, line
    endfor
; L-Z, Q=1.5
    printf, lun, '\cutinhead{\emph{L-Z}, $Q=1.5$}'
    for ii = 0, ntable-1 do begin
       line = strarr(ntags)
       for jj = 0L, ntags-1 do begin
          if (jj lt ntags-1) then suffix = ' & ' else if (ii lt ntable-1) then $
            suffix = ' \\ ' else suffix = ''
          line[jj] = table2[ii].(jj)+suffix
       endfor
       printf, lun, line
    endfor
    printf, lun, '\enddata'
;   printf, lun, '\tablecomments{'+tablecomments+'}'
;   niceprintf, lun, '\tablenotetext'+tablenotetext
    printf, lun, '\end{deluxetable}'
    free_lun, lun

stop
    
return
end


;;; ---------------------------------------------------------------------------    
;;; SDSS MZ relation
;;    texfile = paperpath+'mztable_mzlocal'+filesuffix+'.tex'
;;    splog, 'Writing '+texfile
;;    
;;; read the fitting results and gather the information we need based on
;;; both the reddening-corrected fluxes and the EWs
;;    info = get_mzfit_info(ncalib=ncalib) 
;;    ntags = n_tags(info)
;;
;;;   tablenotetext = [$
;;;     '{a}{Fits were performed using a ``robust'' ordinaryleast-squares bisector fit.}',$
;;;     '{b}{$I$-band magnitudes are on the Vega system, while $g$- and $r$-band magnitudes are AB-relative.}']
;;    
;;    colhead1 = mzget_colhead(['Calibration','',$
;;      '\multicolumn{4}{c}{Polynomial}','',$
;;      '\multicolumn{5}{c}{Double Power-law}','',$
;;      '\multicolumn{5}{c}{Broken Power-law}'])
;;    colhead2 = '\cline{3-6}\cline{8-12}\cline{14-18}'
;;    colhead3 = mzget_colhead(['','',$
;;      '$a_{0}$','$a_{1}$','$a_{2}$','$\sigma$','',$
;;      '\ohstar','\mstar','$\beta_{1}$','$\beta_{2}$','$\sigma$','',$
;;      '\ohstar','\mstar','$\gamma_{1}$','$\gamma_{2}$','$\sigma$'],/nobreak)
;;    texcenter = ['l','c',$
;;      replicate('c',4),'c',$
;;      replicate('c',5),'c',$
;;      replicate('c',5)]
;;
;;    openw, lun, texfile, /get_lun
;;    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;;    printf, lun, '\tablecaption{SDSS \mz{} Relation\label{table:mzlocal}}' 
;;    printf, lun, '\tablewidth{0pt}'
;;    printf, lun, '\tablehead{'
;;    niceprintf, lun, colhead1
;;    niceprintf, lun, colhead2
;;    niceprintf, lun, colhead3
;;    printf, lun, '}'
;;    printf, lun, '\startdata'
;;    printf, lun, '\multicolumn{2}{c}{} & \multicolumn{15}{c}{Reddening-Corrected Fluxes} \\'
;;    printf, lun, '\cline{3-18}'
;;    for ii = 0, ncalib-1 do begin
;;       for jj = 0, ntags-1 do begin
;;          if (jj eq ntags-1) then if (ii eq ncalib-1) then suffix = '' else $
;;            suffix = ' \\' else suffix = ' & '
;;          printf, lun, info[ii].(jj)+suffix
;;       endfor
;;    endfor
;;    
;;    printf, lun, '\enddata'
;;;   printf, lun, '\tablecomments{'+tablecomments+'}'
;;;   niceprintf, lun, '\tablenotetext'+tablenotetext
;;    printf, lun, '\end{deluxetable}'
;;    free_lun, lun
;;
;;function get_mzfit_info, fluxcor=fluxcor, ncalib=ncalib
;;    mzpath = ages_path(/projects)+'mz/'
;;
;;    if keyword_set(fluxcor) then begin
;;       polyfit = mrdfits(mzpath+'mzlocal_sdss_fluxcor_poly.fits.gz',1)
;;       doublepl = mrdfits(mzpath+'mzlocal_sdss_fluxcor_doublepl.fits.gz',1)
;;       brokenpl = mrdfits(mzpath+'mzlocal_sdss_fluxcor_brokenpl.fits.gz',1)
;;    endif else begin
;;       polyfit = mrdfits(mzpath+'mzlocal_sdss_poly.fits.gz',1)
;;       doublepl = mrdfits(mzpath+'mzlocal_sdss_doublepl.fits.gz',1)
;;       brokenpl = mrdfits(mzpath+'mzlocal_sdss_brokenpl.fits.gz',1)
;;    endelse
;;
;;    ncalib = 3
;;    info = replicate({calib: '', empty1: ' ', $
;;      a0: '', a1: '', a2: '', pscatter: '', empty2: ' ', $
;;      ohstar1: '', mstar1: '', beta1: '', beta2: '', dscatter: '', empty3: ' ', $
;;      ohstar2: '', mstar2: '', gamma1: '', gamma2: '', bscatter: ''},ncalib)
;;    info.calib = strupcase(strtrim(polyfit.calib,2))
;;; polynomial
;;    info.a0       = im_sigfigs(polyfit.coeff_bin[0],3)
;;    info.a1       = im_sigfigs(polyfit.coeff_bin[1],4)
;;    info.a2       = im_sigfigs(polyfit.coeff_bin[2],4)
;;    info.pscatter = string(polyfit.scatter_bin,format='(F5.3)')
;;; doublepl
;;    info.ohstar1  = im_sigfigs(doublepl.coeff_bin[0],4)
;;    info.mstar1   = im_sigfigs(doublepl.coeff_bin[1],4)
;;    info.beta1    = im_sigfigs(doublepl.coeff_bin[2],4)
;;    info.beta2    = im_sigfigs(doublepl.coeff_bin[3],4)
;;    info.dscatter = string(doublepl.scatter_bin,format='(F5.3)')
;;; brokenpl
;;    info.ohstar2  = im_sigfigs(brokenpl.coeff_bin[0],4)
;;    info.mstar2   = im_sigfigs(brokenpl.coeff_bin[1],4)
;;    info.gamma1   = im_sigfigs(brokenpl.coeff_bin[2],4)
;;    info.gamma2   = im_sigfigs(brokenpl.coeff_bin[3],4)
;;    info.bscatter = string(brokenpl.scatter_bin,format='(F5.3)')
;;
;;    struct_print, info
;;
;;return, info
;;end
;;




;;
;;; ---------------------------------------------------------------------------    
;;; Table: samples and numbers of galaxies
;;
;;    table = {sample: '', ngal: '', remark: ''}
;;    nages = 8 & nsdss = 7
;;    table = replicate(table,nages+nsdss)
;;    ntable = n_elements(table)
;;    ntags = n_tags(table)
;;
;;    table.sample = [$
;;      'AGES',$
;;      '\hspace{0.2cm}Parent',$
;;      '\hspace{0.2cm}Emission-Line',$
;;      '\hspace{0.5cm}SF\tablenotemark{b}',$
;;      '\hspace{0.5cm}AGN\tablenotemark{d}',$
;;      '\hspace{0.5cm}Unclassified\tablenotemark{e}',$
;;      '\hspace{0.5cm}SF+Unclassified',$
;;      '\hspace{0.2cm}Abundance\tablenotemark{g}',$
;;      'SDSS',$
;;      '\hspace{0.2cm}Parent',$
;;      '\hspace{0.2cm}Emission-Line',$
;;      '\hspace{0.5cm}SF\tablenotemark{b}',$
;;      '\hspace{0.5cm}AGN',$
;;      '\hspace{0.5cm}Unclassified\tablenotemark{f}',$
;;      '\hspace{0.2cm}Abundance\tablenotemark{g}']
;;
;;    table.ngal = [$
;;      '',$
;;      '$11121$',$
;;      '$4151$',$
;;      '$2347$',$
;;      '$884$',$
;;      '$920$',$
;;      '$3267$',$
;;      '$3047$',$
;;      '',$
;;      '$489337$',$
;;      '$119834$',$
;;      '$94619$',$
;;      '$25215$',$
;;      '$848$',$
;;      '$93516$']
;;
;;    table.remark = [$
;;      '',$
;;      '$0.05<z<0.75$; $14<I_{\rm Vega}<19.95$',$
;;      '$\snr(\oii,\hb,\oiii)>4$; $\ewhb>2$~\AA',$
;;      'BPT diagram\tablenotemark{c}',$
;;      'Optical, X-ray, radio, mid-IR',$
;;      'Assumed to be SF; $\lesssim20\%$ may be AGN',$
;;      'Union of SF and Unclassified subsamples',$
;;      'Unambiguous; upper \pagel{} branch',$
;;      '',$
;;      '$0.033<z<0.25$; $14.5<r_{\rm AB}<17.6$',$
;;      '$\snr(\oii,\hb,\oiii)>4$; $\ewhb>2$~\AA',$
;;      'BPT diagram\tablenotemark{c}',$
;;      'BPT diagram\tablenotemark{c}',$
;;      '$\snr(\nii)<4$',$
;;      'Unambiguous; upper \pagel{} branch']
;;    
;;    colhead1 = $
;;      '\colhead{Sample} & '+$
;;      '\colhead{$N_{\rm gal}$\tablenotemark{a}} & '+$
;;      '\colhead{Remarks}'
;;    colhead2 = '\cline{1-3}'
;;    texcenter = ['l','c','c']
;; 
;;;   tablecomments = ['AGES and SDSS sample definitions and sizes.']
;;    caption = 'AGES and SDSS Sample Definitions and Sizes\label{table:samples}'
;;
;;    tablenotetext = [$
;;      '{a}{Number of galaxies in each sample.}',$
;;      '{b}{Star-forming galaxies.}',$
;;      '{c}{The BPT diagram \citep{baldwin81a} refers to the '+$
;;      '\oiii/\hb{} vs.~\nii/\ha{} emission-line diagnostic diagram (see Fig.~\ref{fig:bpt}).}',$
;;      '{d}{AGN in AGES were identified using various multiwavelength '+$
;;      'criteria (see \S\ref{sec:agn}).}',$
;;      '{e}{No classification possible.}',$
;;      '{f}{Unclassifed SDSS galaxies are removed from subsequent analysis.}',$
;;      '{g}{Number of galaxies with well-defined oxygen abundances (see \S\ref{sec:oh}).}']
;;      
;;; write out
;;    texfile = paperpath+'mztable_samples_'+filesuffix+'.tex'
;;    splog, 'Writing '+texfile
;;    openw, lun, texfile, /get_lun
;;    if keyword_set(emulateapj) then begin
;;;      printf, lun, '\LongTables'
;;    endif
;;    printf, lun, '\begin{deluxetable}{'+strjoin(texcenter)+'}'
;;;   printf, lun, '\tabletypesize{\small}'
;;    printf, lun, '\tablecaption{'+caption+'}'
;;    printf, lun, '\tablewidth{0pt}'
;;    printf, lun, '\tablehead{'
;;    niceprintf, lun, colhead1
;;    printf, lun, '}'
;;    printf, lun, '\startdata'
;;    for ii = 0, ntable-1 do begin
;;       line = strarr(ntags)
;;       for jj = 0L, ntags-1 do begin
;;          if (jj lt ntags-1) then suffix = ' & ' else if (ii lt ntable-1) then $
;;            suffix = ' \\ ' else suffix = ''
;;          line[jj] = table[ii].(jj)+suffix
;;       endfor
;;       printf, lun, line
;;       if (ii eq nages-1) then printf, lun, colhead2
;;    endfor
;;    printf, lun, '\enddata'
;;;   printf, lun, '\tablecomments{'+tablecomments+'}'
;;    niceprintf, lun, '\tablenotetext'+tablenotetext
;;    printf, lun, '\end{deluxetable}'
;;    free_lun, lun
;;
;;stop
;;
