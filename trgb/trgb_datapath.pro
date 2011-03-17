function trgb_datapath
; jm00aug5ucb
; datapaths to the TRGB data

	datapaths = strarr(5)
        datapaths[0] = '/deepscr1/ioannis/trgb/'	; HST data
        datapaths[1] = '/deepscr1/ioannis/archive/'	; HST archival data
        datapaths[2] = '/deep3/marc/trgb/data/'		; Keck data
        datapaths[3] = '/deep1/ioannis/trgb/plots/'	; plot subdirectory
        datapaths[4] = '/deep1/ioannis/trgb/results/'	; results subdirectory

return, datapaths
end
