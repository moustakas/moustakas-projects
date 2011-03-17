
pro fig

tstr = 'Velocity field at cz = 500 km s^{-1} and 279 nearby galaxies'

ps_open, 'trgb.fig1', thick=4, /ps_fonts, /encapsulated

plot_smooth, 'setup.dat', 'vfield_iras.dat', $
  vshell=500, nc=17, ti=textoidl(tstr), /no_slabel, /no_zlabel

plot_neighbors, 'neighbor.lis', colu=[1,2,5], /nocolor

ps_close

end
