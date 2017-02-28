for(i in 1:10)
  assign(nam[i], value = matrix(df[,(6+i)], nrow = 985, ncol = 822))

disc = makeBrush(31, "disc")
disc = disc / sum(disc)
offset = 0.05
image_bg = filter2( X191Ir.Ir191Di., disc)
image_th = X191Ir.Ir191Di. > image_bg + offset
display(image_th)


nmask = watershed(distmap(X191Ir.Ir191Di.), 2)

nmask = thresh(X191Ir.Ir191Di., w=3,h=3, offset = 0.05)

nmask = opening(nmask, makeBrush(5, shape = 'disc'))

nmask = fillHull(nmask)

nmask = bwlabel(nmask)
display(colorLabels(nmask))