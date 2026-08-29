#!/usr/bin/env Rscript
## Re-assemble Figure 1 (r2): drop the old panel A (pipeline schematic; the old
## composed file figures_final/"Figure 1-r1.pdf" is kept untouched as backup),
## promote the body map to A and the mirrored bar panel to B, and add the new
## cross-sex fertility volcano (fig1_panelC_volcano.R) as C in the space freed
## on the left.  The r1 figure was composed outside R (no panel sources on
## disk), so panels A/B are cut from the r1 page as 600-dpi raster pieces along
## seams verified pixel-clean: body|bars at x=2019 (y>=1319), schematic|bars at
## x=2220 (y<1319); all coordinates below are 600-dpi pixels of the r1 render.
suppressMessages(library(magick))
D   <- "${PROJ}/figures_final"
r1  <- image_read_pdf(file.path(D, "Figure 1-r1.pdf"), density = 600)
stopifnot(abs(image_info(r1)$width - 4376) < 10)
W <- image_info(r1)$width; H <- image_info(r1)$height

geom <- function(w, h, x, y) sprintf("%dx%d+%d+%d", w, h, x, y)

## The body-map labels and the bar panel's longest disease labels interleave
## around x~2000, so no vertical seam separates the two panels cleanly.  Their
## INK colours do: bar-panel label text is dark achromatic grey (~97), body-map
## labels are either light grey (~181) or saturated colours.  Inside the
## contested band (x 1720-2620) the body-side copy drops dark achromatic ink
## and the bars-side copy drops everything EXCEPT dark achromatic ink; the two
## copies are then pasted with white transparent so neither hides the other.
suppressMessages(library(png))
X0 <- 1720L; X1 <- 2620L; BW <- X1 - X0 + 1L
band <- image_crop(r1, geom(BW, H, X0, 0))
bA <- as.integer(image_data(band))
if (dim(bA)[1] == BW) bA <- aperm(bA, c(2, 1, 3))         # -> [H, BW, 3]
stopifnot(dim(bA)[1] == H, dim(bA)[2] == BW)
mx <- pmax(bA[,,1], bA[,,2], bA[,,3]); mn <- pmin(bA[,,1], bA[,,2], bA[,,3])
chroma <- mx - mn; gray <- (bA[,,1] + bA[,,2] + bA[,,3]) %/% 3L
bL <- bA; bR <- bA
kill <- function(arr, m) { for (ch in 1:3) { p <- arr[,,ch]; p[m] <- 255L; arr[,,ch] <- p }; arr }
mL <- chroma < 45 & gray < 170                            # drop bars text from body side
mR <- chroma >= 45 | gray > 150                           # keep only dark achromatic ink
mL[4184:H, ] <- FALSE                                     # legend rows: body side keeps
mR[4184:H, ] <- TRUE                                      #   all, bars side keeps none
bL <- kill(bL, mL)
bR <- kill(bR, mR)
tmp <- tempfile(fileext = ".png")
mkim <- function(arr) { writePNG(arr / 255, tmp); image_read(tmp) }
imL <- image_composite(r1, mkim(bL), offset = sprintf("+%d+0", X0))
imR <- image_composite(r1, mkim(bR), offset = sprintf("+%d+0", X0))

body   <- image_crop(imL, geom(2050, 2784, 0, 1400))      # body map, old panel B
legend <- image_crop(imL, geom(2560,  225, 0, 4184))      # its category legend
barsT  <- image_crop(r1,  geom(W - 2220, 1318, 2220, 0))  # above body top: plain cut
barsB  <- image_crop(imR, geom(W - X0, H - 1319, X0, 1319))

volc <- image_read_pdf(file.path(D, "Fig1_panelC_volcano_small.pdf"), density = 600)
volc <- image_resize(volc, "1725x")                       # right edge clear of the
                                                          # IPF label (starts x 1832)
cv <- image_blank(W, H, "white")
cv <- image_composite(cv, barsB,  offset = sprintf("+%d+1319", X0))
cv <- image_composite(cv, barsT,  offset = "+2220+0")
cv <- image_composite(cv, image_transparent(body,   "white"), offset = "+0+60")
cv <- image_composite(cv, image_transparent(legend, "white"), offset = "+0+2905")
cv <- image_composite(cv, volc,   offset = "+100+3190")
## erase the baked-in letters (old B at 60-92,136-180 after shift; old C at
## 2256-2340,64-160) and write the new ones
cv <- image_draw(cv)
rect(40, 110, 130, 200, col = "white", border = NA)      # old "B" -> A
rect(2240, 40, 2350, 175, col = "white", border = NA)    # old "C" -> B
text(60,   180, "A", cex = 5, font = 2, family = "Helvetica", adj = c(0, 1))
text(2256, 180, "B", cex = 5, font = 2, family = "Helvetica", adj = c(0, 1))
text(24,  3300, "C", cex = 5, font = 2, family = "Helvetica", adj = c(0, 1))
dev.off()

image_write(cv, file.path(D, "Figure1_r2.png"))
image_write(cv, file.path(D, "Figure1_r2.pdf"), format = "pdf", density = 600)
cat("wrote Figure1_r2.pdf/.png (", W, "x", H, "px @600dpi )\n")
