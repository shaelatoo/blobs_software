pro radon_matrix,coords_x,coords_y,thetas,rs,L,rad

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Calculates the radon transform matrix for the     ;;
;;          given parameters and coordinate spacing described ;;
;;          by the arrays image_coords_x and image_coords_y.  ;;
;;                                                            ;;
;; Inputs: coords_x,coords_y - m,n-length vectors ;;
;;           giving the x,y coordinate positions associated   ;;
;;           with each pixel. For height-time images, y is    ;;
;;           altitude and x is time.                          ;;
;;         thetas,rs - vectors giving the theta,r positions   ;;
;;           desired in transform space                       ;;
;;                                                            ;;
;; Outputs: L - array containing the number of points inter-  ;;
;;            sected by each line                             ;;
;;          rad - array containing the "transform operator" - ;;
;;            see notes at end for interpretation             ;;
;;                                                            ;; 
;; Created: 11/15/07                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; parameters


; input checks
szcoords_x=SIZE(coords_x)
szcoords_y=SIZE(coords_y)
if szcoords_x[0] ne 1 then begin
  print,'Variable coords_x must be a 1D array'
  return
endif
if szcoords_y[0] ne 1 then begin
  print,'Variable coords_y must be a 1D array'
  return
endif
szthetas=SIZE(thetas)
szrs=SIZE(rs)
if szthetas[1] ne 1 then begin
  print,'Variable thetas must be a 1D array'
  return
endif
if szrs[1] ne 1 then begin
  print,'Variable rs must be a 1D array'
  return
endif
avex=MEAN(coords_x)
if (avex/MIN(rs) gt 10) or (avex/MAX(rs) le .01) then begin
  print,'Problem with size of rs'
  return
endif
avey=MEAN(coords_y)
if (avey/MIN(rs) gt 10) or (avey/MAX(rs) le .01) then begin
  print,'Problem with size of rs'
  return
endif
locs=WHERE(coords_x[1:*]-coords_x le 0)
if locs[0] ne -1 then begin
  print,'Problem with x coordinate vector'
  return
endif
locs=WHERE(coords_y[1:*]-coords_y le 0)
if locs[0] ne -1 then begin
  print,'Problem with y coordinate vector'
  return
endif



; initializations
nx=szcoords_x[1]
ny=szcoords_y[1]
nthet=szthetas[1]
nr=szrs[1]
if nx gt ny[1] then nts=nx else nts=ny
if nts MOD 2 ne 0 then nts++
L=INTARR(nthet*nr)
rad=FLTARR(nthet*nr,nx*ny)
t=INTARR(nts+1)-nts/2
minx=MIN(coords_x,max=maxx)
miny=MIN(coords_y,max=maxy)



; calc matrix elements
for thet=0,nthet-1 do begin
  cthet=COS(thetas[thet])
  sthet=SIN(thetas[thet])
  for r=0,nr-1 do begin
    nL=0
    r0=rs[r]
    row=thet*nr+r
    x=-cthet*r0+sthet*t
    y=-sthet*r0-cthet*t
    for point=0,nts-1 do begin
      xp=x[point]
      yp=y[point]
      if (xp le maxx) and (xp ge minx) then begin
        if (yp le maxy) and (yp ge miny) then begin
          garbage=MIN(ABS(xp-coords_x),locx)
          if xp gt coords_x[locx] then begin
            x2=coords_x[locx+1]
            x1=coords_x[locx]
            locx2=locx+1
          endif else begin
            locx2=locx
            locx--
            x2=coords_x[locx2]
            x1=coords_x[locx]
          endelse
          garbage=MIN(ABS(yp-coords_y),locy)
          if yp gt coords_y[locy] then begin
            locy2=locy+1
            y2=coords_y[locy2]
            y1=coords_y[locy]
          endif else begin
            locy2=locy
            locy--
            y2=coords_y[locy2]
            y1=coords_y[locy]
          endelse
          mult=1/(x2-x1)/(y2-y1)
          rad[row,locy*nx+locx]=rad[row,locy*nx+locx]+ $
                      (x2-xp)*(y2-yp)*mult       ;(1,1) CHECK THESE
          rad[row,locy*nx+locx2]=rad[row,locy*nx+locx2]+ $
                      (x2-xp)*(yp-y1)*mult       ;(1,2)
          rad[row,locy2*nx+locx]=rad[row,locy2*nx+locx]+ $
                      (xp-x1)*(y2-yp)*mult       ;(2,1)
          rad[row,locy2*nx+locx2]=rad[row,locy2*nx+locx2]+ $
                      (xp-x1)*(yp-y1)*mult       ;(2,2)
          nL++
        endif
      endif
    endfor
    L[row]=nL
  endfor
endfor



return
end




; notes

; each row of the matrix corresponds to a single (r,theta) pair
; the rows are number in order first of increasing r, then increasing
; theta.

;	r1, theta1
;	r2, theta1
;     .
;     .
;     .
;     r1, theta2
;     r2, theta2
;     .
;	.
;	.
; the image is assumed to be in a column vector that increases in firt 
; x and then y
