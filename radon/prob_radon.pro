pro prob_radon,image,xcoords,ycoords,thetas,rs,probs,rad

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Calculates prob that points along each line in the;;
;;          transform space are from a different distribution ;;
;;          than the image as a whole.                        ;;
;;                                                            ;;
;; Inputs: image - original image                             ;;
;;         xcoords,ycoords - coordinate values associated     ;;
;;           with each column,row of pixels                   ;;
;;         thetas,rs - theta,r values for each column,row in  ;;
;;           transform space                                  ;; 
;;                                                            ;;
;; Outputs: probs - transform space array showing the probab- ;;
;;            ility that the line contains points with a diff-;;
;;            erent mean than that of the original image      ;;
;;          means - the approximate mean along each line      ;;
;;                                                            ;;
;; Created: 01/21/08                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; parameters
nts_factor=2
tmax_factor=1.1
upbnd=20
lwbnd=-20
mintm=0.5



; input checks
szxcoords=SIZE(xcoords)
szycoords=SIZE(ycoords)
if szxcoords[0] ne 1 then begin
  print,'Variable xcoords must be a 1D array'
  trans=-1
  return
endif
if szycoords[0] ne 1 then begin
  print,'Variable ycoords must be a 1D array'
  trans=-1
  return
endif
szthetas=SIZE(thetas)
szrs=SIZE(rs)
if szthetas[0] ne 1 then begin
  print,'Variable thetas must be a 1D array'
  trans=-1
  return
endif
if szrs[0] ne 1 then begin
  print,'Variable rs must be a 1D array'
  trans=-1
  return
endif
locs=WHERE(xcoords[1:*]-xcoords le 0)  ; want strictly increasing x,y
if locs[0] ne -1 then begin
  print,'Problem with x coordinate vector'
  trans=-1
  return
endif
locs=WHERE(ycoords[1:*]-ycoords le 0)
if locs[0] ne -1 then begin
  print,'Problem with y coordinate vector'
  trans=-1
  return
endif



; initializations
nx=szxcoords[1]
ny=szycoords[1]
nthet=szthetas[1]
nr=szrs[1]
trythis=FLTARR(nthet,nr)+1.0
;means=FLTARR(nthet,nr)
rad=FLTARR(nthet,nr)
thetasp=thetas
rsp=rs



; adjust coordinate axes
xcoordsp=xcoords-xcoords[szxcoords[1]/2]
ycoordsp=ycoords-ycoords[szycoords[1]/2]
xcoordsp=xcoordsp/(xcoordsp[nx-1]-xcoordsp[0])*nx
ycoordsp=ycoordsp/(ycoordsp[ny-1]-ycoordsp[0])*ny
minx=MIN(xcoordsp,max=maxx)
miny=MIN(ycoordsp,max=maxy)
imhist=HISTOGRAM(image,locations=imlocs)
imres=GAUSSFIT(imlocs,imhist,est,nterms=3)



; initialize ts
nts=ROUND(SQRT(nx^2+ny^2)*nts_factor)
if nts MOD 2 eq 1 then t=FINDGEN(nts)-nts/2 else $
          t=FINDGEN(nts+1)-nts/2
tmax=tmax_factor*sqrt(nx^2+ny^2)
t=t/nts*tmax*nts_factor



; find points along each line
for thet=0,nthet-1 do begin
  cthet=COS(thetasp[thet])
  sthet=SIN(thetasp[thet])
  for r=0,nr-1 do begin
    pixval=FLTARR(nts)+100.
    r0=rsp[r]
    x=-cthet*r0+sthet*t
    y=-sthet*r0-cthet*t
    for point=0,nts-1 do begin
      xp=x[point]
      yp=y[point]
      if (xp le maxx) and (xp ge minx) then begin
        if (yp le maxy) and (yp ge miny) then begin
          garbage=MIN(ABS(xp-xcoordsp),locx)
          if xp gt xcoordsp[locx] then begin
            locx2=locx++
            x2=xcoordsp[locx2]
            x1=xcoordsp[locx]
          endif else begin
            locx2=locx
            locx--
            x2=xcoordsp[locx2]
            x1=xcoordsp[locx]
          endelse
          garbage=MIN(ABS(yp-ycoordsp),locy)
          if yp gt ycoordsp[locy] then begin
            locy2=locy++
            y2=ycoordsp[locy2]
            y1=ycoordsp[locy]
          endif else begin
            locy2=locy
            locy--
            y2=ycoordsp[locy2]
            y1=ycoordsp[locy]
          endelse
          mult=1/(x2-x1)/(y2-y1)
          pixval[point]=(x2-xp)*(y2-yp)*image[locx,locy]*mult+ $
              (x2-xp)*(yp-y1)*image[locx,locy2]*mult+ $
              (xp-x1)*(y2-yp)*image[locx2,locy]*mult+ $
              (xp-x1)*(yp-y1)*image[locx2,locy2]*mult
        endif
      endif
    endfor
    locs=where(pixval le 99.)
    if n_elements(locs) gt 30 then begin
      pix=pixval[locs]
      tm=TM_TEST(image,pix)
      trythis[thet,r]=tm[1]
;      means[thet,r]=MEAN(pix)
      rad[thet,r]=TOTAL(pix)
    endif
  endfor
endfor
;stop


;stop
probs=1.0-TEMPORARY(trythis)
return
end


; obsolete

;      if tm[1] gt mintm then begin
;        pixhist=HISTOGRAM(pix,locations=pixlocs)
;        pixresl=GAUSSFIT(pixlocs,pixhist,a,nterms=3,estimates=est)
;        if (a[1] lt upbnd) and (a[1] gt lwbnd) then $
;                    means[thet,r]=a[1]
;      endif
;      means[thet,r]=MEAN(pix)
