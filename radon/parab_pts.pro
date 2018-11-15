function parab_pts,a,xs,ys,ts

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Calculates x,y coordinates for points along a     ;;
;;            parabola with vertical symmetry.                ;;
;;                                                            ;;
;; Inputs: a - acceleration of parabola                       ;;
;;         xs,ys - x,y coordinates of parabola at image bounds;;
;;         ts - parametrization of points at which to interp- ;;
;;           olate image                                      ;; 
;;                                                            ;;
;; Returns: 2xn array of x,y coordinates at which to interpol-;;
;;            ate the image for a parabolic projection        ;;
;;                                                            ;;
;;                                                            ;;
;; Created: 04/01/08                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; input checking
if xs[1] le xs[0] or ys[1] le ys[0] then begin
  print,'WTF?'
  stop
endif
;stop


; initializations
num=N_ELEMENTS(ts)
rs=FLTARR(2,num)



; main
b=(a*(xs[0]^2-xs[1]^2)-ys[0]+ys[1])/(xs[1]-xs[0])
c=((ys[0]*xs[1]-ys[1]*xs[0])-a*xs[0]*xs[1]*(xs[0]-xs[1]))/(xs[1]-xs[0])
ts12=2.*a*xs+b
if ts12[0] gt ts12[1] then ts12=REVERSE(ts12)
;tmid=a*(xs[1]+xs[0])+b
;tsp=ts[WHERE(ts ge ts12[0] and ts le ts12[1])]
tsp=ts*(ts12[1]-ts12[0])+ts12[0]
rs[0,*]=1./(2.*a)*(tsp-b)
rs[1,*]=1./(4.*a)*(tsp^2-b^2)+c



return,rs
end



; obsolete

;yxmin=a*xrange[0]^2+b*xrange[0]+c*xrange[0]
;yxmax=a*xrange[1]^2+b*xrange[1]+c*xrange[1]
;if yxmin lt yrange[0] then ycurvemin=yrange[0] else ycurvemin=yxmin
;if yxmax gt yrange[1] then ycurvemax=yrange[1] else ycurvemax=yxmax
;arg=SQRT(b^2-4.*a*(c-ycurvemin))
;xcurvemin=-b/2./a+arg/2./ABS(a)
;arg=SQRT(b^2-4.*a*(c-ycurvemax))
;xcurvemin=-b/2./a+arg/2./ABS(a)
;xhalf=xcurvemin+(xcurvemax-xcurvemin)/2.
;thalf=2.*xhalf*a+b
;newts=ts+thalf
;rs[0,*]=(ts-b)/2./a
;rs[1,*]=(ts^2+4.*a*c-b^2)/4./a

