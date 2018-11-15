function linear,thet,r,ts 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Calculates x,y coordinate points along a straight ;;
;;           line, y=mx+b.                                    ;;
;;                                                            ;;
;; Inputs: thet - slope proxy                                 ;;
;;         r - closest approach of line to origin             ;;
;;         ts - n-elements vector of points at which to calc- ;;
;;           ulate x,y coordinates (distances along line      ;;
;;           from y-intercept)                                ;;
;;                                                            ;;
;; Returns: 2*n array of (x,y) coordinates corresponding to   ;;
;;            each elements in ts array                       ;;
;;                                                            ;;
;; Created: 03/31/08                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


num=N_ELEMENTS(ts)
xs=FLTARR(2,num)
cthet=COS(thet)
sthet=SIN(thet)
xs[0,*]=-cthet*r+sthet*ts
xs[1,*]=-sthet*r-cthet*ts


return,xs
end


