pro cluster_mask3d,trans,mask 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Clusters contiguous groups of non-zero pixels in  ;;
;;          a 3D transform, creating a mask that numbers con- ;;
;;          tiguous clusters sequentially.                    ;;
;;                                                            ;;
;; Inputs: trans - 3D array containing a Radon transform to   ;;
;;           be masked                                        ;;
;;                                                            ;;
;; Outputs: mask - array of same size as trans; entries con-  ;;
;;            tain the number of the contiguous group the     ;;
;;            corresponding entry in trans belongs to         ;;
;;                                                            ;;
;; Dependencies: none                                         ;;
;;                                                            ;;
;; Created: 05/02/08                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; input checking
sztrans=SIZE(trans)
if sztrans[0] ne 3 then begin
  print,'Variable trans must be a 3D array.'
  mask=-1
  return
endif



; initializations
nx=sztrans[1]
ny=sztrans[2]
nplanes=sztrans[3]
mask=INTARR(nx,ny,nplanes)
locs=WHERE(trans ne 0)
if locs[0] eq -1 then begin
  print,'wtf?'
  stop
endif
inds=ARRAY_INDICES(trans,locs)
nclusts=N_ELEMENTS(locs)
mask[inds[0,*],inds[1,*],inds[2,*]]=LINDGEN(nclusts)+1
oldnclusts=nclusts



; aggregate adjacent pixels
adjacent=INTARR(3,3,3)
adjacent[*,*,1]=[[0,1,0],[1,1,1],[0,1,0]]
adjacent[1,1,0]=1
adjacent[1,1,2]=1
unique_vals=mask[UNIQ(mask,SORT(mask))]
repeat begin
  for i=1L,nclusts do begin
    locs=WHERE(mask eq unique_vals[i])
    if locs[0] ne -1 then begin
      szclust=N_ELEMENTS(locs)
      pos=ARRAY_INDICES(mask,locs)
      for j=0L,szclust-1 do begin
        x=pos[0,j]
        y=pos[1,j]
        z=pos[2,j]
        mins=pos[*,j]-1
        maxes=pos[*,j]+1
        temp=adjacent
        if x eq 0 then begin
          temp=temp[1:2,*,*]
          mins[0]=0
        endif else if x eq nx-1 then begin
          temp=temp[0:1,*,*]
          maxes[0]=nx-1
        endif
        if y eq 0 then begin
          temp=temp[*,1:2,*]
          mins[1]=0
        endif else if y eq ny-1 then begin
          temp=temp[*,0:1,*]
          maxes[1]=ny-1
        endif
        if z eq 0 then begin
          temp=temp[*,*,1:2]
          mins[2]=0
        endif else if z eq nplanes-1 then begin
          temp=temp[*,*,0:1]
          maxes[2]=nplanes-1
        endif
        subarray=temp*mask[mins[0]:maxes[0],mins[1]: $
                    maxes[1],mins[2]:maxes[2]]
        subnonzero=WHERE(subarray ne 0)
        submin=MIN(subarray[subnonzero])
        inds=ARRAY_INDICES(subarray,subnonzero)
        mask[mins[0]+inds[0,*],mins[1]+inds[1,*],mins[2]+ $
               inds[2,*]]=submin
      endfor
    endif
  endfor
  oldoldnclusts=oldnclusts
  oldnclusts=nclusts
;  uniques=UNIQ(mask,SORT(mask))
  unique_vals=mask[UNIQ(mask,SORT(mask))]
  nclusts=N_ELEMENTS(unique_vals)-1
endrep until (((nclusts eq oldoldnclusts) and (nclusts eq $
              oldnclusts)) or (nclusts eq 1))



; re-number clusters
;vals=mask[UNIQ(mask,SORT(mask))]
for i=1L,nclusts do begin
  locs=WHERE(mask eq unique_vals[i])
  mask[locs]=i
endfor



end
