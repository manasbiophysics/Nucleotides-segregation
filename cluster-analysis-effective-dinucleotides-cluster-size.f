
c	IMPLICIT CHARACTER*1 (A-Z)
        CHARACTER*80 val1,LINE
        CHARACTER*4 ATOM(30000), RESIDUE(30000),atomn(30000)
        INTEGER   RESNUM(30000)
        character*1 chain(30000)
        character*6 box
        DIMENSION X(30000),Y(30000),Z(30000),dist(3000,3000)
        dimension xp1(1000),yp1(1000),zp1(1000),natom(30000)
        dimension xp2(1000),yp2(1000),zp2(1000)
        DIMENSION xc1(1000),yc1(1000),zc1(1000)
        DIMENSION xn1(1000),yn1(1000),zn1(1000)
        DIMENSION xo1(1000),yo1(1000),zo1(1000)
        DIMENSION xc2(1000),yc2(1000),zc2(1000)
        DIMENSION xn2(1000),yn2(1000),zn2(1000)
        DIMENSION xo2(1000),yo2(1000),zo2(1000)
        dimension xm(30000),ym(30000),zm(30000)
        dimension ncluster(500,500),na(1000),nb(1000),nc(1000)
        dimension ncluster2(500,500),nd(1000),ne(1000),nnc(1000,1000)
        dimension nhist(1000),hist(1000),mcsize(1000)
        
        call getarg(1,val1)


        OPEN(UNIT=2,FILE=val1,STATUS='OLD')
          cut = 10.0       

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Read the coordinate of the phosphate
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           READ(2,*) box,bx,by,bz
c             write(*,*)bx,by,bz
               j = 0
               l = 0
               lc1 = 0
               lc2 = 0
               lc3 = 0
               lc4 = 0
               lt1 = 0
               lt2 = 0
               lt3 = 0
               lt4 = 0

          DO I=1,200000
            READ(2,1,end=99) LINE
            if(LINE(1:4) .EQ. 'ATOM') THEN
              J = J + 1
         READ(LINE,2)atomn(J),natom(J),ATOM(J),RESIDUE(J),chain(j),
     1   RESNUM(j),X(J),Y(J),Z(J)

c         write(*,2)atomn(J),natom(J),ATOM(J),RESIDUE(J),chain(j),
c     1   RESNUM(j),X(J),Y(J),Z(J)
c..............................................................................................
           if(j.le.3100)then

           if(mod(j,62).ne.0)then

              if((ATOM(J)(2:2)).eq.'C')then
                      lc1=lc1+1
                  xc1(lc1) = x(j)
                  yc1(lc1) = y(j)
                  zc1(lc1) = z(j)
               endif
              if((ATOM(J)(2:2)).eq.'N')then
                      lc2=lc2+1
                  xn1(lc2) = x(j)
                  yn1(lc2) = y(j)
                  zn1(lc2) = z(j)
               endif
              if((ATOM(J)(2:2)).eq.'O')then
                      lc3=lc3+1
                  xo1(lc3) = x(j)
                  yo1(lc3) = y(j)
                  zo1(lc3) = z(j)
               endif
              if((ATOM(J)(2:2)).eq.'P')then
                      lc4=lc4+1
                  xp1(lc4) = x(j)
                  yp1(lc4) = y(j)
                  zp1(lc4) = z(j)
               endif
           else
                 l = l +1

                 xcc=0.0
                 ycc=0.0
                 zcc=0.0
               do n=1,lc1
                xcc=xcc+12.0*xc1(n)
                ycc=ycc+12.0*yc1(n)
                zcc=zcc+12.0*zc1(n)
               enddo

                 xnc=0.0
                 ync=0.0
                 znc=0.0
               do n=1,lc2
                xnc=xnc+14.0*xn1(n)
                ync=ync+14.0*yn1(n)
                znc=znc+14.0*zn1(n)
               enddo

                 xoc=0.0
                 yoc=0.0
                 zoc=0.0
              do n=1,lc3
                xoc=xoc+16.0*xo1(n)
                yoc=yoc+16.0*yo1(n)
                zoc=zoc+16.0*zo1(n)
               enddo

                 xpc=0.0
                 ypc=0.0
                 zpc=0.0
               do n=1,lc4
                xpc=xpc+31.0*xp1(n)
                ypc=ypc+31.0*yp1(n)
                zpc=zpc+31.0*zp1(n)
               enddo

            tmass = (12.0*lc1+14.0*lc2+16.0*lc3+31.0*lc4)

                 xm(l) = (xcc+xnc+xoc+xpc)/tmass
                 ym(l) = (ycc+ync+yoc+ypc)/tmass
                 zm(l) = (zcc+znc+zoc+zpc)/tmass

                  lc1 = 0
                  lc2 = 0
                  lc3 = 0
                  lc4 = 0

c              write(*,20)l,l,xm(l),ym(l),zm(l)

             endif
            endif

c..........................................................................

          endif
          enddo
 99     continue
c              write(*,*)'TER'
c.........................................................................


c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
C Generate cm-cm distance matrix
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
                do m=1,l
                  do n=1,l


                Rxmn = xm(n) - xm(m)
                Rymn = ym(n) - ym(m)
                Rzmn = zm(n) - zm(m)
c Consideration of periodic boundary condition
                   Rxmn =  Rxmn - bx*anint(Rxmn/bx)
                   Rymn =  Rymn - by*anint(Rymn/by)
                   Rzmn =  Rzmn - bz*anint(Rzmn/bz)  

             dist(m,n)=sqrt(Rxmn*Rxmn+Rymn*Rymn+Rzmn*Rzmn)

                 if(dist(m,n).le.cut)then
                    dist(m,n)=1.0
                  else
                    dist(m,n)=0.0
                  endif

                  enddo
                 enddo
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ndata = l
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Calculate neighbour list
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      do i = 1,ndata
                         na(i)=0
                      enddo

                  do i =1,ndata
                   do j = 1,ndata
                   
                     if((j.ne.i) .and. (dist(i,j).eq.1.0)) then
                              na(i) = na(i)+1
                              ncluster(i,na(i)) = j
                     endif
                    enddo
c                      write(*,*)i,na(i)
c                      write(10,*)i,(ncluster(i,k),k=1,na(i))
                   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Calculate larger cluster size
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                         do i = 1,ndata
                             nb(i)=na(i)
                            do 100 j = 1,ndata
                             if(i.ne.j)then
                               do k = 1,na(i)
                            if(ncluster(i,k).eq.j) GOTO 100
                                do l = 1,na(j)
                      if(ncluster(i,k).eq.ncluster(j,l)) then
                             do m =1,na(i)
                               if(j.eq.ncluster(i,m))then
                                GOTO 100
                               endif
                             enddo
                                nb(i)=nb(i)+1
                           ncluster(i,nb(i)) = j
                          GOTO 100
                      endif
                                enddo
                               enddo
                             endif
 100                        enddo
c                     write(*,*)i,nb(i)
c                     write(20,*)i,(ncluster(i,k),k=1,nb(i))
                         enddo
c....................................................................
                         do i =1,ndata
                           nnc(1,i) = nb(i)
                         enddo
                   do n = 2,ndata
                        do i = 1,ndata
                             nnc(n,i)=nnc(n-1,i)
                            do 200 j = 1,ndata
                             if(i.ne.j)then
                               do k = 1,nnc(n-1,i)
                            if(ncluster(i,k).eq.j) GOTO 200
                                do l = 1,nnc(n-1,j)
                      if(ncluster(i,k).eq.ncluster(j,l)) then
                             do m =1,nnc(n-1,i)
                               if(j.eq.ncluster(i,m))then
                                GOTO 200
                               endif
                             enddo
                                nnc(n,i)=nnc(n,i)+1
                           ncluster(i,nnc(n,i)) = j
                          GOTO 200
                      endif
                                enddo
                               enddo
                             endif
 200                        enddo
                         enddo
                      enddo
     
                   do i = 1,ndata
c                     write(*,*)i,nc(i)
                     write(30,*)i,(ncluster(i,k),k=1,nnc(ndata,i))
                   enddo
c..........................................................................
                        
                        maxcluster = 0 
                         do i=1,ndata
                           if(nnc(ndata,i).ge.maxcluster)then
                            maxcluster = nnc(ndata,i)
                            ncenter = i
                           endif
                         enddo
                 write(30,*) ncenter,"maxcluster size=", maxcluster
         
!                     if(ncenter.le.25)then
!                        maxcg = 1
!                        maxta = 0
!                     else
!                        maxcg = 0
!                        maxta = 1
!                     endif
                       
!                      do i=1,maxcluster
!                if(ncluster(ncenter,i).le.25)then
!                    maxcg = maxcg+1
!                else
!                    maxta = maxta+1
!                endif
!                      enddo
!                write(30,*) ncenter, maxcg, maxta

c......................................................................
C Calculate 2nd large cluster size
c......................................................................
                         do i =1,ndata
                          nc(i) = nnc(ndata,i)
                         enddo                         


                         do i =1,ndata
                             k=0
                           do j =1,maxcluster
                  if((i.eq.ncenter).or.(i.eq.ncluster(ncenter,j)))then
                                k=k+1
                  endif 
                           enddo
                      if(k.eq.0)then
                          nd(i)=nc(i)
                          do l=1,nc(i)
                      ncluster2(i,l)=ncluster(i,l)
                          enddo
                       else
                          nd(i)=0
                          do l=1,nc(i)
                      ncluster2(i,l)=0
                          enddo
c                  write(*,*)i,(ncluster(i,l),l=1,nc(i))
                      endif
                         enddo       

                         do i=1,ndata
c                     write(*,*)i,(ncluster2(i,l),l=1,nc(i))
                         enddo              
c...................................................................
                         do i = 1,ndata
                             ne(i)=nd(i)
                            do j = 1,ndata
                             if(i.ne.j)then
                               do k = 1,nd(i)
                            if(ncluster2(i,k).eq.j) GOTO 300
                                do l = 1,nd(j)
                      if(ncluster2(i,k).eq.ncluster2(j,l)) then
                             do m =1,nd(i)
                               if(j.eq.ncluster2(i,m))then
                                GOTO 300
                               endif
                             enddo
                                ne(i)=ne(i)+1
                           ncluster2(i,ne(i)) = j
                          GOTO 300
                      endif
                                enddo
                               enddo
                             endif
 300                        enddo
c                     write(*,*)i,ne(i)
                     write(40,*)i,(ncluster2(i,k),k=1,ne(i))
                         enddo
c...........................................................................
                        maxcluster2 = 0
                         do i=1,ndata
                           if(ne(i).ge.maxcluster2)then
                            maxcluster2 = ne(i)
                            ncenter2 = i
                           endif
                         enddo
                   write(40,*) ncenter2, "cluster2=", maxcluster2


c------------------------------------------------------------------ 
C Calculate 3rd large cluster size
C.................................................................
                       maxcluster3 = 0
                         do i=1,ndata
             if((ne(i).lt.maxcluster2).and.ne(i).ge.maxcluster3)then
                       maxcluster3 = ne(i)
                        ncenter3 = i
             endif
                         enddo
                    write(40,*) ncenter3, "cluster3=", maxcluster3
C..........................................................................
C...........................................................................
                          maxcluster = maxcluster+1
                          maxcluster2 = maxcluster2+1
                          maxcluster3 = maxcluster3+1
C..............................................................................
C............................................................................
C Clustering entropy
C...................................................................................
                            mcsize(1) = maxcluster
                            mcsize(2) = maxcluster2
                            mcsize(3) = maxcluster3
               mcsize(4) = ndata - (maxcluster+maxcluster2+maxcluster3)

                              entropy = 0.0
                            do i = 1,3
                               hist(i)=mcsize(i)/50.0
                            enddo
                               hist(4)=1/50.0

                          do i = 1,3
                        if(hist(i).gt.0)then
                    entropy = entropy+(hist(i)*(log(hist(i))))
c               entropy= entropy+(hist(i)*((log(hist(i)))/(log(2.0))))
                        endif
                         enddo
                   entropy = entropy+mcsize(4)*(hist(4)*(log(hist(4))))

                         fentropy = -(entropy)
                         ecluster = exp(fentropy)
C...................................................................................

        write(*,10) val1,maxcluster,maxcg,maxta,maxcluster2,maxcluster3,
     1              fentropy,ecluster
C.............................................................................
C.............................................................................


 1	FORMAT(A80)
 2	FORMAT(a4,i7,1x,2A4,1x,a1,i4,4X,3F8.3)
 10      FORMAT(a40,5i5,2f8.3)
 20     FORMAT('ATOM',i7,1x,' C  ',' MC ',2x,i4,4x,3F8.3)
 	
        stop
        END
