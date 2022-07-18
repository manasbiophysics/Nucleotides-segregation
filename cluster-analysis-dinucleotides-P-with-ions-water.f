
c	IMPLICIT CHARACTER*1 (A-Z)
        CHARACTER*80 val1,LINE
       CHARACTER*4 ATOM(300000), RESIDUE(300000),atomn(300000)
        INTEGER   RESNUM(300000)
        character*1 chain(300000)
        character*6 box
        DIMENSION X(300000),Y(300000),Z(300000),dist(3000,3000)
        dimension xp(30000),yp(30000),zp(30000),natom(30000)
        dimension xion(30000),yion(30000),zion(30000)
        dimension xwat(300000),ywat(300000),zwat(300000)
        dimension ncluster(500,500),na(1000),nb(1000),nc(1000)
        dimension ncluster2(500,500),nd(1000),ne(1000),nnc(1000,1000)
        dimension nhist(1000),hist(1000),mcsize(1000)
        dimension ion(30000),iwat(300000)
        dimension ionf(30000),iwatf(300000)
        
        call getarg(1,val1)


        OPEN(UNIT=2,FILE=val1,STATUS='OLD')
          cut = 12.0       ! P-P distance cutoff 

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Read the coordinate of the phosphate
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           READ(2,*) box,bx,by,bz
c             write(*,*)bx,by,bz
               j = 0
               l = 0
               m1 = 0
               m2 = 0 

          DO I=1,2000000
            READ(2,1,end=99) LINE
            if(LINE(1:4) .EQ. 'ATOM') THEN
              J = J + 1
         READ(LINE,2)atomn(J),natom(J),ATOM(J),RESIDUE(J),chain(j),
     1   RESNUM(j),X(J),Y(J),Z(J)

c         write(*,2)atomn(J),natom(J),ATOM(J),RESIDUE(J),chain(j),
c     1   RESNUM(j),X(J),Y(J),Z(J)

              if(ATOM(J).eq.' P  ')then
                      l=l+1
                  xp(l) = x(j)
                  yp(l) = y(j)
                  zp(l) = z(j)
c              write(*,*) xp(l),yp(l),zp(l)
              endif

              if(ATOM(J).eq.'Na+ ')then                
                      m1=m1+1
                  xion(m1) = x(j)
                  yion(m1) = y(j)
                  zion(m1) = z(j)
c              write(*,*) xion(m1),yion(m1),zion(m1)
              endif

              if((RESIDUE(J).eq.' WAT').and.(ATOM(J).eq.' O  '))then 
                      m2=m2+1
                  xwat(m2) = x(j)
                  ywat(m2) = y(j)
                  zwat(m2) = z(j)
c              write(*,*) xwat(m2),ywat(m2),zwat(m2)
              endif


            endif
          enddo
  99     continue


c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
C Generate P-P distance matrix
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
                do m=1,l
                  do n=1,l

c        dist(m,n)=sqrt((xp(n)-xp(m))**2+(yp(n)-yp(m))**2+
c     1                         (zp(n)-zp(m))**2)

                Rxmn = xp(n) - xp(m)
                Rymn = yp(n) - yp(m)
                Rzmn = zp(n) - zp(m)
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
         

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
                       nion = 0
                       nwat = 0

                      do i=1,maxcluster
                       n = ncluster(ncenter,i)

c                     write(50,*)n

                  do n1=1,m1

                Rx1 = xion(n1) - xp(n)
                Ry1 = yion(n1) - yp(n)
                Rz1 = zion(n1) - zp(n)
c Consideration of periodic boundary condition
                   Rx1 =  Rx1 - bx*anint(Rx1/bx)
                   Ry1 =  Ry1 - by*anint(Ry1/by)
                   Rz1 =  Rz1 - bz*anint(Rz1/bz)

             distion=sqrt(Rx1*Rx1+Ry1*Ry1+Rz1*Rz1)

                 if(distion.le.5.0)then
                      nion = nion + 1
                      ion(nion) = n1
                  endif

                  enddo

 
                 do n2=1,m2

                Rx2 = xwat(n2) - xp(n)
                Ry2 = ywat(n2) - yp(n)
                Rz2 = zwat(n2) - zp(n)
c Consideration of periodic boundary condition
                   Rx2 =  Rx2 - bx*anint(Rx2/bx)
                   Ry2 =  Ry2 - by*anint(Ry2/by)
                   Rz2 =  Rz2 - bz*anint(Rz2/bz)

             distwat=sqrt(Rx2*Rx2+Ry2*Ry2+Rz2*Rz2)

                 if(distwat.le.5.0)then
                      nwat = nwat + 1
                      iwat(nwat) = n2
                  endif

                  enddo


              enddo


                  k1 = 1
                ionf(1) = ion(1)
               do i1=2,nion
                if (any( ionf == ion(i1) )) cycle
                  k1 = k1 + 1
                 ionf(k1) = ion(i1)
                end do

                   k2 = 1
                iwatf(1) = iwat(1)
               do i2=2,nwat
                if (any( iwatf == iwat(i2) )) cycle
                  k2 = k2 + 1
                 iwatf(k2) = iwat(i2)
                end do


                 write(50,*) 'ions= ',k1
                 write(50,*) (ionf(i),i=1,k1)
                 write(50,*) 'water= ',k2
                 write(50,*) (iwatf(i),i=1,k2)

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c                     if(ncenter.le.25)then
c                        maxcg = 1
c                        maxta = 0
c                     else
c                        maxcg = 0
c                        maxta = 1
c                     endif
                       
c                      do i=1,maxcluster
c                if(ncluster(ncenter,i).le.25)then
c                    maxcg = maxcg+1
c                else
c                    maxta = maxta+1
c                endif
c                      enddo
c                write(30,*) ncenter, maxcg, maxta

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
C............................................................................
c............................................................................
C Clustering entropy
C............................................................................
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
        write(*,10) val1,maxcluster,maxcg,maxta,maxcluster2,maxcluster3,
     1              fentropy,ecluster,k1,k2
C.............................................................................
C.............................................................................


 1	FORMAT(A80)
 2	FORMAT(a4,i7,1x,2A4,1x,a1,i4,4X,3F8.3)
 10      FORMAT(a40,5i5,2f8.3,2i7)
 	
        stop
        END
