
c	IMPLICIT CHARACTER*1 (A-Z)
        CHARACTER*80 val1,LINE
        CHARACTER*4 ATOM(300000), RESIDUE(300000),atomn(300000)
        INTEGER   RESNUM(300000)
        character*1 chain(300000)
        character*6 box
        DIMENSION X(300000),Y(300000),Z(300000),dist(3000,3000)
        dimension xp(1000),yp(1000),zp(1000),natom(300000)
        dimension ncluster(500,500),na(1000)
        real      ccf(1000)
        dimension nhist(1000),hist(1000),nc(1000)
C	
        
        call getarg(1,val1)


        OPEN(UNIT=2,FILE=val1,STATUS='OLD')
          cut = 12.0        

c      OPEN(UNIT=10,FILE='graph-cluster-system-CG-50-100-200ns-300K.dat')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Read the coordinate of the phosphate
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           READ(2,*) box,bx,by,bz
c             write(*,*)bx,by,bz
               j = 0
               l = 0
          DO I=1,200000
            READ(2,1,end=99) LINE
            if(LINE(1:4) .EQ. 'ATOM') THEN
              J = J + 1
         READ(LINE,2)atomn(J),natom(J),ATOM(J),RESIDUE(J),chain(j),
     1   RESNUM(j),X(J),Y(J),Z(J)


              if((ATOM(J)(2:2)).eq.'P')then
                      l=l+1
                  xp(l) = x(j)
                  yp(l) = y(j)
                  zp(l) = z(j)
               endif

          endif
          enddo
 99     continue

c             write(*,*)l

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
C Generate P-P distance matrix based on the cutoff distance
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
                do m=1,l
                  do n=1,l


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
c                   do m=1,l
c                write(*,*)(dist(m,n),n=1,l)
c                   enddo
c                 write(*,*) ndata
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Calculate neighbour list for each vartix of graph
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
c                      write(100,*)i,na(i)
c                      write(200,*)i,(ncluster(i,k),k=1,na(i))
                   enddo
C.....................................................................
C calculate average no of neighbour associated with a vartex
C......................................................................
                     sum_na = 0.0
                    do i=1,ndata
                      sum_na = sum_na + na(i)
                    enddo
                     avg_na = sum_na/ndata

C..................................................................
c  Calculate clustering coefficient
C..................................................................
                   do i=1,ndata
                     nc(i) = 0
                     ccf(i) = 0.0
                   enddo
                      sum_ccf = 0.0

                   do i=1,ndata
                     do j=1,na(i)-1
                         ic = ncluster(i,j)
                       do k=j+1,na(i)
                        do l =1,na(ic)
                  if(ncluster(i,k).eq.ncluster(ic,l))then
                          nc(i) = nc(i)+1
                  endif
                        enddo
                       enddo
                     enddo
                        connect = (na(i)*(na(i)-1))/2
                        if(connect.ne.0.0)then
                        ccf(i) = (nc(i)/connect)
                        else
                        ccf(i) = 0.0
                        endif
                        sum_ccf = sum_ccf+ccf(i)
c                      write(300,*)i,na(i),nc(i),connect,ccf(i)
                   enddo
                        asum_ccf = sum_ccf/ndata
C.......................................................................
C Calculation of graph entropy
C.....................................................................
                  do i =0,ndata
                    nhist(i)=0
                  enddo
                    entropy =0.0
                 do i =0,ndata
                   do j =1,ndata
                    if(na(j).eq.i)then
                     nhist(i)=nhist(i)+1
                    endif
                   enddo
                     hist(i)=nhist(i)/real(ndata)
c                   write(*,*)i,hist(i)
                 enddo

                do i =0,ndata
                 if(hist(i).gt.0.0)then
                 entropy= entropy+(hist(i)*log(hist(i)))
                 endif
                enddo
                    fentropy= -(entropy)
c........................................................................


                 write(*,10)(na(i),i=1,ndata),avg_na,asum_ccf,fentropy

C.............................................................................
C.............................................................................


 1	FORMAT(A80)
 2	FORMAT(a4,i7,1x,2A4,1x,a1,i4,4X,3F8.3)
 20     FORMAT('ATOM',i7,1x,' C  ',' MC ',2x,i4,4x,3F8.3)
 10      FORMAT(50i4,3f8.3)
 	
        stop
        END
