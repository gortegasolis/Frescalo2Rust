C A program to calculate similarity of grid squares
      integer neigh,i,ii,m,mm,j,jj,n,nn,is2,i1,i2
      real big,small
      parameter(mm=4000,nn=10000,nndat=5000000)
      parameter(big=1000000,small=0.00005)
      integer itot(mm),iocc(mm),index(mm),jtot(nn),iseqq(mm,mm)
      real simil(mm,mm),sim(mm)
      character*30 ddji(nndat),distii(nndat),dji
      character*20 filein, fileou
      character*10 samp,samp1,spec,any,blnk10/'          '/,word
      character*10 sa(mm),sp(nn),d(3)
      character*80 blurb,blnk80
      character*1 b(80),w(10)
      equivalence (b,blurb),(w,word),(d,dji)
    
      m=0
      n=0
      ndist=0
      ndata=0
      do ii=1,mm
         itot(ii)=0
         enddo
      do jj=1,nn
         jtot(jj)=0
         enddo
C Set up files for reading and writing

      write(*,2000) mm,nn,nndat
 2000 format(//1x,'NEIGHSIM - Neighbourhood similarity based on',
     1       ' training-set species and physical proximity'/
     2       1x,'written by Mark Hill, January-June 2011'//
     3       1x,'NOTE PARAMETER LIMITS: Sites',i5,' Species',i6,
     4       1x,'Training-set number of records',i9//
     5       1x,'Type name of input file with Training-set species ',
     6       'data [sample species] ....')
      call filin(filein,4)
      write(*,2005)
 2005 format(1x,'Type name of input file with physical distances ....')
      call filin(filein,3)
      write(*,2001)
 2001 format(1x,'Type name of Training-set similarity output file ...')
      call filout(fileou,9)
      write(*,2006)
 2006 format(1x,'Type name of weights output file for use',
     1       ' in Frescalo ...')
      call filout(fileou,8)
      write(*,2002)
 2002 format(1x,'Type number of neighbours to include ...')
      read(*,*) neigh

C Read in distance data
   30 continue
      call getd(3,samp,samp1,any,iend,blurb,blnk80,b,word,blnk10,w)
      if(iend.eq.1) goto 50
      ndist=ndist+1
      if(ndist.gt.nndat) then
         write(*,*) ' Too many data items in physical distance file',
     1              ' - limit is',nndat
         call hold
         end if
      if(mod(ndist,20000).eq.0) write(*,*) ndist,' Dist ',samp,samp1,any
      d(1)=samp
      d(2)=samp1
      d(3)=any
      distii(ndist)=dji
      goto 30


   50 continue
      call getd(4,samp,spec,any,iend,blurb,blnk80,b,word,blnk10,w)
      if(iend.eq.1) goto 100
      call addwrd(sa,mm,m,samp)
      if(m.gt.mm) then
         write(*,*) ' Too many samples - limit is',mm
         call hold
         end if
      call addwrd(sp,nn,n,spec)
      if(n.gt.nn) then
         write(*,*) ' Too many species - limit is',nn
         call hold
         end if
      ndata=ndata+1
      if(ndata.gt.nndat) then
         write(*,*) ' Too many data items - limit is',nndat
         call hold
         end if
      if(mod(ndata,20000).eq.0) write(*,*)ndata,' Spdata ',samp,spec,any
      d(1)=spec
      d(2)=samp
      d(3)=any
      ddji(ndata)=dji
      goto 50

  100 continue
      write(*,*) 'Sorting main data ...'
      call sort30(ddji,ndata)
      write(*,*) 'Sort completed'

      do idata=1,ndata
         if(mod(idata,20000).eq.0) write(*,*)
     1                             ' Calculating totals',idata
         dji=ddji(idata)
         spec=d(1)
         samp=d(2)
         call binfnd(sa,m,samp,i)
         call binfnd(sp,n,spec,j)
         itot(i)=itot(i)+1
         jtot(j)=jtot(j)+1
         enddo

C Now start calculating similarity
      do i1=1,m
         do i2=1,m
            simil(i1,i2)=0
            enddo
         enddo
      idata=0
 
      do j=1,n
         do i=1,m
            iocc(i)=0
            enddo
         do iidata=1,jtot(j)
            idata=idata+1
            if(mod(idata,20000).eq.0) write(*,*) ' Similarities',idata
            dji=ddji(idata)
            samp=d(2)
            if(d(1).ne.sp(j)) then
               write(*,*) 'Unequal species',d(1),sp(j)
               call hold
               end if
            call binfnd(sa,m,samp,i)
            iocc(i)=-i
            enddo

         call isort(iocc,m)
         do ic=1,m
            if(iocc(ic).eq.0) then
               mcc=ic-1
               goto 150
               end if
            if(ic.eq.m) then
               mcc=m
               goto 150
               end if
            enddo
  150    continue
C mcc is the length of nonzero items in iocc
         do icc1=1,mcc
            i1=-iocc(icc1)
            do icc2=1,mcc
               i2=-iocc(icc2)
               simil(i1,i2)=simil(i1,i2)+1
C This calculates the number of species in common; later we divide to calc similarity
               enddo
            enddo
         enddo
C We have calculated the similarities, now print them out
C First multiply those within preferred region by big
      neigh1=0
      do idist=1,ndist
         if(mod(idist,20000).eq.0) write(*,*) idist,' Dist ',samp,samp1
         dji=distii(idist)
         samp=d(1)
         samp1=d(2)
         call binfnd(sa,m,samp,i1)
         call binfnd(sa,m,samp1,i2)
         if(i1.ne.0.and.i2.ne.0) then
            simil(i1,i2)=simil(i1,i2)*big
            any=d(3)
            call getnum(any,seq,word,w)
            iseqq(i1,i2)=ifix(seq)
            if(neigh1.lt.ifix(seq)) neigh1=ifix(seq)
            end if
         enddo
      if(neigh1.eq.0) write(8,*) 
     1                'Unrecognized sample names in distance data'

      do i1=1,m
         if(mod(i1,100).eq.0) write(*,*) 'Writing output  ',sa(i1),i1
         do i2=1,m
            sim(i2)=simil(i1,i2)*2.0/(itot(i1)+itot(i2))
            index(i2)=i2
            enddo
         call sort2(sim,index,m)
         do is2=1,neigh
            i2=m-is2+1
            if(i2.lt.1) goto 200
            iis2=index(i2)
            write(9,2030) sa(i1),sa(iis2),is2,sim(i2)/big
 2030 format(2a10,i5,f6.3,i5)
            if(neigh1.eq.0) goto 190
            amult1=(1-((float(is2)-1)/neigh)**2)**4
            amult2=(1-(float(iseqq(i1,iis2)-1)/neigh1)**2)**4
            amult=amult1*amult2
            if(amult.gt.small) write(8,2031) sa(i1),sa(iis2),
     1                         amult,amult1,amult2,neigh,neigh1
 2031 format(2a10,3f7.4,2i6)
  190       continue
            enddo
  200    continue
         enddo
      write(*,*) 'neigh,neigh1',neigh,neigh1

  999 continue
      close(8)
      close(9)
      call hold
      end

      subroutine sort2(dict,type,n)
C To sort two columns at once
      integer n
      real dict(n),djx,djjx,djjjx
      integer type(n),tjx,tjjx,tjjjx,i,j,jj,jjj
      do 10 i=1,n
      j=i
      djx=dict(j)
      tjx=type(j)
    5 if(j.eq.1) goto 8
      jj=j/2
      djjx=dict(jj)
      tjjx=type(jj)
      if(djjx.gt.djx) goto 8
      if(djjx.eq.djx) then
           if(tjjx.ge.tjx) goto 8
           end if
      dict(j)=djjx
      type(j)=tjjx
      j=jj
      goto 5
    8 dict(j)=djx
      type(j)=tjx
   10 continue
      i=n
      goto 14
   12 dict(j)=djx
      type(j)=tjx
   14 if(i.eq.1) return
      djx=dict(i)
      tjx=type(i)
      dict(i)=dict(1)
      type(i)=type(1)
      i=i-1
      j=1
      jj=2
   15 if(i-jj) 12,17,16
   16 djjx=dict(jj)
      tjjx=type(jj)
      jjj=jj+1
      djjjx=dict(jjj)
      tjjjx=type(jjj)
      if(djjx.gt.djjjx) goto 18
      if(djjx.eq.djjjx) then
          if(tjjx.ge.tjjjx) goto 18
          end if
      if(djx.gt.djjjx) goto 12
      if(djx.eq.djjjx) then
          if(tjx.ge.tjjjx) goto 12
          end if
      dict(j)=djjjx
      type(j)=tjjjx
      j=jjj
      jj=j*2
      goto 15
   17 djjx=dict(jj)
      tjjx=type(jj)
   18 if(djx.gt.djjx) goto 12
      if(djx.eq.djjx) then
          if(tjx.ge.tjjx) goto 12
          end if
      dict(j)=djjx
      type(j)=tjjx
      j=jj
      jj=j*2
      goto 15
      end

      subroutine getd(iunit,samp,spec,time,iend,
     1                  blurb,blnk80,b,word,blnk10,w)
C Gets data from iunit, makes first 3 words: samp spec time
      character*10 samp,spec,word,time,blnk10
      character*80 blurb,blnk80
      character*1 b(80),blank/' '/,w(10)
      iend=0
      blurb=blnk80
      read(iunit,1000,end=97) blurb
      goto 98
   97 if(blurb.eq.blnk80) goto 999
   98 continue
 1000 format(a80)
      do k=1,80
         if(b(k).ne.blank) goto 100
         enddo
      return
  100 continue
      k1=k
      do k=k1,80
         if(b(k).eq.blank) goto 200
         enddo
      return
  200 continue
      k2=k-1
      do k=k2+2,80
         if(b(k).ne.blank) goto 300
         enddo
      return
  300 continue
      k3=k
      do k=k3,80
         if(b(k).eq.blank) goto 400
         enddo
      return
  400 continue
      k4=k-1
      do k=k4+2,80
         if(b(k).ne.blank) goto 500
         enddo
      if(k.gt.79) k=80
  500 continue
      k5=k
      do k=k5,80
         if(b(k).eq.blank) goto 600
         enddo
  600 continue
      k6=k-1
      if(k5.eq.80) k6=80
      word=blnk10
      do k=k1,k2
         if(k.lt.k1+10) w(1+k-k1)=b(k)
         enddo
      samp=word
      word=blnk10
      do k=k3,k4
         if(k.lt.k3+10) w(1+k-k3)=b(k)
         enddo
      spec=word
      word=blnk10
      do k=k5,k6
         if(k.lt.k5+9) w(1+k-k5)=b(k)
         enddo
      time=word
      return

  999 continue
      iend=1
      return
      end

      subroutine getnum(weight,wgt,word,w)
C Given a number weight, reads wgt as a real number
      character*10 weight,word
      character*1 w(10),dot/'.'/,blank/' '/
      real wgt
      word=weight
      idot=0
      do k=1,10
         if(w(k).eq.dot) idot=1
         if(w(k).eq.blank) goto 100
         enddo
  100 continue
      if(idot.eq.0) w(k)=dot
      read(word,2991,err=200) wgt
 2991 format(f10.4)
      return
  200 continue
      wgt=0
      return
      end

      subroutine addwrd(sa,mm,m,samp)
C checks samp against a list of words sa(mm), and adds samp to the list if new
      character*10 sa(mm),samp
      if(m.eq.0) then
         i=0
         else
         call binfnd(sa,m,samp,i)
         end if
      if(i.ne.0) return 
      m=m+1
      if(m.gt.mm) return
      sa(m)=samp
      call sort10(sa,m)
      return
      end

        subroutine binfnd(ma,n,na,i)
C Searches for item na in array ma (length n), index of na is i
C i=0 if na cannot be found
        character*10 ma(n),na,iamin,iamid,iamax
        i=0
        imin=1
        iamin=ma(imin)
        imax=n
        iamax=ma(imax)
   10   continue
        if(imax-imin.le.1) then
                if(iamin.eq.na) then
                        i=imin
                        return
                        end if
                if(iamax.eq.na) i=imax
                return
                end if
        imid=(imax+imin)/2
        iamid=ma(imid)
        if (na.le.iamid) then
                imax=imid
                iamax=iamid
                goto 10
                end if
        imin=imid
        iamin=iamid
        goto 10
        end

      subroutine sort10(dict,n)
C To sort character*10 data
      character*10 dict(n),djx,djjx,djjjx
      integer i,j,jj,jjj
      do 10 i=1,n
      j=i
      djx=dict(j)
    5 if(j.eq.1) goto 8
      jj=j/2
      djjx=dict(jj)
      if(djjx.gt.djx) goto 8
      if(djjx.eq.djx) then
           end if
      dict(j)=djjx
      j=jj
      goto 5
    8 dict(j)=djx
   10 continue
      i=n
      goto 14
   12 dict(j)=djx
   14 if(i.eq.1) return
      djx=dict(i)
      dict(i)=dict(1)
      i=i-1
      j=1
      jj=2
   15 if(i-jj) 12,17,16
   16 djjx=dict(jj)
      jjj=jj+1
      djjjx=dict(jjj)
      if(djjx.ge.djjjx) goto 18
      if(djx.ge.djjjx) goto 12
      dict(j)=djjjx
      j=jjj
      jj=j*2
      goto 15
   17 djjx=dict(jj)
   18 if(djx.ge.djjx) goto 12
      dict(j)=djjx
      j=jj
      jj=j*2
      goto 15
      end

      subroutine sort30(dict,n)
C To sort character*30 data
      character*30 dict(n),djx,djjx,djjjx
      integer i,j,jj,jjj
      do 10 i=1,n
      j=i
      djx=dict(j)
    5 if(j.eq.1) goto 8
      jj=j/2
      djjx=dict(jj)
      if(djjx.gt.djx) goto 8
      if(djjx.eq.djx) then
           end if
      dict(j)=djjx
      j=jj
      goto 5
    8 dict(j)=djx
   10 continue
      i=n
      goto 14
   12 dict(j)=djx
   14 if(i.eq.1) return
      djx=dict(i)
      dict(i)=dict(1)
      i=i-1
      j=1
      jj=2
   15 if(i-jj) 12,17,16
   16 djjx=dict(jj)
      jjj=jj+1
      djjjx=dict(jjj)
      if(djjx.ge.djjjx) goto 18
      if(djx.ge.djjjx) goto 12
      dict(j)=djjjx
      j=jjj
      jj=j*2
      goto 15
   17 djjx=dict(jj)
   18 if(djx.ge.djjx) goto 12
      dict(j)=djjx
      j=jj
      jj=j*2
      goto 15
      end

      subroutine isort(dict,n)
C To sort integers
      integer dict(n),djx,djjx,djjjx
      integer i,j,jj,jjj
      do 10 i=1,n
      j=i
      djx=dict(j)
    5 if(j.eq.1) goto 8
      jj=j/2
      djjx=dict(jj)
      if(djjx.gt.djx) goto 8
      if(djjx.eq.djx) then
           end if
      dict(j)=djjx
      j=jj
      goto 5
    8 dict(j)=djx
   10 continue
      i=n
      goto 14
   12 dict(j)=djx
   14 if(i.eq.1) return
      djx=dict(i)
      dict(i)=dict(1)
      i=i-1
      j=1
      jj=2
   15 if(i-jj) 12,17,16
   16 djjx=dict(jj)
      jjj=jj+1
      djjjx=dict(jjj)
      if(djjx.ge.djjjx) goto 18
      if(djx.ge.djjjx) goto 12
      dict(j)=djjjx
      j=jjj
      jj=j*2
      goto 15
   17 djjx=dict(jj)
   18 if(djx.ge.djjx) goto 12
      dict(j)=djjx
      j=jj
      jj=j*2
      goto 15
      end

      subroutine hold
C Asks for user input before exiting, to retain DOS window
      character*1 ch
      write(*,2000)
 2000 format(//'Press <RETURN> to exit'///
     2         '----------------------')
      read(*,1000) ch
 1000 format(a1)
      stop
      end

      subroutine filout(filein,iunit)
      character*20 filein
C Checks whether filein available for writing or already exists
  100 continue
      read(*,1000) filein
 1000 format(a20)
      open(unit=iunit,file=filein,status='new',err=999)
      return
  999 continue
      write(*,2000)
 2000 format(/'  *** ERROR *** File already exists'/
     1        ' Type another name')
      goto 100
      end

      subroutine filin(filein,iunit)
      character*20 filein
C Checks whether filein available for reading
  100 continue
      read(*,1000) filein
 1000 format(a20)
      open(unit=iunit,file=filein,status='old',err=999)
      return
  999 continue
      write(*,2000)
 2000 format(/'  *** ERROR *** File does not exist'/
     1        ' Type another name')
      goto 100
      end
