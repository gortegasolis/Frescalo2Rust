C A program to calculate trend values

      integer m,n,nt
      parameter(mm=4000,nn=2000,nnt=100,nndat=2000000,nnwgt=500000)
C mm is max number of samples; nn of species; nnt of time periods
C nndat is size of main data matrix; nnwgt is size of neighbourhood weight data
      parameter(fmax=0.99999,fmin=1.0E-10,tol=0.0003,irepmx=100)
C fmax is a value less than 1 to ensure that logs of (1-f) can be taken
C fmin is used to ensure that frequencies do not sum to zero in small datasets
C tol is the convergence limit for rescaling;  irepmx is the maximum number of iterations
      parameter(phidef=0.74,blmdef=0.2703)
C phidef is default value of phi, the local frequency
C blmdef is default limit for benchmark spp.
      character*30 ddijt(nndat),ddjti(nndat),dtji,wneigh(nnwgt)
      character*20 filein,fileou,blnk20/'                    '/,filsam
      character*10 samp,samp1,spec,spec1,blnk10/'          '/,word,time
      character*10 sa(mm),sp(nn),tim(nnt),d(3),weight,bnchx(nn)
C bnchx will contain the list of unsuitable bench species
      integer numwgt(mm),idat(mm),jrank(nn),iitot(mm)
      integer lendat(nn,nnt),iocc(mm),jocc(nn),ibench(mm,nn)
      integer idata(mm,nn)
      real timf(nn,nnt),sdtimf(nn,nnt),sampef(mm,nnt),ffij(mm,nn)
      real smpint(mm),fff(mm),wgttot(mm),wgtt2(mm),f(nn),ff(nn),ffff(mm)
      real blim,abtot(mm),bwght(nn)
C bwght are weights, set to 0.001 for unsuitable bench species, 1 for the rest
      character*80 blurb,blnk80
      character*1 b(80),blank/' '/,w(10)
      equivalence (b,blurb),(w,word),(d,dtji),(ddijt,ddjti)

      do k=1,80
         b(k)=blank
         enddo
      blnk80=blurb
      samp1=blnk10
      do iiit=1,nnt
         do ii=1,mm
            sampef(ii,iiit)=1.0E-7
            enddo
         do jj=1,nn
            timf(jj,iiit)=0
            sdtimf(jj,iiit)=0
            lendat(jj,iiit)=0
            enddo
         enddo
      do ii=1,mm
         wgttot(ii)=0
         wgtt2(ii)=0
         numwgt(ii)=0
         idat(ii)=0
         iitot(ii)=0
         do jj=1,nn
            ffij(ii,jj)=0
            ibench(ii,jj)=0
            enddo
         enddo
      do jj=1,nn
         bwght(jj)=1
         enddo

C Set up files for reading and writing
      write(*,1998) mm,nn,nnt,nndat,nnwgt
 1998 format(//1x,'FRESCALO - Trend_analysis using local frequencies'/
     1       1x,'written by Mark Hill, January-June 2011'//
     2       1x,'Input limits:'/
     3       4x,'Number of samples = ', I7 /
     4       4x,'Number of species = ', I7 /
     5       4x,'Number of time periods = ', I7 / 
     6       4x,'Number of observations = ', I7 /
     7       4x,'Number of neighbourhood weights = ', I7 //
     8       1x,'Type name of log file ...')
      call filout(fileou,10)
      write(10,1999) mm,nn,nnt,nndat,nnwgt
 1999 format(1x,'Log file for FRESCALO'//
     2       1x,'Input limits:'/
     3       4x,'Number of samples = ', I7 /
     4       4x,'Number of species = ', I7 /
     5       4x,'Number of time periods = ', I7 / 
     6       4x,'Number of observations = ', I7 /
     7       4x,'Number of neighbourhood weights = ', I7 //
     8       1x,'Type name of log file ...')
      write(10,1000) fileou
 1000 format(a20)
      write(*,2000)
      write(10,2000)
 2000 format(1x,'Type occurrence input file [sample species time] ...')
      call filin(filein,4)
      write(10,1000) filein
      write(*,2500)
      write(10,2500)
 2500 format(1x,'Type neighbourhood weight input file ...')
      call filin(filein,3)
      write(10,1000) filein
      write(*,2501)
      write(10,2501)
   20 continue
 2501 format(1x,'Type file with species to exclude from benchmarks '/
     1          '     or press <RETURN> if no exclusions...')
      read(*,1000) filein
      nbnchx=0
      if(filein.eq.blnk20) then
         write(10,2502)
 2502 format('<No exclusions>')
         goto 30
         end if
      open(unit=2,file=filein,status='old',err=20)
      write(10,1000) filein
      do jj=1,nn
         read(2,1001,end=25) spec
 1001 format(a10)
         nbnchx=jj
         bnchx(nbnchx)=spec
         enddo
   25 continue
   30 continue

      write(*,2001)
      write(10,2001)
 2001 format(1x,'Type name of sample stats output file...')
      call filout(fileou,7)
      write(10,1000) fileou
      write(*,2002)
      write(10,2002)
      filsam=fileou
 2002 format(1x,'Type name of rescaled frequency file...')
      call filout(fileou,8)
      write(10,1000) fileou
      write(*,2004)
      write(10,2004)
 2004 format(1x,'Type name of trend output file...')
      call filout(fileou,9)
      write(10,1000) fileou

      phibig=phidef
   40 continue
      write(*,2009) phidef
      write(10,2009) phidef
 2009 format(1x,'Type target value of local frequency phi ',
     1          '(default=',f4.2,')...')
      read(*,2010,err=45) phibig
 2010 format(f8.4)
   45 continue
      if(phibig.eq.0) phibig=phidef
      if(phibig.gt.0.95.or.phibig.lt.0.50) then
         write(*,2012)
 2012 format(1x,'***ERROR*** Outside range 0.50 to 0.95')
         goto 40
         end if
      write(*,2011) phibig
      write(10,2011) phibig
 2011 format(1x,'Target value is ',f4.2/)

      blim=blmdef
   46 continue
      write(*,2007) blmdef
      write(10,2007) blmdef
 2007 format(1x,'Type value of Benchmark Limit ',
     1        '(default=',f4.2,')...')
      read(*,2010,err=47) blim
   47 continue
      if(blim.eq.0) blim=blmdef
      if(blim.gt.0.5.or.blim.lt.0.08) then
         write(*,2013)
 2013 format(1x,'***ERROR*** Outside range 0.08 to 0.5')
         goto 46
         end if
      write(*,2014) blim
      write(10,2014) blim
 2014 format(1x,'Benchmark limit is ',f4.2/)

C Read in data using subroutine getd
      m=0
      n=0
      nt=0
      ndjti=0
      nwgt=0

   52 continue
C Reads hectad weights(samp,samp1,weight), stores them in wneigh
      if(nwgt.eq.0) write(*,*) 'Reading in smoothing ',
     1                           'weights from samples'
      call getd(3,samp,samp1,weight,iend,blurb,blnk80,b,word,blnk10,w)
      if(iend.eq.1) goto 53
      nwgt=nwgt+1
      if(mod(nwgt,20000).eq.0) write(*,*) 'Weights ',samp,samp1,nwgt
      if(nwgt.gt.nnwgt) then
         write(*,2995) nnwgt
         write(10,2995) nnwgt
 2995 format(1x,'No of neighbourhood weights > maximum which is',i5)
         call hold
         end if
      call addwrd(sa,mm,m,samp)
      call addwrd(sa,mm,m,samp1)
      d(1)=samp1
      d(2)=samp
C This order so that all samp of which samp1 is a neighbour are together
      d(3)=weight
      wneigh(nwgt)=dtji
C This has stored the weight data in wneigh (character*30)
      goto 52

   53 continue
      write(*,*) 'Sorting local frequency weights ...'
      call sort30(wneigh,nwgt)
C We now have weights sorted in wneigh, ready for use
      
  100 continue
      call getd(4,samp,spec,time,iend,blurb,blnk80,b,word,blnk10,w)
      if(iend.eq.1) goto 101
      call binfnd(sa,m,samp,i)
      if(i.eq.0) then
         if(samp.ne.samp1) write(10,2996) samp
 2996 format('*** ',a10,' location ignored - in species data but not ',
     1           'listed in neighbourhood weights')
         samp1=samp
         goto 100
         end if
C This omits all data from samples that do not have weights defined
      call addwrd(sp,nn,n,spec)
      call addwrd(tim,nnt,nt,time)
      if(nt.gt.nnt) then
         write(*,2997) nnt
         write(10,2997) nnt
 2997 format(1x,'Number of time periods exceeds maximum which is',i5)
         call hold
         end if
      ndtji=ndtji+1
      if(mod(ndtji,20000).eq.0) write(*,*) samp, spec, time, ndtji
      d(1)=samp
      d(2)=spec
      d(3)=time
      ddijt(ndtji)=dtji
      goto 100

  101 continue
      write(*,2505) m,n,nt
      write(*,2506) ndtji,nwgt,nbnchx
      write(10,2505) m,n,nt
      write(10,2506) ndtji,nwgt,nbnchx
 2505 format(/1x,'Actual numbers in data'/
     1       4x,'Number of samples      ',i8/
     2       4x,'Number of species      ',i8/
     3       4x,'Number of time periods ',i8)
 2506 format(4x,'Number of observations ',i8/
     1       4x,'Neighbourhood weights  ',i8/
     2       4x,'Benchmark exclusions   ',i8/)
      if(nbnchx.gt.0) write(10,2507)
 2507 format(1x,'Benchmark exclusions')
      if(nbnchx.gt.0) write(10,2508)(bnchx(ibnchx),ibnchx=1,nbnchx)
 2508 format(4x,a10)
      write(*,*) 'Sorting main data ...'
      call sort30(ddijt,ndtji)
C The main dataset is sorted and in ddijt; next to calculate local frequencies
C First calculate sum of weights for each sample
      do iwgt=1,nwgt
         dtji=wneigh(iwgt)
         samp1=d(1)
         samp=d(2)
         weight=d(3)
         call binfnd(sa,m,samp,i)
         call binfnd(sa,m,samp1,ii)
         call getnum(weight,wgt,word,w)
         wgttot(i)=wgttot(i)+wgt
         wgtt2(i)=wgtt2(i)+wgt*wgt
         numwgt(ii)=numwgt(ii)+1
         enddo

C Now calculate the number of data items for each sample in main data
      spec1=blnk10
      do idtji=1,ndtji
         dtji=ddijt(idtji)
         samp=d(1)
         spec=d(2)
         call binfnd(sa,m,samp,i)
         idat(i)=idat(i)+1
         if(spec.ne.spec1) then
            iitot(i)=iitot(i)+1
            spec1=spec
            end if
         enddo

C Now calculate frequencies
      idtji=0
      iwgt=0
      do ii=1,m
         if(mod(ii,100).eq.0) write(*,*) 'frequencies ',sa(ii),ii
         do j=1,n
            jocc(j)=0
            idata(ii,j)=0
            enddo
         do iid=1,idat(ii)
            idtji=idtji+1
            dtji=ddijt(idtji)
            spec=d(2)
            call binfnd(sp,n,spec,j)
            jocc(j)=1
            idata(ii,j)=1
            enddo
         do iiw=1,numwgt(ii)
            iwgt=iwgt+1
            dtji=wneigh(iwgt)
            samp=d(2)
            weight=d(3)
            call binfnd(sa,m,samp,i)
            call getnum(weight,wgt,word,w)
               do j=1,n
               ffij(i,j)=ffij(i,j)+jocc(j)*wgt/(wgttot(i)+1.0E-10)
               enddo
            enddo
         enddo
C This completes the calculation of local frequencies, which are now in ffij(i,j)
C Now downweight species that are unsuitable as benchmarks
      do ibnchx=1,nbnchx
         spec=bnchx(ibnchx)
         call binfnd(sp,n,spec,j)
         if(j.ne.0) bwght(j)=0.001
         enddo

      do i=1,m
         if(mod(ii,100).eq.0) write(*,*) 'rescaling ',sa(ii),ii
         do j=1,n
            f(j)=ffij(i,j)
            jocc(j)=idata(i,j)
            enddo
         samp=sa(i)
         itot=iitot(i)
         wn2=wgttot(i)*wgttot(i)/(wgtt2(i)+1.0E-12)
C wn2 is the Effective Number of weights in the neighbourhood
         call fresca(i,n,itot,jocc,f,ff,jrank,samp,sp,phi1,
     1            phibig,fmax,fmin,wn2,spnum,tol,irepmx)
         ffff(i)=phi1
         abtot(i)=1.0E-7
         do j=1,n
            if(ffij(i,j).ne.0) ffij(i,j)=ff(j)
            jj=jrank(j)
            rank1=j/spnum
            if(rank1.lt.blim.or.j.eq.1)then
C The case j=1 is put in because with small samples rank1 may be greater than 0.2 for j=1
               ibench(i,jj)=1
               abtot(i)=abtot(i)+bwght(jj)
               else
               ibench(i,jj)=0
               end if
            enddo
         enddo

C The next task is to reorder the main data matrix
      do idtji=1,ndtji
         dtji=ddijt(idtji)
         samp=d(1)
         spec=d(2)
         time=d(3)
         if(mod(idtji,20000).eq.0) write(*,*) 'Re-ordering ',
     1                             samp, spec, time, idtji
         d(1)=spec
         d(2)=time
         d(3)=samp
         ddijt(idtji)=dtji
         enddo
      write(*,*) 'Now doing second sort of reordered main data ...'
      call sort30(ddijt,ndtji)
C From now onwards, the sorted form is called ddjti - it is the same vector equivalenced

      do idtji=1,ndtji
         if(mod(idtji,20000).eq.0) write(*,*) 
     1      'Main data to calc sampling effort',ddjti(idtji),idtji
         dtji=ddjti(idtji)
         spec=d(1)
         time=d(2)
         samp=d(3)
         call binfnd(sa,m,samp,i)
         call binfnd(sp,n,spec,j)
         call binfnd(tim,nt,time,iit)
         lendat(j,iit)=lendat(j,iit)+1
         sampef(i,iit)=sampef(i,iit)+ibench(i,j)*bwght(j)/abtot(i)
         enddo

      idtji=0
      jx=1
      iitx=1
C These are markers to ensure that we print out time periods at which species not recorded
      write(9,2060)
 2060 format('Species__  Time______ TFactor St_Dev _Count ___spt',
     1       ' ___est N>0.00 N>0.98')

  130 continue
      do i=1,m
         iocc(i)=0
         enddo     
      idtji=idtji+1
      if(idtji.gt.ndtji) goto 140
      dtji=ddjti(idtji)
      spec=d(1)
      time=d(2)
      samp=d(3)
      call binfnd(sp,n,spec,j)
      call binfnd(tim,nt,time,iit)
      do i=1,m
         smpint(i)= sampef(i,iit)
         fff(i) = ffij(i,j)
         enddo
      call binfnd(sa,m,samp,i)
      iocc(i)=1
      do id=1,lendat(j,iit)-1
         idtji=idtji+1
         if(idtji.gt.ndtji) goto 140
         dtji=ddjti(idtji)
         samp=d(3)
         call binfnd(sa,m,samp,i)
         iocc(i)=1
         enddo
         call tfcalc(tf,sd,spt,jtot,est,iocc,smpint,fff,m,ic1,ic2)
C first we pad out the output with time when species was not recorded
         if(sp(jx).lt.sp(j)) then
            do iiit=iitx,nt
               write(9,2050) sp(jx),tim(iiit),0.0,0.0,0,0.0,0.0,0,0
               enddo
            do jj=jx+1,j-1
               do iiit=1,nt
                  write(9,2050) sp(jj),tim(iiit),0.0,0.0,0,0.0,0.0,0,0
                  enddo
               enddo
            jx=j
            iitx=1                 
            end if
         if(tim(iitx).lt.tim(iit)) then
            do iiit=iitx,iit-1   
               write(9,2050) sp(j),tim(iiit),0.0,0.0,0,0.0,0.0,0,0
               enddo
            iitx=iit
            end if
         write(9,2050) sp(j),tim(iit),tf,sd,jtot,spt,est,ic1,ic2
         if(mod(j,10).eq.0) write(*,2050) sp(j),tim(iit),tf,sd,j
 2050 format(a10,1x,a10,f8.3,f7.3,i7,2f7.1,2i7)
         jx=j
         iitx=iit+1
      
      goto 130

  140 continue
      do iiit=iitx,nt
         write(9,2050) sp(jx),tim(iiit),0.0,0.0,0,0.0,0.0,0,0
         enddo

      close(7)
      close(8)
      close(9)

C Finally test whether given value of phi appears to be unrealistically low
      call sort(ffff,m)
      i985=0.985*m
      phi985=ffff(i985)
      write(*,2513) phi985,phibig
      write(10,2513) phi985,phibig
 2513 format(//1x,'98.5 percentile of input phi ',f5.2/
     1       1x,'Target value of phi          ',f5.2)
      if(phibig.lt.phi985) then
         write(*,2514)
         write(10,2514)
         end if
 2514 format(///1x,'*** BEWARE *** '/
     1          1x,'Target value of phi may be too small'/)
      write(*,2503)
      write(10,2503)
 2503 format(/1x,'Calculation reached completion'/)
      close(10)

      call hold
      end

      subroutine fresca(m,n,itot,jocc,f,ff,jrank,samp1,splist,phi1,
     1                  phibig,fmax,fmin,wn2,spnum,tol,irepmx)
      real f(n),ff(n)
      integer jrank(n),jocc(n)
      character*10 splist(n),samp1
C converts the rank-frequency data f to rescaled values ff
      alpha=1
      do ir=1,irepmx
         do j=1,n
            if(f(j).gt.fmax) f(j)=fmax
            if(f(j).lt.fmin) f(j)=fmin
            ff(j)=-log(1-f(j))
            enddo
         tot=0
         tot2=0
         do j=1,n
            ffij=1-exp(-ff(j)*alpha)
C ffij is the new frequency after recording intensity multiplied by alpha
            tot=tot+ffij
            tot2=tot2+ffij**2
            enddo
         phi=tot2/tot
         an2=tot*tot/tot2
         spnum=tot
         if(ir.lt.20) then
            alpha=alpha*exp(1.86*(log(1-phi)-log(1-phibig)))
C Successive approx based on linear relation - fails with some small datasets
            else
            alpha=alpha*phibig/phi
C Crude successive approximation - slower with big datasets
            end if
         if(ir.eq.1) then
            phi1=phi
            an21=an2
            spnum1=tot
            end if
         if(abs(phi-phibig).lt.tol) goto 200
         enddo
  200 continue
      if(m.eq.1) write(7,2001)
 2001 format('Location',2x,'Loc_no ',' No_spp',' Phi_in','  Alpha',
     1       '  Wgt_n2',' Phi_out','  Spnum_in',' Spnum_out',' Iter')
      alph=alpha
      if(alpha.gt.999.99) alph=999.99
      write(7,2002) samp1,m,itot,phi1,alph,wn2,phi,spnum1,spnum,ir
 2002 format(a10,2i7,f7.3,1x,f6.2,1x,f7.2,1x,f7.3,2f10.1,i5)
      do j=1,n
         ff(j)=-f(j)+j*1.0E-12
         jrank(j)=j
         enddo
      call sort2(ff,jrank,n)
      do j=1,n
         jj=jrank(j)
         fij=f(jj)
         if(fij.gt.fmax) fij=fmax
         ffij=1-exp(alpha*log(1-fij))
         sdfij=sqrt(fij*(1-fij)/wn2)
         ffff=fij+sdfij
         if(1-ffff.lt.1.0E-12) ffff=1-1.0E-12
         fffff=fij-sdfij
         ffsd=1-exp(alpha*log(1-ffff))
         fffsd=1-exp(alpha*log(1-fffff))
         sdij=0.5*(ffsd-fffsd)
C sdij is an estimate of the standard error
         if(j.eq.1.and.m.eq.1) write(8,2003)
 2003 format('Location  ',' Species   ',' Pres','  Freq__','  Freq_1',
     1       ' SD_Frq1','  Rank','  Rank_1')
         ff(jj)=ffij
         if(f(jj).gt. 0.00005) write(8,2004) 
     1      samp1,splist(jj),jocc(jj),f(jj),ffij,sdij,j,j/spnum
 2004 format(a10,1x,a10,1x,i4,3f8.4,1x,i5,f8.3)
         enddo
      return
      end

      subroutine tfcalc(tf,sd,sptot,jtot,esttot,iocc,smpint,
     1                  fff,m,ic1,ic2)
      integer iocc(m)
      real smpint(m),fff(m)
C Calculates a time factor tf, given observed species total sptot for that time
C sd is the standard error
C iocc(i) is 1 if species found at location i at the time, 0 otherwise
C smpint(i) is the sampling intensity at location i and time t
C fff(i) is the smoothed time-independent frequency of species at location i
C ic1 is the number of samples for which there is nonzero probability of occurrence
C ic2 is the number of cases where smptint*fff >0.98
C weights are applied, downweighting cases where smpint<0.1 (no systematic sampling)
      kmax=100
      tf=1
      do k=1,kmax
         esttot=0
         estvar=0
         sptot=0
         jtot=0
         ic1=0
         ic2=0
         do i=1,m
            wgt=1
            if(smpint(i).lt.0.0995) wgt=10*smpint(i)+0.005
            pfac=smpint(i)*fff(i)
C the probability of finding species i is multiplied by the sampling intensity
            if(pfac.gt.0) ic1=ic1+1
            if(pfac.gt.0.98)then
               pfac=0.98
               ic2=ic2+1
               end if                   
C this is a fudge to allow the next step of taking logs; otherwise there is a danger that pfac=1
            plog=-log(1-pfac)
            estval=1-exp(-plog*tf)
            esttot=esttot+wgt*estval
            estvar=estvar+wgt*wgt*estval*(1-estval)
            sptot=sptot+wgt*iocc(i)
            jtot=jtot+iocc(i)
C esttot is the estimated total for species j at time iit
            enddo
         if(abs(sptot-esttot).lt.0.0005) goto 100
         tf=tf*sptot/(esttot+0.0000001)
         enddo

  100 continue
      sptot1=sptot+sqrt(estvar)
C sptot1 is precisely 1 standard deviation bigger than sptot
C recalculate with this target value to obtain standard error of tf

      tf1=tf
      do k=1,kmax
         esttt1=0
         do i=1,m
            wgt=1
            if(smpint(i).lt.0.0995) wgt=10*smpint(i)+0.005
            pfac=smpint(i)*fff(i)
C the probability of finding species i is multiplied by the sampling intensity
            if(pfac.gt.0.98)then
               pfac=0.98
               end if                   
C this is a fudge to allow the next step of taking logs; otherwise there is a danger that pfac=1
            plog=-log(1-pfac)
            estval=1-exp(-plog*tf1)
            esttt1=esttt1+wgt*estval
C esttt1 is the estimated total targeted at 1 sd bigger
            enddo
         if(abs(sptot1-esttt1).lt.0.0005) goto 200
         tf1=tf1*sptot1/(esttt1+0.0000001)
         enddo

  200 continue
      sd=tf1-tf
      return
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
      sa(m)=samp
      call sort10(sa,m)
      return
      end

      subroutine sort(dict,n)
C To sort real numbers
      real dict(n),djx,djjx,djjjx
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
C To sort a character*30 list
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

      subroutine sort10(dict,n)
C To sort a character*10 list
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

      subroutine sort2(dict,type,n)
C To sort two columns at once
C The integers in type(j) are sorted by the real numbers in dict(j)
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
