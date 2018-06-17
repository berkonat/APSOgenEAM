c
c     Copyright (c) 2012 Berk Onat
c     View LICENSE file in this directory for the license of this code.
c
c     EAM File Generation Program in DYNAMO setfl or funcfl format
c
c     Inputs with < : filename
c
c     Number of alloy pairs = n(n-1)/2 
c     Implemented up to n=3 types,  alloy pair num = 1 or 3
c
c
 

      program apsogeneam

      include "mkl_vsl.f77"
      include "mpif.h"

      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      parameter (psopar=1000,psoprm=100,psorand=10000,
     . psornd=500000,rndprm=300,eamprm=100)
      parameter (MASTER = 0)
      parameter (FROM_MASTER = 1)
      parameter (FROM_WORKER = 2)
      character*32 dateh
      character*12 dateread
      character*14 statechar
      character*80 eamfiles(nelmax)
      character*80 eamfileout
      character*80 logfilestr
      character*80 logfilestr2
      character*12 logfile
      character*80 fitvalues
      character*80 expvalues
      character*80 weights
      character*4 myidstr
      integer readpureeams,pairnum
      common /particle/ rv(6,nmax)
      common /forces/ f(3,nmax),fe(nmax),e(nmax),stresst(3,3),
     .                slocal(3,3,nmax),phie(nmax)
      common /density/ rho(nmax),fp(nmax)
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      common /interact/ frho(ngrid,nelmax),z2r(ngrid,nelmax,nelmax),
     .  rhor(ngrid,nelmax),drho,dr,rcutsq,nrho,nr
      common /grids/ frhoar(ngridar,nelmax),frhoar1(ngridar,nelmax),
     .  frhoar2(ngridar,nelmax),frhoar3(ngridar,nelmax),
     .  frhoar4(ngridar,nelmax),frhoar5(ngridar,nelmax),
     .  frhoar6(ngridar,nelmax),rhorar(ngridar,nelmax),
     .  rhorar1(ngridar,nelmax),rhorar2(ngridar,nelmax),
     .  rhorar3(ngridar,nelmax),rhorar4(ngridar,nelmax),
     .  rhorar5(ngridar,nelmax),rhorar6(ngridar,nelmax),
     .  z2rar(ngridar,nelmax,nelmax),z2rar1(ngridar,nelmax,nelmax),
     .  z2rar2(ngridar,nelmax,nelmax),z2rar3(ngridar,nelmax,nelmax),
     .  z2rar4(ngridar,nelmax,nelmax),z2rar5(ngridar,nelmax,nelmax),
     .  z2rar6(ngridar,nelmax,nelmax),drhoar,drar,nrhoar,nrar
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /eamdef/ fheader(80,nelmax),rcut(nelmax),blat(nelmax),
     . lat(nelmax),elmname
      common /inputdef/ nxchain,nyatoms,nzlayer,rnum,rhonum
      common /inputlat/ rhoend,lata0,caratio,latai,latstep,lataf,
     1esub,bmod
      integer nxchain,nyatoms,nzlayer,rnum,rhonum
      double precision lata0,caratio,rhoend
      character elmname(nelmax)*2
      double precision gridpack(ngrid*nelmax)
      integer gridpacknum
      dimension frhoin(ngrid,nelmax),drhoin(nelmax),nrhoin(nelmax)
      dimension rhorin(ngrid,nelmax),zrin(ngrid,nelmax),drin(nelmax),
     1 nrin(nelmax)
      dimension sfrhoin(ngrid,nelmax),sdrhoin(nelmax),snrhoin(nelmax)
      dimension srhorin(ngrid,nelmax),szrin(ngrid,nelmax),sdrin(nelmax),
     1 snrin(nelmax)
      double precision sx(nelmax),gx(nelmax)
      double precision sxup(nelmax),gxup(nelmax)
      double precision sxdw(nelmax),gxdw(nelmax)
      double precision gs,ga,sa,gb,sb,gc,sc
      double precision alpha1,beta1,dm1,rm1,rcutguess1,de,re
      double precision alpha2,beta2,dm2,rm2,rcutguess2
      double precision alpha3,beta3,dm3,rm3,rcutguess3
      double precision latai,lataf,latstep,esub,bmod
      double precision alpha1up,beta1up,dm1up,rm1up,rcutguess1up
      double precision alpha2up,beta2up,dm2up,rm2up,rcutguess2up
      double precision alpha3up,beta3up,dm3up,rm3up,rcutguess3up
      double precision alpha1dw,beta1dw,dm1dw,rm1dw,rcutguess1dw
      double precision alpha2dw,beta2dw,dm2dw,rm2dw,rcutguess2dw
      double precision alpha3dw,beta3dw,dm3dw,rm3dw,rcutguess3dw
      double precision parx(psoprm*psopar),parv(psopar,psoprm)
      double precision bounds(psoprm,2),vmax(psoprm,2)
      double precision wt,evofact,ese(4),esehig,eselow
      double precision c1,c2,c1upd,c2upd,c1c2tot,c1ok,c2ok
      double precision gbest(psoprm),pbest(psopar,psoprm)
      double precision mygbest(psoprm)
      double precision gworse(psoprm),params(eamprm),ppart(psoprm)
      double precision gbestval,mygbestval(2,1),pbestval(psopar)
      double precision testx(psopar)
      double precision gbestoldval,gbestold(psoprm),mygworseval(2,1)
      double precision mygbestoldval,mygbestold(psoprm),gbestpack(2,1)
      double precision partdist(psopar),fval(psopar),gworseval
      double precision gworsepack(2,1),pp
      double precision rand1(psorand),rand2(psorand),dimd
      double precision delta1(psorand),delta2(psorand)
      double precision distmax,distmin,distg,difparm,distpart
      double precision sigma,sigmamax,sigmamin,perturb,gauss,mu
      double precision bounddif,rndsign,funcprm(psoprm),ppartval
      integer deltasign1(psorand),deltasign2(psorand)
      integer*4 nparams
      integer nparts,mynparts
      integer gen,gn,ngens,esestate,evaluate
      integer esenext,esenext2,esenext3,esenext4,eseprev
      integer esehigloc,eselowloc,randdim
      integer gbestid,gworseid,fixcutoff
      integer mygbestid,mygworseid,gbestloc,gworseloc
      integer randchunk,nlast,dims,di(psoprm),npartsin,ngensin
      integer phitype,rhotype,embtype,expnum,ind,inn,elsupd
      integer at1,at2,lattype,weightnum
      integer*4 eamftyp(10),testflag
      integer i,j,k,ii,jj,kk,kkk,seedcount
      integer printon,eamalloy
      integer ftyp,dyntyp,dyntypin
      integer prntflag,prmflag,prmno,logflag
      double precision expdata(100),weightdata(100),rms,tol
      character*8  date
      character*10 time
      character*5  zone
      integer values(8)
      integer lcharnum,dateseed
      double precision testfunc
      integer ceil
      
      integer*8 stream(2)
      integer*8 errcode
      integer*8 brng,seed,method,randnum,randnum1,randnum2
      double precision leftb,rightb,rands(psornd),randn(rndprm)
      double precision seedzigot(10000),randm(rndprm)
      integer seedcells,seedcells1,seedcells2

      integer*4 numprocs,myid,myfileid
      integer*4 mypartsnum(psopar),offsets(psopar)
      integer*4 source,dest,aveparts,extra,ierr
      integer*4 myparts,offset
      integer*4 displs(0:psopar),rcount(0:psopar)
      integer status(MPI_STATUS_SIZE)
      integer*4 mtype,mytype, mytype2
      integer*4 place,blocklen(500),displacement(500)
      double precision buffer(100000),backpack(100)     

     
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr ) 
c     write(*,*)'Myid:',myid
c     write(*,*)myid,'Numprocs:',numprocs
      if (myid .eq. MASTER) then
      write(*,*)'--------------------------------'
      call date_and_time(date,time,zone,values)
      write (*,1)values(3),values(2),values(1),
     &           values(5),values(6),values(7)
      write(dateh,1)INT(values(3)),INT(values(2)),INT(values(1)),
     &              INT(values(5)),INT(values(6)),INT(values(7))
1     format ( '  Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &         i2.2, ':', i2.2, ':', i2.2 )
      write(*,*)'--------------------------------'
      endif
      readpureeams=0

c    electron-volt to joule conversion constant
      ev_joule=0.624150965

   20 format(80a1)
   21 format(2a1)
   30 format(i5,2g15.5,a8)
  150 format(i5)
 9901 format(i5,e24.16,i5,2e24.16)
 9902 format(5e24.16)
      phitype=1
      fixcutoff=1
      rhotype=1
      embtype=1
      lattype=0
      eamalloy=1
      logflag=0
c     -------------------
c      START MASTER ONLY
c     -------------------
      if (myid .eq. MASTER) then
c     ----------------------------------------
c     EAM file parameters     
c
c     For couple Phi cutoff and Rho cutoff use fixcutoff=1
      read(*,*)fixcutoff
c     Phi interaction: Foiles-Baskes-Daw=0 Voter-Chen=1
c     Ref. for F-B-D: PRB,33,7983
c     Ref. for VC: Mater. Res. Proc. Symp. 81 or 
c       Intermatallic Compounds: Vol.e1, Principles Chapter 3
      read(*,*)phitype
c     EAM file format: funcfl=0 setfl=1
      read(*,*)fileformat
c     EAM file header
       if(fileformat.lt.1) then
        eamalloy=0
        write(*,*)' EAM potential generation for PURE element',
     &' with funcfl format'
        if(phitype .ne. 7) then
        write(*,*)' WARNING!: use phitype=7 for funcfl format ',
     &' if you want to use LAMMPS or DYNAMO. LAMMPS/DYNAMO ',
     &' reads only z of Coulomb interation phi=z^2/r from ',
     &' funcfl. This program use phi=z (in other words any) ',
     &' type of pair interaction. phitype=7 converts output ',
     &' to z -> sqrt((z*r)/(27.2*0.529)) which effect the ',
     &' precison. For a better precision use setfl format for ',
     &' pure element files to read in LAMMPS/DYNAMO.'
c     !search for WARNING funcfl words to see the convertion!     
        else
        write(*,*)' WARNING!: using phitype=7 writes corrected ',
     &' LAMMPS or DYNAMO funcfl files. You can not use these ',
     &' files with this program to generate alloy potentials.',
     &' If you want to generate alloy potential use phitype=1 ',
     &' for pure funcfl files first, else the output potential ',
     &' of optimization will be a junk.'
        endif
c
c     Density and Embedding types: 1,2,3...
        read(*,*)rhotype
        read(*,*)embtype
c     Lattice type: 0=FCC 1=BCC 2=HCP 3=DIA(not implemented)
        read(*,*)lattype
        read(*,20)(fheader(j,1),j=1,10)
        lcharnum=1
        do 2 i=1,80
         if(fheader(i,3).ne.' ') then
           lcharnum=i
         endif
2       continue
        if(lcharnum.lt.48) then
      read(dateh,'(32a1)',END=3)(fheader(i,1),i=lcharnum+1,lcharnum+33)
        endif
3       continue
       else
        read(*,*)eamalloy
        if (eamalloy .lt. 1) then
         write(*,*)' EAM potential generation for PURE element',
     &' with setfl format'
c     Density and Embedding types: 1,2,3...
         read(*,*)rhotype
         read(*,*)embtype
c     Lattice type: 0=FCC 1=BCC 2=HCP
         read(*,*)lattype
         read(*,20)(fheader(j,1),j=1,10)
         lcharnum=1
         do 12 i=1,80
          if(fheader(i,3).ne.' ') then
           lcharnum=i
          endif
12       continue
         if(lcharnum.lt.48) then
      read(dateh,'(32a1)',END=13)(fheader(i,1),i=lcharnum+1,lcharnum+33)
         endif
13       continue
        else
         write(*,*)' EAM potential generation for ALLOY setfl format'
c        read(*,20)(fheader(j,1),j=1,80)
c        read(*,20)(fheader(j,2),j=1,80)
         read(*,20)(fheader(j,3),j=1,80)
        endif
       endif
      read(*,*)fitvalues
      read(*,*)expvalues
      open(unit=50,file=expvalues)
      read(50,*)expnum
      read(50,*) (expdata(i),i=1,expnum)
      close(50)
      read(*,*)weights
      open(unit=51,file=weights)
      read(51,*)weightnum
      read(51,*) (weightdata(i),i=1,weightnum)
      close(51)
c     Number of atom types
      if (eamalloy .lt. 1) then
       ntypes=1 
      else
       read(*,*)ntypes
      endif
       if( ntypes .eq. 2 ) then
        read(*,*)at1,at2
        if ((at1 .lt. 1 .or. at2 .lt. 1) .or. (at1 .eq. at2) .or. 
     &  (at1 .gt. at2)) then
          write(*,*)'Wrong Input: You give type1:',at1,' and type2:',
     &at2,'. You can not choose atom types for alloys less',
     &' than 1 or at1 = at2 or at1 > at2. Input must be at1 < at2.'
         call MPI_FINALIZE(ierr)
         stop
        endif
       else if (ntypes .lt. 2 ) then
c     Atoms' element nos and masses
         read(*,*)ielement(1),amass(1)
         at1=0
         at2=0
c      convert mass to DYNAMO/LAMMPS format
       endif
c      at1,at2,phitype,rhotype,embtype,eamalloy,fileformat packing
       eamftyp(1)=ntypes
       eamftyp(2)=at1
       eamftyp(3)=at2
       eamftyp(4)=phitype
       eamftyp(5)=rhotype
       eamftyp(6)=embtype
       eamftyp(7)=eamalloy
       eamftyp(8)=fileformat
       eamftyp(9)=lattype
       eamftyp(10)=fixcutoff
      endif
c     -------------------
c       END MASTER ONLY
c     -------------------
      call MPI_BCAST(eamftyp,10,MPI_INTEGER,MASTER,
     &MPI_COMM_WORLD,ierr)
c      at1,at2,phitype,rhotype,embtype,eamalloy,fileformat unpacking
       ntypes=eamftyp(1)
       at1=eamftyp(2)
       at2=eamftyp(3)
       phitype=eamftyp(4)
       rhotype=eamftyp(5)
       embtype=eamftyp(6)
       eamalloy=eamftyp(7)
       fileformat=eamftyp(8)
       lattype=eamftyp(9)
       fixcutoff=eamftyp(10)
      pairnum=ntypes*(ntypes-1)/2
      if( pairnum .lt. 1 ) pairnum=1


      if( ntypes .gt. 3 ) then
       if (myid .eq. MASTER) then
        write(*,*)' Maximum number of types for alloying elements is 3.'
        write(*,*)' You gave:',ntypes
        write(*,*)' Exiting program...'
       endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_FINALIZE(ierr)
        stop
      elseif (fileformat .eq. 1 .and. ntypes .lt. 2 .and.
     .        eamalloy .eq. 1 ) then
       if (myid .eq. MASTER) then
        write(*,*)' Minimum number of types for alloying elements is 2.'
        write(*,*)' You gave:',ntypes
        write(*,*)' Exiting program...'
       endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_FINALIZE(ierr)
        stop 
      endif
c     -------------------
c      START MASTER ONLY
c     -------------------
      if (myid .eq. MASTER) then
       if (eamalloy .gt. 0) then
c     EAM input file names
        do 9 i=1,ntypes
         read(*,*)eamfiles(i)
    9   continue
       endif
c     EAM output file name
      read(*,*)eamfileout
c     FCC structure input parameters
c     Number of chains in x direction
      read(*,*)nxchain
c     Number of atoms in y direction
      read(*,*)nyatoms
c     Number of layers in z direction
      read(*,*)nzlayer
      vacatom=(nxchain*nyatoms)*((nzlayer/2)-1)+(nxchain/2)*nyatoms
c     print*,'vacatom:',vacatom
      if (eamalloy .lt. 1) then
c     ----------------------------------------
c     Rose Universal Binding Energy lattice
c     constant perturbation parameters
c
c     latai : initial (a) lattice length
      read(*,*)latai
c     latstep : number of steps
      read(*,*)latstep
c     lataf : final (a) lattice length
      read(*,*)lataf
c     choosen atom number for vacancy formation
c     read(*,*)vacatom
c     vacatom=(nxchain*nyatoms)*((nzlayer/2)-1)+(nxchain/2)*nyatoms
c     Sublimation energy
      read(*,*)esub
c     Bulk modulus
      read(*,*)bmod
c     ----------------------------------------
c
      endif
c     Lattice constant
      read(*,*)lata0
      read(*,*)caratio
c     Number of points for pair (alloy) potential
c      number of points in r(A) for phi(r) spline
      read(*,*)rnum
      
      if(eamalloy .lt. 1) then
       if (embtype .lt. 2) then
c      Embedding function TYPE 1 is a spline derived from Rose EOS
c       number of points in rho(A^-3) for F(rho) spline
        read(*,*)rhonum
c       last rho point for rho(A^-3)
        read(*,*)rhoend
       endif
      endif
c     Number of particles and generations for APSO
      read(*,*)npartsin
      read(*,*)ngensin
c     if dyntyp=1, pbests, gbest evaluate in each generation of APSO
c     For dynamic minimum in systems use dyntyp=1
      read(*,*)dyntypin
      read(*,*)prntflag
      read(*,*)logflag
      if(eamalloy .lt. 1) then
       if (rhotype .lt. 2) then
c      Charge Density Function TYPE 1 is Hydrogenic 4s orbital
        read(*,*)beta1dw,beta1up
        if (fixcutoff.lt.1) then
         read(*,*)rcutguess2dw,rcutguess2up
        endif
       elseif (rhotype .lt. 3) then
c      Charge Density Function TYPE 2 is Metal Like Oscilating
        read(*,*)beta1dw,beta1up
        read(*,*)beta2dw,beta2up
       elseif (rhotype .lt. 4) then
c      Charge Density Function TYPE 3 is Thomas-Fermi Screening function
        read(*,*)beta1dw,beta1up
       elseif (rhotype .lt. 5) then
c      Charge Density Function TYPE 1 is Hydrogenic 4s orbital (b1,b2)
        read(*,*)beta1dw,beta1up
        read(*,*)beta2dw,beta2up
       elseif (rhotype .lt. 6) then
c      Charge Density Function TYPE 1 is Hydrogenic 4s orbital (b1,b2,b3)
        read(*,*)beta1dw,beta1up
        read(*,*)beta2dw,beta2up
        read(*,*)beta3dw,beta3up
       elseif (rhotype .lt. 9) then
c      Charge Density Function TYPE 1 is Hydrogenic 4s orbital (b1,b2,b3,b4)
        read(*,*)beta1dw,beta1up
        read(*,*)beta2dw,beta2up
        read(*,*)beta3dw,beta3up
        read(*,*)dm3dw,dm3up
       elseif (rhotype .lt. 10) then
c      Charge Density Function TYPE 1 is Hydrogenic 4s orbital (b1,b2,b3,b4)
        read(*,*)beta1dw,beta1up
        read(*,*)beta2dw,beta2up
        read(*,*)beta3dw,beta3up
        read(*,*)dm3dw,dm3up
       elseif (rhotype .lt. 11) then
c      Charge Density Function TYPE 1 is Hydrogenic 4s orbital (b1,b2,b3,b4)
        read(*,*)beta1dw,beta1up
        read(*,*)beta2dw,beta2up
        read(*,*)beta3dw,beta3up
        read(*,*)dm2dw,dm2up
        read(*,*)dm3dw,dm3up
       elseif (rhotype .lt. 12) then
c      Charge Density Function TYPE 1 is Hydrogenic 4s orbital (b1,b2,b3,b4)
        read(*,*)beta1dw,beta1up
        read(*,*)beta2dw,beta2up
        read(*,*)beta3dw,beta3up
        read(*,*)rm3dw,rm3up
        read(*,*)dm2dw,dm2up
        read(*,*)dm3dw,dm3up
       elseif (rhotype .lt. 13) then
c      Charge Density Function TYPE 1 is Hydrogenic 4s orbital (b1,b2,b3)
        read(*,*)beta1dw,beta1up
        read(*,*)beta2dw,beta2up
        read(*,*)beta3dw,beta3up
       elseif (rhotype .lt. 16) then
c      Charge Density Function TYPE 1 is Hydrogenic 4s orbital (b1,b2,b3)
        read(*,*)beta1dw,beta1up
        read(*,*)beta2dw,beta2up
        read(*,*)beta3dw,beta3up
        read(*,*)dm3dw,dm3up
       endif
      endif
c     ----------------------------------------
c     Fitting parameters for Voter Chen EAM
c 
c     write(*,*)'ntypes:',ntypes,'pairnum:',pairnum
      if (eamalloy .gt. 0) then
c     Transformation parameters
c     alloy scaling parameters
       do 6 i=1,ntypes
        read(*,*)sxdw(i),sxup(i)
        read(*,*)gxdw(i),gxup(i)
    6  continue
      endif
      if(phitype .eq. 1 .or. phitype .eq. 7 .or.
     .   phitype .eq. 8) then
c      Pair Function TYPE 1 is a Morse function
       if( pairnum .lt. 1 ) pairnum=1
c      ii=1
       do 7 i=1,pairnum
        if ( i .eq. 3 ) then
         read(*,*)dm3dw,dm3up
         read(*,*)rm3dw,rm3up
         read(*,*)alpha3dw,alpha3up
         read(*,*)rcutguess3dw,rcutguess3up
        elseif ( i .eq. 2 ) then
         read(*,*)dm2dw,dm2up
         read(*,*)rm2dw,rm2up
         read(*,*)alpha2dw,alpha2up
         read(*,*)rcutguess2dw,rcutguess2up
        else
         read(*,*)dm1dw,dm1up
         read(*,*)rm1dw,rm1up
         read(*,*)alpha1dw,alpha1up
         read(*,*)rcutguess1dw,rcutguess1up
        endif
7      continue
       if (phitype.eq.8) then
         read(*,*)dm2dw,dm2up
         read(*,*)rm2dw,rm2up
         read(*,*)alpha2dw,alpha2up
       endif
      endif
c     ----------------------------------------
c     else
c     ----------------------------------------
c     Fitting parameters for Daw Baskes EAM
c
c     Clementi coef. and Density parameters
c     read(*,*)clementifile
c     read(*,*)nn
c     Z(R) function parameters
c     read(*,*)z0
c     read(*,*)nu
c     read(*,*)alpha1
c     read(*,*)beta1
c     Alloy density parameter
c     read(*,*)ns
c     ----------------------------------------
c     endif
c     Alloy cutoff radius for Phi
c     ----------------------------------------
c     
c     if(ntypes.lt.2) then
c        rcutoff=rcutguess
c     endif
c
       if(eamalloy .gt. 0) then
        do 8 i=1,ntypes
            write(*,*)'Reading pure element ',
     .                 trim(eamfiles(i)),' EAM file...'
            open(unit=57,file=eamfiles(i))
          read(57,20)(fheader(j,i),j=1,80)
          write(elmname(i),21)fheader(1,i),fheader(2,i)
          read(57,30) ielement(i), amass(i), blat(i), lat(i)
          write(*,*)'alat_0:',blat(i)
          write(*,8888)' lattice: ',lat(i)
          read(57,9901) nrhoin(i), drhoin(i), nrin(i), drin(i), rcut(i)
          write(*,*)'Nrho:',nrhoin(i),'Nr:',nrin(i)
          read(57,9902) (frhoin(j,i),j=1,nrhoin(i))
          read(57,9902) (zrin(j,i),j=1,nrin(i))
          read(57,9902) (rhorin(j,i),j=1,nrin(i))
          close(57)
    8   continue
       endif


      endif
c     -------------------
c       END MASTER ONLY
c     -------------------
8888  format(a10,a3)

      call MPI_BCAST(logflag,1,MPI_INTEGER,MASTER,
     &MPI_COMM_WORLD,ierr)

      if(logflag.eq.1) then
       myfileid=myid+1
       write(myidstr,'(I4)') myfileid
       logfilestr="APSOlog."
       if (myfileid .lt. 10) then
        write(logfile,'(a8,a1)')logfilestr,adjustl(myidstr)
       else if (myfileid .lt. 100) then
        write(logfile,'(a8,a2)')logfilestr,adjustl(myidstr)
       else if (myfileid .lt. 1000) then
        write(logfile,'(a8,a3)')logfilestr,adjustl(myidstr)
       else
        write(logfile,'(a8,a4)')logfilestr,adjustl(myidstr)
       endif
        open(unit=11,file=TRIM(logfile))
        write(11,*)'myid: ',myid
        write(11,*)'eamftyp: ',(eamftyp(i),i=1,10)
      endif
      printon=1
      
      if (myid .eq. MASTER) then
c      write(*,'(A,I3,A)')'PROC:',myid,' packing...'
        place=0
        backpack(1)=DBLE(phitype)
        backpack(2)=DBLE(rhotype)
        backpack(3)=DBLE(embtype)
        backpack(4)=DBLE(expnum)
        backpack(5)=DBLE(ntypes)
        backpack(6)=DBLE(ielement(1))
        backpack(7)=DBLE(amass(1))
        backpack(8)=DBLE(at1)
        backpack(9)=DBLE(at2)
        backpack(10)=DBLE(nxchain)
        backpack(11)=DBLE(nyatoms)
        backpack(12)=DBLE(nzlayer)
        backpack(13)=DBLE(latai)
        backpack(14)=DBLE(latstep)
        backpack(15)=DBLE(lataf)
        backpack(16)=DBLE(esub)
        backpack(17)=DBLE(bmod)
        backpack(18)=DBLE(rnum)
        backpack(19)=DBLE(rhonum)
        backpack(20)=DBLE(rhoend)
        backpack(21)=DBLE(npartsin)
        backpack(22)=DBLE(ngensin)
        backpack(23)=DBLE(dyntypin)
        backpack(24)=DBLE(prntflag)
        backpack(25)=lata0
        backpack(26)=caratio
        backpack(27)=beta1dw
        backpack(28)=beta2dw
        backpack(29)=beta3dw
        backpack(30)=beta1up
        backpack(31)=beta2up
        backpack(32)=beta3up
        backpack(33)=sxdw(1)
        backpack(34)=sxdw(2)
        backpack(35)=sxdw(3)
        backpack(36)=sxup(1)
        backpack(37)=sxup(2)
        backpack(38)=sxup(3)
        backpack(39)=gxdw(1)
        backpack(40)=gxdw(2)
        backpack(41)=gxdw(3)
        backpack(42)=gxup(1)
        backpack(43)=gxup(2)
        backpack(44)=gxup(3)
        backpack(45)=dm1dw
        backpack(46)=dm1up
        backpack(47)=rm1dw
        backpack(48)=rm1up
        backpack(49)=alpha1dw
        backpack(50)=alpha1up
        backpack(51)=rcutguess1dw
        backpack(52)=rcutguess1up
        backpack(53)=dm2dw
        backpack(54)=dm2up
        backpack(55)=rm2dw
        backpack(56)=rm2up
        backpack(57)=alpha2dw
        backpack(58)=alpha2up
        backpack(59)=rcutguess2dw
        backpack(60)=rcutguess2up
        backpack(61)=dm3dw
        backpack(62)=dm3up
        backpack(63)=rm3dw
        backpack(64)=rm3up
        backpack(65)=alpha3dw
        backpack(66)=alpha3up
        backpack(67)=rcutguess3dw
        backpack(68)=rcutguess3up
        backpack(69)=nrhoin(1)
        backpack(70)=nrhoin(2)
        backpack(71)=nrhoin(3)
        backpack(72)=drhoin(1)
        backpack(73)=drhoin(2)
        backpack(74)=drhoin(3)
        backpack(75)=DBLE(nrin(1))
        backpack(76)=DBLE(nrin(2))
        backpack(77)=DBLE(nrin(3))
        backpack(78)=drin(1)
        backpack(79)=drin(2)
        backpack(80)=drin(3)
        backpack(81)=rcut(1)
        backpack(82)=rcut(2)
        backpack(83)=rcut(3)
        backpack(84)=fixcutoff
        backpack(85)=DBLE(weightnum)
c      write(*,'(A,I3,A)')'PROC:',myid,' packed!'
      endif

      call MPI_BCAST(backpack,85,MPI_DOUBLE_PRECISION,MASTER,
     &MPI_COMM_WORLD,ierr)
    
      if (myid .ne. MASTER) then
        place=0
c      write(*,'(A,I3,A)')'PROC:',myid,' unpacking...'
        phitype=INT(backpack(1))
        rhotype=INT(backpack(2))
        embtype=INT(backpack(3))
        expnum=INT(backpack(4))
        ntypes=INT(backpack(5))
        ielement(1)=INT(backpack(6))
        amass(1)=backpack(7)
        at1=INT(backpack(8))
        at2=INT(backpack(9))
        nxchain=INT(backpack(10))
        nyatoms=INT(backpack(11))
        nzlayer=INT(backpack(12))
        latai=backpack(13)
        latstep=backpack(14)
        lataf=backpack(15)
        esub=backpack(16)
        bmod=backpack(17)
        rnum=INT(backpack(18))
        rhonum=INT(backpack(19))
        rhoend=backpack(20)
        npartsin=INT(backpack(21))
        ngensin=INT(backpack(22))
        dyntypin=INT(backpack(23))
        prntflag=INT(backpack(24))
        lata0=backpack(25)
        caratio=backpack(26)
        beta1dw=backpack(27)
        beta2dw=backpack(28)
        beta3dw=backpack(29)
        beta1up=backpack(30)
        beta2up=backpack(31)
        beta3up=backpack(32)
        sxdw(1)=backpack(33)
        sxdw(2)=backpack(34)
        sxdw(3)=backpack(35)
        sxup(1)=backpack(36)
        sxup(2)=backpack(37)
        sxup(3)=backpack(38)
        gxdw(1)=backpack(39)
        gxdw(2)=backpack(40)
        gxdw(3)=backpack(41)
        gxup(1)=backpack(42)
        gxup(2)=backpack(43)
        gxup(3)=backpack(44)
        dm1dw=backpack(45)
        dm1up=backpack(46)
        rm1dw=backpack(47)
        rm1up=backpack(48)
        alpha1dw=backpack(49)
        alpha1up=backpack(50)
        rcutguess1dw=backpack(51)
        rcutguess1up=backpack(52)
        dm2dw=backpack(53)
        dm2up=backpack(54)
        rm2dw=backpack(55)
        rm2up=backpack(56)
        alpha2dw=backpack(57)
        alpha2up=backpack(58)
        rcutguess2dw=backpack(59)
        rcutguess2up=backpack(60)
        dm3dw=backpack(61)
        dm3up=backpack(62)
        rm3dw=backpack(63)
        rm3up=backpack(64)
        alpha3dw=backpack(65)
        alpha3up=backpack(66)
        rcutguess3dw=backpack(67)
        rcutguess3up=backpack(68)
        nrhoin(1)=INT(backpack(69))
        nrhoin(2)=INT(backpack(70))
        nrhoin(3)=INT(backpack(71))
        drhoin(1)=backpack(72)
        drhoin(2)=backpack(73)
        drhoin(3)=backpack(74)
        nrin(1)=INT(backpack(75))
        nrin(2)=INT(backpack(76))
        nrin(3)=INT(backpack(77))
        drin(1)=backpack(78)
        drin(2)=backpack(79)
        drin(3)=backpack(80)
        rcut(1)=backpack(81)
        rcut(2)=backpack(82)
        rcut(3)=backpack(83)
        fixcutoff=INT(backpack(84))
        weightnum=INT(backpack(85))
c      write(*,'(A,I3,A)')'PROC:',myid,' unpacked!'
      endif

       if(logflag.eq.1) then
        write(11,*)'fixcutoff: ',fixcutoff
        write(11,*)'phitype: ',phitype
        write(11,*)'rhotype: ',rhotype
        write(11,*)'embtype: ',embtype
        write(11,*)'expnum: ',expnum
        write(11,*)'weightnum: ',weightnum
        write(11,*)'ntypes: ',ntypes
        write(11,*)'ielement(1): ',ielement(1)
        write(11,*)'amass(1): ',amass(1)
        write(11,*)'at1: ',at1
        write(11,*)'at2: ',at2
        write(11,*)'nxchain: ',nxchain
        write(11,*)'nyatoms: ',nyatoms
        write(11,*)'nzlayer: ',nzlayer
        write(11,*)'latai: ',latai
        write(11,*)'latstep: ',latstep
        write(11,*)'lataf: ',lataf
        write(11,*)'esub: ',esub
        write(11,*)'bmod: ',bmod
        write(11,*)'rnum: ',rnum
        write(11,*)'rhonum: ',rhonum
        write(11,*)'rhoend: ',rhoend
        write(11,*)'npartsin: ',npartsin
        write(11,*)'ngensin: ',ngensin
        write(11,*)'dyntypin: ',dyntypin
        write(11,*)'prntflag: ',prntflag
        write(11,*)'lata0: ',lata0
        write(11,*)'caratio: ',caratio
        write(11,*)'beta1dw: ',beta1dw
        write(11,*)'beta2dw: ',beta2dw
        write(11,*)'beta3dw: ',beta3dw
        write(11,*)'beta1up: ',beta1up
        write(11,*)'beta2up: ',beta2up
        write(11,*)'beta3up: ',beta3up
        write(11,*)'sxdw(1): ',sxdw(1)
        write(11,*)'sxdw(2): ',sxdw(2)
        write(11,*)'sxdw(3): ',sxdw(3)
        write(11,*)'sxup(1): ',sxup(1)
        write(11,*)'sxup(2): ',sxup(2)
        write(11,*)'sxup(3): ',sxup(3)
        write(11,*)'gxdw(1): ',gxdw(1)
        write(11,*)'gxdw(2): ',gxdw(2)
        write(11,*)'gxdw(3): ',gxdw(3)
        write(11,*)'gxup(1): ',gxup(1)
        write(11,*)'gxup(2): ',gxup(2)
        write(11,*)'gxup(3): ',gxup(3)
        write(11,*)'dm1dw: ',dm1dw
        write(11,*)'dm1up: ',dm1up
        write(11,*)'rm1dw: ',rm1dw
        write(11,*)'rm1up: ',rm1up
        write(11,*)'alpha1dw: ',alpha1dw
        write(11,*)'alpha1up: ',alpha1up
        write(11,*)'rcutguess1dw: ',rcutguess1dw
        write(11,*)'rcutguess1up: ',rcutguess1up
        write(11,*)'dm2dw: ',dm2dw
        write(11,*)'dm2up: ',dm2up
        write(11,*)'rm2dw: ',rm2dw
        write(11,*)'rm2up: ',rm2up
        write(11,*)'alpha2dw: ',alpha2dw
        write(11,*)'alpha2up: ',alpha2up
        write(11,*)'rcutguess2dw: ',rcutguess2dw
        write(11,*)'rcutguess2up: ',rcutguess2up
        write(11,*)'dm3dw: ',dm3dw
        write(11,*)'dm3up: ',dm3up
        write(11,*)'rm3dw: ',rm3dw
        write(11,*)'rm3up: ',rm3up
        write(11,*)'alpha3dw: ',alpha3dw
        write(11,*)'alpha3up: ',alpha3up
        write(11,*)'rcutguess3dw: ',rcutguess3dw
        write(11,*)'rcutguess3up: ',rcutguess3up
        write(11,*)'nrhoin(1): ',nrhoin(1)
        write(11,*)'nrhoin(2): ',nrhoin(2)
        write(11,*)'nrhoin(3): ',nrhoin(3)
        write(11,*)'drhoin(1): ',drhoin(1)
        write(11,*)'drhoin(2): ',drhoin(2)
        write(11,*)'drhoin(3): ',drhoin(3)
        write(11,*)'nrin(1): ',nrin(1)
        write(11,*)'nrin(2): ',nrin(2)
        write(11,*)'nrin(3): ',nrin(3)
        write(11,*)'drin(1): ',drin(1)
        write(11,*)'drin(2): ',drin(2)
        write(11,*)'drin(3): ',drin(3)
        write(11,*)'rcut(1): ',rcut(1)
        write(11,*)'rcut(2): ',rcut(2)
        write(11,*)'rcut(3): ',rcut(3)
       endif
      
c      write(*,'(A,I3,A,I30)')'PROC:',myid,' Bcast expdata:',expnum
      call MPI_BCAST(expdata,expnum,MPI_DOUBLE_PRECISION,MASTER,
     &MPI_COMM_WORLD,ierr)
      
       if(logflag.eq.1) then
        write(11,*)'expdata: ',(expdata(i),i=1,expnum)
       endif

      call MPI_BCAST(weightdata,weightnum,MPI_DOUBLE_PRECISION,MASTER,
     &MPI_COMM_WORLD,ierr)
      
       if(logflag.eq.1) then
        write(11,*)'weightdata: ',(weightdata(i),i=1,weightnum)
       endif
c     Procs learn pure eam functions if alloy eam generation 
      if (eamalloy .gt. 0) then
c      write(*,'(A,I3,A,I7)')'PROC:',myid,' Bcast frhoin1:',nrhoin(1)
      if (myid .eq. MASTER) then
        k=1
        do 71 i=1,ntypes
           do 72 j=1,nrhoin(i) 
             gridpack(k)=frhoin(j,i)
             k=k+1
72          continue
71       continue
        gridpacknum=k-1
      endif
      call MPI_BCAST(gridpacknum,1,MPI_INTEGER,MASTER,
     &MPI_COMM_WORLD,ierr)
      call MPI_BCAST(gridpack,gridpacknum,MPI_DOUBLE_PRECISION,MASTER,
     &MPI_COMM_WORLD,ierr)
c      write(*,'(A,I3,A)')'PROC:',myid,' Bcast zrin'
      if (myid .ne. MASTER) then
        k=1
        do 144 i=1,ntypes
           do 154 j=1,nrhoin(i) 
             frhoin(j,i)=gridpack(k)
             k=k+1
154          continue
144       continue
      endif
      if (myid .eq. MASTER) then
        k=1
        do 145 i=1,ntypes
           do 155 j=1,nrin(i) 
             gridpack(k)=zrin(j,i)
             k=k+1
155          continue
145       continue
        gridpacknum=k-1
      endif
      call MPI_BCAST(gridpacknum,1,MPI_INTEGER,MASTER,
     &MPI_COMM_WORLD,ierr)
      call MPI_BCAST(gridpack,gridpacknum,MPI_DOUBLE_PRECISION,MASTER,
     &MPI_COMM_WORLD,ierr)
c      write(*,'(A,I3,A)')'PROC:',myid,' Broudcast rhorin'
      if (myid .ne. MASTER) then
        k=1
        do 146 i=1,ntypes
           do 156 j=1,nrhoin(i) 
             zrin(j,i)=gridpack(k)
             k=k+1
156          continue
146       continue
      endif
      if (myid .eq. MASTER) then
        k=1
        do 147 i=1,ntypes
           do 157 j=1,nrin(i) 
             gridpack(k)=rhorin(j,i)
             k=k+1
157          continue
147       continue
        gridpacknum=k-1
      endif
      call MPI_BCAST(gridpacknum,1,MPI_INTEGER,MASTER,
     &MPI_COMM_WORLD,ierr)
      call MPI_BCAST(gridpack,gridpacknum,MPI_DOUBLE_PRECISION,MASTER,
     &MPI_COMM_WORLD,ierr)
      if (myid .ne. MASTER) then
        k=1
        do 148 i=1,ntypes
           do 158 j=1,nrhoin(i) 
             rhorin(j,i)=gridpack(k)
             k=k+1
158          continue
148       continue
      endif

        write(*,'(A,I3,A,3(A,E18.10))')'PROC(',myid,')   ',
     &'frhoin(7,2):',frhoin(7,2),' zrin(7,2):',zrin(7,2),
     &'frhoin(nrhoin-1,1):',frhoin(nrhoin(1)-1,1),' zrin(nrin-1,1):',
     &zrin(nrin(1)-1,1),' rhorin(7,2):',rhorin(7,2),
     &' rhorin(nrin-1,1):',rhorin(nrin(1)-1,1)

      do 78 i=1,ntypes
c         Restart with initial EAM tables
          snrhoin(i)=nrhoin(i)
          sdrhoin(i)=drhoin(i)
          snrin(i)=nrin(i)
          sdrin(i)=drin(i)
          do 73 j=1,nrhoin(i)
             sfrhoin(j,i)=frhoin(j,i)
73         continue
          do 79 j=1,nrin(i)
             szrin(j,i)=zrin(j,i)
             srhorin(j,i)=rhorin(j,i)
79         continue
78    continue
c     end of alloy pure eam reads      
      endif

      rms=200000.0
      pairnum=ntypes*(ntypes-1)/2
      if(pairnum.lt.1) pairnum=1

      if (ntypes.eq.3) then
        if (sxdw(1).eq.sxup(1) .and. sxdw(2).eq.sxup(2) .and. 
     &      sxdw(3).eq.sxup(3) .and. gxdw(1).eq.gxup(1) .and.
     &      gxdw(2).eq.gxup(2) .and. gxdw(3).eq.gxup(3)) then
          if (sxdw(1).eq.1.0) then
            nparams=17
            prmflag=1
            funcprm(1)=gxdw(1)
            funcprm(2)=sxdw(2)
            funcprm(3)=gxdw(2)
            funcprm(4)=sxdw(3)
            funcprm(5)=gxdw(3)
          elseif (sxdw(2).eq.1.0) then
            nparams=17
            prmflag=2
            funcprm(1)=sxdw(1)
            funcprm(2)=gxdw(1)
            funcprm(3)=gxdw(2)
            funcprm(4)=sxdw(3)
            funcprm(5)=gxdw(3)
          elseif (sxdw(3).eq.1.0) then
            nparams=17
            prmflag=3
            funcprm(1)=sxdw(1)
            funcprm(2)=gxdw(1)
            funcprm(3)=sxdw(2)
            funcprm(4)=gxdw(2)
            funcprm(5)=gxdw(3)
          endif
            funcprm(6)=dm1dw
            funcprm(7)=rm1dw
            funcprm(8)=alpha1dw
            funcprm(9)=rcutguess1dw
            funcprm(10)=dm2dw
            funcprm(11)=rm2dw
            funcprm(12)=alpha2dw
            funcprm(13)=rcutguess2dw
            funcprm(14)=dm3dw
            funcprm(15)=rm3dw
            funcprm(16)=alpha3dw
            funcprm(17)=rcutguess3dw
      if (myid .eq. MASTER) then
       write(*,*)'--------------------------------'
       write(*,*)' EAM Main routine call'
       write(*,*)'--------------------------------'
c      prntflag=1
       prmno=0
       call eammain(rms,funcprm,nparams,lcharnum,fitvalues,expdata,
     .eamfileout,snrhoin,sdrhoin,snrin,sdrin,sfrhoin,szrin,srhorin,
     .eamftyp,prmflag,printon,logflag,myid,prmno,weightdata)
       write(*,*)'--------------------------------'
       write(*,*)' RMS:',rms
       write(*,*)'--------------------------------'
      endif
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         call MPI_FINALIZE(ierr)
         stop 
        endif
      else if (ntypes .eq. 2) then
        if (sxdw(1).eq.sxup(1) .and. sxdw(2).eq.sxup(2) .and. 
     &      gxdw(1).eq.gxup(1) .and. gxdw(2).eq.gxup(2)) then
          if (sxdw(1).eq.1.0) then
            nparams=7
            prmflag=1
            funcprm(1)=gxdw(1)
            funcprm(2)=sxdw(2)
            funcprm(3)=gxdw(2)
          elseif (sxdw(2).eq.1.0) then
            nparams=7
            prmflag=2
            funcprm(1)=sxdw(1)
            funcprm(2)=gxdw(1)
            funcprm(3)=gxdw(2)
          endif
            funcprm(4)=dm1dw
            funcprm(5)=rm1dw
            funcprm(6)=alpha1dw
            funcprm(7)=rcutguess1dw
      if (myid .eq. MASTER) then
       write(*,*)'--------------------------------'
       write(*,*)' EAM Main routine call'
       write(*,*)'--------------------------------'
c      prntflag=1
       prmno=0
       call eammain(rms,funcprm,nparams,lcharnum,fitvalues,expdata,
     .eamfileout,snrhoin,sdrhoin,snrin,sdrin,sfrhoin,szrin,srhorin,
     .eamftyp,prmflag,printon,logflag,myid,prmno,weightdata)
       write(*,*)'--------------------------------'
       write(*,*)' RMS:',rms
       write(*,*)'--------------------------------'
      endif
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         call MPI_FINALIZE(ierr)
         stop 
        endif
      else
       if (phitype .eq. 1 .or. phitype .eq. 7) then
         k=1
         funcprm(k)=dm1dw
         funcprm(k+1)=rm1dw
         funcprm(k+2)=alpha1dw
         funcprm(k+3)=rcutguess1dw
         k=k+3
       elseif (phitype .eq. 8) then
         k=1
         funcprm(k)=dm1dw
         funcprm(k+1)=rm1dw
         funcprm(k+2)=alpha1dw
         funcprm(k+3)=rcutguess1dw
         funcprm(k+4)=dm2dw
         funcprm(k+5)=rm2dw
         funcprm(k+6)=alpha2dw
         k=k+6
       elseif (phitype .eq. 2) then
         k=1
         funcprm(k)=dm1dw
         funcprm(k+1)=alpha1dw
         funcprm(k+2)=rcutguess1dw
         k=k+2
       endif
       if (rhotype .eq. 1 .or. rhotype .eq. 3) then
         funcprm(k+1)=beta1dw
         if(fixcutoff.lt.1) then
           funcprm(k+2)=rcutguess2dw
           k=k+1
         endif
         k=k+1
       elseif (rhotype .eq. 2 .or. rhotype .eq. 4) then
         funcprm(k+1)=beta1dw
         funcprm(k+2)=beta2dw
         k=k+2
       elseif (rhotype .eq. 5 .or. rhotype .eq. 12) then
         funcprm(k+1)=beta1dw
         funcprm(k+2)=beta2dw
         funcprm(k+3)=beta3dw
         k=k+3
       elseif (rhotype .eq. 6 .or.
     .         rhotype .eq. 7 .or.
     .         rhotype .eq. 8 .or.
     .         rhotype .eq. 9 .or.
     .         rhotype .eq. 13 .or.
     .         rhotype .eq. 14 .or.
     .         rhotype .eq. 15) then
         funcprm(k+1)=beta1dw
         funcprm(k+2)=beta2dw
         funcprm(k+3)=beta3dw
         funcprm(k+4)=dm3dw
         k=k+4
       elseif (rhotype .eq. 10 ) then
         funcprm(k+1)=beta1dw
         funcprm(k+2)=beta2dw
         funcprm(k+3)=beta3dw
         funcprm(k+4)=dm2dw
         funcprm(k+5)=dm3dw
         k=k+5
       elseif (rhotype .eq. 11 ) then
         funcprm(k+1)=beta1dw
         funcprm(k+2)=beta2dw
         funcprm(k+3)=beta3dw
         funcprm(k+4)=rm3dw
         funcprm(k+5)=dm2dw
         funcprm(k+6)=dm3dw
         k=k+6
       endif
        nparams=k
        prmflag=0
        if (dm1dw.eq.dm1up .and. alpha1dw.eq.alpha1up .and. 
     &  rcutguess1dw.eq.rcutguess1up .and. beta1dw.eq.beta1up) then
      if (myid .eq. MASTER) then
       write(*,*)'--------------------------------'
       write(*,*)' EAM Main routine call for PURE '
       write(*,*)'--------------------------------'
c      prntflag=1
       prmno=0
       call eammain(rms,funcprm,nparams,lcharnum,fitvalues,expdata,
     .eamfileout,snrhoin,sdrhoin,snrin,sdrin,sfrhoin,szrin,srhorin,
     .eamftyp,prmflag,printon,logflag,myid,prmno,weightdata)
       write(*,*)'--------------------------------'
       write(*,*)' RMS:',rms
       write(*,*)'--------------------------------'
      endif
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         call MPI_FINALIZE(ierr)
         stop 
        endif
      endif


c     open(unit=17,file='chist')
c     open(unit=18,file='randhist')
c     --------------------------------------
c       APSO Inputs
c     --------------------------------------
      pairnum=ntypes*(ntypes-1)/2
      if(pairnum.lt.1) pairnum=1
      nparts=npartsin
      if (ntypes .eq. 2) then
        if (sxdw(1).eq.1.0 .and. sxup(1).eq.1.0) then
          nparams=7
          prmflag=1
        elseif (sxdw(2).eq.1.0 .and. sxup(2).eq.1.0) then
          nparams=7
          prmflag=2
        else
          nparams=8
        endif
      elseif (ntypes .eq. 3) then
        if (sxdw(1).eq.1.0 .and. sxup(1).eq.1.0) then
          nparams=17
          prmflag=1
        elseif (sxdw(2).eq.1.0 .and. sxup(2).eq.1.0) then
          nparams=17
          prmflag=2
        elseif (sxdw(3).eq.1.0 .and. sxup(3).eq.1.0) then
          nparams=17
          prmflag=3
        else
          nparams=18
        endif
      else
          nparams=5
          if (fixcutoff.lt.1 .or. rhotype.eq.4) then
          nparams=nparams+1
          else if (rhotype.eq.5 .or. rhotype .eq. 12) then
          nparams=nparams+2
          else if (rhotype.eq.6 .or.
     .             rhotype.eq.7 .or.
     .             rhotype.eq.8 .or.
     .             rhotype.eq.9 .or.
     .             rhotype.eq.13 .or.
     .             rhotype.eq.14 .or.
     .             rhotype.eq.15) then
          nparams=nparams+3
          else if (rhotype.eq.10) then
          nparams=nparams+4
          else if (rhotype.eq.11) then
          nparams=nparams+5
          endif
          if (phitype.eq.8) then
            nparams=nparams+3
          endif
          prmflag=0
      endif
      ngens=ngensin
c     ftyp=0
      dyntyp=dyntypin
      if(logflag.eq.1) then
       write(11,'(A,I3,3(A,I7))')'PROC(',myid,')   Parts:',nparts,
     &              ' Params:',nparams,' Gens:',ngens
      endif
      if (myid.eq.MASTER) then
       write(*,'(A,I3,3(A,I7))')'PROC(',myid,')   Parts:',nparts,
     &              ' Params:',nparams,' Gens:',ngens
      if (dyntyp.eq.1) then
          write(*,'(A,I3,A)')'PROC(',myid,')   Dynamic system: On'
      else
          write(*,'(A,I3,A)')'PROC(',myid,')   Dynamic system: Off'
      endif
      if (prntflag.eq.1) then
          write(*,'(A,I3,A)')'PROC(',myid,')   Display printing: On'
      else
          write(*,'(A,I3,A)')'PROC(',myid,')   Display printing: Off'
      endif
      endif
c      do 11 i=1,nparams
c         bounds(i,1)=-10
c         bounds(i,2)=10
c11    continue
    
      if (myid .eq. MASTER) then 
       if(lattype.eq.0) then
          write(*,*)'Lattice type for EAM calculation:  FCC'
       else
          write(*,*)'Lattice type for EAM calculation:  HCP'
       endif
       write(*,'(A,I4,4(A,I7))')'PROC(',myid,')   Ntypes:',ntypes,
     &     ' expnum:',expnum,'#w:',weightnum,' Type1:',at1,' Type2:',at2
       write(*,'(A,I3,3(A,I4))')'PROC(',myid,')   Supercell:',nxchain,
     &     ' x ',nyatoms,' x ',nzlayer
       write(*,'(A,I3,2(A,E17.10),2(A,I6))')'PROC(',myid,')   Lata0:',
     & lata0,' c/a:',caratio,' rnum:',rnum,' pairnum:',pairnum
       if (eamalloy .gt. 0) then
       do 100 i=1,ntypes 
        write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   sx['
     &,sxdw(i),',',sxup(i),']'
        write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   gx['
     &,gxdw(i),',',gxup(i),']'
100    continue
       endif
       do 101 i=1,pairnum 
        if (i.eq.1) then
          write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   dm1:[',dm1dw,
     &',',dm1up,']'
          write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   rm1:[',rm1dw,
     &',',rm1up,']'
          write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   alpha1:[',
     &alpha1dw,',',alpha1up,']'
          write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,
     &')   rcutguess1:[',rcutguess1dw,',',rcutguess1up,']'
        elseif (i.eq.2) then
          write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   dm2:[',dm2dw,
     &',',dm2up,']'
          write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   rm2:[',rm2dw,
     &',',rm2up,']'
          write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   alpha2:[',
     &alpha2dw,',',alpha2up,']'
          write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,
     &')   rcutguess2:[',rcutguess2dw,',',rcutguess2up,']'
        elseif (i.eq.3) then
          write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   dm3:[',dm3dw,
     &',',dm3up,']'
          write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   rm3:[',rm3dw,
     &',',rm3up,']'
          write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   alpha3:[',
     &alpha3dw,',',alpha3up,']'
          write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,
     &')   rcutguess3:[',rcutguess3dw,',',rcutguess3up,']'
        endif
101    continue
       if (eamalloy .lt. 1) then
          write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   beta1:[',
     &beta1dw,',',beta1up,']'
          if (fixcutoff.lt.1) then
          write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,') rcutguess2:[',
     &rcutguess2dw,',',rcutguess2up,']'
          endif
       endif
      endif
     
 
      if(logflag.eq.1) then
       write(11,'(A,I4,4(A,I7))')'PROC(',myid,')   Ntypes:',ntypes,
     &     ' expnum:',expnum,'#w:',weightnum,' Type1:',at1,' Type2:',at2
       write(11,'(A,I3,3(A,I4))')'PROC(',myid,')   Supercell:',nxchain,
     &     ' x ',nyatoms,' x ',nzlayer
       write(11,'(A,I3,2(A,E17.10),2(A,I6))')'PROC(',myid,')   Lata0:',
     & lata0,' c/a:',caratio,' rnum:',rnum,' pairnum:',pairnum
       if (eamalloy .gt. 0) then
       do 107 i=1,ntypes 
        write(11,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   sx['
     &,sxdw(i),',',sxup(i),']'
        write(11,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   gx['
     &,gxdw(i),',',gxup(i),']'
107    continue
       endif
       do 108 i=1,pairnum 
        if (i.eq.1) then
        write(11,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   dm1:[',dm1dw,
     &',',dm1up,']'
        write(11,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   rm1:[',rm1dw,
     &',',rm1up,']'
        write(11,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   alpha1:[',
     &alpha1dw,',',alpha1up,']'
        write(11,'(A,I3,A,2(E17.10,A))')'PROC(',myid,
     &')   rcutguess1:[',rcutguess1dw,',',rcutguess1up,']'
        elseif (i.eq.2) then
        write(11,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   dm2:[',dm2dw,
     &',',dm2up,']'
        write(11,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   rm2:[',rm2dw,
     &',',rm2up,']'
        write(11,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   alpha2:[',
     &alpha2dw,',',alpha2up,']'
        write(11,'(A,I3,A,2(E17.10,A))')'PROC(',myid,
     &')   rcutguess2:[',rcutguess2dw,',',rcutguess2up,']'
        elseif (i.eq.3) then
        write(11,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   dm3:[',dm3dw,
     &',',dm3up,']'
        write(11,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   rm3:[',rm3dw,
     &',',rm3up,']'
        write(11,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   alpha3:[',
     &alpha3dw,',',alpha3up,']'
        write(11,'(A,I3,A,2(E17.10,A))')'PROC(',myid,
     &')   rcutguess3:[',rcutguess3dw,',',rcutguess3up,']'
        endif
108    continue
       if (eamalloy .lt. 1) then
          write(11,'(A,I3,A,2(E17.10,A))')'PROC(',myid,')   beta1:[',
     &beta1dw,',',beta1up,']'
          if (fixcutoff.lt.1) then
          write(*,'(A,I3,A,2(E17.10,A))')'PROC(',myid,') rcutguess2:[',
     &rcutguess2dw,',',rcutguess2up,']'
          endif
       endif
      endif


      ii=1
      if (eamalloy .gt. 0) then
       if (sxdw(1).ne.1.0 .and. sxup(1).ne.1.0) then
          bounds(ii,1)=sxdw(1)
          bounds(ii,2)=sxup(1)
          ii=ii+1
       endif
          bounds(ii,1)=gxdw(1)
          bounds(ii,2)=gxup(1)
          ii=ii+1
       if (sxdw(2).ne.1.0 .and. sxup(2).ne.1.0) then
          bounds(ii,1)=sxdw(2)
          bounds(ii,2)=sxup(2)
          ii=ii+1
       endif
          bounds(ii,1)=gxdw(2)
          bounds(ii,2)=gxup(2)
          ii=ii+1
       if (nparams.gt.8) then
        if (sxdw(3).ne.1.0 .and. sxup(3).ne.1.0) then
          bounds(ii,1)=sxdw(3)
          bounds(ii,2)=sxup(3)
          ii=ii+1
        endif
          bounds(ii,1)=gxdw(3)
          bounds(ii,2)=gxup(3)
          ii=ii+1
       endif
          bounds(ii,1)=dm1dw
          bounds(ii,2)=dm1up
          bounds(ii+1,1)=rm1dw
          bounds(ii+1,2)=rm1up
          bounds(ii+2,1)=alpha1dw
          bounds(ii+2,2)=alpha1up
          bounds(ii+3,1)=rcutguess1dw
          bounds(ii+3,2)=rcutguess1up
          ii=ii+4
       if (nparams.gt.8) then
          bounds(ii,1)=dm2dw
          bounds(ii,2)=dm2up
          bounds(ii+1,1)=rm2dw
          bounds(ii+1,2)=rm2up
          bounds(ii+2,1)=alpha2dw
          bounds(ii+2,2)=alpha2up
          bounds(ii+3,1)=rcutguess2dw
          bounds(ii+3,2)=rcutguess2up
          bounds(ii+4,1)=dm3dw
          bounds(ii+4,2)=dm3up
          bounds(ii+5,1)=rm3dw
          bounds(ii+5,2)=rm3up
          bounds(ii+6,1)=alpha3dw
          bounds(ii+6,2)=alpha3up
          bounds(ii+7,1)=rcutguess3dw
          bounds(ii+7,2)=rcutguess3up
          ii=ii+8
       endif
      else
          bounds(ii,1)=dm1dw
          bounds(ii,2)=dm1up
          bounds(ii+1,1)=rm1dw
          bounds(ii+1,2)=rm1up
          bounds(ii+2,1)=alpha1dw
          bounds(ii+2,2)=alpha1up
          bounds(ii+3,1)=rcutguess1dw
          bounds(ii+3,2)=rcutguess1up
          ii=ii+4
          if(phitype.eq.8) then
          bounds(ii,1)=dm2dw
          bounds(ii,2)=dm2up
          bounds(ii+1,1)=rm2dw
          bounds(ii+1,2)=rm2up
          bounds(ii+2,1)=alpha2dw
          bounds(ii+2,2)=alpha2up
          ii=ii+3
          endif 
          bounds(ii,1)=beta1dw
          bounds(ii,2)=beta1up
          ii=ii+1
          if (fixcutoff.lt.1) then
            bounds(ii,1)=rcutguess2dw
            bounds(ii,2)=rcutguess2up
            ii=ii+1
          endif
          if (rhotype.eq.4) then
            bounds(ii,1)=beta2dw
            bounds(ii,2)=beta2up
            ii=ii+1
          else if (rhotype.eq.5 .or.
     .             rhotype.eq.12) then
            bounds(ii,1)=beta2dw
            bounds(ii,2)=beta2up
            bounds(ii+1,1)=beta3dw
            bounds(ii+1,2)=beta3up
            ii=ii+2
          else if (rhotype.eq.6 .or.
     .             rhotype.eq.7 .or.
     .             rhotype.eq.8 .or.
     .             rhotype.eq.9 .or.
     .             rhotype.eq.13 .or.
     .             rhotype.eq.14 .or.
     .             rhotype.eq.15) then
            bounds(ii,1)=beta2dw
            bounds(ii,2)=beta2up
            bounds(ii+1,1)=beta3dw
            bounds(ii+1,2)=beta3up
            bounds(ii+2,1)=dm3dw
            bounds(ii+2,2)=dm3up
            ii=ii+3
          else if (rhotype.eq.10) then
            bounds(ii,1)=beta2dw
            bounds(ii,2)=beta2up
            bounds(ii+1,1)=beta3dw
            bounds(ii+1,2)=beta3up
            bounds(ii+2,1)=dm2dw
            bounds(ii+2,2)=dm2up
            bounds(ii+3,1)=dm3dw
            bounds(ii+3,2)=dm3up
            ii=ii+4
          else if (rhotype.eq.11) then
            bounds(ii,1)=beta2dw
            bounds(ii,2)=beta2up
            bounds(ii+1,1)=beta3dw
            bounds(ii+1,2)=beta3up
            bounds(ii+2,1)=rm3dw
            bounds(ii+2,2)=rm3up
            bounds(ii+3,1)=dm2dw
            bounds(ii+3,2)=dm2up
            bounds(ii+4,1)=dm3dw
            bounds(ii+4,2)=dm3up
            ii=ii+5
          endif
      endif
      ii=ii-1
      if (myid .eq. MASTER) then
       write(*,*)'Params in bounds:',ii
       write(*,*)'Fixed Sx Param :',prmflag
      endif
      if(logflag.eq.1) then
       write(11,*)'Params in bounds:',ii
       write(11,*)'Fixed Sx Param :',prmflag
      endif

c     --------------------------------------
      aveparts=nparts/INT(numprocs)
      extra=MOD(nparts,numprocs)
      offset=1
      print*,'ok! numprocs:',numprocs,myid
      do 102 dest=1,numprocs
        if (dest .le. extra) then
           mypartsnum(dest)=aveparts + 1
           rcount(dest-1)=(aveparts + 1)
        else
           mypartsnum(dest)=aveparts
           rcount(dest-1)=aveparts
        endif
        offsets(dest)=offset
        displs(dest-1)=(offset-1)
        offset=offset+mypartsnum(dest)
c       write(*,*)'dest: ',dest,myid
102   continue
      myparts=mypartsnum(myid+1)
      offset=offsets(myid+1)
      write(*,'(A,I3,A,I5,A,I5)')'PROC(',myid,
     &')   I have ',myparts,' parts starting from offset:',offset
      write(*,'(A,I3,A)')'PROC(',myid,')   Initializing first swarm.'
      if(logflag.eq.1) then
      write(11,'(A,I3,A,I5,A,I5)')'PROC(',myid,
     &')   I have ',myparts,' parts starting from offset:',offset
      write(11,'(A,I3,A)')'PROC(',myid,')   Initializing first swarm.'
      endif
c     --------------------------------------
      
      call MPI_TYPE_VECTOR(nparams,1,nparams,
     &             MPI_DOUBLE_PRECISION,mytype,ierr)
      call MPI_TYPE_COMMIT(mytype,ierr)
      call MPI_TYPE_CONTIGUOUS(nparams,MPI_DOUBLE_PRECISION,
     &                         mytype2,ierr)
      call MPI_TYPE_COMMIT(mytype2,ierr)
  
      print*,'PROC(',myid,') Testing MPI!'
      do 111 i=offset,myparts+offset-1
            testx(i)=myid+1
111    continue 
       call MPI_ALLGATHERV(testx(offset),myparts,
     & MPI_DOUBLE_PRECISION,testx,rcount,displs,
     & MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      if(logflag.eq.1) then
      endif
c     print*,'PROC(',myid,') Test array:',(testx(j),j=1,nparts)
      if (myid .eq. MASTER) then
       testflag=0
       do 114 j=1,nparts
        if (testx(j).lt.1) then
           testflag=1
        endif
114     continue
      endif
      if (testflag.eq.1) then
        write(*,*)'Can not go on with this MPI program. '
        write(*,*)'MPI_ALLGATHERV function is not working! '
        write(*,*)'You must recompile the code with a '
        write(*,*)'proper MPI Library! Exiting...'
        call MPI_FINALIZE(ierr)
        stop
      endif
      print*,'MPI test is ok! Good to go! ',myid
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
 
c     --------------------------------------
c       Initialize Default APSO Parameters 
c     --------------------------------------
      gen=0
      wt=0.9
      c1=2.0
      c2=2.0
      c1ok=c1
      c2ok=c2
      sigmamin=0.1
      sigmamax=1.0
      mu=0.0
      gbestid=1
      mygbestid=1
      mygbestval(1,1)=100000.0
      mygbestoldval=100000.0
      gbestval=100000.0
      gbestoldval=100000.0
      gworseid=1
      gworseval=0.0
      mygworseval(1,1)=0.0
      esestate=1
      do 22 i=1,nparts
         pbestval(i)=100000.0
         partdist(i)=0.0
         do 23 j=1,nparams
           gbest(j)=0.0
           mygbest(j)=0.0
           pbest(i,j)=0.0
23       continue
22    continue
      do 24 i=1,4
         ese(i)=0
24    continue
       if (myid .eq. MASTER) then
         write(*,*)'Bounds:'
       endif
       if(logflag.eq.1) then
         write(11,*)'Bounds:'
       endif
      do 25 j=1,nparams
         vmax(j,1) = -0.2d0 * (bounds(j,2)-bounds(j,1))
         vmax(j,2) =  0.2d0 * (bounds(j,2)-bounds(j,1))
        if (myid .eq. MASTER) then
         write(*,*)bounds(j,1),bounds(j,2),
     &     ' min:',vmax(j,1),vmax(j,2)
        endif
        if(logflag.eq.1) then
         write(11,*)bounds(j,1),bounds(j,2),
     &     ' min:',vmax(j,1),vmax(j,2)
        endif
25    continue   
        
      randchunk=2
      randnum=10000
      leftb=-10000000
      rightb=10000000
      brng=VSL_BRNG_MT19937
      method=VSL_METHOD_DUNIFORM_STD_ACCURATE
      seedcells1=loc(getseed())
c     print*,seedcells1
      seedcells2=(values(5)+1)*(values(6)+2)*(values(7)+7)
      seed=seedcells1*seedcells2
      if (seed.eq.0) then 
        seedcells1=loc(getseed())
        if (seedcells1.eq.0) then 
            seedcells1=3
        endif
        seed=seedcells1*seedcells2
      endif
      seed=seed+myid
c     call MPI_BCAST(seed,1,MPI_DOUBLE_PRECISION,MASTER,
c    &MPI_COMM_WORLD,ierr)
c     if (MOD(seed,2).eq.0) then
c        seed=seed-1
c     endif
      write(*,'(A,I3,A,I)')'PROC(',myid,')  Seed:',seed
      if(logflag.eq.1) then
        write(11,'(A,I3,A,I)')'PROC(',myid,')  Seed:',seed
      endif
      errcode=vslnewstream(stream,brng,seed)
      errcode=vdrnguniform(method,stream,randnum,seedzigot,leftb,rightb)
      errcode=vsldeletestream(stream)
     
      seedcount=2 
      leftb=0
      rightb=1
      randnum=3*nparts*nparams*(randchunk)+(randchunk+1)*6
      if (randnum.gt.psornd) then
        if (myid .eq. MASTER) then
         write(*,*)'Can not go on with this random number size ',randnum
         write(*,*)'The limit for one time random numbers is',psornd
         write(*,*)'Rise the psornd parameter in program. Exiting...'
        endif
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         call MPI_FINALIZE(ierr)
         stop
      endif
      brng=VSL_BRNG_MT19937
      method=VSL_METHOD_DUNIFORM_STD_ACCURATE
      seedcells=seedzigot(seedcount)
c     seedcells=loc(getseed())
c     print*,seedcells
      seedcells2=(values(5)+1)*(values(6)+2)*(values(7)+7)
      seed=seedcells1*seedcells2
      if (seed.eq.0) then 
        seedcells1=loc(getseed())
        if (seedcells1.eq.0) then 
            seedcells1=3
        endif
        seed=seedcells1*seedcells2
      endif
      seed=seed+myid
c     if (MOD(seed,2).eq.0) then
c        seed=seed-1
c     endif
      write(*,'(A,I)')'Seed:',seed
      if(logflag.eq.1) then
        write(11,'(A,I)')'Seed:',seed
      endif
      errcode=vslnewstream(stream,brng,seed)
      errcode=vdrnguniform(method,stream,randnum,rands,leftb,rightb)
      errcode=vsldeletestream(stream)
c     call MPI_BCAST(rands,randnum,MPI_DOUBLE_PRECISION,MASTER,
c    &MPI_COMM_WORLD,ierr)

      k=1
      nlast=randchunk+1
      do 26 i=1,nlast
         if (rands(k+1).lt.0.5) then
            deltasign1(i)=-1
         else
            deltasign1(i)=1
         endif
         if (rands(k+2).lt.0.5) then
            deltasign2(i)=-1
         else
            deltasign2(i)=1
         endif
         k=k+2
26    continue
      nlast=randchunk+1
      do 27 i=1,nlast
         delta1(i)=(0.05+0.05*rands(k+1))
         delta2(i)=(0.05+0.05*rands(k+2))
         k=k+2
27    continue
      nlast=nparts*nparams*(randchunk)
      do 57 i=1,nlast
         rand1(i)=rands(k+1)
         rand2(i)=rands(k+2)
         k=k+2
57    continue

c     --------------------------------------
c       Generate Initial Swarm Randomly
c     --------------------------------------
      k=k+((offset-1)*nparams)*2
      do 28 i=offset,myparts+offset-1
         do 29 j=1,nparams
            ind=j+(i-1)*nparams
            parx(ind)=(bounds(j,2)-bounds(j,1))*rands(k+1)+bounds(j,1)
            parv(i,j)=(vmax(j,2)-vmax(j,1))*rands(k+2)+vmax(j,1)
            funcprm(j)=parx(ind)
            pbest(i,j)=parx(ind)
            k=k+2
29       continue
c       write(*,*)'PROC(',myid,') part:',i,':',(funcprm(j),j=1,nparams)
      call eammain(rms,funcprm,nparams,lcharnum,fitvalues,expdata,
     .eamfileout,snrhoin,sdrhoin,snrin,sdrin,sfrhoin,szrin,srhorin,
     .eamftyp,prmflag,prntflag,logflag,myid,i,weightdata)
        write(*,*)'PROC(',myid,') part: ',i,'  rms: ',rms
         fval(i)=rms
         pbestval(i)=fval(i)
c     ------------------------------------------
c      Find Particle with Globally Best Fitness
c     ------------------------------------------
          if (fval(i).lt.gbestval) then
             mygbestval(1,1)=fval(i)
             gbestoldval=fval(i)
             gbestval=fval(i)
             gbestid=i
             mygbestval(2,1)=myid
             do 33 j=1,nparams
                ind=j+(i-1)*nparams
                gbest(j)=parx(ind)
                gbestold(j)=parx(ind)
33           continue
          endif
          if (fval(i).gt.gworseval) then
             mygworseval(1,1)=fval(i)
             gworseval=fval(i)
             gworseid=i
             mygworseval(2,1)=myid
             do 34 j=1,nparams
                ind=j+(i-1)*nparams
                gworse(j)=parx(ind)
34           continue
          endif
28    continue
       
c     --------------------------------------
c       APSO main loop starts
c     --------------------------------------
      do 77 gen=0,ngens
       if (myid .eq. MASTER) then
        write(*,'(A,I3,A,I6)')'PROC(',myid,')  Generation:',gen
       endif
       gn=MOD(gen,randchunk)
      
c      ---------------------------------------------
c       Rebuilding Random Number Arrays Starts Here
c      ---------------------------------------------
c       Rebuild every randchunk (default=50) generations
       if (gn.eq.0) then
         leftb=0
         rightb=1
         randnum=3*nparts*nparams*(randchunk)+(randchunk+1)*6
         call date_and_time(date,time,zone,values)
         brng=VSL_BRNG_MT19937
         method=VSL_METHOD_DUNIFORM_STD_ACCURATE
         seedcount=seedcount+1
         if (MOD(seedcount,10000).eq.0) then
             seedcount=2
         endif
         seedcells1=seedzigot(seedcount)
         seedcells2=(values(5)+1)*(values(6)+2)*(values(7)+7)
c        if (myid .eq. MASTER) then
           write(*,*)'PROC(',myid,')  timeseed:',seedcells2,
     &' memseed:',seedcells1
c        endif
         seed=seedcells1*seedcells2
         if (seed.eq.0) then 
            seedcount=seedcount+1
            if (MOD(seedcount,10000).eq.0) then
                seedcount=2
            endif
            seedcells1=seedzigot(seedcount)
            if (seedcells1.eq.0) then 
                seedcells1=3
            endif
            seed=seedcells1*seedcells2
         endif
         seed=seed+myid
c        if (MOD(seed,2).eq.0) then
c           seed=seed-1
c        endif
        write(*,'(A,I,A,E17.10,A,E17.10)')'   Random rebuilding. Seed:',
     .  seed,' rand1:',rands(10),' rand2',rands(11)
         errcode=vslnewstream(stream,brng,seed)
         errcode=vdrnguniform(method,stream,randnum,rands,leftb,rightb)
         errcode=vsldeletestream(stream)
         k=1
         nlast=randchunk+1
         do 35 i=1,nlast
            if (rands(k+1).lt.0.5) then
               deltasign1(i)=-1
            else
               deltasign1(i)=1
            endif
            if (rands(k+2).lt.0.5) then
               deltasign2(i)=-1
            else
               deltasign2(i)=1
            endif
            k=k+2
35       continue
         nlast=randchunk+1
         do 36 i=1,nlast
            delta1(i)=(0.05+0.05*rands(k+1))
            delta2(i)=(0.05+0.05*rands(k+2))
            k=k+2
36       continue
         nlast=nparts*nparams*(randchunk)
         do 37 i=1,nlast
            rand1(i)=rands(k+1)
            rand2(i)=rands(k+2)
            k=k+2
37       continue
       endif
c      ----------------------------------------------
c       Rebuilding Random Number Arrays are Finished
c      ----------------------------------------------

       if (gen.eq.0) then
c     --------------------------------------
c      Generate Ramdom Swarm for parx again
c     --------------------------------------
      k=1
      k=k+((offset-1)*nparams)*2
      do 58 i=offset,myparts+offset-1
        if (fval(i).ge.1000.0) then
         do 59 j=1,nparams
            ind=j+(i-1)*nparams
            parx(ind)=(bounds(j,2)-bounds(j,1))*rands(k+1)+bounds(j,1)
            parv(i,j)=(vmax(j,2)-vmax(j,1))*rands(k+2)+vmax(j,1)
            funcprm(j)=parx(ind)
            k=k+2
59       continue
c     write(*,'(A,I3,A,I3,A,7(E18.10))')'PROC(',myid,')   parx(',
c    &i,')',(parx(j+(i-1)*nparams),j=1,nparams)
      call eammain(rms,funcprm,nparams,lcharnum,fitvalues,expdata,
     .eamfileout,snrhoin,sdrhoin,snrin,sdrin,sfrhoin,szrin,srhorin,
     .eamftyp,prmflag,prntflag,logflag,myid,i,weightdata)
      write(*,'(A,I3,A,E18.10)')'PROC(',myid,')   rms:',rms
         fval(i)=rms
        endif
c     ------------------------------------------
c      Find Particle with Globally Best Fitness
c     ------------------------------------------
          if (fval(i).lt.pbestval(i)) then
            pbestval(i)=fval(i)
            do 55 j=1,nparams
               ind=j+(i-1)*nparams
               pbest(i,j)=parx(ind)
55          continue
          endif
          if (fval(i).lt.gbestval) then
             tol=abs(mygbestval(1,1)-fval(i))
             mygbestval(1,1)=fval(i)
             gbestoldval=gbestval
             gbestval=fval(i)
             gbestid=i
             mygbestval(2,1)=myid
             do 53 j=1,nparams
                ind=j+(i-1)*nparams
                gbestold(j)=gbest(j)
                gbest(j)=parx(ind)
53           continue
          endif
          if (fval(i).gt.gworseval) then
             mygworseval(1,1)=fval(i)
             gworseval=fval(i)
             gworseid=i
             mygworseval(2,1)=myid
             do 54 j=1,nparams
                ind=j+(i-1)*nparams
                gworse(j)=parx(ind)
54           continue
          endif
58    continue
c     --------------------------------------
c       Find global gbest through all procs
c     --------------------------------------
      call MPI_ALLREDUCE( mygbestval, gbestpack, 1,
     &MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_WORLD,ierr)
       gbestval=gbestpack(1,1)
       gbestloc=gbestpack(2,1)
c     write(*,*)'PROC(',myid,') ALLREDUCE1 PASSED'
      call MPI_BCAST(gbest,nparams,MPI_DOUBLE_PRECISION,gbestloc,
     &MPI_COMM_WORLD,ierr)
c     write(*,*)'PROC(',myid,') BCAST1 PASSED'
      call MPI_BCAST(gbestid,1,MPI_INTEGER,gbestloc,
     &MPI_COMM_WORLD,ierr)
c     write(*,*)'PROC(',myid,') BCAST2 PASSED'
c     --------------------------------------
c       Find global gworse through all procs
c     --------------------------------------
      call MPI_ALLREDUCE( mygworseval, gworsepack, 1,
     &MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_WORLD,ierr)
       gworseval=gworsepack(1,1)
       gworseloc=gworsepack(2,1)
c     write(*,*)'PROC(',myid,') ALLREDUCE2 PASSED'
      call MPI_BCAST(gworse,nparams,MPI_DOUBLE_PRECISION,gworseloc,
     &MPI_COMM_WORLD,ierr)
c     write(*,*)'PROC(',myid,') BCAST3 PASSED'
      call MPI_BCAST(gworseid,1,MPI_INTEGER,gworseloc,
     &MPI_COMM_WORLD,ierr)
c     write(*,*)'PROC(',myid,') BCAST4 PASSED'
c     if (myid .eq. MASTER) then 
       tol=abs(gbestval-gbestoldval)
      if (gbestid.ge.offset.and.gbestid.le.(offset+myparts-1)) then
         fval(gbestid)=gbestval
         do 62 j=1,nparams
            parx(j+(gbestid-1)*nparams)=gbest(j)
62       continue
      endif
      if (gworseid.ge.offset.and.gworseid.le.(offset+myparts-1)) then
         fval(gworseid)=gworseval
         do 63 j=1,nparams
            parx(j+(gworseid-1)*nparams)=gworse(j)
63       continue
      endif
       if (nparams.eq.5) then
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,5(E17.10))')'PROC(',myid,
     &')   gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       else if (nparams.eq.6) then
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,6(E18.10))')'PROC(',myid,
     &')   gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       else if (nparams.eq.7) then
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,7(E17.10))')'PROC(',myid,
     &')   gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       else if (nparams.eq.8) then
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,8(E17.10))')'PROC(',myid,
     &')   gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       else if (nparams.eq.9) then
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,9(E17.10))')'PROC(',myid,
     &')   gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       else
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,17(E17.10))')'PROC(',myid,
     &')   gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       endif
c     --------------------------------------
c      All Gather parx from all to all procs
c     --------------------------------------
      call MPI_ALLGATHERV(parx((offset-1)*nparams+1),myparts,
     &mytype2,parx,rcount,displs,
     &mytype2,MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHERV(pbestval(offset),myparts,
     &MPI_DOUBLE_PRECISION,pbestval,rcount,displs,
     &MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHERV(fval(offset),myparts,
     &MPI_DOUBLE_PRECISION,fval,rcount,displs,
     &MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
c     --------------------------------------
c      End Generation of Ramdom Swarm
c     --------------------------------------
      endif

c     if(gen.gt.0) then
       if (myid.eq.MASTER) then
        write(*,*)'PROC:',myid,'pbestval:',(pbestval(j),j=1,nparts)
        write(*,*)'PROC:',myid,'partsval:',(fval(j),j=1,nparts)
        do i=1,nparts
         write(*,*)'PROC:',myid,'part:',i,'->',
     &(parx(j+(i-1)*nparams),j=1,nparams)
        enddo
      endif
       if (gbestval.ge.1000.0) then
         write(*,*)'Bad swarm condition: GBEST<1000 condition exceeded!'
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         call MPI_FINALIZE(ierr)
         stop 
       endif
c     endif

c     Evolutionary State Estimation (ESE) stage

      distmin=1000000.0
      distmax=-1.0
      distg=10000.0
      do 38 i=1,nparts
       distpart=0.0
       do 39 j=1,nparts
         difparm=0.0
         if (i.ne.j) then
           do 40 k=1,nparams
             ind=k+(i-1)*nparams
             inn=k+(j-1)*nparams
             difparm=difparm+(parx(ind)-parx(inn))*(parx(ind)-parx(inn))
40         continue
         endif  
         distpart=distpart+sqrt(difparm)
39     continue
       distpart=distpart/(nparts-1)
       if (distpart.lt.distmin) then
           distmin=distpart
       endif
       if (distpart.gt.distmax) then
           distmax=distpart
       endif
       if (i.eq.gbestid) then
           distg=distpart
       endif
38    continue
     
c     Calculate Evolutionary Factor
      if (gen .eq. 0) then 
         evofact=1.0
      else
         evofact=(distg-distmin)/(distmax-distmin)
      endif
    
      if (evofact .gt. 1.0) then
          evofact=1.0
      endif 
      if (evofact .lt. 0.0) then
          evofact=0.0
      endif 

c     if (myid .eq. MASTER) then
c       write(*,'(A,I3,A,E17.10)')'PROC(',myid,')  evofact:',evofact
c     endif

c     Evolutionary State Calculation
      if (evofact.ge.0.0 .and. evofact.le.0.4) then
        ese(1)=0.0
        ese(4)=0.0
        if (evofact.le.0.1) then
            ese(2)=0.0
            ese(3)=1.0
        elseif (evofact.le.0.3) then
            ese(3)=-(5.0*evofact)+1.5
            if (evofact.le.0.2) then
                ese(2)=0.0
            else
                ese(2)=(10.0*evofact)-2.0
            endif
        else
            ese(2)=1.0
            ese(3)=0.0
        endif
      elseif (evofact.gt.0.4 .and. evofact.le.0.6) then
        ese(1)=(5.0*evofact)-2.0
        ese(2)=-(5.0*evofact)+3.0
        ese(3)=0.0
        ese(4)=0.0
      elseif (evofact.gt.0.6 .and. evofact.le.1.0) then
        ese(2)=0.0
        ese(3)=0.0
        if (evofact.le.0.7) then
            ese(1)=1.0
            ese(4)=0.0
        elseif (evofact.le.0.9) then
            ese(4)=(5.0*evofact)-3.5
            if (evofact.le.0.8) then
                ese(1)=-(10.0*evofact)+8.0
            else
                ese(1)=0.0
            endif
        else
            ese(4)=1.0
        endif
      endif

      esehig=0.0
      eselow=0.0
      esehigloc=1
      eselowloc=1
      do 41 i=1,4
          if(ese(i).gt.esehig) then
             esehig=ese(i)
             esehigloc=i
          endif
          if(ese(i).lt.esehig .and. eselow.lt.esehig
     &       .and. ese(i).ne.0.0) then
             eselow=ese(i)
             eselowloc=i
          endif
41    continue

      if (gen .eq. 0) then
          esestate=1
      else
          if (esestate .ne. esehigloc) then
              esenext=MOD(esestate+1,4)
              esenext2=MOD(esestate+2,4)
              esenext3=MOD(esestate+3,4)
              eseprev=MOD(esestate-1,4)
              if(esenext.eq.0) then
                 esenext=4
              endif
              if(esenext2.eq.0) then
                 esenext2=4
              endif
              if(esenext3.eq.0) then
                 esenext3=4
              endif
              if(eseprev.eq.0) then
                 eseprev=4
              endif
              if (esehigloc.eq.esenext) then
                  esestate=esenext
              elseif (esehigloc.eq.esenext2) then
                  esestate=esenext
              elseif (esehigloc.eq.esenext3) then
                  esestate=esenext
              else
                  esestate=esestate
              endif
          endif
      endif


c     Adaptive Parameter Control Stage for PSO

c     Control of Acceleration Coefficients

c        State       Strategy        c1                  c2
c     Exploration   esestate=1    Increase            Decrease
c     Exploitation  esestate=2    Increase slightly   Decrease slightly
c     Convergence   esestate=3    Increase slightly   Increase slightly
c     Jumping-out   esestate=4    Decrease            Increase
      if (myid .eq. MASTER) then
      if (esestate.eq.1) then
          statechar=' Exploration  '
      elseif (esestate.eq.2) then
          statechar=' Exploitation '
      elseif (esestate.eq.3) then
          statechar=' Convergence  '
      elseif (esestate.eq.4) then
          statechar=' Jumping-out  '
      else
          statechar='WARN! No state'
      endif
          write(*,'(A,I3,A,A)')'PROC(',myid,') ',statechar
      endif

      if(gen .eq. 0) then
         c1 = 2.0
         c2 = 2.0
         c1ok = c1
         c2ok = c2
         c1c2tot = c1ok + c2ok
      else
         if (esestate.eq.2) then
             c2upd = c2 - ( 0.5 * delta2(gn+1) )
             c1upd = c1 + ( 0.5 * delta1(gn+1) )
         elseif (esestate.eq.3) then
             c2upd = c2 + ( 0.5 * delta2(gn+1) )
             c1upd = c1 + ( 0.5 * delta1(gn+1) )
         elseif (esestate.eq.4) then
             c2upd = c2 + delta2(gn+1)
             c1upd = c1 - delta1(gn+1)
         else
             c2upd = c2 - delta2(gn+1)
             c1upd = c1 + delta1(gn+1)
         endif
         if (c1upd .ge. 1.5) then
             c1ok = c1upd
c        else
c            c1ok = c1
         endif
         if (c2upd .ge. 1.5) then
             c2ok = c2upd
c        else
c            c2ok = c2
         endif
         c1c2tot = c1ok + c2ok
         if(c1c2tot .gt. 4.0) then
            c2 = c2ok*(4.0/c1c2tot)
            c1 = c1ok*(4.0/c1c2tot)
         else
            c2 = c2ok
            c1 = c1ok
         endif
      endif
      
c     write(*,'(A,I3,2(A,E17.10))')'PROC(',myid,')  c1:',c1,'  c2:',c2


c     Elitist Learning Strategy (ELS) stage
c     Only for Convergence state of ESE
      if (esestate.eq.3) then
        elsupd=0
        leftb=0
        rightb=1
        randnum1=2*nparams+2
        call date_and_time(date,time,zone,values)
        brng=VSL_BRNG_MT19937
        method=VSL_METHOD_DUNIFORM_STD_ACCURATE
        seedcount=seedcount+1
        if (MOD(seedcount,10000).eq.0) then
            seedcount=2
        endif
        seedcells1=seedzigot(seedcount)
        seedcells2=(values(5)+1)*(values(6)+2)*(values(7)+7)
        seed=seedcells1*seedcells2
        if (seed.eq.0) then 
           seedcount=seedcount+1
           seedcells1=seedzigot(seedcount)
           if (seedcells1.eq.0) then 
                seedcells1=3
           endif
           seed=seedcells1*seedcells2
        endif
        seed=seed+myid
        errcode=vslnewstream(stream,brng,seed)
        errcode=vdrnguniform(method,stream,randnum1,randm,leftb,rightb)
        errcode=vsldeletestream(stream)
c       write(*,'(A,E17.10,A,E17.10,$)')'U Seed: ',seed,
c    &'randm: ',randm(2)
        ndim=nparams
        perturb=randm(2)
        perturb=perturb*ndim
        dims=ceil(perturb,ndim)
c       write(*,'(A,I3,A,E17.10)')'PROC(',myid,')  randm: ',randm(2)
        do 42 j=1,dims
           dimd=randm(j+2)*nparams
           di(j)=ceil(dimd,nparams)
42      continue
c       write(*,'(A,I3,2(A,I3))')'PROC(',myid,')  di(1): ',
c    &di(1),' di(2): ',di(2)
        sigma=sigmamax-(sigmamax-sigmamin)*(gen/ngens)
        sigma=sigma*sigma
        randnum2=2*dims+2
        brng=VSL_BRNG_MT19937
        method=VSL_METHOD_DGAUSSIAN_ICDF
        seedcount=seedcount+1
        if (MOD(seedcount,10000).eq.0) then
            seedcount=2
        endif
        seedcells1=seedzigot(seedcount)
        seedcells2=(values(5)+1)*(values(6)+2)*(values(7)+7)
        seed=seedcells1*seedcells2
        if (seed.eq.0) then 
           seedcount=seedcount+1
           seedcells1=seedzigot(seedcount)
           if (seedcells1.eq.0) then 
                seedcells1=3
           endif
           seed=seedcells1*seedcells2
        endif
        seed=seed+myid
c       if (MOD(seed,2).eq.0) then
c          seed=seed-1
c       endif
c       write(*,'(A,E17.10)')'G Seed:',seed
        errcode=vslnewstream(stream,brng,seed)
        errcode=vdrnggaussian(method,stream,randnum2,randn,mu,sigma)
        errcode=vsldeletestream(stream)
c     write(*,'(A,I3,A,I3,2(A,E17.10))')'PROC(',myid,
c    .')  Random Gaussian, Seed: ',seed,' dims: ',
c    .dims,' rand1: ',randn(3),' rand2 ',randn(4)
        do 43 j=1,nparams
           ppart(j)=gbest(j)
43      continue
        do 44 j=1,dims
           gauss=randn(j+2)
           bounddif=bounds(di(j),2)-bounds(di(j),1)
           ppart(di(j))=ppart(di(j))+bounddif*gauss
           if (ppart(di(j)).gt.bounds(di(j),2)) then
               ppart(di(j))=bounds(di(j),2)
           endif
           if (ppart(di(j)).lt.bounds(di(j),1)) then
               ppart(di(j))=bounds(di(j),1)
           endif
44      continue
        prmno=0
      call eammain(rms,ppart,nparams,lcharnum,fitvalues,expdata,
     .eamfileout,snrhoin,sdrhoin,snrin,sdrin,sfrhoin,szrin,srhorin,
     .eamftyp,prmflag,prntflag,logflag,myid,prmno,weightdata)
        ppartval=rms
        if (dyntyp.gt.0) then
        write(*,*)'PROC:',myid,'ELS dynamic gbest update'
        prmno=0
      call eammain(rms,gbest,nparams,lcharnum,fitvalues,expdata,
     .eamfileout,snrhoin,sdrhoin,snrin,sdrin,sfrhoin,szrin,srhorin,
     .eamftyp,prmflag,prntflag,logflag,myid,prmno,weightdata)
        gbestval=rms
        endif
        if (ppartval.lt.gbestval) then
             elsupd=1
             gbestoldval=gbestval
             gbestval=ppartval
             mygbestval(1,1)=ppartval
             mygbestval(2,1)=myid
             do 45 j=1,nparams
                gbestold(j)=gbest(j)
                gbest(j)=ppart(j)
45           continue
         else
             gworseval=ppartval
             mygworseval(1,1)=ppartval
             mygworseval(2,1)=myid
             do 46 j=1,nparams
                gworse(j)=ppart(j)
46           continue
         endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c     --------------------------------------
c       Find global gbest through all procs
c     --------------------------------------
      call MPI_ALLREDUCE( mygbestval, gbestpack, 1,
     &MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_WORLD,ierr)
       gbestval=gbestpack(1,1)
       gbestloc=gbestpack(2,1)
c     write(*,*)'PROC(',myid,') ALLREDUCE1 PASSED'
      call MPI_BCAST(gbest,nparams,MPI_DOUBLE_PRECISION,gbestloc,
     &MPI_COMM_WORLD,ierr)
c     write(*,*)'PROC(',myid,') BCAST1 PASSED'
      call MPI_BCAST(gbestid,1,MPI_INTEGER,gbestloc,
     &MPI_COMM_WORLD,ierr)
c     write(*,*)'PROC(',myid,') BCAST2 PASSED'
c     --------------------------------------
c       Find global gworse through all procs
c     --------------------------------------
      call MPI_ALLREDUCE( mygworseval, gworsepack, 1,
     &MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_WORLD,ierr)
       gworseval=gworsepack(1,1)
       gworseloc=gworsepack(2,1)
c     write(*,*)'PROC(',myid,') ALLREDUCE2 PASSED'
      call MPI_BCAST(gworse,nparams,MPI_DOUBLE_PRECISION,gworseloc,
     &MPI_COMM_WORLD,ierr)
c     write(*,*)'PROC(',myid,') BCAST3 PASSED'
      call MPI_BCAST(gworseid,1,MPI_INTEGER,gworseloc,
     &MPI_COMM_WORLD,ierr)
c     write(*,*)'PROC(',myid,') BCAST4 PASSED'

      if (gbestid.ge.offset.and.gbestid.le.(offset+myparts-1)) then
         fval(gbestid)=gbestval
         do 60 j=1,nparams
            parx(j+(gbestid-1)*nparams)=gbest(j)
60       continue
      endif
      if (gworseid.ge.offset.and.gworseid.le.(offset+myparts-1)) then
         fval(gworseid)=gworseval
         do 61 j=1,nparams
            parx(j+(gworseid-1)*nparams)=gworse(j)
61       continue
      endif
      tol=abs(gbestval-gbestoldval)
      if (myid .eq. MASTER .and. elsupd .eq. 1) then 
       elsupd=0
       if (nparams.eq.5) then
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,5(E17.10))')'PROC(',myid,
     &') ELS UPDATED  gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       else if (nparams.eq.6) then
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,6(E18.10))')'PROC(',myid,
     &') ELS UPDATED  gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       else if (nparams.eq.7) then
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,7(E17.10))')'PROC(',myid,
     &') ELS UPDATED  gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       else if (nparams.eq.8) then
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,8(E17.10))')'PROC(',myid,
     &') ELS UPDATED  gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       else if (nparams.eq.9) then
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,9(E17.10))')'PROC(',myid,
     &') ELS UPDATED  gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       else
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,17(E17.10))')'PROC(',myid,
     &') ELS UPDATED  gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       endif
      endif
      call MPI_ALLGATHERV(parx((offset-1)*nparams+1),myparts,
     &mytype2,parx,rcount,displs,
     &mytype2,MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHERV(pbestval(offset),myparts,
     &MPI_DOUBLE_PRECISION,pbestval,rcount,displs,
     &MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHERV(fval(offset),myparts,
     &MPI_DOUBLE_PRECISION,fval,rcount,displs,
     &MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      endif
c     ELS ends here

c     Adaptive Control of the Inertia Weight
      if(gen.eq.0) then
         wt=0.9
      else
         wt=(1.0)/(1.0+(1.5*exp(-2.6*evofact)))
      endif

c     --------------------------------------
c         Update Swarm Particle X and V 
c     --------------------------------------
c      do 47 i=1,nparts
      do 47 i=offset,myparts+offset-1
       evaluate=1
c      write(*,'(3(A,I3))')'PROC(',myid,
c    &')   evaluate:',i,'  ',evaluate
       do 48 j=1,nparams
        ind=j+(i-1)*nparams
        randid=(gn*nparts*nparams)+ind
        parv(i,j)=wt*parv(i,j)+c1*rand1(randid)*(pbest(i,j)-parx(ind))
     &                        +c2*rand2(randid)*(gbest(j)-parx(ind))
        if ( parv(i,j) .gt. vmax(j,2) ) then
             parv(i,j) = vmax(j,2)
        endif
        if ( parv(i,j) .lt. vmax(j,1) ) then
             parv(i,j) = vmax(j,1)
        endif
        parx(ind)=parx(ind)+parv(i,j)
        if ( parx(ind) .lt. bounds(j,1) ) then 
             evaluate=0
c       else
c            evaluate=1
        endif
        if ( parx(ind) .gt. bounds(j,2) ) then
             evaluate=0
c       else
c            evaluate=1
        endif
        funcprm(j)=parx(ind)
48     continue
c      write(*,'(3(A,I3))')'PROC(',myid,
c    &')   evaluate:',i,'  ',evaluate
       if (evaluate.eq.1) then
      call eammain(rms,funcprm,nparams,lcharnum,fitvalues,expdata,
     .eamfileout,snrhoin,sdrhoin,snrin,sdrin,sfrhoin,szrin,srhorin,
     .eamftyp,prmflag,prntflag,logflag,myid,i,weightdata)
        write(*,*)'PROC(',myid,') part: ',i,'  rms: ',rms
          fval(i)=rms
c      Check if pBest is still the pBest
        if (dyntyp.gt.0) then
          write(*,*)'PROC:',myid,'APSO MAIN dynamic gbest update'
          do 49 j=1,nparams
                funcprm(j)=pbest(i,j)
49        continue
      call eammain(rms,funcprm,nparams,lcharnum,fitvalues,expdata,
     .eamfileout,snrhoin,sdrhoin,snrin,sdrin,sfrhoin,szrin,srhorin,
     .eamftyp,prmflag,prntflag,logflag,myid,i,weightdata)
          pbestval(i)=rms
        endif
c      Find pBest and gBest
          if (fval(i).lt.pbestval(i)) then
            pbestval(i)=fval(i)
            do 50 j=1,nparams
               ind=j+(i-1)*nparams
               pbest(i,j)=parx(ind)
50          continue
          endif
c      Check if gBest is still the gBest
          if (dyntyp.gt.0) then
      call eammain(rms,gbest,nparams,lcharnum,fitvalues,expdata,
     .eamfileout,snrhoin,sdrhoin,snrin,sdrin,sfrhoin,szrin,srhorin,
     .eamftyp,prmflag,prntflag,logflag,myid,i,weightdata)
          gbestval=rms
          endif
c     ------------------------------------------
c      Find Particle with Globally Best Fitness
c     ------------------------------------------
          if (fval(i).lt.gbestval) then
            gbestoldval=gbestval
            gbestval=fval(i)
            gbestid=i
            mygbestval(1,1)=fval(i)
            mygbestval(2,1)=myid
            do 51 j=1,nparams
               ind=j+(i-1)*nparams
               gbestold(j)=gbest(j)
               gbest(j)=parx(ind)
51          continue
          endif
          if (fval(i).gt.gworseval) then
            gworseid=i
            gworseval=fval(i)
            mygworseval(1,1)=fval(i)
            mygworseval(2,1)=myid
            do 52 j=1,nparams
               ind=j+(i-1)*nparams
               gworse(j)=parx(ind)
52          continue
          endif
       endif
47    continue
c     --------------------------------------
c       Find global gbest through all procs
c     --------------------------------------
      call MPI_ALLREDUCE( mygbestval, gbestpack, 1,
     &MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_WORLD,ierr)
       gbestval=gbestpack(1,1)
       gbestloc=gbestpack(2,1)
c     write(*,*)'PROC(',myid,') ALLREDUCE1 PASSED'
      call MPI_BCAST(gbest,nparams,MPI_DOUBLE_PRECISION,gbestloc,
     &MPI_COMM_WORLD,ierr)
c     write(*,*)'PROC(',myid,') BCAST1 PASSED'
      call MPI_BCAST(gbestid,1,MPI_INTEGER,gbestloc,
     &MPI_COMM_WORLD,ierr)
c     write(*,*)'PROC(',myid,') BCAST2 PASSED'
c     --------------------------------------
c       Find global gworse through all procs
c     --------------------------------------
      call MPI_ALLREDUCE( mygworseval, gworsepack, 1,
     &MPI_2DOUBLE_PRECISION,MPI_MINLOC,MPI_COMM_WORLD,ierr)
       gworseval=gworsepack(1,1)
       gworseloc=gworsepack(2,1)
c     write(*,*)'PROC(',myid,') ALLREDUCE2 PASSED'
      call MPI_BCAST(gworse,nparams,MPI_DOUBLE_PRECISION,gworseloc,
     &MPI_COMM_WORLD,ierr)
c     write(*,*)'PROC(',myid,') BCAST3 PASSED'
      call MPI_BCAST(gworseid,1,MPI_INTEGER,gworseloc,
     &MPI_COMM_WORLD,ierr)
c     write(*,*)'PROC(',myid,') BCAST4 PASSED'
      tol=abs(gbestval-gbestoldval)
      call MPI_ALLGATHERV(parx((offset-1)*nparams+1),myparts,
     &mytype2,parx,rcount,displs,
     &mytype2,MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHERV(pbestval(offset),myparts,
     &MPI_DOUBLE_PRECISION,pbestval,rcount,displs,
     &MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call MPI_ALLGATHERV(fval(offset),myparts,
     &MPI_DOUBLE_PRECISION,fval,rcount,displs,
     &MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      if (myid .eq. MASTER) then 
       if (nparams.eq.5) then
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,5(E17.10))')'PROC(',myid,
     &')   gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       write(*,'(A,I5,A,A,E17.10)')'Generation:',gen,
     &statechar,' GBEST:',gbestval
       else if (nparams.eq.6) then
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,6(E18.10))')'PROC(',myid,
     &')   gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       write(*,'(A,I5,A,A,E18.10)')'Generation:',gen,
     &statechar,' GBEST:',gbestval
       else if (nparams.eq.7) then
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,7(E17.10))')'PROC(',myid,
     &')   gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       write(*,'(A,I5,A,A,E17.10)')'Generation:',gen,
     &statechar,' GBEST:',gbestval
       else if (nparams.eq.8) then
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,8(E17.10))')'PROC(',myid,
     &')   gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       write(*,'(A,I5,A,A,E17.10)')'Generation:',gen,
     &statechar,' GBEST:',gbestval
       else if (nparams.eq.9) then
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,9(E17.10))')'PROC(',myid,
     &')   gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       write(*,'(A,I5,A,A,E17.10)')'Generation:',gen,
     &statechar,' GBEST:',gbestval
       else
       write(*,'(A,I3,A,I3,A,I3,A,E17.10,A,17(E17.10))')'PROC(',myid,
     &')   gbestid:',gbestid,'gbestloc:',gbestloc,
     &'gbestval:',gbestval,' gbest:',(gbest(j),j=1,nparams)
       write(*,'(A,I5,A,A,E17.10)')'Generation:',gen,statechar,
     &' GBEST:',gbestval
       endif
      endif

77    continue
c     --------------------------------------
c       APSO main loop ends
c     --------------------------------------

      if (myid .eq. MASTER) then 
      write(*,*)' '
      write(*,*)'--------------------------------'
      write(*,*)'tolerance:',tol
      write(*,*)'gbest:',gbestval,'params:'
      write(*,*)(gbest(j),j=1,nparams)

 
      write(*,*)'--------------------------------'
      call date_and_time(date,time,zone,values)
      write (*,1)values(3),values(2),values(1),
     &           values(5),values(6),values(7)
      write(dateh,1)INT(values(3)),INT(values(2)),INT(values(1)),
     &              INT(values(5)),INT(values(6)),INT(values(7))
      write(*,*)'--------------------------------'

      write(*,*)'--------------------------------'
      write(*,*)' EAM Main routine call'
      write(*,*)'--------------------------------'
      prmno=0
c     prntflag=1
      call eammain(rms,gbest,nparams,lcharnum,fitvalues,expdata,
     .eamfileout,snrhoin,sdrhoin,snrin,sdrin,sfrhoin,szrin,srhorin,
     .eamftyp,prmflag,printon,logflag,myid,prmno,weightdata)

      write(*,*)'--------------------------------'
      write(*,*)' RMS:',rms
      write(*,*)'--------------------------------'
      call date_and_time(date,time,zone,values)
      write (*,1)values(3),values(2),values(1),
     &           values(5),values(6),values(7)
      write(dateh,1)INT(values(3)),INT(values(2)),INT(values(1)),
     &              INT(values(5)),INT(values(6)),INT(values(7))
      write(*,*)'--------------------------------'
      else
      write(*,*)'PROC(',myid,') WAITING gbest PRINT CALC FOR FINALIZE.'
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
      end


      double precision function getseed()
      double precision x
      getseed=x
      end

      integer function ceil(xi,ndims)
      double precision xi
      integer i,ndims,xf
      xf=ndims
      do 10 i=1,ndims
         if (xi.gt.(i-1)) then
             xf=i
         endif
10    continue
      ceil=xf
      end
      
      double precision function testfunc(u,n,g,functype)
      parameter (psopar=1000,psoprm=100,psorand=10000,
     . psornd=500000,rndprm=300,eamprm=100)
      double precision u(psoprm)
      integer n,g,functype
      double precision rmssum,r
      integer j
      rmssum=0.0
      r=0.0
      if (g.lt.50) then
         r=-5.0
      else
         r=5.0
      endif
      if (functype.eq.1) then
        do 10 j=1,n
         rmssum=rmssum-u(j)*sin(sqrt(abs(u(j))))
10      continue
      else
        do 11 j=1,n
         rmssum=rmssum+(u(j)-r)*(u(j)-r)
11      continue
      endif
      testfunc=rmssum
      end



      
      subroutine eammain(rmssum,params,nparams,lcharnum,fitvalues,
     .expdata,eamfileout,snrhoin,sdrhoin,snrin,sdrin,sfrhoin,
     .szrin,srhorin,eamftyp,prmflag,prntflag,logflag,procnum,
     .prmid,weightdata)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      parameter (psopar=1000,psoprm=100,psorand=10000,
     . psornd=500000,rndprm=300,eamprm=100)
      PARAMETER (pnts=81)
      double precision rmssum,params(eamprm)
      integer lattype
      integer*4 nparams
      integer*4 eamftyp(10)
      character*32 dateh
      character*80 eamfiles(nelmax)
      character*80 eamfileout
      character*80 fitvalues
      character*80 vmdfile
      character*1 elmtyp
      integer readpureeams,pairnum,prntflag,prmflag
      integer logflag
      integer*4 procnum,prmid
      double precision rrr,rrrho
      common /particle/ rv(6,nmax)
      common /unitcell/ ucell(nmax)
      integer ua1,ua2,ua3,ua4
      data conmas/1.0365e-4/
      double precision pi,sqrt2,sqrt3
      double precision rhoinei(nmax,neimax),rhojnei(nmax,neimax)
      common /forces/ f(3,nmax),fe(nmax),e(nmax),stresst(3,3),
     .                slocal(3,3,nmax),phie(nmax)
      common /density/ rho(nmax),fp(nmax)
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      common /interact/ frho(ngrid,nelmax),z2r(ngrid,nelmax,nelmax),
     .  rhor(ngrid,nelmax),drho,dr,rcutsq,nrho,nr
      common /grids/ frhoar(ngridar,nelmax),frhoar1(ngridar,nelmax),
     .  frhoar2(ngridar,nelmax),frhoar3(ngridar,nelmax),
     .  frhoar4(ngridar,nelmax),frhoar5(ngridar,nelmax),
     .  frhoar6(ngridar,nelmax),rhorar(ngridar,nelmax),
     .  rhorar1(ngridar,nelmax),rhorar2(ngridar,nelmax),
     .  rhorar3(ngridar,nelmax),rhorar4(ngridar,nelmax),
     .  rhorar5(ngridar,nelmax),rhorar6(ngridar,nelmax),
     .  z2rar(ngridar,nelmax,nelmax),z2rar1(ngridar,nelmax,nelmax),
     .  z2rar2(ngridar,nelmax,nelmax),z2rar3(ngridar,nelmax,nelmax),
     .  z2rar4(ngridar,nelmax,nelmax),z2rar5(ngridar,nelmax,nelmax),
     .  z2rar6(ngridar,nelmax,nelmax),drhoar,drar,nrhoar,nrar
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      dimension sfrhoin(ngrid,nelmax),sdrhoin(nelmax),snrhoin(nelmax)
      dimension srhorin(ngrid,nelmax),szrin(ngrid,nelmax),sdrin(nelmax),
     1 snrin(nelmax)
      dimension frhoin(ngrid,nelmax),drhoin(nelmax),nrhoin(nelmax)
      dimension rhorin(ngrid,nelmax),zrin(ngrid,nelmax),drin(nelmax),
     1 nrin(nelmax)
      dimension zrtemp(ngrid,nelmax)
      common /eamdef/ fheader(80,nelmax),rcut(nelmax),blat(nelmax),
     . lat(nelmax),elmname
      character elmname(nelmax)*2
      common /inputdef/ nxchain,nyatoms,nzlayer,rnum,rhonum
      common /inputlat/ rhoend,lata0,caratio,latai,latstep,lataf,
     1esub,bmod
      integer nxchain,nyatoms,nzlayer,rnum,rhonum
      double precision lata0,caratio,rhoend,saveper(2,3),saveperlen(3)
      double precision rr,zrr,zrcut,zdrrcut,rhorr,rhorcut,rhodrrcut
      double precision demin1,demin2,demax,dedif
      double precision decal,dephi,defe1,defe2,derho1,derho2,dez2ij
      double precision recal,rei,ref,res,hcpc,hcpc1,hcpc2,c
      dimension rhoeq(nelmax),arhor(nelmax)
      double precision adrin1,adrin2,adrin3
      real*4 aadrin1,aadrin2,aadrin3
      double precision azrin1(ngrid),azrin2(ngrid),azrin3(ngrid)
      double precision azrtemp1(ngrid),azrtemp2(ngrid),azrtemp3(ngrid)
      double precision rdrhoar,rdrar
      real pij,rij,pp
      integer kij,i1,i2,anrin1,anrin2,anrin3,ren
      integer nu,factor,kk
      double precision sx(nelmax),gx(nelmax)
      double precision gs,ga,sa,gb,sb,gc,sc,rmaxguess
      double precision alpha1,beta1,dm1,rm1,rcutguess1,de,re
      double precision alpha2,beta2,dm2,rm2,rcutguess2
      double precision alpha3,beta3,dm3,rm3,rcutguess3
      double precision alphax(nelmax),betax(nelmax)
      double precision newq,apart,astarcut,qstar
      double precision dmx(nelmax),rmx(nelmax),rcutguess(nelmax)
      double precision esub,bmod,astar,binde,bindf,bindfmod
      double precision lata,latai,lataf,latstep,epslon,fqn,latacut
      double precision cof1,cof2,cof3,cof4,latano,rcutoff,rcutoffsq
      double precision phi,phi2,eqrho,eqphi,eqphi2,anorm,rho,rhoij
      double precision kai,fa(ngrid),rhoa(ngrid),rhoeam(ngrid),ream
      double precision rhomod(ngrid),phimod(ngrid),rcutall
      double precision rdis(nmax,nmax,2),rdist,xa,disnei,zr(ngrid)
      double precision rdisd(nmax,nmax,3),phiequil(ngrid),phidr
      integer latno,vacatom,vacatom1,vacatom2,is,id,i,j,k,neigh(nmax)
      double precision diski(3),rhos,rhod,fa2(ngrid)
      integer iat,jat,jj,ii,nii,vacnum1(10),vacnum2(10),migatoms(2)
      double precision fspline(ngrid),vacmig1(10),vacmig2(10)
      double precision yp1,ypn,y2p(ngrid),finter,rhofirst,fafirst
      double precision c11(10),c12(10),c44(10),bbb(10),bbb2(10)
      double precision c13(10),c33(10),c55(10),c66(10),bbb3(10),bbb4(10)
      double precision embed1(10),embed2(10),vacform1(10),vacform2(10)
      double precision embed(5),datome(10),datomr(10),ca(nmax),c66ca(10)
      double precision omega0,omega,ev_joule,shifts(3),xshift,yshift
      double precision z2ij,z2pij,phiij,lrho,lfrho,lphi,mev_to_mjoule
      double precision totenergy,minenergy,latfrho,latphi,latis
      double precision etot0,etot1,etot2,twin111(4),sf111(4)
      double precision sisf111(4),apb100(4),apb111(4),fmax
      double precision b2alat(10*nelmax),b2ecoh(10*nelmax)
      double precision l10alat(10*nelmax),l10ecoh(10*nelmax)
      double precision l11alat(10*nelmax),l11ecoh(10*nelmax)
      double precision l12alat(10*nelmax),l12ecoh(10*nelmax)
      double precision l13alat(10*nelmax),l13ecoh(10*nelmax)
      double precision hexalat(10*nelmax),hexecoh(10*nelmax)
      double precision hexca(10*nelmax),expvals(100),fitvals(200)
      double precision errors(100),expdata(100),weightdata(200)
      double precision fpp(nmax),econ(7,nmax),uu(7),ww(7),vv(7),
     .  bulkmod(7,nmax),esublim(nmax),evacf(nmax),dcon(3,3)
      integer nknot,fileformat,phitype,l1,b2,fixcutoff,ntyps
      integer l1fcc,stacks(3),zplanes,zplns(2)
      integer xplanes,yplanes,t1,t2,iter,typs(10),tn,tnum(nelmax)
      integer at1,at2,fitnums
      character*8  date
      character*10 time
      character*5  zone
      integer values(8)
      integer lcharnum,istop,iclock
      double precision atobohr,kxr,kxr1,kxr2,kxr3,m1,m2,m3,m4,mm
      double precision kpi,kpi1,kpi2,kpi4,kpnt(3)
      double precision latak,kpoints(200,3)
      integer points,kx,ky,kz,kcount,brillouinzone,a,b
      double precision pcell(3,32),pvec(3)
      integer ptyp(32),patm,allphonon
      complex*16 self1(3,3),self2(3,3),ikr
      complex*16 dmat11(3,3,pnts),dmat12(3,3,pnts)
      complex*16 dmat21(3,3,pnts),dmat22(3,3,pnts)
      integer slf1,slf2,d11,d12,d21,d22,kr,kc,stck(3)
      double precision phg,phxt,phxl,phlt,phll,phkt1,phkt2,phkl
c      LAPACK Parameters
      INTEGER NN,NNN,LDA,LDB,LWMAX,LWMAXB
      PARAMETER ( NN = 6 , NNN = 3 )
      PARAMETER ( LDA = NN , LDB = NNN)
      PARAMETER ( LWMAX = 1000, LWMAXB = 1000 )
c      LAPACK Local Scalars
      INTEGER INFO,LWORK,INFOB,LWORKB
c      LAPACK Local Arrays
      DOUBLE PRECISION W(NN),RWORK(3*NN-2)
      DOUBLE PRECISION WB(NNN),RWORKB(3*NNN-2)
      COMPLEX*16 dmic(LDA,NN),WORK(LWMAX)
      COMPLEX*16 dmicb(LDB,NNN),WORKB(LWMAXB)
c      External Subroutines
      EXTERNAL ZHEEV


c     electron-volt to joule conversion constant
      ev_joule=0.624150965
      mev_to_mjoule=16.0217646
      pi=3.1415926535897932384626433832795028841971693993751058209749445
      sqrt3=1.7320508075688772935274463415058723669428052538103806280558
      sqrt2=1.4142135623730950488016887242096980785696718753769480731767

c      at1,at2,phitype,rhotype,embtype,eamalloy,fileformat unpacking
c      ntypes=eamftyp(1)
       at1=eamftyp(2)
       at2=eamftyp(3)
       phitype=eamftyp(4)
       rhotype=eamftyp(5)
       embtype=eamftyp(6)
       eamalloy=eamftyp(7)
       fileformat=eamftyp(8)
       lattype=eamftyp(9)
       fixcutoff=eamftyp(10)
     
      if (eamalloy .gt. 0) then
c     Read EAM tables for alloy generation
       do 8 i=1,ntypes
c         Restart with initial EAM tables
          nrhoin(i)=snrhoin(i)
          drhoin(i)=sdrhoin(i)
          nrin(i)=snrin(i)
          drin(i)=sdrin(i)
          do 3 j=1,snrhoin(i)
             frhoin(j,i)=sfrhoin(j,i)
3         continue
          do 9 j=1,snrin(i)
             zrin(j,i)=szrin(j,i)
             rhorin(j,i)=srhorin(j,i)
9         continue
8      continue
      endif

      pairnum=ntypes*(ntypes-1)/2
      if (ntypes .le. 1) pairnum=1       

      ii=1
      if (eamalloy .gt. 0) then
c      if (nparams.eq.7) then
       if (ntypes.eq.2) then
        if (prmflag.eq.1) then 
         sx(1)=1.0
         gx(1)=params(ii)
         sx(2)=params(ii+1)
         gx(2)=params(ii+2)
         ii=ii+3
        elseif (prmflag.eq.2) then
         sx(1)=params(ii)
         gx(1)=params(ii+1)
         sx(2)=1.0
         gx(2)=params(ii+2)
         ii=ii+3
        else
         sx(1)=params(ii)
         gx(1)=params(ii+1)
         sx(2)=params(ii+2)
         gx(2)=params(ii+3)
         ii=ii+4
        endif
c      elseif (nparams.eq.8) then
c        sx(1)=params(ii)
c        gx(1)=params(ii+1)
c        sx(2)=params(ii+2)
c        gx(2)=params(ii+3)
c        ii=ii+4
c      elseif (nparams.eq.17) then
       elseif (ntypes.eq.3) then
        if (prmflag.eq.1) then
         sx(1)=1.0
         gx(1)=params(ii)
         sx(2)=params(ii+1)
         gx(2)=params(ii+2)
         sx(3)=params(ii+3)
         gx(3)=params(ii+4)
         ii=ii+5
        elseif (prmflag.eq.2) then
         sx(1)=params(ii)
         gx(1)=params(ii+1)
         sx(2)=1.0
         gx(2)=params(ii+2)
         sx(3)=params(ii+3)
         gx(3)=params(ii+4)
         ii=ii+5
        elseif (prmflag.eq.3) then
         sx(1)=params(ii)
         gx(1)=params(ii+1)
         sx(2)=params(ii+2)
         gx(2)=params(ii+3)
         sx(3)=1.0
         gx(3)=params(ii+4)
         ii=ii+5
        else
         sx(1)=params(ii)
         gx(1)=params(ii+1)
         sx(2)=params(ii+2)
         gx(2)=params(ii+3)
         sx(3)=params(ii+4)
         gx(3)=params(ii+5)
         ii=ii+6
        endif
c     elseif (nparams.eq.18) then
c       sx(1)=params(ii)
c       gx(1)=params(ii+1)
c       sx(2)=params(ii+2)
c       gx(2)=params(ii+3)
c       sx(3)=params(ii+4)
c       gx(3)=params(ii+5)
c       ii=ii+6
       endif
      endif
c     do 6 i=1,ntypes
c       sx(i)=params(ii)
c       gx(i)=params(ii+1)
c       ii=ii+2
c   6 continue
       do 7 i=1,pairnum
        if (phitype .eq. 1 .or. phitype .eq. 7) then
         if ( i .eq. 3 ) then
          dm3=params(ii)
          rm3=params(ii+1)
          alpha3=params(ii+2)
          rcutguess3=params(ii+3)
          ii=ii+4
         elseif ( i .eq. 2 ) then
          dm2=params(ii)
          rm2=params(ii+1)
          alpha2=params(ii+2)
          rcutguess2=params(ii+3)
          ii=ii+4
         else
          dm1=params(ii)
          rm1=params(ii+1)
          alpha1=params(ii+2)
          rcutguess1=params(ii+3)
          ii=ii+4
         endif
        elseif (phitype .eq. 2) then
         if ( i .eq. 3 ) then
          dm3=params(ii)
          alpha3=params(ii+1)
          rcutguess3=params(ii+2)
          ii=ii+3
         elseif ( i .eq. 2 ) then
          dm2=params(ii)
          alpha2=params(ii+1)
          rcutguess2=params(ii+2)
          ii=ii+3
         else
          dm1=params(ii)
          alpha1=params(ii+1)
          rcutguess1=params(ii+2)
          ii=ii+3
         endif
        elseif (phitype .eq. 8 .and.
     .          i .eq. 1) then
          dm1=params(ii)
          rm1=params(ii+1)
          alpha1=params(ii+2)
          rcutguess1=params(ii+3)
          dm2=params(ii+4)
          rm2=params(ii+5)
          alpha2=params(ii+6)
          ii=ii+7
        endif
7      continue
      if (eamalloy .lt. 1) then
       if (rhotype .eq. 1) then
         beta1=params(ii)
         if (fixcutoff.lt.1) then
           rcutguess2=params(ii+1)
           ii=ii+1
         endif
         ii=ii+1
       elseif (rhotype .eq. 2) then
         beta1=params(ii)
         beta2=params(ii+1)
         ii=ii+2
       elseif (rhotype .eq. 3) then
         beta1=params(ii)
         ii=ii+1
       elseif (rhotype .eq. 4) then
         beta1=params(ii)
         beta2=params(ii+1)
         ii=ii+2
       elseif (rhotype .eq. 5 .or.
     .         rhotype .eq. 12) then
         beta1=params(ii)
         beta2=params(ii+1)
         beta3=params(ii+2)
         ii=ii+3
       elseif (rhotype .eq. 6 .or.
     .         rhotype .eq. 7 .or.
     .         rhotype .eq. 8 .or.
     .         rhotype .eq. 9 .or.
     .         rhotype .eq. 13 .or.
     .         rhotype .eq. 14 .or.
     .         rhotype .eq. 15) then
         beta1=params(ii)
         beta2=params(ii+1)
         beta3=params(ii+2)
         dm3=params(ii+3)
         ii=ii+4
       elseif (rhotype .eq. 10) then
         beta1=params(ii)
         beta2=params(ii+1)
         beta3=params(ii+2)
         dm2=params(ii+3)
         dm3=params(ii+4)
         ii=ii+5
       elseif (rhotype .eq. 11) then
         beta1=params(ii)
         beta2=params(ii+1)
         beta3=params(ii+2)
         rm3=params(ii+3)
         dm2=params(ii+4)
         dm3=params(ii+5)
         ii=ii+6
       endif
      endif

c      write(*,*)'PROC(',procnum,') part:',prmid,
c    .           ' PARAMS: ',(params(i),i=1,nparams)
c     if(phitype) then
      if (prntflag.eq.1.or.logflag.eq.1) then
       if (phitype .eq. 1 .or. phitype .eq. 7 .or.
     .     phitype .eq. 8) then
        if(logflag.eq.1) then
        write(11,*)'(Voter-Chen) Type EAM'
        write(11,*)'Using Morse Potential for alloy ', 
     .    ' Phi(Rij) pair interaction.'
        else
        write(*,*)'(Voter-Chen) Type EAM'
        write(*,*)'Using Morse Potential for alloy ', 
     .    ' Phi(Rij) pair interaction.'
        endif
       elseif (phitype .eq. 2) then
        if(logflag.eq.1) then
        write(11,*)'(S. Chen Comp. Mate. Sci. 2004) Type EAM'
        write(11,*)'Using Yukawa Type Potential for',
     .    ' Phi(Rij) pair interaction.'
        else
        write(*,*)'(S. Chen Comp. Mate. Sci. 2004) Type EAM'
        write(*,*)'Using Yukawa Type Potential for',
     .    ' Phi(Rij) pair interaction.'
        endif
       endif
       do 4 i=1,pairnum
        if (i.eq.3) then
         if(logflag.eq.1) then
         write(11,*)'DM',i,':',dm3
         else
         write(*,*)'DM',i,':',dm3
         endif
         if (phitype .eq. 1 .or. phitype .eq. 7) then
          if(logflag.eq.1) then
          write(11,*)'RM',i,':',rm3
          else
          write(*,*)'RM',i,':',rm3
          endif
         endif
         if(logflag.eq.1) then
         write(11,*)'Alpha',i,':',alpha3
         write(11,*)'Rcut Guess',i,':',rcutguess3
         else
         write(*,*)'Alpha',i,':',alpha3
         write(*,*)'Rcut Guess',i,':',rcutguess3
         endif
        else if (i.eq.2) then
         if(logflag.eq.1) then
         write(11,*)'DM',i,':',dm2
         else
         write(*,*)'DM',i,':',dm2
         endif
         if (phitype .eq. 1 .or. phitype .eq. 7) then
          if(logflag.eq.1) then
          write(11,*)'RM',i,':',rm2
          else
          write(*,*)'RM',i,':',rm2
          endif
         endif
         if(logflag.eq.1) then
         write(11,*)'Alpha',i,':',alpha2
         write(11,*)'Rcut Guess',i,':',rcutguess2
         else
         write(*,*)'Alpha',i,':',alpha2
         write(*,*)'Rcut Guess',i,':',rcutguess2
         endif
        else
         if(logflag.eq.1) then
         write(11,*)'DM',i,':',dm1
         else
         write(*,*)'DM',i,':',dm1
         endif
         if (phitype .eq. 1 .or. phitype .eq. 7 .or.
     .       phitype .eq. 8   ) then
          if(logflag.eq.1) then
          write(11,*)'RM',i,':',rm1
          else
          write(*,*)'RM',i,':',rm1
          endif
         endif
         if(logflag.eq.1) then
         write(11,*)'Alpha',i,':',alpha1
         write(11,*)'Rcut Guess',i,':',rcutguess1
         else
         write(*,*)'Alpha',i,':',alpha1
         write(*,*)'Rcut Guess',i,':',rcutguess1
         endif
         if (phitype .eq. 8) then
           if(logflag.eq.1) then
            write(11,*)'DM2',i,':',dm2
            write(11,*)'RM2',i,':',rm2
            write(11,*)'Aplha2',i,':',alpha2
           else
            write(*,*)'DM2',i,':',dm2
            write(*,*)'RM2',i,':',rm2
            write(*,*)'Aplha2',i,':',alpha2
           endif
         endif
        endif
4      continue
       if(eamalloy .lt. 1) then
        if (rhotype .eq. 1) then
         if(logflag.eq.1) then
        write(11,*)'Using Hydrogenic 4s orbital for', 
     .    ' Rho(Rij) density function.'
        write(11,*)'Beta: ',beta1
         if(fixcutoff.lt.1) then
         write(11,*)'Rcut Guess',i,':',rcutguess2
         endif
         else
        write(*,*)'Using Hydrogenic 4s orbital for', 
     .    ' Rho(Rij) density function.'
        write(*,*)'Beta: ',beta1
         if(fixcutoff.lt.1) then
         write(11,*)'Rcut Guess',i,':',rcutguess2
         endif
         endif
        else if (rhotype .eq. 4) then
         if(logflag.eq.1) then
        write(11,*)'Using Hydrogenic 4s orbital for', 
     .    ' Rho(Rij) density function.'
        write(11,*)'Beta1: ',beta1
        write(11,*)'Beta2: ',beta2
         else
        write(*,*)'Using Hydrogenic 4s orbital for', 
     .    ' Rho(Rij) density function.'
        write(*,*)'Beta1: ',beta1
        write(*,*)'Beta2: ',beta2
         endif
        else if (rhotype .eq. 5 .or. rhotype .eq. 12) then
         if(logflag.eq.1) then
        write(11,*)'Using Hydrogenic 4s orbital for', 
     .    ' Rho(Rij) density function.'
        write(11,*)'Beta1: ',beta1
        write(11,*)'Beta2: ',beta2
        write(11,*)'Beta3: ',beta3
         else
        write(*,*)'Using Hydrogenic 4s orbital for', 
     .    ' Rho(Rij) density function.'
        write(*,*)'Beta1: ',beta1
        write(*,*)'Beta2: ',beta2
        write(*,*)'Beta3: ',beta3
         endif
        else if (rhotype .eq. 6 .or.
     .           rhotype .eq. 7 .or.
     .           rhotype .eq. 8 .or.
     .           rhotype .eq. 9 .or.
     .           rhotype .eq. 10.or.
     .           rhotype .eq. 11.or.
     .           rhotype .eq. 13.or.
     .           rhotype .eq. 14.or.
     .           rhotype .eq. 15) then
         if(logflag.eq.1) then
        write(11,*)'Using Hydrogenic 4s orbital for', 
     .    ' Rho(Rij) density function.'
        write(11,*)'Beta1: ',beta1
        write(11,*)'Beta2: ',beta2
        write(11,*)'Beta3: ',beta3
        write(11,*)'Beta4: ',dm3
         else
        write(*,*)'Using Hydrogenic 4s orbital for', 
     .    ' Rho(Rij) density function.'
        write(*,*)'Beta1: ',beta1
        write(*,*)'Beta2: ',beta2
        write(*,*)'Beta3: ',beta3
        write(*,*)'Beta4: ',dm3
         endif
        elseif (rhotype .eq. 2) then
         if(logflag.eq.1) then
        write(11,*)'Using Metal Valance e^- and Oscillation for', 
     .    ' Rho(Rij) density function.'
        write(11,*)'Beta: ',beta1
        write(11,*)'nfactor: ',beta2
         else
        write(*,*)'Using Metal Valance e^- and Oscillation for', 
     .    ' Rho(Rij) density function.'
        write(*,*)'Beta: ',beta1
         endif
        elseif (rhotype .eq. 3) then
         if(logflag.eq.1) then
        write(11,*)'Using Thomas-Fermi Screening Type',
     .    ' Rho(Rij) density function.'
        write(11,*)'Beta: ',beta1
         else
        write(*,*)'Using Thomas-Fermi Screening Type',
     .    ' Rho(Rij) density function.'
        write(*,*)'Beta: ',beta1
         endif
        endif
       endif
      endif

      if (prntflag.eq.1) then
       if (eamalloy .gt. 0) then
        do 5 i=1,ntypes
c     if(i.lt.2)then
c             elmtyp='A'
c     else if(i.lt.3) then
c             elmtyp='B'
c     else
c             elmtyp='C'
c     endif
         write(*,*)'Rcut for ',elmname(i),': ',rcut(i)
         write(*,*)'s',elmname(i),': ',sx(i)
         write(*,*)'g',elmname(i),': ',gx(i)
   5    continue
       endif
      endif
c     write(*,*)fileformat,ntypes,amass(1),ielement(1),lata0
c     write(*,*)eamfileout,initf,rnum,rhonum,rhoend
       stacks(1)=nxchain
       stacks(2)=nyatoms
       stacks(3)=nzlayer
       natoms=nxchain*nyatoms*nzlayer
       vacatom=(nxchain*nyatoms)*((nzlayer/2)-1)+(nxchain/2)*nyatoms
      if (prntflag.eq.1) then
       write(*,7373)nxchain,nyatoms,nzlayer,natoms
      endif
      if(logflag.eq.1) then
       write(11,7373)nxchain,nyatoms,nzlayer,natoms
      endif
7373   format('Using ',i3,' x',i3,' x',i3,' structure,   natoms:',i6)

c
c     if pairnum=1    pair  is  12
c     if pairnum=3    pairs are 12, 13, 23
c
         anrin1=rnum
         adrin1=rcutguess1/(rnum-1)
        if (fixcutoff.lt.1) then
          anrin2=rnum
          adrin2=rcutguess2/(rnum-1)
        endif
        if (eamalloy .lt. 1) then
         nrin(1)=anrin1
         nrhoin(1)=rhonum
         drhoin(1)=rhoend/(rhonum-1)
         drin(1)=adrin1 
c        if (fixcutoff.lt.1) then
c         drin(1) = amax1(drin(1),adrin2)
c        endif
        endif
        if( eamalloy .gt. 0 .and. pairnum .eq. 3 ) then
          anrin2=rnum
          adrin2=rcutguess2/(rnum-1)
          anrin3=rnum
          adrin3=rcutguess3/(rnum-1)
        endif

c     ------------------------------------------------------------
      if (eamalloy .lt. 1) then
c      Find Pure Rho Function
       if(rhotype .eq. 1) then
        if(fixcutoff.lt.1) then
        rhorcut=(rcutguess2**6) * ( exp(-beta1*rcutguess2) +
     .                     (512)*exp(-2.0*beta1*rcutguess2) )
        rhodrrcut=6*(rcutguess2**5) * ( exp(-beta1*rcutguess2) +
     .                   (512)*exp(-2.0*beta1*rcutguess2) ) +
     .       (rcutguess2**6) * ( -beta1*exp(-beta1*rcutguess2) -
     .              (1024)*beta1*exp(-2.0*beta1*rcutguess2) )
        else
        rhorcut=(rcutguess1**6) * ( exp(-beta1*rcutguess1) +
     .                     (512)*exp(-2.0*beta1*rcutguess1) )
        rhodrrcut=6*(rcutguess1**5) * ( exp(-beta1*rcutguess1) +
     .                   (512)*exp(-2.0*beta1*rcutguess1) ) +
     .       (rcutguess1**6) * ( -beta1*exp(-beta1*rcutguess1) -
     .              (1024)*beta1*exp(-2.0*beta1*rcutguess1) )
        endif
       else if (rhotype.eq.4) then
        rhorcut=(rcutguess1**6) * ( exp(-beta1*rcutguess1) +
     .                     (512)*exp(-2.0*beta2*rcutguess1) )
        rhodrrcut=6*(rcutguess1**5) * ( exp(-beta1*rcutguess1) +
     .                   (512)*exp(-2.0*beta2*rcutguess1) ) +
     .       (rcutguess1**6) * ( -beta1*exp(-beta1*rcutguess1) -
     .              (1024)*beta2*exp(-2.0*beta2*rcutguess1) )
       else if (rhotype.eq.5) then
        rhorcut=(rcutguess1**6) * ( exp(-beta1*rcutguess1) +
     .                     (512)*exp(-2.0*beta2*rcutguess1) ) * 
     .                          exp(-((rcutguess1-1)**2)/beta3)
        rhodrrcut=(6*(rcutguess1**5) * ( exp(-beta1*rcutguess1) +
     .                       (512)*exp(-2.0*beta2*rcutguess1) ) +
     .        (rcutguess1**6) * ( -beta1*exp(-beta1*rcutguess1) -
     .                (1024)*beta2*exp(-2.0*beta2*rcutguess1) ) ) *
     .                            exp(-((rcutguess1-1)**2)/beta3) -
     .              ((rcutguess1**6) * ( exp(-beta1*rcutguess1) +
     .                        (512)*exp(-2.0*beta2*rcutguess1) )) * 
     .   2*((rcutguess1-1.0)/beta3)*exp(-((rcutguess1-1)**2)/beta3)
       else if (rhotype.eq.6) then
        rhorcut=(rcutguess1**6) * ( exp(-beta1*rcutguess1) +
     .                     (512)*exp(-2.0*beta2*rcutguess1) ) * 
     .                      exp(-((rcutguess1-1.0)**2)/beta3) *
     .                          exp(-((rcutguess1-1.0)**4)/dm3)
        rhodrrcut=(6*(rcutguess1**5) * ( exp(-beta1*rcutguess1) +
     .                       (512)*exp(-2.0*beta2*rcutguess1) ) +
     .        (rcutguess1**6) * ( -beta1*exp(-beta1*rcutguess1) -
     .                (1024)*beta2*exp(-2.0*beta2*rcutguess1) ) ) *
     .                         (exp(-((rcutguess1-1.0)**2)/beta3) *
     .                          exp(-((rcutguess1-1.0)**4)/dm3) ) +
     .              ((rcutguess1**6) * ( exp(-beta1*rcutguess1) +
     .                       (512)*exp(-2.0*beta2*rcutguess1) )) * 
     .(-2*((rcutguess1-1.0)/beta3)*exp(-((rcutguess1-1.0)**2)/beta3) *
     .                               exp(-((rcutguess1-1.0)**4)/dm3) -
     .      4*((rcutguess1-1.0)/dm3)*exp(-((rcutguess1-1.0)**4)/dm3) *
     .                            exp(-((rcutguess1-1.0)**2)/beta3)  )
       else if (rhotype.eq.7) then
        rhorcut=(rcutguess1**6) * ( exp(-beta1*rcutguess1) +
     .                     (512)*exp(-2.0*beta2*rcutguess1) ) * 
     .                      exp(((rcutguess1-1.0)**2)/beta3) *
     .                          exp(-((rcutguess1-1.0)**4)/dm3)
        rhodrrcut=(6*(rcutguess1**5) * ( exp(-beta1*rcutguess1) +
     .                       (512)*exp(-2.0*beta2*rcutguess1) ) +
     .        (rcutguess1**6) * ( -beta1*exp(-beta1*rcutguess1) -
     .                (1024)*beta2*exp(-2.0*beta2*rcutguess1) ) ) *
     .                         (exp(((rcutguess1-1.0)**2)/beta3) *
     .                          exp(-((rcutguess1-1.0)**4)/dm3) ) +
     .              ((rcutguess1**6) * ( exp(-beta1*rcutguess1) +
     .                       (512)*exp(-2.0*beta2*rcutguess1) )) * 
     .(2*((rcutguess1-1.0)/beta3)*exp(((rcutguess1-1.0)**2)/beta3) *
     .                               exp(-((rcutguess1-1.0)**4)/dm3) -
     .      4*((rcutguess1-1.0)/dm3)*exp(-((rcutguess1-1.0)**4)/dm3) *
     .                            exp(((rcutguess1-1.0)**2)/beta3)  )
       else if (rhotype.eq.8) then
        rhorcut=((rcutguess1-beta1)**5) * ( beta2*
     .               exp(-beta3*(rcutguess1-dm3))) 
        rhodrrcut= 5*((rcutguess1-beta1)**5) *
     .  (beta2*exp(-beta3*(rcutguess1-dm3))) + 
     .  ((rcutguess1-beta1)**5) * (-beta2*beta3 *
     .               exp(-beta3*(rcutguess1-dm3))) 
       else if (rhotype.eq.9) then
        rhorcut=beta1*(rcutguess1**3) * 
     .           ( exp(-beta2*rcutguess1-dm3)) *
     .             exp(-beta3*(rcutguess1)**2)
        rhodrrcut= 3*beta1*(rcutguess1**2) *
     .           ( exp(-beta2*rcutguess1-dm3) ) * 
     .             exp(-beta3*(rcutguess1)**2) +
     .             beta1*(rcutguess1**3)*(-beta2 *
     .              exp(-beta2*(rcutguess1-dm3)) * 
     .               exp(-beta3*(rcutguess1)**2) - 
     .               exp(-beta2*(rcutguess1-dm3)) * 
     .        2*beta3*exp(-beta3*(rcutguess1)**2) )
       else if (rhotype.eq.10) then
        rhorcut=beta1*(rcutguess1**8) * 
     .           (exp(-beta2*(rcutguess1-dm2)) +
     .            exp(-beta3*(rcutguess1-dm3)**2))
        rhodrrcut= 8*beta1*(rcutguess1**7) *
     .         ( exp(-beta2*(rcutguess1-dm2))  + 
     .           exp(-beta3*(rcutguess1-dm3)**2) ) +
     .             beta1*(rcutguess1**3)*(-beta2 *
     .              exp(-beta2*(rcutguess1-dm2)) -
     .        beta3*exp(-beta3*(rcutguess1-dm3)**2) )
       else if (rhotype.eq.11) then
        rhorcut=(rcutguess1**6) * ( exp(-beta1*rcutguess1) +
     .                     (512)*exp(-2.0*beta2*rcutguess1) ) * 
     .                       exp(beta3*((rcutguess1-rm3)**2)) *
     .                          exp(dm2*((rcutguess1-dm3)**4))
        rhodrrcut=(6*(rcutguess1**5) * ( exp(-beta1*rcutguess1) +
     .                       (512)*exp(-2.0*beta2*rcutguess1) ) +
     .        (rcutguess1**6) * ( -beta1*exp(-beta1*rcutguess1) -
     .                (1024)*beta2*exp(-2.0*beta2*rcutguess1) ) ) *
     .                         (exp(beta3*((rcutguess1-rm3)**2)) *
     .                          exp(dm2*((rcutguess1-dm3)**4)) ) +
     .              ((rcutguess1**6) * ( exp(-beta1*rcutguess1) +
     .                       (512)*exp(-2.0*beta2*rcutguess1) )) * 
     . ( 2*beta3*(rcutguess1-rm3)*exp(beta3*((rcutguess1-rm3)**2)) *
     .                              exp(dm2*((rcutguess1-dm3)**4)) +
     .  4*dm2*((rcutguess1-dm3)**3)*exp(dm2*((rcutguess1-dm3)**4)) *
     .                            exp(beta3*((rcutguess1-dm1)**2))  )
       else if (rhotype.eq.12) then
        rhorcut=((rcutguess1)**6) * 
     .                ( exp(-beta1*(rcutguess1-beta3)) +
     .            (512)*exp(-2.0*beta2*(rcutguess1-beta3)) )
        rhodrrcut=6*(rcutguess1**5) * 
     .                ( exp(-beta1*(rcutguess1-beta3)) +
     .            (512)*exp(-2.0*beta2*(rcutguess1-beta3)) ) +
     .   (rcutguess1**6) * ( -beta1*exp(-beta1*(rcutguess1-beta3)) -
     .            (1024)*beta2*exp(-2.0*beta2*(rcutguess1-beta3)) )
       else if (rhotype.eq.13) then
        rhorcut=((rcutguess1)**6) * 
     .        ( exp(-beta1*rcutguess1) +
     .        (512)*exp(-2.0*beta1*rcutguess1) ) +
     .        beta2*(1.0/beta3)*(1.0/dm3)*
     .        exp(-0.5*((rcutguess1-dm3)*beta3)**2)+
     .       (beta2/10.0)*exp(-0.5*((rcutguess1-(dm3+0.5))*beta3)**2)
        rhodrrcut=6*(rcutguess1**5) * ( exp(-beta1*rcutguess1) +
     .                      (512)*exp(-2.0*beta1*rcutguess1) ) +
     .   (rcutguess1**6) * ( -beta1*exp(-beta1*rcutguess1) -
     .            (1024)*beta1*exp(-2.0*beta1*rcutguess1) ) -
     .beta2*(1.0/beta3)*(1.0/dm3)*0.5*((2.0*beta3)*(rcutguess1-dm3))*
     .                        exp(-0.5*((rcutguess1-dm3)*beta3)**2) -
     .    (beta2/10.0)*0.5*((2.0*beta3)*(rcutguess1-(dm3+0.5)))*
     .                    exp(-0.5*((rcutguess1-(dm3+0.5))*beta3)**2) 
       else if (rhotype.eq.14) then
        rhorcut=dtanh(20.*rcutguess1*rcutguess1)*(((rcutguess1)**6) * 
     .        ( exp(-beta1*rcutguess1) +
     .        (512)*exp(-2.0*beta1*rcutguess1) ) +
     .        beta2*(1.0/beta3)*(1.0/dm3)*
     .        exp(-0.5*((rcutguess1-dm3)*beta3)**2)+
     .       (beta2/10.0)*exp(-0.5*((rcutguess1-(dm3+0.5))*beta3)**2))
        rhodrrcut=(1.-(dtanh(20.*rcutguess1*rcutguess1)*
     .                 dtanh(20.*rcutguess1*rcutguess1)) ) *
     .   (6*(rcutguess1**5) * ( exp(-beta1*rcutguess1) +
     .                      (512)*exp(-2.0*beta1*rcutguess1) ) +
     .   (rcutguess1**6) * ( -beta1*exp(-beta1*rcutguess1) -
     .            (1024)*beta1*exp(-2.0*beta1*rcutguess1) ) -
     .beta2*(1.0/beta3)*(1.0/dm3)*0.5*((2.0*beta3)*(rcutguess1-dm3))*
     .                        exp(-0.5*((rcutguess1-dm3)*beta3)**2) -
     .    (beta2/10.0)*0.5*((2.0*beta3)*(rcutguess1-(dm3+0.5)))*
     .                    exp(-0.5*((rcutguess1-(dm3+0.5))*beta3)**2)) +
     .  dtanh(20.*(rcutguess1*rcutguess1))*(6*(rcutguess1**5) * 
     .          ( exp(-beta1*rcutguess1) +
     .                      (512)*exp(-2.0*beta1*rcutguess1) ) +
     .   (rcutguess1**6) * ( -beta1*exp(-beta1*rcutguess1) -
     .            (1024)*beta1*exp(-2.0*beta1*rcutguess1) ) -
     .beta2*(1.0/beta3)*(1.0/dm3)*0.5*((2.0*beta3)*(rcutguess1-dm3))*
     .                        exp(-0.5*((rcutguess1-dm3)*beta3)**2) -
     .    (beta2/10.0)*0.5*((2.0*beta3)*(rcutguess1-(dm3+0.5)))*
     .                    exp(-0.5*((rcutguess1-(dm3+0.5))*beta3)**2)) 
       else if (rhotype.eq.15) then
        rhorcut=dtanh(20.*rcutguess1*rcutguess1)*(((rcutguess1)**6) * 
     .        ( exp(-beta1*rcutguess1) +
     .        (512)*exp(-2.0*beta1*rcutguess1) ) +
     .        beta2*(1.0/beta3)*(1.0/dm3)*
     .        exp(-0.5*((rcutguess1-dm3)*beta3)**2)-
     .       (beta2/10.0)*exp(-0.5*((rcutguess1-(dm3+0.5))*beta3)**2) +
     .                          0.4*exp(-0.5*10.0*(rcutguess1-1.0)**2) )
        rhodrrcut=(1.-(dtanh(20.*rcutguess1*rcutguess1)*
     .                 dtanh(20.*rcutguess1*rcutguess1)) ) *
     .   (6*(rcutguess1**5) * ( exp(-beta1*rcutguess1) +
     .                      (512)*exp(-2.0*beta1*rcutguess1) ) +
     .   (rcutguess1**6) * ( -beta1*exp(-beta1*rcutguess1) -
     .            (1024)*beta1*exp(-2.0*beta1*rcutguess1) ) -
     .beta2*(1.0/beta3)*(1.0/dm3)*0.5*((2.0*beta3)*(rcutguess1-dm3))*
     .                        exp(-0.5*((rcutguess1-dm3)*beta3)**2) +
     .    (beta2/10.0)*0.5*((2.0*beta3)*(rcutguess1-(dm3+0.5)))*
     .                    exp(-0.5*((rcutguess1-(dm3+0.5))*beta3)**2)+
     .                         0.4*exp(-0.5*10.0*(rcutguess1-1.0)**2))+
     .  dtanh(20.*(rcutguess1*rcutguess1))*(6*(rcutguess1**5) * 
     .          ( exp(-beta1*rcutguess1) +
     .                      (512)*exp(-2.0*beta1*rcutguess1) ) +
     .   (rcutguess1**6) * ( -beta1*exp(-beta1*rcutguess1) -
     .            (1024)*beta1*exp(-2.0*beta1*rcutguess1) ) -
     .beta2*(1.0/beta3)*(1.0/dm3)*0.5*((2.0*beta3)*(rcutguess1-dm3))*
     .                        exp(-0.5*((rcutguess1-dm3)*beta3)**2) +
     .    (beta2/10.0)*0.5*((2.0*beta3)*(rcutguess1-(dm3+0.5)))*
     .                    exp(-0.5*((rcutguess1-(dm3+0.5))*beta3)**2)+
     .                   0.4*0.5*10.0*2.0*(rcutguess1-1.0)*
     .                                exp(-0.5*10*(rcutguess1-1.0)**2)) 
       else if (rhotype.eq.2) then
        rhorcut= rm1*(rcutguess1**nfactor) *
     .    (1.0 + sin(rcutguess1*pi/(2.0*lata0)) ) * 
     .    exp(-beta1*rcutguess1)
        rhodrrcut=rm1*nfactor*(rcutguess1**(nfactor-1)) *
     .    (1.0 + sin(rcutguess1*pi/(2.0*lata0)) ) * 
     .    exp(-beta1*rcutguess1)+rm1*(rcutguess1**nfactor) * 
     .    (-beta1*exp(-beta1*rcutguess1) +
     .    pi/(2.0*lata0)*cos(rcutguess1*pi/(2.0*lata0)) *
     .    exp(-beta1*rcutguess1)-beta1*exp(-beta1*rcutguess1)*
     .    sin(rcutguess1*pi/(2.0*pi*lata0)))
       else if (rhotype.eq.3) then
        rhorcut= exp(-beta1*rcutguess1/lata0) / 
     .              (rcutguess1/lata0)
        rhodrrcut= (-(beta1/lata0)*exp(-beta1*rcutguess1/lata0) * 
     .  (rcutguess1/lata0) + (exp(-beta1*rcutguess1/lata0)*(1./lata0)))
     .  / (rcutguess1/lata0)**2
       endif
       do 10 j=1,nrin(1) 
c      if(fixcutoff.lt.1) then      
c       rhorr=(j-1)*adrin2
c      else
        rhorr=(j-1)*drin(1)
c      endif
        if (rhotype .eq. 1) then
         if (fixcutoff.lt.1) then
          if (rhorr.lt.rcutguess2) then
          rhorin(j,1)=(rhorr**6)*( exp(-beta1*rhorr) +
     .                      (512)*exp(-2.0*beta1*rhorr) )
          else
          rhorin(j,1)=0.d0
          endif
         else
          rhorin(j,1)=(rhorr**6)*( exp(-beta1*rhorr) +
     .                      (512)*exp(-2.0*beta1*rhorr) )
         endif
        else if (rhotype .eq. 4) then
          rhorin(j,1)=(rhorr**6)*( exp(-beta1*rhorr) +
     .                      (512)*exp(-2.0*beta2*rhorr) )
        else if (rhotype .eq. 5) then
          rhorin(j,1)=(rhorr**6)*( exp(-beta1*rhorr) +
     .                 (512)*exp(-2.0*beta2*rhorr) ) *
     .                    exp(-((rhorr-1.0)**2)/beta3) 
        else if (rhotype .eq. 6) then
          rhorin(j,1)=(rhorr**6)*( exp(-beta1*rhorr) +
     .                 (512)*exp(-2.0*beta2*rhorr) ) *
     .                  exp(-((rhorr-1.0)**2)/beta3) * 
     .                    exp(-((rhorr-1.0)**4)/dm3) 
        else if (rhotype .eq. 7) then
          rhorin(j,1)=(rhorr**6)*( exp(-beta1*rhorr) +
     .                 (512)*exp(-2.0*beta2*rhorr) ) *
     .                  exp(((rhorr-1.0)**2)/beta3) * 
     .                    exp(-((rhorr-1.0)**4)/dm3) 
        else if (rhotype .eq. 8) then
          rhorin(j,1)=((rhorr-beta1)**5) * ( beta2*
     .                  exp(-beta3*(rhorr-dm3))) 
        else if (rhotype .eq. 9) then
          rhorin(j,1)=beta1*(rhorr**3) *
     .           ( exp(-beta2*rhorr-dm3) *
     .             exp(-beta3*(rhorr)**2) )
        else if (rhotype .eq. 10) then
          rhorin(j,1)=beta1*(rhorr**8) *
     .           ( exp(-beta2*(rhorr-dm2)) +
     .             exp(-beta3*(rhorr-dm3)**2) )
        else if (rhotype .eq. 11) then
          rhorin(j,1)=(rhorr**6)*( exp(-beta1*rhorr) +
     .                 (512)*exp(-2.0*beta2*rhorr) ) *
     .                  exp(beta3*((rhorr-dm1)**2)) * 
     .                    exp(dm2*((rhorr-dm3)**4)) 
        else if (rhotype .eq. 12) then
          rhorin(j,1)=((rhorr)**6) *
     .           ( exp(-beta1*(rhorr-beta3)) +(512)*
     .             exp(-2.0*beta2*(rhorr-beta3)) )
        else if (rhotype .eq. 13) then
          rhorin(j,1)=((rhorr)**6) * ( exp(-beta1*rhorr) +(512)*
     .    exp(-2.0*beta1*rhorr) ) + beta2*(1.0/beta3)*(1.0/dm3)*
     .    exp(-0.5*((rhorr-dm3)*beta3)**2) +
     .    (beta2/10.0)*exp(-0.5*((rhorr-(dm3+0.5))*beta3)**2)
        else if (rhotype .eq. 14) then
          rhorin(j,1)=dtanh(20.*(rhorr**2))*(((rhorr)**6) *
     .    ( exp(-beta1*rhorr) +(512)*
     .    exp(-2.0*beta1*rhorr) ) + beta2*(1.0/beta3)*(1.0/dm3)*
     .    exp(-0.5*((rhorr-dm3)*beta3)**2) +
     .    (beta2/10.0)*exp(-0.5*((rhorr-(dm3+0.5))*beta3)**2))
        else if (rhotype .eq. 15) then
          rhorin(j,1)=dtanh(20.*(rhorr**2))*(((rhorr)**6) *
     .    ( exp(-beta1*rhorr) +(512)*
     .    exp(-2.0*beta1*rhorr) ) + beta2*(1.0/beta3)*(1.0/dm3)*
     .    exp(-0.5*((rhorr-dm3)*beta3)**2) -
     .    (beta2/10.0)*(1.0/beta3)*(1.0/(dm3+0.5))*
     .    exp(-0.5*((rhorr-(dm3+0.5))*beta3)**2) +
     .           0.4*exp(-0.5*10.0*(rhorr-1.0)**2) )
        else if (rhotype .eq. 2) then
         rhorr=rhorr/0.5291
         rhorin(j,1)= rm1*(rhorr**nfactor) *
     .     (1.0 + sin(pi*rhorr/(2.0*lata0)) ) * 
     .     exp(-beta1*rhorr)
        elseif (rhotype .eq. 3) then
         rhorr=rhorr/lata0
         rhorin(j,1)= exp(-beta1*rhorr) / rhorr
        endif
        if(rhotype .eq. 1 .or. rhotype .ge. 3 ) then
         if (fixcutoff.lt.1) then
         rhorin(j,1) = rhorin(j,1) - rhorcut +
     .       (rcutguess2/20.0)*(1.0-(rhorr/rcutguess2)**20)*
     .        rhodrrcut
         else
         rhorin(j,1) = rhorin(j,1) - rhorcut +
     .       (rcutguess1/20.0)*(1.0-(rhorr/rcutguess1)**20)*
     .        rhodrrcut
         endif 
        endif
        if(rhorin(j,1).gt.1000000) then
           rhorin(j,1)=1000000
        elseif(rhorin(j,1).lt.-1000000) then
           rhorin(j,1)=-1000000
        endif
10    continue
       rhorin(nrin(1),1)=0.d0
      endif
c     ------------------------------------------------------------

c     ------------------------------------------------------------
c      Finding Alloy Phi from Morse pair potential
c     ------------------------------------------------------------
      if (phitype .eq. 1 .or. phitype .eq. 7) then
c      phi(r)=D_m(1-exp(ALPHA_m*(r-R_m)))^2-D_m
c      D_m, R_m, and ALPHA_m is depth and position of minimum
        zrcut=dm1*(1.0-exp(-alpha1*
     .               (rcutguess1-rm1)))**2-dm1
        zdrrcut=2.0*dm1*(alpha1*exp(-alpha1*
     .              (rcutguess1-rm1)))*(1.0-exp(-alpha1*
     .              (rcutguess1-rm1)))
      elseif (phitype .eq. 8) then
        zrcut=dm1*(1.0-exp(-alpha1*
     .               (rcutguess1-rm1)))**2-dm1 +
     .        dm2*(1.0-exp(-alpha2*
     .               (rcutguess1-rm2)))**2-dm2 
        zdrrcut=2.0*dm1*(alpha1*exp(-alpha1*
     .              (rcutguess1-rm1)))*(1.0-exp(-alpha1*
     .                                (rcutguess1-rm1))) +
     .          2.0*dm2*(alpha2*exp(-alpha2*
     .              (rcutguess1-rm2)))*(1.0-exp(-alpha2*
     .                                (rcutguess1-rm2))) 
      elseif (phitype .eq. 2) then
        zrcut=dm1*exp(-alpha1*rcutguess1)/rcutguess1
        zdrrcut=-dm1*exp(-alpha1*rcutguess1)/(rcutguess1**2) -
     .           dm1*alpha1*exp(-alpha1*rcutguess1)/rcutguess1
c     elseif(phitype.eq.3) then
c       zrcut=ak(1)*(rk(1)-rcutguess1/lata0)**3
c       zdrrcut=-3*ak(1)*(rk(1)-rcutguess1/lata0)**2
c       k=pasn
      endif

        if (eamalloy .gt. 0) then
c
c change zrin to new calculated z values
        do 9002 j=1,anrin1
         zrr=(j-1)*adrin1
         azrin1(j)=dm1*(1.0-exp(-alpha1*(zrr-rm1)))**2-dm1
c       Modify Z(R) to satisfy Z(Rcutoff)=0 by using
c        smoothing function with m=20
         azrin1(j) = azrin1(j) - zrcut
     .     + (rcutguess1/20.0)*(1.0-(zrr/rcutguess1)**20)*zdrrcut
         if( azrin1(j) .gt. 1000 ) then
            azrin1(j)=1000
         elseif( azrin1(j) .lt. -1000 ) then
            azrin1(j)=-1000
         endif
9002    continue
        azrin1(anrin1)=0.d0

        else
c change zrin to new calculated z values
        do 9001 j=1,nrin(1)
         zrr=(j-1)*drin(1)
         if(phitype.eq.1 .or. phitype .eq. 7) then
          if (fixcutoff.lt.1) then
           if (zrr.lt.rcutguess1) then
          zrin(j,1)=dm1*(1.0-exp(-alpha1*(zrr-rm1)))**2-dm1
           else
          zrin(j,1)=0.d0
           endif
          else
          zrin(j,1)=dm1*(1.0-exp(-alpha1*(zrr-rm1)))**2-dm1
          endif
         elseif(phitype.eq.8) then
          zrin(j,1)=dm1*(1.0-exp(-alpha1*(zrr-rm1)))**2-dm1 +
     .              dm2*(1.0-exp(-alpha2*(zrr-rm2)))**2-dm2 
         elseif(phitype.eq.2) then
          zrr=zrr/0.5291
c         dm1=dm1/0.62415
          zrin(j,1)=(dm1/zrr)*exp(-alpha1*zrr)
         endif
c       Modify Z(R) to satisfy Z(Rcutoff)=0 by using
c        smoothing function with m=20
         if(phitype .eq. 1 .or. phitype .ge. 3 .or. 
     .      phitype .eq. 7 .or. phitype .eq. 8) then
         zrin(j,1) = zrin(j,1) - zrcut
     .     + (rcutguess1/20.0)*(1.0-(zrr/rcutguess1)**20)*zdrrcut
         endif
         if(zrin(j,1).gt.10000) then
            zrin(j,1)=10000
         elseif(zrin(j,1).lt.-10000) then
            zrin(j,1)=-10000
         endif
9001    continue
        zrin(nrin(1),1)=0.d0

        endif

        if( eamalloy .gt. 0 .and. pairnum .eq. 3 ) then
c     ------------------------------------------------------------
c      Finding Alloy Phi from Morse pair potential
c     ------------------------------------------------------------
        if (phitype .eq. 1 .or. phitype .eq. 7) then
        zrcut=dm2*(1.0-exp(-alpha2*
     .               (rcutguess2-rm2)))**2-dm2
        zdrrcut=2.0*dm2*(alpha2*exp(-alpha2*
     .              (rcutguess2-rm2)))*(1.0-exp(-alpha2*
     .              (rcutguess2-rm2)))
        elseif (phitype .eq. 2) then
        zrcut=dm2*exp(-alpha2*rcutguess2)/rcutguess2
        zdrrcut=-dm2*exp(-alpha2*rcutguess2)/(rcutguess2**2) -
     .           dm2*alpha2*exp(-alpha2*rcutguess2)/rcutguess2
        endif
        
c
c change zrin to new calculated z values
        do 9003 j=1,anrin2
         zrr=(j-1)*adrin2
         azrin2(j)=dm2*(1.0-exp(-alpha2*(zrr-rm2)))**2-dm2
c       Modify Z(R) to satisfy Z(Rcutoff)=0 by using
c        smoothing function with m=20
         azrin2(j) = azrin2(j) - zrcut
     .     + (rcutguess2/20.0)*(1.0-(zrr/rcutguess2)**20)*zdrrcut
         if( azrin2(j) .gt. 1000 ) then
            azrin2(j)=1000
         elseif( azrin2(j) .lt. -1000 ) then
            azrin2(j)=-1000
         endif
9003    continue
        azrin2(anrin2)=0.d0

c     ------------------------------------------------------------
c      Finding Alloy Phi from Morse pair potential
c     ------------------------------------------------------------

        if (phitype .eq. 1 .or. phitype .eq. 7) then
        zrcut=dm3*(1.0-exp(-alpha3*
     .               (rcutguess3-rm3)))**2-dm3
        zdrrcut=2.0*dm3*(alpha3*exp(-alpha3*
     .              (rcutguess3-rm3)))*(1.0-exp(-alpha3*
     .              (rcutguess3-rm3)))
        elseif (phitype .eq. 2) then
        zrcut=dm3*exp(-alpha3*rcutguess3)/rcutguess3
        zdrrcut=-dm3*exp(-alpha3*rcutguess3)/(rcutguess3**2) -
     .           dm3*alpha3*exp(-alpha3*rcutguess3)/rcutguess3
        endif
c
c change zrin to new calculated z values
        do 9004 j=1,anrin3
         zrr=(j-1)*adrin3
         azrin3(j)=dm3*(1.0-exp(-alpha3*(zrr-rm3)))**2-dm3
c       Modify Z(R) to satisfy Z(Rcutoff)=0 by using
c        smoothing function with m=20
         azrin3(j) = azrin3(j) - zrcut
     .     + (rcutguess3/20.0)*(1.0-(zrr/rcutguess3)**20)*zdrrcut
         if( azrin3(j) .gt. 1000 ) then
            azrin3(j)=1000
         elseif( azrin3(j) .lt. -1000 ) then
            azrin3(j)=-1000
         endif
9004    continue
        azrin3(anrin3)=0.d0
        
        endif

c95    continue
c
c
c     ------------------------------------------------------
c      Applying invariant transformations for pure elements
c                 but not invariant for alloys
c     ------------------------------------------------------

      if (eamalloy .gt. 0) then
      do 96 i1=1,ntypes
       do 97 j=1,nrin(i1)
        zrin(j,i1)=zrin(j,i1)-2.0*gx(i1)*rhorin(j,i1)
97    continue
96    continue

      do 98 i1=1,ntypes
       do 99 j=1,nrhoin(i1)
        rhorr = (j-1)*drhoin(i1)
        frhoin(j,i1)=frhoin(j,i1)+gx(i1)*rhorr
99    continue
98    continue

      do 103 i1=1,ntypes
c     only relative sx can change alloy energy
c     we can take sx(1)=1.0
      drhoin(i1)=drhoin(i1)*sx(i1)
      do 104 j=1,nrin(i1)
        rhorin(j,i1) = sx(i1)*rhorin(j,i1)
104    continue
103    continue
c     -----------------Transformation end------------------
c     -----------------------------------------------------
      else
        rcut(1)=rcutguess1
        if (fixcutoff.lt.1) then
         rmax = rcut(1)
         rmaxi= rcutguess2
         rmax = amax1(rmax,rmaxi)
         rcut(1) = rmax
        endif
        blat(1)=lata0
        if(lattype .eq. 0) then
          lat(1)='FCC'
        else if(lattype .eq. 1) then
          lat(1)='BCC'
        else if(lattype .eq. 2) then
          lat(1)='HCP'
        else if(lattype .eq. 3) then
          lat(1)='DIA'
        else
          lat(1)='FCC'
        endif
        ielement(2)=ielement(1)
        amass(2)=amass(1)
        blat(2)=blat(1)
        lat(2)=lat(1)
        nrhoin(2)=nrhoin(1)
        drhoin(2)=drhoin(1)
        nrin(2)=nrin(1)
        drin(2)=drin(1)
        rcut(2)=rcut(1)
c        do 9000 j=1,nrhoin(1)
c        frhoin(j,i)=frhoin(j,1)
c9000    continue
        do 9000 j=1,nrin(1)
        zrin(j,2)=zrin(j,1)
        rhorin(j,2)=rhorin(j,1)
9000    continue
      endif

c
c
c determine common grid spacings and number
c
c
c take largest grid spacing
c take smallest maximum
c
        if(eamalloy .lt. 1) then
          ntypes=2
        endif
        if(eamalloy .lt. 1) then
          dr = drin(1)
        else
          dr = adrin1
        endif
        drho = drhoin(1)
        rmax = rcutguess1
        rhomax = (nrhoin(1)-1)*drhoin(1)
c       do 77 i1=1,pairnum
c         write(*,*)'Guess',i1,'Rcut:',rcutguess(i1)
c         write(*,*)'Guess',i1,'nr:',nint(rcutguess(i1)/adrin(i1))+1
c         write(*,*)'Guess',i1,'dr:',adrin(i1)
        if (prntflag.eq.1) then
         if(eamalloy .lt. 1) then
          if(logflag.eq.1) then
          write(11,*)'Guess 1 Rcut:',rcutguess1
          write(11,*)'Guess 1 nr:',nint(rcutguess1/drin(1))+1
          write(11,*)'Guess 1 dr:',drin(1)
           if (fixcutoff.lt.1) then
            write(11,*)'Guess 2 Rcut:',rcutguess2
           endif
          else
          write(*,*)'Guess 1 Rcut:',rcutguess1
          write(*,*)'Guess 1 nr:',nint(rcutguess1/drin(1))+1
          write(*,*)'Guess 1 dr:',drin(1)
           if (fixcutoff.lt.1) then
            write(*,*)'Guess 2 Rcut:',rcutguess2
           endif
          endif
         else
          write(*,*)'Guess 1 Rcut:',rcutguess1
          write(*,*)'Guess 1 nr:',nint(rcutguess1/adrin1)+1
          write(*,*)'Guess 1 dr:',adrin1
         endif
         if( pairnum .eq. 3 .and. eamalloy .gt. 0) then
          write(*,*)'Guess 2 Rcut:',rcutguess2
          write(*,*)'Guess 2 nr:',nint(rcutguess2/adrin2)+1
          write(*,*)'Guess 2 dr:',adrin2
          write(*,*)'Guess 3 Rcut:',rcutguess3
          write(*,*)'Guess 3 nr:',nint(rcutguess3/adrin3)+1
          write(*,*)'Guess 3 dr:',adrin3
         endif
        endif
c77      continue
c       do 78 i1=1,pairnum
c       dr = amax1(dr,adrin(i1))
c       rmaxi = rcutguess(i1)
c       rmax = amax1(rmax,rmaxi)
        if(eamalloy .lt. 1) then
        aadrin1=drin(1)
        dr = amax1(dr,aadrin1)
        else
        aadrin1=adrin1
        dr = amax1(dr,aadrin1)
        endif
        rmaxi = rcutguess1
        rmax = amax1(rmax,rmaxi)
        if (fixcutoff.lt.1) then
         aadrin2=adrin2
         dr = amax1(dr,aadrin2)
         rmaxi = rcutguess2
         rmax = amax1(rmax,rmaxi)
        endif
        if( pairnum .eq. 3 .and. eamalloy .gt. 0) then
         aadrin2=adrin2
         aadrin3=adrin3
         dr = amax1(dr,aadrin2)
         rmaxi = rcutguess2
         rmax = amax1(rmax,rmaxi)
         dr = amax1(dr,aadrin3)
         rmaxi = rcutguess3
         rmax = amax1(rmax,rmaxi)
        endif
c78      continue
7979    format(a6,i1,a6,g24.16)
       if (prntflag.eq.1) then
        do 79 i1=1,ntypes
          write(*,7979)' Type ',i1,' Rcut:',rcut(i1)
          write(*,7979)' Type ',i1,'   nr:',nrin(i1)
          write(*,7979)' Type ',i1,'   dr:',drin(i1)
          write(*,7979)' Type ',i1,' nrho:',nrhoin(i1)
          write(*,7979)' Type ',i1,' drho:',drhoin(i1)
79      continue
       endif
        do 80 i1=1,ntypes
        dr = amax1(dr,drin(i1))
        drho = amax1(drho,drhoin(i1))
        rmaxi = (nrin(i1)-1)*drin(i1)
        rhomaxi = (nrhoin(i1)-1)*drhoin(i1)
        rmax = amax1(rmax,rmaxi)
        rhomax = amax1(rhomax,rhomaxi)
80      continue
        rcutoff = rmax
        nr = nint(rmax/dr)+1
        nrho = nint(rhomax/drho)+1
       if (prntflag.eq.1) then
        write(*,*)'All Rcut:',rcutoff
        write(*,*)'All nr:',nr
        write(*,*)'All dr:',dr
        write(*,*)'All rmax:',rmax
        write(*,*)'All nrho:',nrho
        write(*,*)'All drho:',drho
        write(*,*)'All rhomax:',rhomax
       endif
       if(logflag.eq.1) then
        write(11,*)'All Rcut:',rcutoff
        write(11,*)'All nr:',nr
        write(11,*)'All dr:',dr
        write(11,*)'All rmax:',rmax
        write(11,*)'All nrho:',nrho
        write(11,*)'All drho:',drho
        write(11,*)'All rhomax:',rhomax
       endif
c
c set up the z(r) and rho(r) grids
c
        do 90 i1=1,ntypes
        do 85 j=1,nr
          rr = (j-1)*dr
c
c  do four-point lagrange interpolation
c
          pp = rr/drin(i1) + 1.0
          kk = pp
          kk = min0(kk,nrin(i1)-2)
          kk = max0(kk,2)
          pp = pp - kk
c       make sure that p is less than 2.0
c       then if r is out of range, p = 2.0 and rhor = last value of rhorin
          pp = amin1(pp,2.)
          cof1 = -0.166666667*pp*(pp-1.)*(pp-2.)
          cof2 = 0.5*(pp**2-1.)*(pp-2.)
          cof3 = -0.5*pp*(pp+1.)*(pp-2.)
          cof4 = 0.166666667*pp*(pp**2-1.)
         rhor(j,i1) = cof1*rhorin(kk-1,i1)
     1      + cof2*rhorin(kk,i1)
     2      + cof3*rhorin(kk+1,i1)
     3      + cof4*rhorin(kk+2,i1)
         zrtemp(j,i1) = cof1*zrin(kk-1,i1)
     1      + cof2*zrin(kk,i1)
     2      + cof3*zrin(kk+1,i1)
     3      + cof4*zrin(kk+2,i1)
85      continue
90      continue
       if (eamalloy.lt.1 .and. fixcutoff.lt.1) then
        do 91 i1=1,ntypes
        do 84 j=1,nr
          rr = (j-1)*dr
c
c  do four-point lagrange interpolation
c
          pp = rr/adrin2 + 1.0
          kk = pp
          kk = min0(kk,nrin(i1)-2)
          kk = max0(kk,2)
          pp = pp - kk
c       make sure that p is less than 2.0
c       then if r is out of range, p = 2.0 and rhor = last value of rhorin
          pp = amin1(pp,2.)
          cof1 = -0.166666667*pp*(pp-1.)*(pp-2.)
          cof2 = 0.5*(pp**2-1.)*(pp-2.)
          cof3 = -0.5*pp*(pp+1.)*(pp-2.)
          cof4 = 0.166666667*pp*(pp**2-1.)
         rhor(j,i1) = cof1*rhorin(kk-1,i1)
     1      + cof2*rhorin(kk,i1)
     2      + cof3*rhorin(kk+1,i1)
     3      + cof4*rhorin(kk+2,i1)
84      continue
91      continue
       endif
  
        if (eamalloy .gt. 0) then
c        do 87 i1=1,pairnum
        do 86 j=1,nr
          rr = (j-1)*dr
c
c  do four-point lagrange interpolation
c
          pp = rr/adrin1 + 1.0
          kk = pp
          kk = min0(kk,anrin1-2)
          kk = max0(kk,2)
          pp = pp - kk
c       make sure that p is less than 2.0
c       then if r is out of range, p = 2.0 and zr = last value of zrin
          pp = amin1(pp,2.)
          cof1 = -0.166666667*pp*(pp-1.)*(pp-2.)
          cof2 = 0.5*(pp**2-1.)*(pp-2.)
          cof3 = -0.5*pp*(pp+1.)*(pp-2.)
          cof4 = 0.166666667*pp*(pp**2-1.)
         azrtemp1(j) = cof1*azrin1(kk-1)
     1      + cof2*azrin1(kk)
     2      + cof3*azrin1(kk+1)
     3      + cof4*azrin1(kk+2)
86      continue

        if (pairnum.eq.3) then

        do 87 j=1,nr
          rr = (j-1)*dr
c
c  do four-point lagrange interpolation
c
          pp = rr/adrin2 + 1.0
          kk = pp
          kk = min0(kk,anrin2-2)
          kk = max0(kk,2)
          pp = pp - kk
c       make sure that p is less than 2.0
c       then if r is out of range, p = 2.0 and zr = last value of zrin
          pp = amin1(pp,2.)
          cof1 = -0.166666667*pp*(pp-1.)*(pp-2.)
          cof2 = 0.5*(pp**2-1.)*(pp-2.)
          cof3 = -0.5*pp*(pp+1.)*(pp-2.)
          cof4 = 0.166666667*pp*(pp**2-1.)
         azrtemp2(j) = cof1*azrin2(kk-1)
     1      + cof2*azrin2(kk)
     2      + cof3*azrin2(kk+1)
     3      + cof4*azrin2(kk+2)
87      continue

        do 88 j=1,nr
          rr = (j-1)*dr
c
c  do four-point lagrange interpolation
c
          pp = rr/adrin3 + 1.0
          kk = pp
          kk = min0(kk,anrin3-2)
          kk = max0(kk,2)
          pp = pp - kk
c       make sure that p is less than 2.0
c       then if r is out of range, p = 2.0 and zr = last value of zrin
          pp = amin1(pp,2.)
          cof1 = -0.166666667*pp*(pp-1.)*(pp-2.)
          cof2 = 0.5*(pp**2-1.)*(pp-2.)
          cof3 = -0.5*pp*(pp+1.)*(pp-2.)
          cof4 = 0.166666667*pp*(pp**2-1.)
         azrtemp3(j) = cof1*azrin3(kk-1)
     1      + cof2*azrin3(kk)
     2      + cof3*azrin3(kk+1)
     3      + cof4*azrin3(kk+2)
88      continue

        endif

c87      continue
c
c set up the f(rho) grid
c
        do 101 i1=1,ntypes
        do 102 j=1,nrho
        rr = (j-1)*drho
c
c  do four-point lagrange interpolation
c
        pp = rr/drhoin(i1) + 1.0
        kk = pp
        kk = min0(kk,nrhoin(i1)-2)
        kk = max0(kk,2)
        pp = pp - kk
c       make sure that p is less than 2.0
c       then if r is out of range, p = 2.0 and frho = last value of frhoin
        pp = amin1(pp,2.)
        cof1 = -0.166666667*pp*(pp-1.)*(pp-2.)
        cof2 = 0.5*(pp**2-1.)*(pp-2.)
        cof3 = -0.5*pp*(pp+1.)*(pp-2.)
        cof4 = 0.166666667*pp*(pp**2-1.)
        frho(j,i1) = cof1*frhoin(kk-1,i1)
     1      + cof2*frhoin(kk,i1)
     2      + cof3*frhoin(kk+1,i1)
     3      + cof4*frhoin(kk+2,i1)
102     continue
101     continue

        endif
c
c set up the z2 grid (zi*zj)
c
        do 110 i1=1,ntypes
        do 110 i2=1,ntypes
        do 105 j=1,nr
          if( i1 .eq. i2 ) then
             z2r(j,i1,i2) = zrtemp(j,i1)
          else
           if( ntypes .eq. 3 ) then
             if( ( i1 .eq. 2 .and. i2 .eq. 3 ) .or. 
     $           ( i1 .eq. 3 .and. i2 .eq. 2 ) ) then
                 z2r(j,i1,i2) = azrtemp3(j)
             elseif( ( i1 .eq. 1 .and. i2 .eq. 3 ) .or. 
     $               ( i1 .eq. 3 .and. i2 .eq. 1 ) ) then
                 z2r(j,i1,i2) = azrtemp2(j)
             else
                 z2r(j,i1,i2) = azrtemp1(j)
             endif
           else
             if (eamalloy .gt. 0) then
                 z2r(j,i1,i2) = azrtemp1(j)
             else
              if(phitype.eq.1.or.phitype.ge.3 .or. phitype .ge. 7) then 
                 z2r(j,i1,i2) = zrtemp(j,i1)
              elseif(phitype.eq.2) then 
                 z2r(j,i1,i2) = zrtemp(j,i1)*zrtemp(j,i2)
c             elseif(phitype.eq.4) then 
c                z2r(j,i1,i2) = 27.2*0.529*zrtemp(j,i1)*zrtemp(j,i2)
              endif
             endif
           endif
          endif
105     continue
110     continue

c
c       set the lattice constant to that for type 1 by default
      if (eamalloy .lt. 1) then
        alat = lata0
      else
        alat = blat(1)
      endif
      latty = lat(1)
c
c
c  now set up the dense grids
c
      nrhoar = nrho
      drhoar = drho
      nrar = nr
      drar = dr
      do 500 i = 1,ntypes
       if (eamalloy .gt. 0) then
         do 515 j = 1,nrhoar
            frhoar(j,i) = frho(j,i)
  515    continue
         frhoar1(1,i) = frhoar(2,i)-frhoar(1,i)
         frhoar1(2,i) = 0.5*(frhoar(3,i)-frhoar(1,i))
         frhoar1(nrhoar-1,i) = 0.5*(frhoar(nrhoar,i)-frhoar(nrhoar-2,i))
         frhoar1(nrhoar,i) = frhoar(nrhoar,i)-frhoar(nrhoar-1,i)
         do 520 j = 3,nrhoar-2
            frhoar1(j,i) = ((frhoar(j-2,i)-frhoar(j+2,i)) +
     $                  8.*(frhoar(j+1,i)-frhoar(j-1,i)))/12.
  520       continue
         do 525 j = 1,nrhoar-1
            frhoar2(j,i) = 3.*(frhoar(j+1,i)-frhoar(j,i))
     $                   - 2.*frhoar1(j,i) - frhoar1(j+1,i)
            frhoar3(j,i) = frhoar1(j,i) + frhoar1(j+1,i)
     $                   - 2.*(frhoar(j+1,i)-frhoar(j,i))
  525       continue
         frhoar2(nrhoar,i) = 0.
         frhoar3(nrhoar,i) = 0.
         do 528 j = 1,nrhoar
            frhoar4(j,i) = frhoar1(j,i)/drhoar
            frhoar5(j,i) = 2.*frhoar2(j,i)/drhoar
            frhoar6(j,i) = 3.*frhoar3(j,i)/drhoar
  528       continue
       endif
         do 530 j = 1,nrar
            rhorar(j,i) = rhor(j,i)
  530       continue
         rhorar1(1,i) = rhorar(2,i)-rhorar(1,i)
         rhorar1(2,i) = 0.5*(rhorar(3,i)-rhorar(1,i))
         rhorar1(nrar-1,i) = 0.5*(rhorar(nrar,i)-rhorar(nrar-2,i))
         rhorar1(nrar,i) = 0.
         do 532 j = 3,nrar-2
            rhorar1(j,i) = ((rhorar(j-2,i)-rhorar(j+2,i)) +
     $                  8.*(rhorar(j+1,i)-rhorar(j-1,i)))/12.
  532       continue
         do 535 j = 1,nrar-1
            rhorar2(j,i) = 3.*(rhorar(j+1,i)-rhorar(j,i))
     $                   - 2.*rhorar1(j,i) - rhorar1(j+1,i)
            rhorar3(j,i) = rhorar1(j,i) + rhorar1(j+1,i)
     $                   - 2.*(rhorar(j+1,i)-rhorar(j,i))
  535       continue
         rhorar2(nrar,i) = 0.
         rhorar3(nrar,i) = 0.
         do 538 j = 1,nrar
            rhorar4(j,i) = rhorar1(j,i)/drar
            rhorar5(j,i) = 2.*rhorar2(j,i)/drar
            rhorar6(j,i) = 3.*rhorar3(j,i)/drar
  538       continue
      i1 = i
      do 540 i2 = 1,ntypes
         do 550 j = 1,nrar
            z2rar(j,i1,i2) = z2r(j,i1,i2)
  550       continue
         z2rar1(1,i1,i2) = z2rar(2,i1,i2)-z2rar(1,i1,i2)
         z2rar1(2,i1,i2) = 0.5*(z2rar(3,i1,i2)-z2rar(1,i1,i2))
         z2rar1(nrar-1,i1,i2) =
     $        0.5*(z2rar(nrar,i1,i2)-z2rar(nrar-2,i1,i2))
         z2rar1(nrar,i1,i2) = 0.
         do 560 j = 3,nrar-2
            z2rar1(j,i1,i2) = ((z2rar(j-2,i1,i2)-z2rar(j+2,i1,i2)) +
     $                  8.*(z2rar(j+1,i1,i2)-z2rar(j-1,i1,i2)))/12.
  560       continue
         do 565 j = 1,nrar-1
            z2rar2(j,i1,i2) = 3.*(z2rar(j+1,i1,i2)-z2rar(j,i1,i2))
     $                   - 2.*z2rar1(j,i1,i2) - z2rar1(j+1,i1,i2)
            z2rar3(j,i1,i2) = z2rar1(j,i1,i2) + z2rar1(j+1,i1,i2)
     $                   - 2.*(z2rar(j+1,i1,i2)-z2rar(j,i1,i2))
  565       continue
         z2rar2(nrar,i1,i2) = 0.
         z2rar3(nrar,i1,i2) = 0.
         do 568 j = 1,nrar
            z2rar4(j,i1,i2) = z2rar1(j,i1,i2)/drar
            z2rar5(j,i1,i2) = 2.*z2rar2(j,i1,i2)/drar
            z2rar6(j,i1,i2) = 3.*z2rar3(j,i1,i2)/drar
  568       continue
  540    continue
  500    continue
c     --------- End of EAM splines initialization ----------
c     ------------------------------------------------------
c
      rmssum=1000000.0

c     ******************************************************
c     Generate PURE element F(rho) function
      if (eamalloy .lt. 1) then
c     Calculating Rose et al. Universal Binding Energy for given lattice
c
c     For PURE l1fcc=0
      l1fcc=0
      t1=1
      t2=1
      call loadfcc100(t1,t2,l1fcc,lata0,stacks)
      call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms4:',rmssum
        endif
         return
       endif
      if(lattype.eq.0) then
         omega=lata0*lata0*lata0/4.d0
      else if(lattype.eq.1) then
         omega=lata0*lata0*lata0/2.d0
      else if(lattype.eq.2) then
         omega=2.*sqrt3*hcpc*lata0*lata0/4.
      endif
c     angstrom^3/atom -> cm^3/mol 
c     when multiplied with bmod, bmod will be eV
      omega=omega*0.62415
c
c     latacut = a for nnd=rcutoff
      latacut=sqrt2*rcutoff
      xa=(latacut/lata0-1.d0)
      apart=sqrt((9.0*bmod*omega)/esub)
      astarcut=apart*xa
c     a*=a*cut
      latano=(lataf-latai)/latstep
      latno=int(latano)
      iat=vacatom
c     print*,'iat=',iat
c     print*,'neigh(iat)=',neigh(iat)
      do 40 i=1,latno+1
      lata=latai+(i-1)*latstep
c    Calculating Rho(r) and Phi(r) over neighbors
      eqrho=0.d0
      eqphi=0.d0
      rhoij=0.d0
      do 41 j=1,neigh(iat)
       rdist=(lata/lata0)*rdis(iat,j,1)
       jat=rdis(iat,j,2)
c       if(lata.eq.lata0)then
c         write(72,7272)iat,jat,phi
c       endif
        rdrar = 1.0/drar
        pij = rdist*rdrar + 1.0
        kij = pij
        kij = min0(kij,nrar-1)
        pij = pij - kij
c          make sure that p is less than 1.0
c          then if r is out of range, p = 1.0 and rho = last value of rhor
        pij = amin1(pij,1.0)
        rhoij = rhoij +
     $  ((rhorar3(kij,1)*pij+rhorar2(kij,1))*pij+
     $    rhorar1(kij,1))*pij+rhorar(kij,1)
        phi = ((z2rar3(kij,1,1)*pij+z2rar2(kij,1,1))*pij+
     $    z2rar1(kij,1,1))*pij+z2rar(kij,1,1)
        if(phitype.eq.1.or.phitype.ge.3 .or. phitype.ge.7) then
          eqphi=eqphi+0.5*phi
        else
          eqphi=eqphi+0.5*(phi/rdist)
        endif
c       print*,'rhoij:',rhoij
41      continue
        phiequil(i)=eqphi
        rhoeam(i)=rhoij
c
c
c     ----------------------------------------------------
c      The crystal energy should be zero when nearest
c       neighbor distance equals to rcutoff. To accomplish
c       this Foiles (1985b) suggested to modify universal
c       function as following.
c      
      xa=(lata/lata0-1.d0)
      apart=sqrt((9.0*bmod*omega)/esub)
      astar=apart*xa
      qstar=astarcut
c     -----------------------------
c     do-until F77 implemented loop
c     do
42    continue
       epslon=(1.d0+qstar)*exp(-qstar)
       newq=sqrt(1.d0-epslon)*apart*xa
       bindf=(1.d0+newq)*exp(-newq)
       bindfmod=(bindf-epslon)/(1.d0-epslon)
       qstar=0.5*(qstar+sqrt(1.d0-bindfmod)*astarcut)
       bindnewf=-esub*(1.d0+qstar)*exp(-qstar)
      if (bindnewf.gt.1.0E-14) goto 42
c     until bindnewf < Tolerance
c     -----------------------------
      epslon=(1.d0+qstar)*exp(-qstar)
      newq=sqrt(1.d0-epslon)*apart*xa
      bindf=(1.d0+newq)*exp(-newq)
      bindfmod=(bindf-epslon)/(1.d0-epslon)
c     ----------------------------------------------------
      binde=-esub*bindfmod
      fa2(i)=binde-phiequil(i)
c     write(*,'(5(E25.16))')lata,binde,fa2(i),phiequil(i),rhoeam(i)
40    continue
7272  format(2(i5,1x),E25.16)


      ream=0.
      j=1
      do 5000 i=1,latno+1
        if(rhoeam(latno+2-i).gt.0.d0) then
          if(ream.lt.rhoeam(latno+2-i)) then
             ream=rhoeam(latno+2-i)
          endif
          if((rhoa(j).ne.rhoeam(latno+2-i)).and.
     .          fa(j).ne.fa2(latno+2-i)) then
            j=j+1
            rhoa(j)=rhoeam(latno+2-i)
            fa(j)=fa2(latno+2-i)
          endif
        else
          rhofirst=rhoeam(latno+2-i)
          fafirst=fa2(latno+2-i)
        endif
5000  continue
      rhoa(1)=0.
      fa(1)=0.
      if(rhoa(j).lt.rhoend) then
        if (logflag.eq.1) then
        write(11,*)'WARNING! PROC(',procnum,') ream:',ream,
     .'rhoa:',rhoa(j),' rhoend:',rhoend
        endif
        if (prntflag.eq.1) then
          write(*,*)'rho:',rhoa(j),' rhoeam:',(rhoeam(i),i=1,latno+1),
     .    ' <rhoend is out off scope for F spline'
        endif
        if (logflag.eq.1) then
          write(11,*)'rhoa:',(rhoa(i),i=1,j),
     .    ' rhoeam:',(rhoeam(i),i=1,latno+1),
     .    ' <rhoend is out off scope for F spline'
        endif
        rmssum=3000000.0
        return
      end if
      nknot=j
      if(prntflag .gt. 0) then
        write(*,*)nknot,'-> spline ->',nrhoin(1)
      endif
      if (logflag.eq.1) then
        write(11,*)nknot,'-> spline ->',nrhoin(1)
      endif
c      do 5500 i=1,nknot
c      write(*,*)rhoa(i),fa(i)
c5500  continue
      call splineinit(rhoa,fa,nknot,y2p)
      do 300 i=1,nrhoin(1)
        rrrho=(i-1)*drhoin(1)
        call splineinter(rhoa,fa,y2p,nknot,rrrho,finter)
        frhoin(i,1)=finter
        frhoin(i,2)=finter
300   continue

c
c set up the f(rho) dense grid
c
        do 2700 i1=1,ntypes
        do 2795 j=1,nrho
        rr = (j-1)*drho
c
c  do four-point lagrange interpolation
c
        pp = rr/drhoin(i1) + 1.0
        kk = pp
        kk = min0(kk,nrhoin(i1)-2)
        kk = max0(kk,2)
        pp = pp - kk
c       make sure that p is less than 2.0
c       then if r is out of range, p = 2.0 and rhor = last value of rhorin
        pp = amin1(pp,2.)
        cof1 = -0.166666667*pp*(pp-1.)*(pp-2.)
        cof2 = 0.5*(pp**2-1.)*(pp-2.)
        cof3 = -0.5*pp*(pp+1.)*(pp-2.)
        cof4 = 0.166666667*pp*(pp**2-1.)
        frho(j,i1) = cof1*frhoin(kk-1,i1)
     1      + cof2*frhoin(kk,i1)
     2      + cof3*frhoin(kk+1,i1)
     3      + cof4*frhoin(kk+2,i1)
2795      continue
2700     continue

c
c set up the f(rho) spline
c
      do 1700 i = 1,ntypes
         do 1710 j = 1,nrhoar
            frhoar(j,i) = frho(j,i)
1710       continue
         frhoar1(1,i) = frhoar(2,i)-frhoar(1,i)
         frhoar1(2,i) = 0.5*(frhoar(3,i)-frhoar(1,i))
         frhoar1(nrhoar-1,i) = 0.5*(frhoar(nrhoar,i)-frhoar(nrhoar-2,i))
         frhoar1(nrhoar,i) = frhoar(nrhoar,i)-frhoar(nrhoar-1,i)
         do 1720 j = 3,nrhoar-2
            frhoar1(j,i) = ((frhoar(j-2,i)-frhoar(j+2,i)) +
     $                  8.*(frhoar(j+1,i)-frhoar(j-1,i)))/12.
1720       continue
         do 1725 j = 1,nrhoar-1
            frhoar2(j,i) = 3.*(frhoar(j+1,i)-frhoar(j,i))
     $                   - 2.*frhoar1(j,i) - frhoar1(j+1,i)
            frhoar3(j,i) = frhoar1(j,i) + frhoar1(j+1,i)
     $                   - 2.*(frhoar(j+1,i)-frhoar(j,i))
1725       continue
         frhoar2(nrhoar,i) = 0.
         frhoar3(nrhoar,i) = 0.
         do 1730 j = 1,nrhoar
            frhoar4(j,i) = frhoar1(j,i)/drhoar
            frhoar5(j,i) = 2.*frhoar2(j,i)/drhoar
            frhoar6(j,i) = 3.*frhoar3(j,i)/drhoar
1730       continue
1700    continue

      endif
c     Generating of PURE F(rho) ends here
c     ******************************************************

c        
c       Writing EAM function tables to given eamfileout file
c
  20    format(80a1)
  30    format(i5,2g15.5,a8)
 150    format(i5)
9901    format(i5,e24.16,i5,2e24.16)
9902    format(5e24.16)
2     format ( '  Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &         i2.2, ':', i2.2, ':', i2.2 )

      if (prntflag.eq.1) then
c       write(*,*)'Writing EAM tables to file: ',eamfileout
        open(unit=75,file=eamfileout)
        open(unit=77,file='eamout.setfl')
c     if(fileformat.lt.1) then
      call date_and_time(date,time,zone,values)
      write(dateh,2)INT(values(3)),INT(values(2)),INT(values(1)),
     &              INT(values(5)),INT(values(6)),INT(values(7))
      if(lcharnum.lt.48) then
      read(dateh,'(32a1)',END=23)(fheader(i,1),i=lcharnum+1,lcharnum+33)
23    continue
      endif
       if(fileformat.lt.1) then
        write(75,20)(fheader(j,1),j=1,80)
        write(75,30) ielement(1), amass(1), blat(1), lat(1)
        write(75,9901) nrhoin(1), drhoin(1), nrin(1), drin(1), rcut(1)
        write(75,9902) (frhoin(j,1),j=1,nrhoin(1))
        if (phitype .eq. 7) then
c        WARNING funcfl
         write(75,9902) (sqrt(zrin(j,1)*(j-1)*drin(1) /
     &                                       27.2*0.529),j=1,nrin(1))
c       In DYNAMO funcfl z2r(j,i1,i2) = 27.2*0.529*zrtemp(j,i1)*zrtemp(j,i2)
        else
          write(75,9902) (zrin(j,1),j=1,nrin(1))
        endif
        write(75,9902) (rhorin(j,1),j=1,nrin(1))
        write(75,*)''
        read(dateh,'(32a1)',END=22)(fheader(i,2),i=1,33)
22      continue
        if (ntypes.lt.2) then
          ntyps=2
        else
          ntyps=ntypes
        endif
        write(77,20)(fheader(j,1),j=1,80)
        write(77,20)(fheader(j,2),j=1,80)
        write(77,*)'eamout.setfl file: Setfl EAM function files'
        write(77,150) ntyps
        write(77,9901) nrho, drho, nr, dr, rcutoff
        do 161 i1=1,ntyps
        write(77,30) ielement(i1),amass(i1),blat(i1),lat(i1)
        write(77,9902) (frho(j,i1),j=1,nrho)
        write(77,9902) (rhor(j,i1),j=1,nr)
161     continue
        do 162 i1=1,ntyps
        do 162 i2=1,i1
162     write(77,9902) (z2r(j,i1,i2)*((j-1)*dr),j=1,nr)
        write(77,*)''
       else
        write(75,20)(fheader(j,1),j=1,80)
        write(75,20)(fheader(j,2),j=1,80)
        write(75,20)(fheader(j,3),j=1,80)
        write(75,150) ntypes
        if (ntypes.gt.nelmax) then
          write(6,*)'error: number of types is greater than nelmax'
          stop
        endif
        rcutall=rcutoff
        write(75,9901) nrho, drho, nr, dr, rcutall
        do 160 i1=1,ntypes
        write(75,30) ielement(i1),amass(i1),blat(i1),lat(i1)
        write(75,9902) (frho(j,i1),j=1,nrho)
        write(75,9902) (rhor(j,i1),j=1,nr)
160     continue
        do 170 i1=1,ntypes
        do 170 i2=1,i1
170     write(75,9902) (z2r(j,i1,i2)*((j-1)*dr),j=1,nr)
        write(75,*)''
       endif
       close(75)
      else
        rcutall=rcutoff
      endif


c****************************************************
c*   Calculation of Equilibrium Lattice Constant    *
c*                for B2 BCC Alloys                 *
c****************************************************
       rmssum=1000000.0
       
c      ---------------------------------------------               
      if(eamalloy .lt. 1) then
c     PURE CALCULATION GOES HERE AFTER
c      ---------------------------------------------               
       if(latype .eq. 0) then
        if (prntflag.eq.1) then
         write(*,*)'------- PURE FCC -------'
        endif
        if (prntflag.eq.1) then
         write(*,*)'PROC:',procnum,' part:',prmid,' PURE FCC '
        endif
        l1fcc=0
        t1=1
        t2=2
        call loadfcc100(t1,t2,l1fcc,lata0,stacks)
        call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
        if(istop .eq. 1) then
         if (prntflag.eq.1) then
          write(*,*)'rms4-2:',rmssum
         endif
          return
        endif
        if (prntflag.eq.1) then
         write(*,*)'FCC Test a_0:',lata0
        endif
        call calclatcon(minenergy,latis,lata0,neigh,rdis,iter,prntflag)
        l12alat(1)=latis
        if (prntflag.eq.1) then
         write(*,*)'FCC a_0:',latis
         write(*,*)'FCC Ecoh:',minenergy/natoms
        endif
        l12ecoh(1)=minenergy/natoms
        call loadfcc100(t1,t2,l1fcc,latis,stacks)
        call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
        if(istop .eq. 1) then
         if (prntflag.eq.1) then
          write(*,*)'rms5-2:',rmssum
         endif
        tn=1
        typs(1)=t1
        call calctotetype(tn,typs,tnum,embed,etot1,neigh,rdis)
        etot1=etot1/natoms
        if (prntflag.eq.1) then
         print*,'FCC Ecoh:',etot1
         print*,'FCC Embed1:',embed(1),'num:',tnum(1)
        endif
          return
        endif


       else if(latype .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'------- PURE BCC -------'
        endif
        if (prmid.gt.0) then
         write(*,*)'PROC:',procnum,' part:',prmid,' PURE BCC '
        endif
        b2=0
        t1=1
        t2=2
       call loadbcc100(t1,t2,b2,lata0,stacks)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms2-2:',rmssum
        endif
         return
       endif
       if (prntflag.eq.1) then
        write(*,*)'BCC Test a_0:',lata0
       endif
       call calclatcon(minenergy,latis,lata0,neigh,rdis,iter,prntflag)
       l12alat(2)=latis
       if (prntflag.eq.1) then
        write(*,*)'BCC a_0:',latis
       endif
       call loadbcc100(t1,t2,b2,latis,stacks)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms12-2:',rmssum
        endif
         return
       endif
       call calctote(etot1,neigh,rdis)
       etot1=etot1/natoms
       l12ecoh(1)=etot1
       if (prntflag.eq.1) then
        print*,'BCC Ecoh:',etot1
       endif
       tn=2
       typs(1)=t1
       typs(2)=t2
       call calctotetype(tn,typs,tnum,embed,etot1,neigh,rdis)
       etot1=etot1/natoms
       if (prntflag.eq.1) then
        print*,'BCC Ecoh:',etot1
        print*,'BCC Embed1:',embed(1),'num:',tnum(1)
        print*,'BCC Embed2:',embed(2),'num:',tnum(2)
       endif



       else if(latype .eq. 2) then
        if (prntflag.eq.1) then
         write(*,*)'------- PURE HCP -------'
        endif
        if (prmid.gt.0) then
         write(*,*)'PROC:',procnum,' part:',prmid,' HCP'
        endif
         l1=1
         t1=1
         t2=2
         hcpc=caratio*lata0
        if (prntflag.eq.1) then
         write(*,*)'HCP Test c/a:',caratio
         write(*,*)'HCP Test a_0:',lata0
        endif
        call loadhcp111(t1,t2,l1,lata0,caratio,stacks)
        call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
        if(istop.eq.1) then 
         if (prntflag.eq.1) then
          write(*,*)'rms13-2:',rmssum
         endif
          return
        endif
c      if (prntflag.eq.1) then
c       write(*,*)'CoNi Test a_0:',lata0
c      endif
        rmssum=10000.0
       call calclatca(minenergy,latis,lata0,c,caratio,neigh,rdis,rdisd,
     &rmssum,prntflag)
        if (rmssum.gt.99999.0) then
          rmssum=100000.0
          write(*,*)'calclatca: no min latis,c/a'
          latis=2.0
          c=1.0
        endif
        l12alat(1)=latis
        hexca(1)=c
        if (prntflag.eq.1) then
         write(*,*)'HCP rmssum:',rmssum
         write(*,*)'HCP a_0:',latis
         write(*,*)'HCP c/a:',c
        endif
        hcpc=c*latis
        call loadhcp111(t1,t2,l1,latis,c,stacks)
        call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
        if(istop.eq.1) then 
         if (prntflag.eq.1) then
          write(*,*)'rms14-2:',rmssum
         endif
          return
        endif
c       call calcminehex(etot1,latis,lata0,c,c0,neigh,rdis,rdisd)
        call calctote(etot1,neigh,rdis)
        etot1=etot1/natoms
        l12ecoh(1)=etot1
        if (prntflag.eq.1) then
         print*,'CoNi HEX Ecoh:',etot1
        endif
        tn=2
        typs(1)=t1
        typs(2)=t2
        call calctotetype(tn,typs,tnum,embed,etot1,neigh,rdis)
        etot1=etot1/natoms
        if (prntflag.eq.1) then
         print*,'HCP Ecoh:',etot1
         print*,'HCP Embed1:',embed(1),'num:',tnum(1)
         print*,'HCP Embed2:',embed(2),'num:',tnum(2)
        endif


       endif

c   *****************************************************
c   * EAM fcn file ready to use in Sublimation Energy,  *
c   *     Vacancy Formation Energy, Bulk Modulus,       *
c   *     and Elastic Constant Calculations for FCC     *
c   *****************************************************
    
      do 1845 iat=1,natoms
        e(iat)=0.0
        fe(iat)=0.0
        fp(iat)=0.0
        fpp(iat)=0.0
        phie(iat)=0.0
        rho(iat)=0.0
        evacf(iat)=0.0
1845  continue
      
      rdrar = 1.0/drar 
      do 1875 iat=1,natoms
c
c     compute the contribution to rho(i) from particle j
c     Calculating the embedding energy F(Rho) in fe(i)
c      and Pair Energy contribution Phi(Rij) in phie(i)
c
ccdir$ novector
      do 1855 j=1,neigh(iat)
      rdist=rdis(iat,j,1)
      jat=rdis(iat,j,2)
      pij = rdist*rdrar + 1.0
      kij = pij
      kij = min0(kij,nrar-1)
      pij = pij - kij
c       make sure that p is less than 1.00
c       then if r is out of range, p = 1.00 and rho = last value of rhor
      pij = amin1(pij,1.0)
      ity=itype(iat)
      jty=itype(jat)
      rhojnei(iat,j) = ((rhorar3(kij,jty)*pij+
     $   rhorar2(kij,jty))*pij+
     $   rhorar1(kij,jty))*pij+
     $   rhorar(kij,jty)
      rhoinei(iat,j) = ((rhorar3(kij,ity)*pij+
     $   rhorar2(kij,ity))*pij+
     $   rhorar1(kij,ity))*pij+
     $   rhorar(kij,ity)
      rho(iat) = rho(iat) + rhojnei(iat,j)
      z2ij = ((z2rar3(kij,ity,jty)*pij+
     $          z2rar2(kij,ity,jty))*pij+
     $          z2rar1(kij,ity,jty))*pij+
     $          z2rar(kij,ity,jty)
      z2pij = (z2rar6(kij,ity,jty)*pij+z2rar5(kij,ity,jty))*pij+
     $         z2rar4(kij,ity,jty)
c     if(phitype) then
        phiij = z2ij
c     else
c       pij=1.0/rdist
c       phiij = z2ij * pij
c     endif
      phie(iat) = phie(iat) + 0.5 * phiij
1855  continue
ccdir$ vector
c        
1875  continue
c
c       store fsubi in e(i), the derivative of fsubi in fp
c       and the second derivative of fsubi in fpp
c
      rdrhoar = 1.0/drhoar
      do 2005 iat=1,natoms
      pij = rho(iat)*rdrhoar + 1.0
      kij = pij
      kij = max(1,min(kij,nrhoar-1))
      pij = pij - kij
      ity=itype(iat)
      fe(iat) = ((frhoar3(kij,ity)*pij+
     $        frhoar2(kij,ity))*pij+
     $        frhoar1(kij,ity))*pij+
     $        frhoar(kij,ity)
      fp(iat) = (frhoar6(kij,ity)*pij+
     $                  frhoar5(kij,ity))*pij+
     $                  frhoar4(kij,ity)
      fpp(iat) = (2.0*frhoar6(kij,ity)*pij+
     $                  frhoar5(kij,ity))/drhoar
2005  continue
c
c       store vacancy formation energy in evacf(i)
c
      rdrhoar = 1.0/drhoar
      do 3105 iat=1,natoms
        do 3115 j=1,neigh(iat)
        jat=rdis(iat,j,2)
        pij = (rho(jat)-rhoinei(iat,j))*rdrhoar + 1.0
        kij = pij
        kij = max(1,min(kij,nrhoar-1))
        pij = pij - kij
        ity=itype(iat)
        jty=itype(jat)
        evacf(iat) = evacf(iat) + ( fe(jat) - 
     $      (((frhoar3(kij,jty)*pij+
     $        frhoar2(kij,jty))*pij+
     $        frhoar1(kij,jty))*pij+
     $        frhoar(kij,jty))  )
3115    continue
      e(iat) = phie(iat) + fe(iat)
3105  continue

      do 3125 iat=1,natoms
      evacf(iat) = evacf(iat) + phie(iat)
3125  continue
      

      do 3156 iat=1,natoms
        ca(iat) = 0.0
        do 3157 cij=1,7
         econ(cij,iat) = 0.0
         bulkmod(cij,iat) = 0.0
         uu(cij) = 0.0
         ww(cij) = 0.0
         vv(cij) = 0.0
3157    continue
         do 3505 j=1,neigh(iat)
           jat=rdis(iat,j,2)
           do 3515 kc=1,3
3515       diski(kc) = rdisd(iat,j,kc)
c          call elascon(phitype,diski,jat,iat,uu,ww,vv,bulkmod) 
           call elastic(phitype,diski,jat,iat,uu,ww,vv,bulkmod) 
3505     continue
        econ(1,iat) = 0.5*uu(1) + fp(iat) *
     .                   ww(1) + fpp(iat) * vv(1)*vv(1)
        econ(2,iat) = 0.5*uu(2) + fp(iat) *
     .                   ww(2) + fpp(iat) * vv(1)*vv(2)
        econ(3,iat) = 0.5*uu(3) + fp(iat) *
     .                   ww(3) + fpp(iat) * vv(3)*vv(3)
        econ(4,iat) = 0.5*uu(4) + fp(iat) *
     .                   ww(4) + fpp(iat) * vv(1)*vv(4)
        econ(5,iat) = 0.5*uu(5) + fp(iat) *
     .                   ww(5) + fpp(iat) * vv(4)*vv(4)
        econ(6,iat) = 0.5*uu(6) + fp(iat) *
     .                   ww(6) + fpp(iat) * vv(5)*vv(6)
        econ(7,iat) = 0.5*uu(7) + fp(iat) *
     .                   ww(7) + fpp(iat) * vv(7)*vv(7)
        ca(iat)= fpp(iat)*(vv(1)*vv(2)-vv(7)*vv(7))
        ca(iat)= ca(iat)/ev_joule
        econ(1,iat) = econ(1,iat)/ev_joule
        econ(2,iat) = econ(2,iat)/ev_joule
        econ(3,iat) = econ(3,iat)/ev_joule
        econ(4,iat) = econ(4,iat)/ev_joule
        econ(5,iat) = econ(5,iat)/ev_joule
        econ(6,iat) = econ(6,iat)/ev_joule
        econ(7,iat) = econ(7,iat)/ev_joule
        bulkmod(1,iat) = 0.5*bulkmod(1,iat) + fp(iat) * 
     .    bulkmod(2,iat) + fpp(iat) * bulkmod(3,iat)**2
        bulkmod(4,iat) = 0.5*bulkmod(4,iat) + fp(iat) * 
     .    bulkmod(5,iat) + fpp(iat) * bulkmod(6,iat)**2
        bulkmod(1,iat) = bulkmod(1,iat)/ev_joule
        bulkmod(4,iat) = bulkmod(4,iat)/ev_joule
3156     continue

        if(lattype.eq.0) then
          omega0=latis*latis*latis/4.d0
        else if(lattype.eq.1) then
          omega0=latis*latis*latis/2.d0
        else if(lattype.eq.2) then
          omega0=2.*sqrt3*hcpc*latis*latis/4.
        endif
        c11(1)=0.0
        c12(1)=0.0
        c13(1)=0.0
        c33(1)=0.0
        c44(1)=0.0
        c55(1)=0.0
        c66(1)=0.0
        c66ca(1)=0.0
        bbb(1)=0.0
        bbb2(1)=0.0
        bbb3(1)=0.0
        bbb4(1)=0.0
        embed1(1)=0.0
        embed2(1)=0.0
        vacform1(1)=0.0
        vacform2(1)=0.0
        vacnum1(1)=0
        vacnum2(1)=0
        do 3205 iat=1,natoms
        bbb2(1)=bbb2(1)+(econ(1,iat)+2*econ(2,iat))/(3*omega0)
        c11(1)=c11(1)+econ(1,iat)/omega0
        c12(1)=c12(1)+econ(2,iat)/omega0
        c44(1)=c44(1)+econ(3,iat)/omega0
        c13(1)=c13(1)+econ(4,iat)/omega0
        c33(1)=c33(1)+econ(5,iat)/omega0
        c55(1)=c55(1)+econ(6,iat)/omega0
        c66(1)=c66(1)+econ(7,iat)/omega0
        c66ca(1)=c66ca(1)+(econ(2,iat)-ca(iat))/omega0
        bbb(1)=bbb(1)+bulkmod(1,iat)/(9*omega0)
        bbb3(1)=bbb3(1)+bulkmod(4,iat)/(9.*omega0)
        bbb4(1)=bbb4(1)+(2.*econ(1,iat)+2.*econ(2,iat)+
     .             4.*econ(4,iat)+econ(5,iat))/(9.*omega0)
        if( itype(iat) .eq. 1 ) then
          vacform1(1)=vacform1(1)+evacf(iat)
          embed1(1)=embed1(1)+e(iat)
          vacnum1(1)=vacnum1(1)+1
        else
          vacform2(1)=vacform2(1)+evacf(iat)
          embed2(1)=embed2(1)+e(iat)
          vacnum2(1)=vacnum2(1)+1
        endif
c       write(77,*)iat
c       write(77,*)econ(1,iat)/omega0,econ(2,iat)/omega0,
c    .econ(3,iat)/omega0
c       write(77,*)bulkmod(1,iat)/(9*omega0),bbb2
3205     continue
        bbb(1)=bbb(1)/natoms
        bbb2(1)=bbb2(1)/natoms
        bbb3(1)=bbb3(1)/natoms
        bbb4(1)=bbb4(1)/natoms
        c11(1)=c11(1)/natoms
        c12(1)=c12(1)/natoms
        c44(1)=c44(1)/natoms
        c13(1)=c13(1)/natoms
        c33(1)=c33(1)/natoms
        c55(1)=c55(1)/natoms
        c66(1)=c66(1)/natoms
        c66ca(1)=c66ca(1)/natoms
        if( vacnum1(1) .gt. 0 ) then
           vacform1(1)=vacform1(1)/vacnum1(1)
           embed1(1)=embed1(1)/vacnum1(1)
        endif
        if( vacnum2(1) .gt. 0 ) then
           vacform2(1)=vacform2(1)/vacnum2(1)
           embed2(1)=embed2(1)/vacnum2(1)
        endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         Calculation of Diatomic Strength and         c
c          Diatomic Length for Alloys                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        t1=1
        t2=1
      call calcdiatom(t1,t2,decal,recal)
      if (prntflag.eq.1) then
       write(*,*)'de:',decal,'re:',recal
      endif
      datome(1)=decal
      datomr(1)=recal
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Cu PHONON CALCULATION 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        slf1=0
        d11=0
      do 179 j=1,NNN
179      WB(j)=0.
      do 184 i2=1,NNN
      do 184 j2=1,NNN
184      dmicb(i2,j2)=DCMPLX(0.,0.)
      do 185 i2=1,3
      do 185 j2=1,3
        self1(i2,j2)=DCMPLX(0.,0.)
185    continue 
      do 186 ii=1,pnts
      do 191 i2=1,3
      do 191 j2=1,3
        dmat11(i2,j2,ii)=DCMPLX(0.,0.)
191    continue 
186    continue
c
c       Brillouin Zone:
c       Simple Cubic = 1
c       Face Centered Cubic = 2 
c       Face Centered Cubic = 3 G,X,L,K
c       brillouinzone=1
c       brillouinzone=2
        brillouinzone=3
        latak=latis
        kpi=(2.0d0*pi/latak)
        kpi1=(1.0d0*pi/latak)
        kpi2=(0.5d0*pi/latak)
        kpi4=(4.0d0*pi/latak)
        points=20
        kcount=0
        ua1=1
        ua2=2
        ua3=3
        ua4=4
        allphonon=0
        if (allphonon.eq.1) then
        do 318 kx=0,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=0.0d0
           kpoints(kcount,3)=0.0d0+(kx/DBLE(points))*kpi
318     continue
        do 319 kx=1,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=0.0d0+(kx/DBLE(points))*kpi
           kpoints(kcount,3)=1.0d0*kpi
319     continue
        do 320 kx=1,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=1.0d0*kpi-(kx/DBLE(points))*kpi
           kpoints(kcount,3)=1.0d0*kpi-(kx/DBLE(points))*kpi
320     continue
        do 321 kx=1,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0+(kx/DBLE(points))*kpi1
           kpoints(kcount,2)=0.0d0+(kx/DBLE(points))*kpi1
           kpoints(kcount,3)=0.0d0+(kx/DBLE(points))*kpi1
321     continue
        endif
c          GAMMA (000) point
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=0.0d0
           kpoints(kcount,3)=0.0d0
c          X (001) point
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=0.0d0
           kpoints(kcount,3)=kpi
c          K=(011) x 0.75 point  X=(011)=(001)
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=0.75*kpi
           kpoints(kcount,3)=0.75*kpi
c          L (111) point
           kcount=kcount+1
           kpoints(kcount,1)=kpi1
           kpoints(kcount,2)=kpi1
           kpoints(kcount,3)=kpi1
          if(prntflag.eq.1) then
           write(*,*)'  kx:      ky:      kz:'
           do 1412 kx=1,kcount
             kpnt(1)=kpoints(kx,1)
             kpnt(2)=kpoints(kx,2)
             kpnt(3)=kpoints(kx,3)
             write(*,3232)kpnt(1),kpnt(2),kpnt(3)
1412       continue
          endif
          m1=amass(itype(ua1))*conmas
          m2=amass(itype(ua2))*conmas
          m3=amass(itype(ua3))*conmas
          m4=amass(itype(ua4))*conmas
          if(prntflag.eq.1) then
           if (brillouinzone.eq.1) write(*,*)' Brillouin Zone : SC'
           if (brillouinzone.eq.2) write(*,*)' Brillouin Zone : FCC'
           write(*,*)' Lattice constant:',latak
           write(*,*)' Unit cell types for 1-4 unit cell atoms:'
           write(*,7557)IDINT(ucell(ua1)),IDINT(ucell(ua2)),
     .                  IDINT(ucell(ua3)),IDINT(ucell(ua4))
           write(*,*)' Atom numbers for 1-4 unit cell atoms:'
           write(*,7557)ua1,ua2,ua3,ua4
           write(*,*)' Mass for 1-4 unit cell atoms:'
           write(*,*)'m1=',m1
           write(*,*)'m2=',m2
           write(*,*)'m3=',m3
           write(*,*)'m4=',m4
       write(*,*)' Force Constants calculation for k-points:',kcount
          endif

        do 152 iat=1,natoms
           do 1501 j=1,neigh(iat)
           jat=rdis(iat,j,2)
c       ATOMS SELF CONTRIBUTION
          if (iat .eq. ua1 ) then
           call forcon(neigh,rdis,rdisd,j,iat,fp,fpp,dcon)
c          if( iat .eq. ua1 .or. jat .eq. ua1) then
           do 1320 kc=1,3 
           do 1320 kr=1,3 
1320       self1(kr,kc)=self1(kr,kc)+DCMPLX(dcon(kr,kc),0.)
             slf1=slf1+1
          endif
c       ATOMS SELF CONTRIBUTION END
c       ---------------------------
c        ATOM TYPE 1 CONTRIBUTIONS
          if (iat .eq. ua1) then
           call forcon(neigh,rdis,rdisd,j,iat,fp,fpp,dcon)
           do 1311 kc=1,3
1311       diski(kc) = rdisd(iat,j,kc)
          do 1401 kx=1,kcount
        kxr1=kpoints(kx,1)*diski(1)
        kxr2=kpoints(kx,2)*diski(2)
        kxr3=kpoints(kx,3)*diski(3)
        kxr=kxr1+kxr2+kxr3
           do 1322 kc=1,3 
           do 1322 kr=1,3 
                dmat11(kr,kc,kx)=dmat11(kr,kc,kx)+
     .    DCMPLX(dcon(kr,kc)*DCOS(kxr),dcon(kr,kc)*DSIN(kxr))
1322       continue
                d11=d11+1
1401      continue
          endif
c        ATOM TYPE 1 CONTRIBUTIONS END
c       ------------------------------
1501     continue
152     continue
      if (prntflag.eq.1) then
       write(*,*)' # of Self Matrix Elements:'
       write(*,*)slf1
       write(*,*)' # of Dynamic Matrix Elements:'
       write(*,*)d11/kcount
       write(*,*)' Dynamical Matrix calculation for k-points:',kcount
      endif
      do 1350 kx=1,kcount
        do 1351 a=1,3
        do 1351 b=1,3
        dmat11(a,b,kx)=(dmat11(a,b,kx)-self1(a,b))/
     .                (amass(itype(ua1))*conmas)
1351    continue
1350  continue
      if (prntflag.eq.1) then
       write(*,*)' Solving Eigenvalue Problem for k-points:',kcount
      endif
      do 1355 kx=1,kcount
        do 1366 a=1,3
        do 1366 b=1,3
1366       dmicb(a,b)=dmat11(a,b,kx)
c       write(*,*)dmic
c       Query Optimal Workspace (LWORK=-1)
        LWORKB=-1
        CALL ZHEEV('Vector','Lower',NNN,dmicb,LDB,WB,WORKB,LWORKB,
     .  RWORKB,INFOB)
c       Solve Problem with LWORK workspace
        LWORKB = MIN( LWMAXB, INT( WORKB( 1 ) ) )
        CALL ZHEEV('Vector','Lower',NNN,dmicb,LDB,WB,WORKB,LWORKB,
     .  RWORKB,INFOB)
        kpnt(1)=kpoints(kx,1)
        kpnt(2)=kpoints(kx,2)
        kpnt(3)=kpoints(kx,3)
c       Check Convergence
        if (INFOB.eq.0) then
        do 1466 j=1,NNN
        WB(j)=DABS(DSQRT(WB(j)))/(2.*pi)
c       if eigenvalue is NaN then append it to 0.
1466    continue
        do 1368 j=1,NNN
1368    if (WB(j) .ne. WB(j)) WB(j)=0.
        if (allphonon.eq.1) then
         if (kx.eq.82) then
             phg=WB(3)
         else if (kx.eq.83) then
             phxl=WB(3)
             phxt=WB(2)
         else if (kx.eq.84) then
             phkl=WB(3)
             phkt1=WB(2)
             phkt2=WB(1)
         else if (kx.eq.85) then
             phll=WB(3)
             phlt=WB(2)
         endif
        else
         if (kx.lt.2) then
             phg=WB(3)
         else if (kx.lt.3) then
             phxl=WB(3)
             phxt=WB(2)
         else if (kx.lt.4) then
             phkl=WB(3)
             phkt1=WB(2)
             phkt2=WB(1)
         else if (kx.lt.5) then
             phll=WB(3)
             phlt=WB(2)
         endif
       endif
       if (prntflag.eq.1) then
          write(*,5555)kpnt(1),kpnt(2),kpnt(3),' : ',(WB(j),j=1,NNN)
       endif
c         write(77,5775)kpnt(1),kpnt(2),kpnt(3),(W(j),j=1,NN)
        else if (INFOB .gt. 0) then
       if (prntflag.eq.1) then
          write(*,*)'Eigenvalue ',INFOB,' is not converged for k=',
     .               kpnt(1),kpnt(2),kpnt(3)
       endif
        else
       if (prntflag.eq.1) then
          write(*,*)INFOB,' th argument had an illegal value for k=',
     .               kpnt(1),kpnt(2),kpnt(3)
       endif
        endif
1355  continue
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     END OF PHONON CALCULATION 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c*****************************************************
c*          Calculation PURE Stacking Fault and      *
c*           Antiphase Boundary Energies             *
c*            (100), (111) SISF and APB              *
c*****************************************************
       l1=0
       t1=1
       t2=1
       stacks(1)=10
       stacks(2)=10
       stacks(3)=18
       call loadfcc111(t1,t2,l1,latis,stacks)
       vmdfile='fcc111-00.lammpstrj'
       call writevmd(vmdfile)
       saveper(1,3)=perlb(3)
       saveper(2,3)=perub(3)
       perlb(3)=-10000
       perub(3)=10000
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms9:',rmssum
        endif
         return
       endif
       call calcper()
       call calctote(etot1,neigh,rdis)
       shifts(1)=2.0*(latis*sqrt3/(2.0*sqrt2))/3.0
       shifts(2)=0.
       shifts(3)=0.
       zplanes=stacks(3)/2
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms10:',rmssum
        endif
         return
       endif
       call calcper()
       call calctote(etot2,neigh,rdis)
       sisf111(1)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       sisf111(1)=sisf111(1)*mev_to_mjoule*10
       if (prntflag.eq.1) then
        print*,'etot1:',etot1,'etot2:',etot2
        print*,'SISF(111):',sisf111(1)
       endif
       if (prntflag.eq.1) then
       vmdfile='fcc111-2.lammpstrj'
       stacks(1)=10
       stacks(2)=10
       stacks(3)=18
       call writevmd(vmdfile)
       call loadfcc111(t1,t2,l1,latis,stacks)
       vmdfile='fcc111-0.lammpstrj'
       call writevmd(vmdfile)
       saveper(1,3)=perlb(3)
       saveper(2,3)=perub(3)
       perlb(3)=-10000
       perub(3)=10000
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
        if(istop .eq. 1) then
         if (prntflag.eq.1) then
          write(*,*)'rms4-2:',rmssum
         endif
          return
        endif
       call calcper()
       call calctote(etot1,neigh,rdis)
       zplns(1)=(stacks(3)/2)
       zplns(2)=stacks(3)-(stacks(3)/2)+1
       call loadtwin111(t1,t2,l1fcc,latis,stacks,zplns)
       saveper(1,3)=perlb(3)
       saveper(2,3)=perub(3)
       perlb(3)=-10000
       perub(3)=10000
       shifts(1)=-2.0*(latis*sqrt3/(2.0*sqrt2))/3.0
       shifts(2)=0.
       shifts(3)=0.
       zplanes=(stacks(3)/2)+1
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms10:',rmssum
        endif
         return
       endif
       call calctote(etot2,neigh,rdis)
       vmdfile='twin111.lammpstrj'
       call writevmd(vmdfile)
       twin111(1)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       twin111(1)=twin111(1)*0.5
       twin111(1)=twin111(1)*mev_to_mjoule*10
       if (prntflag.eq.1) then
        print*,'etot1:',etot1,'etot2:',etot2
        print*,'TWIN(111):',twin111(1)
       endif
       endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         Calculation of Vacancy Migration Energy      c
c                    on Static Path                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       stacks(1)=10
       stacks(2)=10
       stacks(3)=12
       l1fcc=0
       t1=1
       t2=1
       migatoms(1)=55
       migatoms(2)=55
       call loadfcc100(t1,t2,l1fcc,latis,stacks)
       call removeatoms(migatoms)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
        if(istop .eq. 1) then
         if (prntflag.eq.1) then
          write(*,*)'rms14-2:',rmssum
         endif
          return
        endif
       if (prntflag.eq.1) then
        vmdfile='vacmig1.lammpstrj'
        call writevmd(vmdfile)
       endif
       call calcper()
       call calctote(etot1,neigh,rdis)
       shifts(1)=latis/4.
       shifts(2)=latis/4.
       shifts(3)=0.0
       migatoms(1)=45
       migatoms(2)=45
       call shiftatoms(migatoms,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
        if(istop .eq. 1) then
         if (prntflag.eq.1) then
          write(*,*)'rms14-3:',rmssum
         endif
          return
        endif
       call calcper()
       call calctote(etot2,neigh,rdis)
       vacmig1(1)=dabs(etot2-etot1)
       if (prntflag.eq.1) then
        print*,'etot1:',etot1,'etot2:',etot2
        print*,'Emig:',vacmig1(1)
        vmdfile='vacmig2.lammpstrj'
        call writevmd(vmdfile)
       endif
         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         fitnums=0      
         if (weightdata(1).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(1)*weightdata(1)
         fitvals(fitnums)=l12alat(1)*weightdata(1)
         endif
         if (weightdata(2).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(2)*weightdata(2)
         fitvals(fitnums)=-l12ecoh(1)*weightdata(2)
         endif
         if (weightdata(3).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(3)*weightdata(3)
         fitvals(fitnums)=-embed1(1)*weightdata(3)
         endif
         if (weightdata(4).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(4)*weightdata(4)
         fitvals(fitnums)=bbb(1)*weightdata(4)
         endif
         if (weightdata(5).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(5)*weightdata(5)
         fitvals(fitnums)=c11(1)*weightdata(5)
         endif
         if (weightdata(6).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(6)*weightdata(6)
         fitvals(fitnums)=c12(1)*weightdata(6)
         endif
         if (weightdata(7).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(7)*weightdata(7)
         fitvals(fitnums)=c44(1)*weightdata(7)
         endif
         if (weightdata(8).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(8)*weightdata(8)
         fitvals(fitnums)=-vacform1(1)*weightdata(8)
         endif
         if (weightdata(9).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(9)*weightdata(9)
         fitvals(fitnums)=-datome(1)*weightdata(9)
         endif
         if (weightdata(10).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(10)*weightdata(10)
         fitvals(fitnums)=datomr(1)*weightdata(10)
         endif
         if (weightdata(11).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(11)*weightdata(11)
         fitvals(fitnums)=phxt*weightdata(11)
         endif
         if (weightdata(12).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(12)*weightdata(12)
         fitvals(fitnums)=phxl*weightdata(12)
         endif
         if (weightdata(13).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(13)*weightdata(13)
         fitvals(fitnums)=vacmig1(1)*weightdata(13)
         endif
         if (weightdata(14).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(14)*weightdata(14)
         fitvals(fitnums)=sisf111(1)*weightdata(14)
         endif
c     ----------------------------------------------
       rmssum=0.0
       if (prntflag.eq.1) then
        write(*,*)'fitnmus:',fitnums
       endif
       do 4005 i=1,fitnums
        if (fitvals(i).eq.0.0) then
            fitvals(i)=1000000
        endif
        if (fitvals(i).lt.0.01) then
            fitvals(i)=fitvals(i)*1000000
        endif
c       errors(i)=abs(expvals(i)-fitvals(i))
        errors(i)=expvals(i)-fitvals(i)
        errors(i)=errors(i)*errors(i)
        rmssum=rmssum+errors(i)
        if (datomr(1).lt.-datome(1)) then
            rmssum=rmssum+25.
        endif
c       write(*,*)expvals(i),fitvals(i),errors(i)
4005  continue
      rmssum=sqrt(abs(rmssum))
c     ----------------------------------------------
       if (prntflag.eq.1) then
       open(unit=70,file=fitvalues)
      write(*,*)'--------------------'
      write(70,7871)l12alat(1),hexca(1),-l12ecoh(1),-embed1(1),
     .-embed2(1),bbb(1),bbb3(1),c11(1),c12(1),c13(1),c33(1),
     .c44(1),c66(1),-datome(1),datomr(1),-vacform1(1),-vacform2(1),
     .sisf111(1),apb100(1),apb111(1)
      write(*,7770)l12alat(1),hexca(1),-l12ecoh(1),-embed1(1),
     .-embed2(1),bbb(1),c11(1),c12(1),c44(1),c13(1),c33(1),c55(1),
     .c66(1),-vacform1(1),-vacform2(1),bbb3(1),-datome(1),datomr(1)
      write(*,7568)sisf111(1),apb100(1),apb111(1),vacmig1(1),
     .twin111(1),sf111(1)
      write(*,7578)phg,phxt,phxl,phlt,phll,phkt1,phkt2,phkl
7770    format((' PURE',/,' a_0:',f9.5,/,' c/a:',f9.5,/,
     .' Ecoh:',f9.5,/,' Ecoh1:',f9.5,/,' Ecoh2:',f9.5,/,
     .' bmod:',f9.5,/,' c11:',f9.5,/,' c12:',f9.5,/,' c44:',f9.5,/,
     .' c13:',f9.5,/,' c33:',f9.5,/,' c55:',f9.5,/,
     .' c66:',f9.5,/,' Evf1:',f9.5,/,' Evf2:',f9.5,/,' Bmod:',f9.5,/,
     .' De:',f9.5,/,' Re',f9.5,/))
7871    format(20(E25.16,1x))
7568    format((' SISF(111):',f11.6,/,' APB(100):',f11.6,/,
     .' APB(111):',f11.6,/,' Evm1:',f11.6,/,
     .' Twin(111):',f11.6,/,' SF(111):',f11.6,/))
7578    format((' v(G):',f10.5,/,' vT(X):',f10.5,/,
     .          ' vL(X):',f10.5,/,' vT(L):',f10.5,/,
     .          ' vL(L):',f10.5,/,' vT1(K):',f10.5,/,
     .          ' vT2(K):',f10.5,/,' vL(K):',f10.5,/))
c      write(*,*)'--------------------'
c      write(*,*)'RMS:',rmssum
        close(70)
       endif
c     ---------------------------------------------               
c     PURE CALCULATION ENDS HERE
      else
c     ALLOY CALCULATION GOES HERE AFTER
c     ---------------------------------------------               

c      ---------------------------------------------               
c      FCC B2 NiCu - CuNi 
       if( ntypes .lt. 3 .and. ( at1 .eq. 2 .and. at2 .eq. 4 ) ) then
c****************************************************
c*   Calculation of Equilibrium Lattice Constant    *
c*           for L1_0 and L1_2 FCC Alloy            *
c****************************************************
c     
c     L1_0 Cu3Ni Calculations
c
      if (prntflag.eq.1) then
       write(*,*)'------- Cu3Ni FCC L1_0 -------'
      endif
      if (prntflag.eq.1) then
        write(*,*)'PROC:',procnum,' part:',prmid,' Cu3Ni FCC L1_0 '
      endif
      l1fcc=1
      if( ntypes .eq. 2 ) then
       t1=1
       t2=2
      else
       t1=2
       t2=4
      endif
      call loadfcc100(t1,t2,l1fcc,lata0,stacks)
      call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms1-Cu3Ni-L10:',rmssum
        endif
         return
       endif
      if (prntflag.eq.1) then
       write(*,*)'Test a_0:',lata0
      endif
      call calclatcon(minenergy,latis,lata0,neigh,rdis,iter,prntflag)
      l10alat(1)=latis
       if (prntflag.eq.1) then
        write(*,*)'L1_0 FCC Cu3Ni a_0:',latis
       endif
      call loadfcc100(t1,t2,l1fcc,latis,stacks)
      call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms2-Cu3Ni-L10:',rmssum
        endif
         return
       endif
       tn=2
       typs(1)=t1
       typs(2)=t2
       call calctotetype(tn,typs,tnum,embed,etot1,neigh,rdis)
       l10ecoh(2)=etot1
       etot1=etot1/natoms
       l10ecoh(1)=etot1
       embed1(3)=embed(1)
       embed2(3)=embed(2)
       if (prntflag.eq.1) then
        write(*,*)'L1_0 FCC Cu3Ni Ecoh:',etot1
        write(*,*)'L1_0 FCC Cu3Ni Embed1:',embed(1),'num:',tnum(1)
        write(*,*)'L1_0 FCC Cu3Ni Embed2:',embed(2),'num:',tnum(2)
       endif
c     
c     L1_2 Cu3Ni Calculations
c
      if (prntflag.eq.1) then
       write(*,*)'------- Cu3Ni FCC L1_2 -------'
      endif
      if (prntflag.eq.1) then
        write(*,*)'PROC:',procnum,' part:',prmid,' Cu3Ni FCC L1_2 '
      endif
      l1fcc=3
      if( ntypes .eq. 2 ) then
       t1=1
       t2=2
      else
       t1=2
       t2=4
      endif
      call loadfcc100(t1,t2,l1fcc,lata0,stacks)
      call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms4:',rmssum
        endif
         return
       endif
      if (prntflag.eq.1) then
       write(*,*)'Test a_0:',lata0
      endif
      call calclatcon(minenergy,latis,lata0,neigh,rdis,iter,prntflag)
      l12alat(1)=latis
      if (prntflag.eq.1) then
       write(*,*)'L1_2 FCC Cu3Ni a_0:',latis
       write(*,*)'L1_2 FCC Cu3Ni Ecoh:',minenergy/natoms
      endif
      l12ecoh(1)=minenergy/natoms
      
      l12ecoh(2)=l10ecoh(2)-minenergy
      if (prntflag.eq.1) then
      write(*,*)'CuNi L10-L12 Total eV:',l12ecoh(2)
      endif
c     l1fcc=3
c     t1=2
c     t2=1
c     call loadfcc111(t1,t2,l1fcc,lata0,stacks)
c     call loadneigh(rcutoff,neigh,rdis,rdisd)
c     write(*,*)'Test a_0:',lata0
c     call calclatcon(minenergy,latis,lata0,neigh,rdis)
c     l12alat2=latis
c     write(*,*)'L12 111 FCC a_0:',latis
c     call loadfcc111(t1,t2,l1fcc,latis,stacks)
c     call loadneigh(rcutoff,neigh,rdis,rdisd)
c     call calctote(etot1,neigh,rdis)
c     etot1=etot1/natoms
c     l12ecoh2=etot1
c     write(*,*)'L12 111 FCC Ecoh:',etot1
c     vmdfile='cuni-l12-111.lammpstrj'
c     call writevmd(vmdfile)
c     l1fcc=4
c     call loadfcc111(t1,t2,l1fcc,lata0,stacks)
c     call loadneigh(rcutoff,neigh,rdis,rdisd)
c     write(*,*)'Test a_0:',lata0
c     call calclatcon(minenergy,latis,lata0,neigh,rdis)
c     l13alat2=latis
c     write(*,*)'L13 111 FCC a_0:',latis
c     call loadfcc111(t1,t2,l1fcc,latis,stacks)
c     call loadneigh(rcutoff,neigh,rdis,rdisd)
c     call calctote(etot1,neigh,rdis)
c     etot1=etot1/natoms
c     l13ecoh2=etot1
c     write(*,*)'L13 111 FCC Ecoh:',etot1
c     vmdfile='cuni-l13-111.lammpstrj'
c     call writevmd(vmdfile)
c     l1fcc=2
c     call loadfcc111(t1,t2,l1fcc,lata0,stacks)
c     call loadneigh(rcutoff,neigh,rdis,rdisd)
c     write(*,*)'Test a_0:',lata0
c     call calclatcon(minenergy,latis,lata0,neigh,rdis)
c     l11alat2=latis
c     write(*,*)'L11 111 FCC a_0:',latis
c     call loadfcc111(t1,t2,l1fcc,latis,stacks)
c     call loadneigh(rcutoff,neigh,rdis,rdisd)
c     call calctote(etot1,neigh,rdis)
c     etot1=etot1/natoms
c     l11ecoh2=etot1
c     write(*,*)'L11 111 FCC Ecoh:',etot1
c     vmdfile='cuni-l11-111.lammpstrj'
c     call writevmd(vmdfile)
c     l1fcc=1
c     call loadfcc111(t1,t2,l1fcc,lata0,stacks)
c     call loadneigh(rcutoff,neigh,rdis,rdisd)
c     write(*,*)'Test a_0:',lata0
c     call calclatcon(minenergy,latis,lata0,neigh,rdis)
c     l10alat2=latis
c     write(*,*)'L10 111 FCC a_0:',latis
c     call loadfcc111(t1,t2,l1fcc,latis,stacks)
c     call loadneigh(rcutoff,neigh,rdis,rdisd)
c     call calctote(etot1,neigh,rdis)
c     etot1=etot1/natoms
c     l10ecoh2=etot1
c     write(*,*)'L10 111 FCC Ecoh:',etot1
c     vmdfile='cuni-l10-111.lammpstrj'
c     call writevmd(vmdfile)
c     call loadfcc100(t1,t2,l1fcc,lata0,stacks)
c     call loadneigh(rcutoff,neigh,rdis,rdisd)
c     write(*,*)'Test a_0:',lata0
c     call calclatcon(minenergy,latis,lata0,neigh,rdis)
c     lata0=latis
c     write(*,*)'100 L12 FCC a_0:',lata0
c     write(*,*)'100 L12 FCC Ecoh:',minenergy/natoms
c     stop

c*****************************************************
      l1fcc=3
      if( ntypes .eq. 2 ) then
       t1=1
       t2=2
      else
       t1=2
       t2=4
      endif
      call loadfcc100(t1,t2,l1fcc,latis,stacks)
      call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms4-L1_2:',rmssum
        endif
         return
       endif

      omegafcc=lata0*lata0*lata0/4.d0
      omegabcc=lata0*lata0*lata0/2.d0
      omegahcp=2.*sqrt3*hcpc*lata0*lata0/4.d0
c     angstrom^3/atom -> cm^3/mol 
c     when multiplied with bmod, bmod will be eV
      omega=latis*latis*latis/4.d0
      omega=omega*0.62415
c


c   *****************************************************
c   * EAM fcn file ready to use in Sublimation Energy,  *
c   *     Vacancy Formation Energy, Bulk Modulus,       *
c   *     and Elastic Constant Calculations for FCC     *
c   *****************************************************
    
      do 1451 iat=1,natoms
        e(iat)=0.0
        fe(iat)=0.0
        fp(iat)=0.0
        fpp(iat)=0.0
        phie(iat)=0.0
        rho(iat)=0.0
        evacf(iat)=0.0
1451  continue
      
      rdrar = 1.0/drar 
      do 1481 iat=1,natoms
c
c     compute the contribution to rho(i) from particle j
c     Calculating the embedding energy F(Rho) in fe(i)
c      and Pair Energy contribution Phi(Rij) in phie(i)
c
ccdir$ novector
      do 1461 j=1,neigh(iat)
      rdist=rdis(iat,j,1)
      jat=rdis(iat,j,2)
      pij = rdist*rdrar + 1.0
      kij = pij
      kij = min0(kij,nrar-1)
      pij = pij - kij
c       make sure that p is less than 1.00
c       then if r is out of range, p = 1.00 and rho = last value of rhor
      pij = amin1(pij,1.0)
      ity=itype(iat)
      jty=itype(jat)
      rhojnei(iat,j) = ((rhorar3(kij,jty)*pij+
     $   rhorar2(kij,jty))*pij+
     $   rhorar1(kij,jty))*pij+
     $   rhorar(kij,jty)
      rhoinei(iat,j) = ((rhorar3(kij,ity)*pij+
     $   rhorar2(kij,ity))*pij+
     $   rhorar1(kij,ity))*pij+
     $   rhorar(kij,ity)
      rho(iat) = rho(iat) + rhojnei(iat,j)
      z2ij = ((z2rar3(kij,ity,jty)*pij+
     $          z2rar2(kij,ity,jty))*pij+
     $          z2rar1(kij,ity,jty))*pij+
     $          z2rar(kij,ity,jty)
      z2pij = (z2rar6(kij,ity,jty)*pij+z2rar5(kij,ity,jty))*pij+
     $         z2rar4(kij,ity,jty)
c     if(phitype) then
        phiij = z2ij
c     else
c       pij=1.0/rdist
c       phiij = z2ij * pij
c     endif
      phie(iat) = phie(iat) + 0.5 * phiij
1461  continue
ccdir$ vector
c        
1481  continue
c
c       store fsubi in e(i), the derivative of fsubi in fp
c       and the second derivative of fsubi in fpp
c
      rdrhoar = 1.0/drhoar
      do 2412 iat=1,natoms
      pij = rho(iat)*rdrhoar + 1.0
      kij = pij
      kij = max(1,min(kij,nrhoar-1))
      pij = pij - kij
      ity=itype(iat)
      fe(iat) = ((frhoar3(kij,ity)*pij+
     $        frhoar2(kij,ity))*pij+
     $        frhoar1(kij,ity))*pij+
     $        frhoar(kij,ity)
      fp(iat) = (frhoar6(kij,ity)*pij+
     $                  frhoar5(kij,ity))*pij+
     $                  frhoar4(kij,ity)
      fpp(iat) = (2.0*frhoar6(kij,ity)*pij+
     $                  frhoar5(kij,ity))/drhoar
2412  continue
c
c       store vacancy formation energy in evacf(i)
c
      rdrhoar = 1.0/drhoar
      do 3401 iat=1,natoms
        do 3411 j=1,neigh(iat)
        jat=rdis(iat,j,2)
        pij = (rho(jat)-rhoinei(iat,j))*rdrhoar + 1.0
        kij = pij
        kij = max(1,min(kij,nrhoar-1))
        pij = pij - kij
        ity=itype(iat)
        jty=itype(jat)
        evacf(iat) = evacf(iat) + ( fe(jat) - 
     $      (((frhoar3(kij,jty)*pij+
     $        frhoar2(kij,jty))*pij+
     $        frhoar1(kij,jty))*pij+
     $        frhoar(kij,jty))  )
3411    continue
      e(iat) = phie(iat) + fe(iat)
3401  continue

      do 3421 iat=1,natoms
      evacf(iat) = evacf(iat) + phie(iat)
3421  continue
      

      do 3453 iat=1,natoms
        ca(iat) = 0.0
        do 3454 cij=1,7
         econ(cij,iat) = 0.0
         bulkmod(cij,iat) = 0.0
         uu(cij) = 0.0
         ww(cij) = 0.0
         vv(cij) = 0.0
3454    continue
         do 3402 j=1,neigh(iat)
           jat=rdis(iat,j,2)
           do 3412 kc=1,3
3412       diski(kc) = rdisd(iat,j,kc)
c          call elascon(phitype,diski,jat,iat,uu,ww,vv,bulkmod) 
           call elastic(phitype,diski,jat,iat,uu,ww,vv,bulkmod) 
3402     continue
        econ(1,iat) = 0.5*uu(1) + fp(iat) *
     .                   ww(1) + fpp(iat) * vv(1)*vv(1)
        econ(2,iat) = 0.5*uu(2) + fp(iat) *
     .                   ww(2) + fpp(iat) * vv(1)*vv(2)
        econ(3,iat) = 0.5*uu(3) + fp(iat) *
     .                   ww(3) + fpp(iat) * vv(3)*vv(3)
        econ(4,iat) = 0.5*uu(4) + fp(iat) *
     .                   ww(4) + fpp(iat) * vv(1)*vv(4)
        econ(5,iat) = 0.5*uu(5) + fp(iat) *
     .                   ww(5) + fpp(iat) * vv(4)*vv(4)
        econ(6,iat) = 0.5*uu(6) + fp(iat) *
     .                   ww(6) + fpp(iat) * vv(5)*vv(6)
        econ(7,iat) = 0.5*uu(7) + fp(iat) *
     .                   ww(7) + fpp(iat) * vv(7)*vv(7)
        ca(iat)= fpp(iat)*(vv(1)*vv(2)-vv(7)*vv(7))
        ca(iat)= ca(iat)/ev_joule
        econ(1,iat) = econ(1,iat)/ev_joule
        econ(2,iat) = econ(2,iat)/ev_joule
        econ(3,iat) = econ(3,iat)/ev_joule
        econ(4,iat) = econ(4,iat)/ev_joule
        econ(5,iat) = econ(5,iat)/ev_joule
        econ(6,iat) = econ(6,iat)/ev_joule
        econ(7,iat) = econ(7,iat)/ev_joule
        bulkmod(1,iat) = 0.5*bulkmod(1,iat) + fp(iat) * 
     .    bulkmod(2,iat) + fpp(iat) * bulkmod(3,iat)**2
        bulkmod(4,iat) = 0.5*bulkmod(4,iat) + fp(iat) * 
     .    bulkmod(5,iat) + fpp(iat) * bulkmod(6,iat)**2
        bulkmod(1,iat) = bulkmod(1,iat)/ev_joule
        bulkmod(4,iat) = bulkmod(4,iat)/ev_joule
3453     continue

        omega0fcc = lata0*lata0*lata0/4.d0
        omega0bcc = lata0*lata0*lata0/2.d0
        omega0hcp = 2.*sqrt3*hcpc*lata0*lata0/4.d0
        omega0 = latis*latis*latis/4.d0
        c11(1)=0.0
        c12(1)=0.0
        c13(1)=0.0
        c33(1)=0.0
        c44(1)=0.0
        c55(1)=0.0
        c66(1)=0.0
        c66ca(1)=0.0
        bbb(1)=0.0
        bbb2(1)=0.0
        bbb3(1)=0.0
        bbb4(1)=0.0
        embed1(1)=0.0
        embed2(1)=0.0
        vacform1(1)=0.0
        vacform2(1)=0.0
        vacnum1(1)=0
        vacnum2(1)=0
        do 3221 iat=1,natoms
        bbb2(1)=bbb2(1)+(econ(1,iat)+2*econ(2,iat))/(3*omega0)
        c11(1)=c11(1)+econ(1,iat)/omega0
        c12(1)=c12(1)+econ(2,iat)/omega0
        c44(1)=c44(1)+econ(3,iat)/omega0
        c13(1)=c13(1)+econ(4,iat)/omega0
        c33(1)=c33(1)+econ(5,iat)/omega0
        c55(1)=c55(1)+econ(6,iat)/omega0
        c66(1)=c66(1)+econ(7,iat)/omega0
        c66ca(1)=c66ca(1)+(econ(2,iat)-ca(iat))/omega0
        bbb(1)=bbb(1)+bulkmod(1,iat)/(9*omega0)
        bbb3(1)=bbb3(1)+bulkmod(4,iat)/(9.*omega0)
        bbb4(1)=bbb4(1)+(2.*econ(1,iat)+2.*econ(2,iat)+
     .             4.*econ(4,iat)+econ(5,iat))/(9.*omega0)
        if( itype(iat) .eq. 1 ) then
          vacform1(1)=vacform1(1)+evacf(iat)
          embed1(1)=embed1(1)+e(iat)
          vacnum1(1)=vacnum1(1)+1
        else
          vacform2(1)=vacform2(1)+evacf(iat)
          embed2(1)=embed2(1)+e(iat)
          vacnum2(1)=vacnum2(1)+1
        endif
c       write(77,*)iat
c       write(77,*)econ(1,iat)/omega0,econ(2,iat)/omega0,
c    .econ(3,iat)/omega0
c       write(77,*)bulkmod(1,iat)/(9*omega0),bbb2
3221     continue
        bbb(1)=bbb(1)/natoms
        bbb2(1)=bbb2(1)/natoms
        bbb3(1)=bbb3(1)/natoms
        bbb4(1)=bbb4(1)/natoms
        c11(1)=c11(1)/natoms
        c12(1)=c12(1)/natoms
        c44(1)=c44(1)/natoms
        c13(1)=c13(1)/natoms
        c33(1)=c33(1)/natoms
        c55(1)=c55(1)/natoms
        c66(1)=c66(1)/natoms
        c66ca(1)=c66ca(1)/natoms
        if( vacnum1(1) .gt. 0 ) then
           vacform1(1)=vacform1(1)/vacnum1(1)
           embed1(1)=embed1(1)/vacnum1(1)
        endif
        if( vacnum2(1) .gt. 0 ) then
           vacform2(1)=vacform2(1)/vacnum2(1)
           embed2(1)=embed2(1)/vacnum2(1)
        endif
      if (prntflag.eq.1) then
        write(*,*)'L1_2 FCC Cu3Ni Embed1:',embed1(1),'num:',vacnum1(1)
        write(*,*)'L1_2 FCC Cu3Ni Embed2:',embed2(1),'num:',vacnum2(1)
      endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         Calculation of Diatomic Strength and         c
c          Diatomic Length for Alloys                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if( ntypes .eq. 2 ) then
        t1=1
        t2=2
      else
        t1=1
        t2=2
      endif
      call calcdiatom(t1,t2,decal,recal)
      if (prntflag.eq.1) then
       write(*,*)'de:',decal,'re:',recal
      endif
      datome(1)=decal
      datomr(1)=recal
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       if (prntflag.eq.1) then
c*****************************************************
c*   Calculation Cu3Ni Alloy Stacking Fault and      *
c*           Antiphase Boundary Energies             *
c*            (100), (111) SISF and APB              *
c*****************************************************
       call calcper()
       call calctote(etot1,neigh,rdis)
       shifts(1)=abs(rv(1,stacks(2)+1)-rv(1,1))
       shifts(2)=abs(rv(2,stacks(2)+1)-rv(2,1))
       shifts(3)=0.0
       zplanes=stacks(3)/2
       if (prntflag.eq.1) then
        print*,'Shifting',zplanes,'planes of FCC100 to x:',
     .shifts(1),'and y:',shifts(2)
        print*,'area:',perlen(1)*perlen(2)
       endif
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms6:',rmssum
        endif
         return
       endif
       call calctote(etot2,neigh,rdis)
c      vmdfile='ni3al-100.lammpstrj'
c      call writevmd(vmdfile)
       apb100(1)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       apb100(1)=apb100(1)/2.0
       apb100(1)=10000*apb100(1)/ev_joule
       if (prntflag.eq.1) then
        print*,'etot1:',etot1,'etot2:',etot2
        print*,'APB(100):',apb100(1)
       endif
       l1=3
       if( ntypes .eq. 2 ) then
         t1=1
         t2=2
       else
         t1=1
         t2=2
       endif
       call loadfcc111(t1,t2,l1,latis,stacks)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms7:',rmssum
        endif
         return
       endif
c      call calcper()
       call calctote(etot1,neigh,rdis)
c      print*,'L1 Ecoh for 100 stacking:',etot1/natoms
c      vmdfile='Ni3Al-L12-111.lammpstrj'
c      call writevmd(vmdfile)
c      vmdfile='ni3al.lammpstrj'
c      call writevmd(vmdfile)
       shifts(1)=0.0
       shifts(2)=abs(rv(2,1)-rv(2,2))
       shifts(3)=0.0
       zplanes=4
       if (prntflag.eq.1) then
        print*,'Shifting',zplanes,'planes of FCC111 to x:',
     .shifts(1),'and y:',shifts(2)
        print*,'area:',perlen(1)*perlen(2)
       endif
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms8:',rmssum
        endif
         return
       endif
       call calctote(etot2,neigh,rdis)
       apb111(1)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       apb111(1)=apb111(1)/2.0
       apb111(1)=10000*apb111(1)/ev_joule
       if (prntflag.eq.1) then
        print*,'etot1:',etot1,'etot2:',etot2
        print*,'APB(111):',apb111(1)
       endif
c      vmdfile='ni3al-apb.lammpstrj'
c      call writevmd(vmdfile)
       stacks(3)=stacks(3)-1
       call loadfcc111(t1,t2,l1,latis,stacks)
c      call calcper()
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms9:',rmssum
        endif
         return
       endif
       call calctote(etot1,neigh,rdis)
       xshift = rv(1,1) - rv(1,stacks(2)+2)
       yshift = rv(2,1) - rv(2,stacks(2)+2)
       shifts(1)=2.0*abs(xshift)/3.0
       shifts(2)=2.0*abs(yshift)/3.0
       shifts(3)=0.0
        zplanes=4
       if (prntflag.eq.1) then
        print*,'Shifting',zplanes,'planes of FCC111 to x:',
     .shifts(1),'and y:',shifts(2)
        print*,'area:',perlen(1)*perlen(2)
       endif
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms10:',rmssum
        endif
         return
       endif
       call calctote(etot2,neigh,rdis)
       sisf111(1)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       sisf111(1)=sisf111(1)/2.0
       sisf111(1)=10000*sisf111(1)/ev_joule
       if (prntflag.eq.1) then
        print*,'etot1:',etot1,'etot2:',etot2
        print*,'SISF(111):',sisf111(1)
       endif
       stacks(3)=stacks(3)+1
c      vmdfile='ni3al-sisf.lammpstrj'
c      call writevmd(vmdfile)
       endif

c****************************************************
c*   Calculation of Equilibrium Lattice Constant    *
c*           for L1_3 and L1_1 FCC Alloy            *
c****************************************************
c     
c     L1_3 Cu3Ni Calculations
c
      if (prntflag.eq.1) then
       write(*,*)'------- Cu3Ni FCC L1_3 -------'
      endif
      if (prntflag.eq.1) then
        write(*,*)'PROC:',procnum,' part:',prmid,' Cu3Ni FCC L1_3 '
      endif
      l1fcc=4
      if( ntypes .eq. 2 ) then
       t1=1
       t2=2
      else
       t1=2
       t2=4
      endif
      call loadfcc100(t1,t2,l1fcc,lata0,stacks)
      call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms1-Cu3Ni-L13:',rmssum
        endif
         return
       endif
      if (prntflag.eq.1) then
       write(*,*)'Test a_0:',lata0
      endif
      call calclatcon(minenergy,latis,lata0,neigh,rdis,iter,prntflag)
      l13alat(1)=latis
       if (prntflag.eq.1) then
        write(*,*)'L1_3 FCC Cu3Ni a_0:',latis
       endif
      call loadfcc100(t1,t2,l1fcc,latis,stacks)
      call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms2-Cu3Ni-L13:',rmssum
        endif
         return
       endif
       tn=2
       typs(1)=t1
       typs(2)=t2
       call calctotetype(tn,typs,tnum,embed,etot1,neigh,rdis)
       l13ecoh(2)=etot1
       etot1=etot1/natoms
       l13ecoh(1)=etot1
       embed1(4)=embed(1)
       embed2(4)=embed(2)
       if (prntflag.eq.1) then
        write(*,*)'L1_3 FCC Cu3Ni Ecoh:',etot1
        write(*,*)'L1_3 FCC Cu3Ni Embed1:',embed(1),'num:',tnum(1)
        write(*,*)'L1_3 FCC Cu3Ni Embed2:',embed(2),'num:',tnum(2)
       endif
c     
c     L1_1 Cu3Ni Calculations
c
      if (prntflag.eq.1) then
       write(*,*)'------- Cu3Ni FCC L1_1 -------'
      endif
      if (prntflag.eq.1) then
        write(*,*)'PROC:',procnum,' part:',prmid,' Cu3Ni FCC L1_1 '
      endif
      l1fcc=2
      if( ntypes .eq. 2 ) then
       t1=1
       t2=2
      else
       t1=2
       t2=4
      endif
      call loadfcc100(t1,t2,l1fcc,lata0,stacks)
      call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms4:',rmssum
        endif
         return
       endif
      if (prntflag.eq.1) then
       write(*,*)'Test a_0:',lata0
      endif
      call calclatcon(minenergy,latis,lata0,neigh,rdis,iter,prntflag)
      l11alat(1)=latis
      if (prntflag.eq.1) then
       write(*,*)'L1_1 FCC Cu3Ni a_0:',latis
       write(*,*)'L1_1 FCC Cu3Ni Ecoh:',minenergy/natoms
      endif
      l11ecoh(1)=minenergy/natoms
      
      l11ecoh(2)=l13ecoh(2)-minenergy
      if (prntflag.eq.1) then
      write(*,*)'CuNi L13-L11 Total eV:',l11ecoh(2)
      endif

      l1fcc=2
      if( ntypes .eq. 2 ) then
       t1=1
       t2=2
      else
       t1=2
       t2=4
      endif
      stck(1)=stacks(1)
      stck(2)=stacks(2)
      stck(3)=stacks(3)
      stacks(1)=4
      stacks(2)=4
      stacks(3)=4
      patm=32
      pvec(1)=2.*latis
      pvec(2)=2.*latis
      pvec(3)=2.*latis
      ptyp(1)=1
      ptyp(2)=2
      ptyp(3)=2
      ptyp(4)=1
      ptyp(5)=2
      ptyp(6)=1
      ptyp(7)=1
      ptyp(8)=2
      ptyp(9)=2
      ptyp(10)=1
      ptyp(11)=2
      ptyp(12)=1
      ptyp(13)=1
      ptyp(14)=2
      ptyp(15)=1
      ptyp(16)=2
      ptyp(17)=2
      ptyp(18)=1
      ptyp(19)=1
      ptyp(20)=2
      ptyp(21)=1
      ptyp(22)=2
      ptyp(23)=2
      ptyp(24)=1
      ptyp(25)=1
      ptyp(26)=2
      ptyp(27)=1
      ptyp(28)=2
      ptyp(29)=2
      ptyp(30)=1
      ptyp(31)=2
      ptyp(32)=1
      pcell(1,1)=0.
      pcell(2,1)=0.      
      pcell(3,1)=0.
      pcell(1,2)=0.
      pcell(2,2)=0.5
      pcell(3,2)=0.
      pcell(1,3)=0.25
      pcell(2,3)=0.25
      pcell(3,3)=0.
      pcell(1,4)=0.25
      pcell(2,4)=0.75
      pcell(3,4)=0.
      pcell(1,5)=0.5
      pcell(2,5)=0.
      pcell(3,5)=0.
      pcell(1,6)=0.5
      pcell(2,6)=0.5
      pcell(3,6)=0.
      pcell(1,7)=0.75
      pcell(2,7)=0.25
      pcell(3,7)=0.
      pcell(1,8)=0.75
      pcell(2,8)=0.75
      pcell(3,8)=0.
      pcell(1,9)=0.
      pcell(2,9)=0.25
      pcell(3,9)=0.25
      pcell(1,10)=0.
      pcell(2,10)=0.75
      pcell(3,10)=0.25
      pcell(1,11)=0.25
      pcell(2,11)=0.
      pcell(3,11)=0.25
      pcell(1,12)=0.25
      pcell(2,12)=0.5
      pcell(3,12)=0.25
      pcell(1,13)=0.5
      pcell(2,13)=0.25
      pcell(3,13)=0.25
      pcell(1,14)=0.5
      pcell(2,14)=0.75
      pcell(3,14)=0.25
      pcell(1,15)=0.75
      pcell(2,15)=0.
      pcell(3,15)=0.25
      pcell(1,16)=0.75
      pcell(2,16)=0.5
      pcell(3,16)=0.25
      pcell(1,17)=0.
      pcell(2,17)=0.
      pcell(3,17)=0.5
      pcell(1,18)=0.
      pcell(2,18)=0.5
      pcell(3,18)=0.5
      pcell(1,19)=0.25
      pcell(2,19)=0.25
      pcell(3,19)=0.5
      pcell(1,20)=0.25
      pcell(2,20)=0.75
      pcell(3,20)=0.5
      pcell(1,21)=0.5
      pcell(2,21)=0.
      pcell(3,21)=0.5
      pcell(1,22)=0.5
      pcell(2,22)=0.5
      pcell(3,22)=0.5
      pcell(1,23)=0.75
      pcell(2,23)=0.25
      pcell(3,23)=0.5
      pcell(1,24)=0.75
      pcell(2,24)=0.75
      pcell(3,24)=0.5
      pcell(1,25)=0.
      pcell(2,25)=0.25
      pcell(3,25)=0.75
      pcell(1,26)=0.
      pcell(2,26)=0.75
      pcell(3,26)=0.75
      pcell(1,27)=0.25
      pcell(2,27)=0.
      pcell(3,27)=0.75
      pcell(1,28)=0.25
      pcell(2,28)=0.5
      pcell(3,28)=0.75
      pcell(1,29)=0.5
      pcell(2,29)=0.25
      pcell(3,29)=0.75
      pcell(1,30)=0.5
      pcell(2,30)=0.75
      pcell(3,30)=0.75
      pcell(1,31)=0.75
      pcell(2,31)=0.
      pcell(3,31)=0.75
      pcell(1,32)=0.75
      pcell(2,32)=0.5
      pcell(3,32)=0.75
      call loadprimcell(pcell,ptyp,patm,pvec,stacks)
c     print*,'natoms:',natoms
c     vmdfile='cuni-l11.input'
c     call writedyn86(vmdfile) 
c     vmdfile='cuni-l11.lammpstrj'
c     call writevmd(vmdfile) 
      call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms4-L1_1:',rmssum
        endif
         return
       endif
c*****************************************************
      omegafcc=lata0*lata0*lata0/4.d0
      omegabcc=lata0*lata0*lata0/2.d0
      omegahcp=2.*sqrt3*hcpc*lata0*lata0/4.d0
c     angstrom^3/atom -> cm^3/mol 
c     when multiplied with bmod, bmod will be eV
      omega=latis*latis*latis/4.d0
      omega=omega*0.62415
c

c   *****************************************************
c   * EAM fcn file ready to use in Sublimation Energy,  *
c   *     Vacancy Formation Energy, Bulk Modulus,       *
c   *     and Elastic Constant Calculations for FCC     *
c   *****************************************************
    
      do 1950 iat=1,natoms
        e(iat)=0.0
        fe(iat)=0.0
        fp(iat)=0.0
        fpp(iat)=0.0
        phie(iat)=0.0
        rho(iat)=0.0
        evacf(iat)=0.0
1950  continue
      
      rdrar = 1.0/drar 
      do 1980 iat=1,natoms
c
c     compute the contribution to rho(i) from particle j
c     Calculating the embedding energy F(Rho) in fe(i)
c      and Pair Energy contribution Phi(Rij) in phie(i)
c
ccdir$ novector
      do 1960 j=1,neigh(iat)
      rdist=rdis(iat,j,1)
      jat=rdis(iat,j,2)
      pij = rdist*rdrar + 1.0
      kij = pij
      kij = min0(kij,nrar-1)
      pij = pij - kij
c       make sure that p is less than 1.00
c       then if r is out of range, p = 1.00 and rho = last value of rhor
      pij = amin1(pij,1.0)
      ity=itype(iat)
      jty=itype(jat)
      rhojnei(iat,j) = ((rhorar3(kij,jty)*pij+
     $   rhorar2(kij,jty))*pij+
     $   rhorar1(kij,jty))*pij+
     $   rhorar(kij,jty)
      rhoinei(iat,j) = ((rhorar3(kij,ity)*pij+
     $   rhorar2(kij,ity))*pij+
     $   rhorar1(kij,ity))*pij+
     $   rhorar(kij,ity)
      rho(iat) = rho(iat) + rhojnei(iat,j)
      z2ij = ((z2rar3(kij,ity,jty)*pij+
     $          z2rar2(kij,ity,jty))*pij+
     $          z2rar1(kij,ity,jty))*pij+
     $          z2rar(kij,ity,jty)
      z2pij = (z2rar6(kij,ity,jty)*pij+z2rar5(kij,ity,jty))*pij+
     $         z2rar4(kij,ity,jty)
c     if(phitype) then
        phiij = z2ij
c     else
c       pij=1.0/rdist
c       phiij = z2ij * pij
c     endif
      phie(iat) = phie(iat) + 0.5 * phiij
1960  continue
ccdir$ vector
c        
1980  continue
c
c       store fsubi in e(i), the derivative of fsubi in fp
c       and the second derivative of fsubi in fpp
c
      rdrhoar = 1.0/drhoar
      do 2111 iat=1,natoms
      pij = rho(iat)*rdrhoar + 1.0
      kij = pij
      kij = max(1,min(kij,nrhoar-1))
      pij = pij - kij
      ity=itype(iat)
      fe(iat) = ((frhoar3(kij,ity)*pij+
     $        frhoar2(kij,ity))*pij+
     $        frhoar1(kij,ity))*pij+
     $        frhoar(kij,ity)
      fp(iat) = (frhoar6(kij,ity)*pij+
     $                  frhoar5(kij,ity))*pij+
     $                  frhoar4(kij,ity)
      fpp(iat) = (2.0*frhoar6(kij,ity)*pij+
     $                  frhoar5(kij,ity))/drhoar
2111  continue
c
c       store vacancy formation energy in evacf(i)
c
      rdrhoar = 1.0/drhoar
      do 3410 iat=1,natoms
        do 3420 j=1,neigh(iat)
        jat=rdis(iat,j,2)
        pij = (rho(jat)-rhoinei(iat,j))*rdrhoar + 1.0
        kij = pij
        kij = max(1,min(kij,nrhoar-1))
        pij = pij - kij
        ity=itype(iat)
        jty=itype(jat)
        evacf(iat) = evacf(iat) + ( fe(jat) - 
     $      (((frhoar3(kij,jty)*pij+
     $        frhoar2(kij,jty))*pij+
     $        frhoar1(kij,jty))*pij+
     $        frhoar(kij,jty))  )
3420    continue
      e(iat) = phie(iat) + fe(iat)
3410  continue

      do 3430 iat=1,natoms
      evacf(iat) = evacf(iat) + phie(iat)
3430  continue
      

      do 3460 iat=1,natoms
        ca(iat) = 0.0
        do 3461 cij=1,7
         econ(cij,iat) = 0.0
         bulkmod(cij,iat) = 0.0
         uu(cij) = 0.0
         ww(cij) = 0.0
         vv(cij) = 0.0
3461    continue
         do 3810 j=1,neigh(iat)
           jat=rdis(iat,j,2)
           do 3820 kc=1,3
3820       diski(kc) = rdisd(iat,j,kc)
c          call elascon(phitype,diski,jat,iat,uu,ww,vv,bulkmod) 
           call elastic(phitype,diski,jat,iat,uu,ww,vv,bulkmod) 
3810     continue
        econ(1,iat) = 0.5*uu(1) + fp(iat) *
     .                   ww(1) + fpp(iat) * vv(1)*vv(1)
        econ(2,iat) = 0.5*uu(2) + fp(iat) *
     .                   ww(2) + fpp(iat) * vv(1)*vv(2)
        econ(3,iat) = 0.5*uu(3) + fp(iat) *
     .                   ww(3) + fpp(iat) * vv(3)*vv(3)
        econ(4,iat) = 0.5*uu(4) + fp(iat) *
     .                   ww(4) + fpp(iat) * vv(1)*vv(4)
        econ(5,iat) = 0.5*uu(5) + fp(iat) *
     .                   ww(5) + fpp(iat) * vv(4)*vv(4)
        econ(6,iat) = 0.5*uu(6) + fp(iat) *
     .                   ww(6) + fpp(iat) * vv(5)*vv(6)
        econ(7,iat) = 0.5*uu(7) + fp(iat) *
     .                   ww(7) + fpp(iat) * vv(7)*vv(7)
        ca(iat)= fpp(iat)*(vv(1)*vv(2)-vv(7)*vv(7))
        ca(iat)= ca(iat)/ev_joule
        econ(1,iat) = econ(1,iat)/ev_joule
        econ(2,iat) = econ(2,iat)/ev_joule
        econ(3,iat) = econ(3,iat)/ev_joule
        econ(4,iat) = econ(4,iat)/ev_joule
        econ(5,iat) = econ(5,iat)/ev_joule
        econ(6,iat) = econ(6,iat)/ev_joule
        econ(7,iat) = econ(7,iat)/ev_joule
        bulkmod(1,iat) = 0.5*bulkmod(1,iat) + fp(iat) * 
     .    bulkmod(2,iat) + fpp(iat) * bulkmod(3,iat)**2
        bulkmod(4,iat) = 0.5*bulkmod(4,iat) + fp(iat) * 
     .    bulkmod(5,iat) + fpp(iat) * bulkmod(6,iat)**2
        bulkmod(1,iat) = bulkmod(1,iat)/ev_joule
        bulkmod(4,iat) = bulkmod(4,iat)/ev_joule
3460     continue

        omega0fcc = lata0*lata0*lata0/4.d0
        omega0bcc = lata0*lata0*lata0/2.d0
        omega0hcp = 2.*sqrt3*hcpc*lata0*lata0/4.d0
        omega0 = latis*latis*latis/4.d0
        c11(2)=0.0
        c12(2)=0.0
        c13(2)=0.0
        c33(2)=0.0
        c44(2)=0.0
        c55(2)=0.0
        c66(2)=0.0
        c66ca(2)=0.0
        bbb(2)=0.0
        bbb2(2)=0.0
        bbb3(2)=0.0
        bbb4(2)=0.0
        embed1(2)=0.0
        embed2(2)=0.0
        vacform1(2)=0.0
        vacform2(2)=0.0
        vacnum1(2)=0
        vacnum2(2)=0
        do 3311 iat=1,natoms
        bbb2(2)=bbb2(2)+(econ(1,iat)+2*econ(2,iat))/(3*omega0)
        c11(2)=c11(2)+econ(1,iat)/omega0
        c12(2)=c12(2)+econ(2,iat)/omega0
        c44(2)=c44(2)+econ(3,iat)/omega0
        c13(2)=c13(2)+econ(4,iat)/omega0
        c33(2)=c33(2)+econ(5,iat)/omega0
        c55(2)=c55(2)+econ(6,iat)/omega0
        c66(2)=c66(2)+econ(7,iat)/omega0
        c66ca(2)=c66ca(2)+(econ(2,iat)-ca(iat))/omega0
        bbb(2)=bbb(2)+bulkmod(1,iat)/(9*omega0)
        bbb3(2)=bbb3(2)+bulkmod(4,iat)/(9.*omega0)
        bbb4(2)=bbb4(2)+(2.*econ(1,iat)+2.*econ(2,iat)+
     .             4.*econ(4,iat)+econ(5,iat))/(9.*omega0)
        if( itype(iat) .eq. 1 ) then
          vacform1(2)=vacform1(2)+evacf(iat)
          embed1(2)=embed1(2)+e(iat)
          vacnum1(2)=vacnum1(2)+1
        else
          vacform2(2)=vacform2(2)+evacf(iat)
          embed2(2)=embed2(2)+e(iat)
          vacnum2(2)=vacnum2(2)+1
        endif
3311     continue
        bbb(2)=bbb(2)/natoms
        bbb2(2)=bbb2(2)/natoms
        bbb3(2)=bbb3(2)/natoms
        bbb4(2)=bbb4(2)/natoms
        c11(2)=c11(2)/natoms
        c12(2)=c12(2)/natoms
        c44(2)=c44(2)/natoms
        c13(2)=c13(2)/natoms
        c33(2)=c33(2)/natoms
        c55(2)=c55(2)/natoms
        c66(2)=c66(2)/natoms
        c66ca(2)=c66ca(2)/natoms
        if( vacnum1(2) .gt. 0 ) then
           vacform1(2)=vacform1(2)/vacnum1(2)
           embed1(2)=embed1(2)/vacnum1(2)
        endif
        if( vacnum2(2) .gt. 0 ) then
           vacform2(2)=vacform2(2)/vacnum2(2)
           embed2(2)=embed2(2)/vacnum2(2)
        endif
c      -------------------------------------------------
      if (prntflag.eq.1) then
        write(*,*)'L1_1 FCC CuNi Embed1:',embed1(2),'num:',vacnum1(2)
        write(*,*)'L1_1 FCC CuNi Embed2:',embed2(2),'num:',vacnum2(2)
      endif
       
       if (prntflag.eq.1) then
c*****************************************************
c*   Calculation Cu3Ni Alloy Stacking Fault and      *
c*           Antiphase Boundary Energies             *
c*            (100), (111) SISF and APB              *
c*****************************************************
      l1fcc=2
      if( ntypes .eq. 2 ) then
       t1=1
       t2=2
      else
       t1=2
       t2=4
      endif
      if (prntflag.eq.1) then
       write(*,*)'L1_1 FCC Cu3Ni a_0:',latis
      endif
      call loadfcc100(t1,t2,l1fcc,latis,stacks)
      call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rrms6:',rmssum
        endif
         return
       endif
       call calcper()
       call calctote(etot1,neigh,rdis)
       shifts(1)=abs(rv(1,stacks(2)+1)-rv(1,1))
       shifts(2)=abs(rv(2,stacks(2)+1)-rv(2,1))
       shifts(3)=0.0
       zplanes=stacks(3)/2
       if (prntflag.eq.1) then
        print*,'Shifting',zplanes,'planes of FCC100 to x:',
     .shifts(1),'and y:',shifts(2)
        print*,'area:',perlen(1)*perlen(2)
       endif
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms6:',rmssum
        endif
         return
       endif
       call calctote(etot2,neigh,rdis)
c      vmdfile='ni3al-100.lammpstrj'
c      call writevmd(stacks,vmdfile)
       apb100(2)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       apb100(2)=apb100(2)/2.0
       apb100(2)=10000*apb100(2)/ev_joule
       if (prntflag.eq.1) then
        print*,'etot1:',etot1,'etot2:',etot2
        print*,'APB(100):',apb100(2)
       endif
       l1=2
       if( ntypes .eq. 2 ) then
         t1=1
         t2=2
       else
         t1=1
         t2=2
       endif
       call loadfcc111(t1,t2,l1,latis,stacks)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms7:',rmssum
        endif
         return
       endif
       call calctote(etot1,neigh,rdis)
c      print*,'L1 Ecoh for 100 stacking:',etot1/natoms
c      vmdfile='Ni3Al-L12-111.lammpstrj'
c      call writevmd(stacks,vmdfile)
c      vmdfile='ni3al.lammpstrj'
c      call writevmd(stacks,vmdfile)
       shifts(1)=0.0
       shifts(2)=abs(rv(2,1)-rv(2,2))
       shifts(3)=0.0
       zplanes=4
       if (prntflag.eq.1) then
        print*,'Shifting',zplanes,'planes of FCC111 to x:',
     .shifts(1),'and y:',shifts(2)
        print*,'area:',perlen(1)*perlen(2)
       endif
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms8:',rmssum
        endif
         return
       endif
       call calctote(etot2,neigh,rdis)
       apb111(2)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       apb111(2)=apb111(2)/2.0
       apb111(2)=10000*apb111(2)/ev_joule
       if (prntflag.eq.1) then
        print*,'etot1:',etot1,'etot2:',etot2
        print*,'APB(111):',apb111(2)
       endif
c      vmdfile='ni3al-apb.lammpstrj'
c      call writevmd(stacks,vmdfile)
       stacks(3)=stacks(3)-1
       call loadfcc111(t1,t2,l1,latis,stacks)
c      call calcper()
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms9:',rmssum
        endif
         return
       endif
       call calctote(etot1,neigh,rdis)
       xshift = rv(1,1) - rv(1,stacks(2)+2)
       yshift = rv(2,1) - rv(2,stacks(2)+2)
       shifts(1)=2.0*abs(xshift)/3.0
       shifts(2)=2.0*abs(yshift)/3.0
       shifts(3)=0.0
        zplanes=4
       if (prntflag.eq.1) then
        print*,'Shifting',zplanes,'planes of FCC111 to x:',
     .shifts(1),'and y:',shifts(2)
        print*,'area:',perlen(1)*perlen(2)
       endif
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms10:',rmssum
        endif
         return
       endif
       call calctote(etot2,neigh,rdis)
       sisf111(2)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       sisf111(2)=sisf111(2)/2.0
       sisf111(2)=10000*sisf111(2)/ev_joule
       if (prntflag.eq.1) then
        print*,'etot1:',etot1,'etot2:',etot2
        print*,'SISF(111):',sisf111(2)
       endif
       stacks(3)=stacks(3)+1
c      vmdfile='ni3al-sisf.lammpstrj'
c      call writevmd(stacks,vmdfile)
       endif

        allphonon=0
        if (allphonon.eq.1) then
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     PHONON CALCULATION 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        slf1=0
        slf2=0
        d11=0
        d12=0
        d21=0
        d22=0
      do 69 j=1,NN
69      W(j)=0.
      do 74 i2=1,NN
      do 74 j2=1,NN
74      dmic(i2,j2)=DCMPLX(0.,0.)
      do 75 i2=1,3
      do 75 j2=1,3
        self1(i2,j2)=DCMPLX(0.,0.)
        self2(i2,j2)=DCMPLX(0.,0.)
75    continue 
      do 76 ii=1,pnts
      do 81 i2=1,3
      do 81 j2=1,3
        dmat11(i2,j2,ii)=DCMPLX(0.,0.)
        dmat12(i2,j2,ii)=DCMPLX(0.,0.)
        dmat21(i2,j2,ii)=DCMPLX(0.,0.)
        dmat22(i2,j2,ii)=DCMPLX(0.,0.)
81    continue 
76    continue
c
c       Brillouin Zone:
c       Simple Cubic = 1
c       Face Centered Cubic = 2 
c       brillouinzone=1
        brillouinzone=2
c       latak=3.62D+00
        latak=latis
        kpi=(2.0d0*pi/latak)
        kpi1=(1.0d0*pi/latak)
        kpi2=(0.5d0*pi/latak)
        kpi4=(4.0d0*pi/latak)
        points=20
        kcount=0
        ua1=1
        ua2=2
        ua3=stacks(2)+1
        ua4=stacks(2)+2
        if(brillouinzone.eq.1) then
        do 114 kx=0,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=0.0d0
           kpoints(kcount,3)=0.0d0+(kx/DBLE(points))*kpi1
114     continue
        do 115 kx=1,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=0.0d0+(kx/DBLE(points))*kpi1
           kpoints(kcount,3)=1.0d0*kpi1
115     continue
        do 116 kx=1,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=1.0d0*kpi1-(kx/DBLE(points))*kpi1
           kpoints(kcount,3)=1.0d0*kpi1-(kx/DBLE(points))*kpi1
116     continue
        do 117 kx=1,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0+(kx/DBLE(points))*kpi1
           kpoints(kcount,2)=0.0d0+(kx/DBLE(points))*kpi1
           kpoints(kcount,3)=0.0d0+(kx/DBLE(points))*kpi1
117     continue
        else if (brillouinzone.eq.2) then
        do 118 kx=0,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=0.0d0
           kpoints(kcount,3)=0.0d0+(kx/DBLE(points))*kpi
118     continue
        do 119 kx=1,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=0.0d0+(kx/DBLE(points))*kpi
           kpoints(kcount,3)=1.0d0*kpi
119     continue
        do 120 kx=1,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=1.0d0*kpi-(kx/DBLE(points))*kpi
           kpoints(kcount,3)=1.0d0*kpi-(kx/DBLE(points))*kpi
120     continue
        do 121 kx=1,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0+(kx/DBLE(points))*kpi1
           kpoints(kcount,2)=0.0d0+(kx/DBLE(points))*kpi1
           kpoints(kcount,3)=0.0d0+(kx/DBLE(points))*kpi1
121     continue
        endif
           write(*,*)'  kx:      ky:      kz:'
           do 1112 kx=1,kcount
             kpnt(1)=kpoints(kx,1)
             kpnt(2)=kpoints(kx,2)
             kpnt(3)=kpoints(kx,3)
             write(*,3232)kpnt(1),kpnt(2),kpnt(3)
1112       continue
          m1=amass(itype(ua1))*conmas
          m2=amass(itype(ua2))*conmas
          m3=amass(itype(ua3))*conmas
          m4=amass(itype(ua4))*conmas
           if (brillouinzone.eq.1) write(*,*)' Brillouin Zone : SC'
           if (brillouinzone.eq.2) write(*,*)' Brillouin Zone : FCC'
           write(*,*)' Lattice constant:',latak
           write(*,*)' Unit cell types for 1-4 unit cell atoms:'
           write(*,7557)IDINT(ucell(ua1)),IDINT(ucell(ua2)),
     .                  IDINT(ucell(ua3)),IDINT(ucell(ua4))
           write(*,*)' Atom numbers for 1-4 unit cell atoms:'
           write(*,7557)ua1,ua2,ua3,ua4
           write(*,*)' Mass for 1-4 unit cell atoms:'
           write(*,*)'m1=',m1
           write(*,*)'m2=',m2
           write(*,*)'m3=',m3
           write(*,*)'m4=',m4
      write(*,*)' Force Constants calculation for k-points:',kcount
7557    format(4(1x,i3))

        do 122 iat=1,natoms
           do 1001 j=1,neigh(iat)
           jat=rdis(iat,j,2)
c       ATOMS SELF CONTRIBUTION
          if (iat .eq. ua1 .or. jat .eq. ua1 .or.
     .        iat .eq. ua2 .or. jat .eq. ua2 ) then
           call forcon(neigh,rdis,rdisd,j,iat,fp,fpp,dcon)
           if( iat .eq. ua1 ) then
           do 1020 kc=1,3 
           do 1020 kr=1,3 
1020       self1(kr,kc)=self1(kr,kc)+DCMPLX(dcon(kr,kc), 0.)
             slf1=slf1+1
           endif
           if( iat .eq. ua2 ) then
           do 1021 kc=1,3 
           do 1021 kr=1,3 
1021       self2(kr,kc)=self2(kr,kc)+DCMPLX(dcon(kr,kc), 0.)
             slf2=slf2+1
           endif
          endif
c       ATOMS SELF CONTRIBUTION END
c       ---------------------------
c        ATOM TYPE 1 CONTRIBUTIONS
          if (iat .eq. ua1 .or. jat .eq. ua1) then
           call forcon(neigh,rdis,rdisd,j,iat,fp,fpp,dcon)
           do 1011 kc=1,3
1011       diski(kc) = rdisd(iat,j,kc)
          do 1101 kx=1,kcount
        kxr1=kpoints(kx,1)*diski(1)
        kxr2=kpoints(kx,2)*diski(2)
        kxr3=kpoints(kx,3)*diski(3)
        kxr=kxr1+kxr2+kxr3
          if ( iat .eq. ua1 ) then
             if (IDINT(ucell(jat)) .eq. 1) then
           do 1022 kc=1,3 
           do 1022 kr=1,3 
                dmat11(kr,kc,kx)=dmat11(kr,kc,kx)+
     .    DCMPLX(dcon(kr,kc)*DCOS(kxr),dcon(kr,kc)*DSIN(kxr))
1022       continue
                d11=d11+1
             endif
             if (IDINT(ucell(jat)) .eq. 2) then
           do 1023 kc=1,3 
           do 1023 kr=1,3 
                dmat12(kr,kc,kx)=dmat12(kr,kc,kx)+
     .    DCMPLX(dcon(kr,kc)*DCOS(kxr),dcon(kr,kc)*DSIN(kxr))
1023       continue
                d12=d12+1
             endif
          endif
1101      continue
          endif
c        ATOM TYPE 1 CONTRIBUTIONS END
c       ------------------------------
c        ATOM TYPE 2 CONTRIBUTIONS
          if (iat .eq. ua2 .or. jat .eq. ua2) then
           call forcon(neigh,rdis,rdisd,j,iat,fp,fpp,dcon)
           do 1012 kc=1,3
1012       diski(kc) = rdisd(iat,j,kc)
          do 1102 kx=1,kcount
        kxr1=kpoints(kx,1)*diski(1)
        kxr2=kpoints(kx,2)*diski(2)
        kxr3=kpoints(kx,3)*diski(3)
        kxr=kxr1+kxr2+kxr3
          if ( iat .eq. ua2) then
             if (IDINT(ucell(jat)) .eq. 1) then
           do 1027 kc=1,3 
           do 1027 kr=1,3 
                dmat21(kr,kc,kx)=dmat21(kr,kc,kx)+
     .    DCMPLX(dcon(kr,kc)*DCOS(kxr),dcon(kr,kc)*DSIN(kxr))
1027       continue
                d21=d21+1
             endif
             if (IDINT(ucell(jat)) .eq. 2) then
           do 1028 kc=1,3 
           do 1028 kr=1,3 
                dmat22(kr,kc,kx)=dmat22(kr,kc,kx)+
     .    DCMPLX(dcon(kr,kc)*DCOS(kxr),dcon(kr,kc)*DSIN(kxr))
1028       continue
                d22=d22+1
             endif
          endif
1102     continue
           endif
c        ATOM TYPE 2 CONTRIBUTIONS END
c       ------------------------------
1001     continue
122     continue

       write(*,*)' # of Self Matrix Elements:'
       write(*,*)slf1,slf2
       write(*,*)' # of Dynamic Matrix Elements:'
       write(*,*)d11/kcount,d12/kcount,d21/kcount,d22/kcount
       write(*,*)' Dynamical Matrix calculation for k-points:',kcount
     
      do 1050 kx=1,kcount
       mm=1./(amass(itype(ua1))*amass(itype(ua1))*conmas*conmas)
       do 1051 a=1,3
       do 1051 b=1,3
1051   dmat11(a,b,kx)=(dmat11(a,b,kx)-self1(a,b))*dsqrt(mm)
       mm=1./(amass(itype(ua1))*amass(itype(ua2))*conmas*conmas)
       do 1052 a=1,3
       do 1052 b=1,3
1052   dmat12(a,b,kx)=dmat12(a,b,kx)*dsqrt(mm)
       mm=1./(amass(itype(ua2))*amass(itype(ua1))*conmas*conmas)
       do 1053 a=1,3
       do 1053 b=1,3
1053   dmat21(a,b,kx)=dmat21(a,b,kx)*dsqrt(mm)
       mm=1./(amass(itype(ua2))*amass(itype(ua2))*conmas*conmas)
       do 1054 a=1,3
       do 1054 b=1,3
1054   dmat22(a,b,kx)=(dmat22(a,b,kx)-self1(a,b))*dsqrt(mm)
1050  continue
       write(*,*)' Solving Eigenvalue Problem for k-points:',kcount
      do 1055 kx=1,kcount
        do 1056 a=1,3
        do 1056 b=1,3
1056       dmic(a,b)=dmat11(a,b,kx)
        do 1057 a=1,3
        do 1057 b=1,3
1057       dmic(a,b+3)=dmat12(a,b,kx)
        do 1058 a=1,3
        do 1058 b=1,3
1058       dmic(a+3,b)=dmat21(a,b,kx)
        do 1059 a=1,3
        do 1059 b=1,3
1059       dmic(a+3,b+3)=dmat22(a,b,kx)
c       write(*,*)dmic
c       Query Optimal Workspace (LWORK=-1)
        LWORK=-1
        CALL ZHEEV('Vector','Lower',NN,dmic,LDA,W,WORK,LWORK,
     .  RWORK,INFO)
c       Solve Problem with LWORK workspace
        LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
        CALL ZHEEV('Vector','Lower',NN,dmic,LDA,W,WORK,LWORK,
     .  RWORK,INFO)
        kpnt(1)=kpoints(kx,1)
        kpnt(2)=kpoints(kx,2)
        kpnt(3)=kpoints(kx,3)
c       Check Convergence
        if (INFO.eq.0) then
        do 1066 j=1,NN
        W(j)=DABS(DSQRT(W(j)))/(2.*pi)
c       if eigenvalue is NaN then append it to 0.
1066    continue
        do 1068 j=1,NN
1068    if (W(j) .ne. W(j)) W(j)=0.
          write(*,5555)kpnt(1),kpnt(2),kpnt(3),' : ',(W(j),j=1,NN)
c         write(77,5775)kpnt(1),kpnt(2),kpnt(3),(W(j),j=1,NN)
        else if (INFO .gt. 0) then
          write(*,*)'Eigenvalue ',INFO,' is not converged for k=',
     .               kpnt(1),kpnt(2),kpnt(3)
        else
          write(*,*)INFO,' th argument had an illegal value for k=',
     .               kpnt(1),kpnt(2),kpnt(3)
        endif
1055  continue
5555  format(3(f7.4,1x),a,6(E13.6,1x))
5557  format(3(f7.4,1x),a,3(E13.6,1x))
3232  format(3f9.4)
c5775  format(3(f25.16,1x),6(E25.16,1x))
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     END OF PHONON CALCULATION 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if (prntflag.eq.1) then
c*****************************************************
c*   Calculation Cu3Ni Alloy Stacking Fault and      *
c*           Antiphase Boundary Energies             *
c*            (100), (111) SISF and APB              *
c*****************************************************
      l1fcc=2
       t1=2
       t2=2
      latak=3.615
      call loadfcc100(t1,t2,l1fcc,latak,stacks)
      call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms4:',rmssum
        endif
         return
       endif
       call calcper()
       call calctote(etot1,neigh,rdis)
       shifts(1)=abs(rv(1,stacks(2)+1)-rv(1,1))
       shifts(2)=abs(rv(2,stacks(2)+1)-rv(2,1))
       shifts(3)=0.0
       zplanes=stacks(3)/2
       if (prntflag.eq.1) then
        print*,'Shifting',zplanes,'planes of FCC100 to x:',
     .shifts(1),'and y:',shifts(2)
        print*,'area:',perlen(1)*perlen(2)
       endif
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms6:',rmssum
        endif
         return
       endif
       call calctote(etot2,neigh,rdis)
c      vmdfile='ni3al-100.lammpstrj'
c      call writevmd(stacks,vmdfile)
       apb100(3)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       apb100(3)=apb100(3)/2.0
       apb100(3)=10000*apb100(3)/ev_joule
       if (prntflag.eq.1) then
        print*,'etot1:',etot1,'etot2:',etot2
        print*,'APB(100):',apb100(3)
       endif
       l1=2
         t1=2
         t2=2
       call loadfcc111(t1,t2,l1,latak,stacks)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms7:',rmssum
        endif
         return
       endif
c      call calcper()
       call calctote(etot1,neigh,rdis)
c      print*,'L1 Ecoh for 100 stacking:',etot1/natoms
c      vmdfile='Ni3Al-L12-111.lammpstrj'
c      call writevmd(stacks,vmdfile)
c      vmdfile='ni3al.lammpstrj'
c      call writevmd(stacks,vmdfile)
       shifts(1)=0.0
       shifts(2)=abs(rv(2,1)-rv(2,2))
       shifts(3)=0.0
       zplanes=4
       if (prntflag.eq.1) then
        print*,'Shifting',zplanes,'planes of FCC111 to x:',
     .shifts(1),'and y:',shifts(2)
        print*,'area:',perlen(1)*perlen(2)
       endif
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms8:',rmssum
        endif
         return
       endif
       call calctote(etot2,neigh,rdis)
       apb111(3)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       apb111(3)=apb111(3)/2.0
       apb111(3)=10000*apb111(3)/ev_joule
       if (prntflag.eq.1) then
        print*,'etot1:',etot1,'etot2:',etot2
        print*,'APB(111):',apb111(3)
       endif
c      vmdfile='ni3al-apb.lammpstrj'
c      call writevmd(stacks,vmdfile)
       stacks(3)=stacks(3)-1
       call loadfcc111(t1,t2,l1,latak,stacks)
c      call calcper()
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms9:',rmssum
        endif
         return
       endif
       call calctote(etot1,neigh,rdis)
       xshift = rv(1,1) - rv(1,stacks(2)+2)
       yshift = rv(2,1) - rv(2,stacks(2)+2)
       shifts(1)=2.0*abs(xshift)/3.0
       shifts(2)=2.0*abs(yshift)/3.0
       shifts(3)=0.0
        zplanes=4
       if (prntflag.eq.1) then
        print*,'Shifting',zplanes,'planes of FCC111 to x:',
     .shifts(1),'and y:',shifts(2)
        print*,'area:',perlen(1)*perlen(2)
       endif
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms10:',rmssum
        endif
         return
       endif
       call calctote(etot2,neigh,rdis)
       sisf111(3)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       sisf111(3)=sisf111(3)/2.0
       sisf111(3)=10000*sisf111(3)/ev_joule
       if (prntflag.eq.1) then
        print*,'etot1:',etot1,'etot2:',etot2
        print*,'SISF(111):',sisf111(3)
       endif
       stacks(3)=stacks(3)+1
c      vmdfile='ni3al-sisf.lammpstrj'
c      call writevmd(stacks,vmdfile)
       endif
       if (prntflag.eq.1) then
c*****************************************************
c*   Calculation Cu3Ni Alloy Stacking Fault and      *
c*           Antiphase Boundary Energies             *
c*            (100), (111) SISF and APB              *
c*****************************************************
      l1fcc=2
       t1=1
       t2=1
       latak=3.52308
      call loadfcc100(t1,t2,l1fcc,latak,stacks)
      call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms4:',rmssum
        endif
         return
       endif
       call calcper()
       call calctote(etot1,neigh,rdis)
       shifts(1)=abs(rv(1,stacks(2)+1)-rv(1,1))
       shifts(2)=abs(rv(2,stacks(2)+1)-rv(2,1))
       shifts(3)=0.0
       zplanes=stacks(3)/2
       if (prntflag.eq.1) then
        print*,'Shifting',zplanes,'planes of FCC100 to x:',
     .shifts(1),'and y:',shifts(2)
        print*,'area:',perlen(1)*perlen(2)
       endif
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms6:',rmssum
        endif
         return
       endif
       call calctote(etot2,neigh,rdis)
c      vmdfile='ni3al-100.lammpstrj'
c      call writevmd(stacks,vmdfile)
       apb100(4)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       apb100(4)=apb100(4)/2.0
       apb100(4)=10000*apb100(4)/ev_joule
       if (prntflag.eq.1) then
        print*,'etot1:',etot1,'etot2:',etot2
        print*,'APB(100):',apb100(4)
       endif
       l1=2
         t1=1
         t2=1
       call loadfcc111(t1,t2,l1,latak,stacks)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms7:',rmssum
        endif
         return
       endif
c      call calcper()
       call calctote(etot1,neigh,rdis)
c      print*,'L1 Ecoh for 100 stacking:',etot1/natoms
c      vmdfile='Ni3Al-L12-111.lammpstrj'
c      call writevmd(stacks,vmdfile)
c      vmdfile='ni3al.lammpstrj'
c      call writevmd(stacks,vmdfile)
       shifts(1)=0.0
       shifts(2)=abs(rv(2,1)-rv(2,2))
       shifts(3)=0.0
       zplanes=4
       if (prntflag.eq.1) then
        print*,'Shifting',zplanes,'planes of FCC111 to x:',
     .shifts(1),'and y:',shifts(2)
        print*,'area:',perlen(1)*perlen(2)
       endif
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms8:',rmssum
        endif
         return
       endif
       call calctote(etot2,neigh,rdis)
       apb111(4)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       apb111(4)=apb111(4)/2.0
       apb111(4)=10000*apb111(4)/ev_joule
       if (prntflag.eq.1) then
        print*,'etot1:',etot1,'etot2:',etot2
        print*,'APB(111):',apb111(4)
       endif
c      vmdfile='ni3al-apb.lammpstrj'
c      call writevmd(stacks,vmdfile)
       stacks(3)=stacks(3)-1
       call loadfcc111(t1,t2,l1,latak,stacks)
c      call calcper()
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms9:',rmssum
        endif
         return
       endif
       call calctote(etot1,neigh,rdis)
       xshift = rv(1,1) - rv(1,stacks(2)+2)
       yshift = rv(2,1) - rv(2,stacks(2)+2)
       shifts(1)=2.0*abs(xshift)/3.0
       shifts(2)=2.0*abs(yshift)/3.0
       shifts(3)=0.0
        zplanes=4
       if (prntflag.eq.1) then
        print*,'Shifting',zplanes,'planes of FCC111 to x:',
     .shifts(1),'and y:',shifts(2)
        print*,'area:',perlen(1)*perlen(2)
       endif
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms10:',rmssum
        endif
         return
       endif
       call calctote(etot2,neigh,rdis)
       sisf111(4)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       sisf111(4)=sisf111(4)/2.0
       sisf111(4)=10000*sisf111(4)/ev_joule
       if (prntflag.eq.1) then
        print*,'etot1:',etot1,'etot2:',etot2
        print*,'SISF(111):',sisf111(4)
       endif
       stacks(3)=stacks(3)+1
c      vmdfile='ni3al-sisf.lammpstrj'
c      call writevmd(stacks,vmdfile)
       endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Cu PHONON CALCULATION 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      stacks(1)=8
      stacks(2)=8
      stacks(3)=8
      patm=4
      pvec(1)=3.52
      pvec(2)=3.52
      pvec(3)=3.52
      ptyp(1)=1
      ptyp(2)=1
      ptyp(3)=1
      ptyp(4)=1
      pcell(1,1)=0.
      pcell(2,1)=0.      
      pcell(3,1)=0.
      pcell(1,2)=0.5
      pcell(2,2)=0.5
      pcell(3,2)=0.
      pcell(1,3)=0.
      pcell(2,3)=0.5
      pcell(3,3)=0.5
      pcell(1,4)=0.5
      pcell(2,4)=0.
      pcell(3,4)=0.5
      call loadprimcell(pcell,ptyp,patm,pvec,stacks)
c     vmdfile='ni-fcc.input'
c     call writedyn86(vmdfile) 
c     vmdfile='ni-fcc.lammpstrj'
c     call writevmd(vmdfile) 
      call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms44-NI:',rmssum
        endif
         return
       endif
        slf1=0
        d11=0
      do 169 j=1,NNN
169      WB(j)=0.
      do 174 i2=1,NNN
      do 174 j2=1,NNN
174      dmicb(i2,j2)=DCMPLX(0.,0.)
      do 175 i2=1,3
      do 175 j2=1,3
        self1(i2,j2)=DCMPLX(0.,0.)
175    continue 
      do 176 ii=1,pnts
      do 181 i2=1,3
      do 181 j2=1,3
        dmat11(i2,j2,ii)=DCMPLX(0.,0.)
181    continue 
176    continue
c
c       Brillouin Zone:
c       Simple Cubic = 1
c       Face Centered Cubic = 2 
c       brillouinzone=1
        brillouinzone=2
c       latak=3.62D+00
        latak=3.52
        kpi=(2.0d0*pi/latak)
        kpi1=(1.0d0*pi/latak)
        kpi2=(0.5d0*pi/latak)
        kpi4=(4.0d0*pi/latak)
        points=20
        kcount=0
        ua1=1
        ua2=2
        ua3=3
        ua4=4
        if(brillouinzone.eq.1) then
        do 214 kx=0,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=0.0d0
           kpoints(kcount,3)=0.0d0+(kx/DBLE(points))*kpi1
214     continue
        do 215 kx=1,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=0.0d0+(kx/DBLE(points))*kpi1
           kpoints(kcount,3)=1.0d0*kpi1
215     continue
        do 216 kx=1,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=1.0d0*kpi1-(kx/DBLE(points))*kpi1
           kpoints(kcount,3)=1.0d0*kpi1-(kx/DBLE(points))*kpi1
216     continue
        do 217 kx=1,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0+(kx/DBLE(points))*kpi1
           kpoints(kcount,2)=0.0d0+(kx/DBLE(points))*kpi1
           kpoints(kcount,3)=0.0d0+(kx/DBLE(points))*kpi1
217     continue
        else if (brillouinzone.eq.2) then
        do 218 kx=0,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=0.0d0
           kpoints(kcount,3)=0.0d0+(kx/DBLE(points))*kpi
218     continue
        do 219 kx=1,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=0.0d0+(kx/DBLE(points))*kpi
           kpoints(kcount,3)=1.0d0*kpi
219     continue
        do 220 kx=1,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0
           kpoints(kcount,2)=1.0d0*kpi-(kx/DBLE(points))*kpi
           kpoints(kcount,3)=1.0d0*kpi-(kx/DBLE(points))*kpi
220     continue
        do 221 kx=1,points
           kcount=kcount+1
           kpoints(kcount,1)=0.0d0+(kx/DBLE(points))*kpi1
           kpoints(kcount,2)=0.0d0+(kx/DBLE(points))*kpi1
           kpoints(kcount,3)=0.0d0+(kx/DBLE(points))*kpi1
221     continue
        endif
           write(*,*)'  kx:      ky:      kz:'
           do 1312 kx=1,kcount
             kpnt(1)=kpoints(kx,1)
             kpnt(2)=kpoints(kx,2)
             kpnt(3)=kpoints(kx,3)
             write(*,3232)kpnt(1),kpnt(2),kpnt(3)
1312       continue
          m1=amass(itype(ua1))*conmas
          m2=amass(itype(ua2))*conmas
          m3=amass(itype(ua3))*conmas
          m4=amass(itype(ua4))*conmas
           if (brillouinzone.eq.1) write(*,*)' Brillouin Zone : SC'
           if (brillouinzone.eq.2) write(*,*)' Brillouin Zone : FCC'
           write(*,*)' Lattice constant:',latak
           write(*,*)' Unit cell types for 1-4 unit cell atoms:'
           write(*,7557)IDINT(ucell(ua1)),IDINT(ucell(ua2)),
     .                  IDINT(ucell(ua3)),IDINT(ucell(ua4))
           write(*,*)' Atom numbers for 1-4 unit cell atoms:'
           write(*,7557)ua1,ua2,ua3,ua4
           write(*,*)' Mass for 1-4 unit cell atoms:'
           write(*,*)'m1=',m1
           write(*,*)'m2=',m2
           write(*,*)'m3=',m3
           write(*,*)'m4=',m4
      write(*,*)' Force Constants calculation for k-points:',kcount

        do 142 iat=1,natoms
           do 1201 j=1,neigh(iat)
           jat=rdis(iat,j,2)
c       ATOMS SELF CONTRIBUTION
          if (iat .eq. ua1 ) then
           call forcon(neigh,rdis,rdisd,j,iat,fp,fpp,dcon)
           do 1220 kc=1,3 
           do 1220 kr=1,3 
1220       self1(kr,kc)=self1(kr,kc)+dcon(kr,kc)
             slf1=slf1+1
          endif
c       ATOMS SELF CONTRIBUTION END
c       ---------------------------
c        ATOM TYPE 1 CONTRIBUTIONS
          if (iat .eq. ua1) then
           call forcon(neigh,rdis,rdisd,j,iat,fp,fpp,dcon)
           do 1211 kc=1,3
1211       diski(kc) = rdisd(iat,j,kc)
          do 1301 kx=1,kcount
        kxr1=kpoints(kx,1)*diski(1)
        kxr2=kpoints(kx,2)*diski(2)
        kxr3=kpoints(kx,3)*diski(3)
        kxr=kxr1+kxr2+kxr3
           do 1222 kc=1,3 
           do 1222 kr=1,3 
                dmat11(kr,kc,kx)=dmat11(kr,kc,kx)+
     .    DCMPLX(dcon(kr,kc)*DCOS(kxr),dcon(kr,kc)*DSIN(kxr))
1222       continue
                d11=d11+1
1301      continue
          endif
c        ATOM TYPE 1 CONTRIBUTIONS END
c       ------------------------------
1201     continue
142     continue

       write(*,*)' # of Self Matrix Elements:'
       write(*,*)slf1
       write(*,*)' # of Dynamic Matrix Elements:'
       write(*,*)d11/kcount
       write(*,*)' Dynamical Matrix calculation for k-points:',kcount
     
      do 1250 kx=1,kcount
        do 1251 a=1,3
        do 1251 b=1,3
        dmat11(a,b,kx)=(dmat11(a,b,kx)-self1(a,b))/
     .                (amass(itype(ua1))*conmas)
1251    continue
1250  continue
       write(*,*)' Solving Eigenvalue Problem for k-points:',kcount
      do 1255 kx=1,kcount
        do 1256 a=1,3
        do 1256 b=1,3
1256       dmicb(a,b)=dmat11(a,b,kx)
c       write(*,*)dmic
c       Query Optimal Workspace (LWORK=-1)
        LWORKB=-1
        CALL ZHEEV('Vector','Lower',NNN,dmicb,LDB,WB,WORKB,LWORKB,
     .  RWORKB,INFOB)
c       Solve Problem with LWORK workspace
        LWORKB = MIN( LWMAXB, INT( WORKB( 1 ) ) )
        CALL ZHEEV('Vector','Lower',NNN,dmicb,LDB,WB,WORKB,LWORKB,
     .  RWORKB,INFOB)
        kpnt(1)=kpoints(kx,1)
        kpnt(2)=kpoints(kx,2)
        kpnt(3)=kpoints(kx,3)
c       Check Convergence
        if (INFOB.eq.0) then
        do 1266 j=1,NNN
        WB(j)=DABS(DSQRT(WB(j)))/(2.*pi)
c       if eigenvalue is NaN then append it to 0.
1266    continue
        do 1268 j=1,NNN
1268    if (WB(j) .ne. WB(j)) WB(j)=0.
          write(*,5555)kpnt(1),kpnt(2),kpnt(3),' : ',(WB(j),j=1,NNN)
c         write(77,5775)kpnt(1),kpnt(2),kpnt(3),(W(j),j=1,NN)
        else if (INFOB .gt. 0) then
          write(*,*)'Eigenvalue ',INFOB,' is not converged for k=',
     .               kpnt(1),kpnt(2),kpnt(3)
        else
          write(*,*)INFOB,' th argument had an illegal value for k=',
     .               kpnt(1),kpnt(2),kpnt(3)
        endif
1255  continue
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     END OF PHONON CALCULATION 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif
       
       stacks(1)=stck(1)
       stacks(2)=stck(2)
       stacks(3)=stck(3)
       endif
c      END of CuNi NiCu
c      ---------------------------------------------               


       if( ( ntypes .lt. 3 .and. ( at1 .eq. 2 .and. at2 .eq. 3 ) ) 
     $ .or. ( ntypes .gt. 2 ) ) then
c      ---------------------------------------------               
c      AlNi BCC B2
       if (prntflag.eq.1) then
        write(*,*)'------- AlNi BCC B2 -------'
       endif
       if (prmid.gt.0) then
        write(*,*)'PROC:',procnum,' part:',prmid,' AlNi BCC B2 '
       endif
       b2=1
       if( ntypes .eq. 2 ) then
        t1=1
        t2=2
       else
        t1=2
        t2=3
       endif
       call loadbcc100(t1,t2,b2,lata0,stacks)
c      vmdfile='AlNi-B2.lammpstrj'
c      call writevmd(stacks,vmdfile)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms2:',rmssum
        endif
         return
       endif
       if (prntflag.eq.1) then
        write(*,*)'AlNi Test a_0:',lata0
       endif
       call calclatcon(minenergy,latis,lata0,neigh,rdis,iter,prntflag)
       b2alat(1)=latis
       if (prntflag.eq.1) then
        write(*,*)'AlNi  B2 a_0:',latis
       endif
       call loadbcc100(t1,t2,b2,latis,stacks)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms3:',rmssum
        endif
         return
       endif
       tn=2
       typs(1)=t1
       typs(2)=t2
       call calctotetype(tn,typs,tnum,embed,etot1,neigh,rdis)
       etot1=etot1/natoms
       b2ecoh(1)=etot1
       embed1(4)=embed(1)
       embed2(4)=embed(2)
       if (prntflag.eq.1) then
        print*,'AlNi B2 Ecoh:',etot1
        print*,'AlNi B2 Embed1:',embed(1),'num:',tnum(1)
        print*,'AlNi B2 Embed2:',embed(2),'num:',tnum(2)
       endif

c****************************************************
c*   Calculation of Equilibrium Lattice Constant    *
c*                for L1 FCC Alloy                  *
c****************************************************
     
c     Ni3Al FCC L1_2
      if (prntflag.eq.1) then
       write(*,*)'------- Ni3Al FCC L1_2 -------'
      endif
      if (prmid.gt.0) then
        write(*,*)'PROC:',procnum,' part:',prmid,' Ni3Al FCC L1_2 '
      endif
      l1fcc=3
      if( ntypes .eq. 2 ) then
       t1=2
       t2=1
      else
       t1=3
       t2=2
      endif
      call loadfcc100(t1,t2,l1fcc,lata0,stacks)
      call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms4:',rmssum
        endif
         return
       endif
      if (prntflag.eq.1) then
       write(*,*)'Test a_0:',lata0
      endif
      call calclatcon(minenergy,latis,lata0,neigh,rdis,iter,prntflag)
      l12alat(1)=latis
      if (prntflag.eq.1) then
       write(*,*)'L1 FCC a_0:',latis
       write(*,*)'L1 FCC Ecoh:',minenergy/natoms
      endif
      l12ecoh(1)=minenergy/natoms
c      vmdfile='Ni3Al-L12.lammpstrj'
c      call writevmd(stacks,vmdfile)
      
c*****************************************************

c     ------------------------------------------------------------
c       Generating bulk positions for FCC structure
c     ------------------------------------------------------------
        l1fcc=3
        if( ntypes .eq. 2 ) then
         t1=2 
         t2=1
        else
         t1=3 
         t2=2
        endif
        call loadfcc100(t1,t2,l1fcc,latis,stacks)
        call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop .eq. 1) then
        if (prntflag.eq.1) then
         write(*,*)'rms5:',rmssum
        endif
         return
       endif
c     ------------------------------------------------------------

      omegafcc=lata0*lata0*lata0/4.d0
      omegabcc=lata0*lata0*lata0/2.d0
      omegahcp=2.*sqrt3*hcpc*lata0*lata0/4.d0
c     angstrom^3/atom -> cm^3/mol 
c     when multiplied with bmod, bmod will be eV
      omega=latis*latis*latis/4.d0
      omega=omega*0.62415
c


c   *****************************************************
c   * EAM fcn file ready to use in Sublimation Energy,  *
c   *     Vacancy Formation Energy, Bulk Modulus,       *
c   *     and Elastic Constant Calculations for FCC     *
c   *****************************************************
    
      do 1840 iat=1,natoms
        e(iat)=0.0
        fe(iat)=0.0
        fp(iat)=0.0
        fpp(iat)=0.0
        phie(iat)=0.0
        rho(iat)=0.0
        evacf(iat)=0.0
1840  continue
      
      rdrar = 1.0/drar 
      do 1870 iat=1,natoms
c
c     compute the contribution to rho(i) from particle j
c     Calculating the embedding energy F(Rho) in fe(i)
c      and Pair Energy contribution Phi(Rij) in phie(i)
c
ccdir$ novector
      do 1850 j=1,neigh(iat)
      rdist=rdis(iat,j,1)
      jat=rdis(iat,j,2)
      pij = rdist*rdrar + 1.0
      kij = pij
      kij = min0(kij,nrar-1)
      pij = pij - kij
c       make sure that p is less than 1.00
c       then if r is out of range, p = 1.00 and rho = last value of rhor
      pij = amin1(pij,1.0)
      ity=itype(iat)
      jty=itype(jat)
      rhojnei(iat,j) = ((rhorar3(kij,jty)*pij+
     $   rhorar2(kij,jty))*pij+
     $   rhorar1(kij,jty))*pij+
     $   rhorar(kij,jty)
      rhoinei(iat,j) = ((rhorar3(kij,ity)*pij+
     $   rhorar2(kij,ity))*pij+
     $   rhorar1(kij,ity))*pij+
     $   rhorar(kij,ity)
      rho(iat) = rho(iat) + rhojnei(iat,j)
      z2ij = ((z2rar3(kij,ity,jty)*pij+
     $          z2rar2(kij,ity,jty))*pij+
     $          z2rar1(kij,ity,jty))*pij+
     $          z2rar(kij,ity,jty)
      z2pij = (z2rar6(kij,ity,jty)*pij+z2rar5(kij,ity,jty))*pij+
     $         z2rar4(kij,ity,jty)
c     if(phitype) then
        phiij = z2ij
c     else
c       pij=1.0/rdist
c       phiij = z2ij * pij
c     endif
      phie(iat) = phie(iat) + 0.5 * phiij
1850  continue
ccdir$ vector
c        
1870  continue
c
c       store fsubi in e(i), the derivative of fsubi in fp
c       and the second derivative of fsubi in fpp
c
      rdrhoar = 1.0/drhoar
      do 2000 iat=1,natoms
      pij = rho(iat)*rdrhoar + 1.0
      kij = pij
      kij = max(1,min(kij,nrhoar-1))
      pij = pij - kij
      ity=itype(iat)
      fe(iat) = ((frhoar3(kij,ity)*pij+
     $        frhoar2(kij,ity))*pij+
     $        frhoar1(kij,ity))*pij+
     $        frhoar(kij,ity)
      fp(iat) = (frhoar6(kij,ity)*pij+
     $                  frhoar5(kij,ity))*pij+
     $                  frhoar4(kij,ity)
      fpp(iat) = (2.0*frhoar6(kij,ity)*pij+
     $                  frhoar5(kij,ity))/drhoar
2000  continue
c
c       store vacancy formation energy in evacf(i)
c
      rdrhoar = 1.0/drhoar
      do 3100 iat=1,natoms
        do 3110 j=1,neigh(iat)
        jat=rdis(iat,j,2)
        pij = (rho(jat)-rhoinei(iat,j))*rdrhoar + 1.0
        kij = pij
        kij = max(1,min(kij,nrhoar-1))
        pij = pij - kij
        ity=itype(iat)
        jty=itype(jat)
        evacf(iat) = evacf(iat) + ( fe(jat) - 
     $      (((frhoar3(kij,jty)*pij+
     $        frhoar2(kij,jty))*pij+
     $        frhoar1(kij,jty))*pij+
     $        frhoar(kij,jty))  )
3110    continue
      e(iat) = phie(iat) + fe(iat)
3100  continue

      do 3120 iat=1,natoms
      evacf(iat) = evacf(iat) + phie(iat)
3120  continue
      

      do 3150 iat=1,natoms
        ca(iat) = 0.0
        do 3151 cij=1,7
         econ(cij,iat) = 0.0
         bulkmod(cij,iat) = 0.0
         uu(cij) = 0.0
         ww(cij) = 0.0
         vv(cij) = 0.0
3151    continue
         do 3500 j=1,neigh(iat)
           jat=rdis(iat,j,2)
           do 3510 kc=1,3
3510       diski(kc) = rdisd(iat,j,kc)
c          call elascon(phitype,diski,jat,iat,uu,ww,vv,bulkmod) 
           call elastic(phitype,diski,jat,iat,uu,ww,vv,bulkmod) 
3500     continue
        econ(1,iat) = 0.5*uu(1) + fp(iat) *
     .                   ww(1) + fpp(iat) * vv(1)*vv(1)
        econ(2,iat) = 0.5*uu(2) + fp(iat) *
     .                   ww(2) + fpp(iat) * vv(1)*vv(2)
        econ(3,iat) = 0.5*uu(3) + fp(iat) *
     .                   ww(3) + fpp(iat) * vv(3)*vv(3)
        econ(4,iat) = 0.5*uu(4) + fp(iat) *
     .                   ww(4) + fpp(iat) * vv(1)*vv(4)
        econ(5,iat) = 0.5*uu(5) + fp(iat) *
     .                   ww(5) + fpp(iat) * vv(4)*vv(4)
        econ(6,iat) = 0.5*uu(6) + fp(iat) *
     .                   ww(6) + fpp(iat) * vv(5)*vv(6)
        econ(7,iat) = 0.5*uu(7) + fp(iat) *
     .                   ww(7) + fpp(iat) * vv(7)*vv(7)
        ca(iat)= fpp(iat)*(vv(1)*vv(2)-vv(7)*vv(7))
        ca(iat)= ca(iat)/ev_joule
        econ(1,iat) = econ(1,iat)/ev_joule
        econ(2,iat) = econ(2,iat)/ev_joule
        econ(3,iat) = econ(3,iat)/ev_joule
        econ(4,iat) = econ(4,iat)/ev_joule
        econ(5,iat) = econ(5,iat)/ev_joule
        econ(6,iat) = econ(6,iat)/ev_joule
        econ(7,iat) = econ(7,iat)/ev_joule
        bulkmod(1,iat) = 0.5*bulkmod(1,iat) + fp(iat) * 
     .    bulkmod(2,iat) + fpp(iat) * bulkmod(3,iat)**2
        bulkmod(4,iat) = 0.5*bulkmod(4,iat) + fp(iat) * 
     .    bulkmod(5,iat) + fpp(iat) * bulkmod(6,iat)**2
        bulkmod(1,iat) = bulkmod(1,iat)/ev_joule
        bulkmod(4,iat) = bulkmod(4,iat)/ev_joule
3150     continue

        omega0fcc = lata0*lata0*lata0/4.d0
        omega0bcc = lata0*lata0*lata0/2.d0
        omega0hcp = 2.*sqrt3*hcpc*lata0*lata0/4.d0
        omega0 = latis*latis*latis/4.d0
        c11(1)=0.0
        c12(1)=0.0
        c13(1)=0.0
        c33(1)=0.0
        c44(1)=0.0
        c55(1)=0.0
        c66(1)=0.0
        c66ca(1)=0.0
        bbb(1)=0.0
        bbb2(1)=0.0
        bbb3(1)=0.0
        bbb4(1)=0.0
        embed1(1)=0.0
        embed2(1)=0.0
        vacform1(1)=0.0
        vacform2(1)=0.0
        vacnum1(1)=0
        vacnum2(1)=0
        do 3200 iat=1,natoms
        bbb2(1)=bbb2(1)+(econ(1,iat)+2*econ(2,iat))/(3*omega0)
        c11(1)=c11(1)+econ(1,iat)/omega0
        c12(1)=c12(1)+econ(2,iat)/omega0
        c44(1)=c44(1)+econ(3,iat)/omega0
        c13(1)=c13(1)+econ(4,iat)/omega0
        c33(1)=c33(1)+econ(5,iat)/omega0
        c55(1)=c55(1)+econ(6,iat)/omega0
        c66(1)=c66(1)+econ(7,iat)/omega0
        c66ca(1)=c66ca(1)+(econ(2,iat)-ca(iat))/omega0
        bbb(1)=bbb(1)+bulkmod(1,iat)/(9*omega0)
        bbb3(1)=bbb3(1)+bulkmod(4,iat)/(9.*omega0)
        bbb4(1)=bbb4(1)+(2.*econ(1,iat)+2.*econ(2,iat)+
     .             4.*econ(4,iat)+econ(5,iat))/(9.*omega0)
        if( itype(iat) .eq. 2 ) then
          vacform1(1)=vacform1(1)+evacf(iat)
          embed1(1)=embed1(1)+e(iat)
          vacnum1(1)=vacnum1(1)+1
        else
          vacform2(1)=vacform2(1)+evacf(iat)
          embed2(1)=embed2(1)+e(iat)
          vacnum2(1)=vacnum2(1)+1
        endif
c       write(77,*)iat
c       write(77,*)econ(1,iat)/omega0,econ(2,iat)/omega0,
c    .econ(3,iat)/omega0
c       write(77,*)bulkmod(1,iat)/(9*omega0),bbb2
3200     continue
        bbb(1)=bbb(1)/natoms
        bbb2(1)=bbb2(1)/natoms
        bbb3(1)=bbb3(1)/natoms
        bbb4(1)=bbb4(1)/natoms
        c11(1)=c11(1)/natoms
        c12(1)=c12(1)/natoms
        c44(1)=c44(1)/natoms
        c13(1)=c13(1)/natoms
        c33(1)=c33(1)/natoms
        c55(1)=c55(1)/natoms
        c66(1)=c66(1)/natoms
        c66ca(1)=c66ca(1)/natoms
        if( vacnum1(1) .gt. 0 ) then
           vacform1(1)=vacform1(1)/vacnum1(1)
           embed1(1)=embed1(1)/vacnum1(1)
        endif
        if( vacnum2(1) .gt. 0 ) then
           vacform2(1)=vacform2(1)/vacnum2(1)
           embed2(1)=embed2(1)/vacnum2(1)
        endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         Calculation of Diatomic Strength and         c
c          Diatomic Length for Alloys                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if( ntypes .eq. 2 ) then
        t1=1
        t2=2
      else
        t1=2
        t2=3
      endif
      call calcdiatom(t1,t2,decal,recal)
      if (prntflag.eq.1) then
       write(*,*)'de:',decal,'re:',recal
      endif
      datome(1)=decal
      datomr(1)=recal
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       if (prntflag.eq.1) then
c*****************************************************
c*   Calculation Ni3Al Alloy Stacking Fault and      *
c*           Antiphase Boundary Energies             *
c*            (100), (111) SISF and APB              *
c*****************************************************
       call calcper()
       call calctote(etot1,neigh,rdis)
       shifts(1)=abs(rv(1,stacks(2)+1)-rv(1,1))
       shifts(2)=abs(rv(2,stacks(2)+1)-rv(2,1))
       shifts(3)=0.0
       zplanes=stacks(3)/2
       if (prntflag.eq.1) then
        print*,'Shifting',zplanes,'planes of FCC100 to x:',
     .shifts(1),'and y:',shifts(2)
        print*,'area:',perlen(1)*perlen(2)
       endif
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms6:',rmssum
        endif
         return
       endif
       call calctote(etot2,neigh,rdis)
c      vmdfile='ni3al-100.lammpstrj'
c      call writevmd(stacks,vmdfile)
       apb100(1)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       apb100(1)=apb100(1)/2.0
       apb100(1)=10000*apb100(1)/ev_joule
       if (prntflag.eq.1) then
        print*,'APB(100):',apb100(1)
       endif
       l1=3
       if( ntypes .eq. 2 ) then
         t1=2
         t2=1
       else
         t1=3
         t2=2
       endif
       call loadfcc111(t1,t2,l1,latis,stacks)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms7:',rmssum
        endif
         return
       endif
c      call calcper()
       call calctote(etot1,neigh,rdis)
c      print*,'L1 Ecoh for 100 stacking:',etot1/natoms
c      vmdfile='Ni3Al-L12-111.lammpstrj'
c      call writevmd(stacks,vmdfile)
c      vmdfile='ni3al.lammpstrj'
c      call writevmd(stacks,vmdfile)
       shifts(1)=0.0
       shifts(2)=abs(rv(2,1)-rv(2,2))
       shifts(3)=0.0
       zplanes=4
       if (prntflag.eq.1) then
        print*,'Shifting',zplanes,'planes of FCC111 to x:',
     .shifts(1),'and y:',shifts(2)
        print*,'area:',perlen(1)*perlen(2)
       endif
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms8:',rmssum
        endif
         return
       endif
       call calctote(etot2,neigh,rdis)
       apb111(1)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       apb111(1)=apb111(1)/2.0
       apb111(1)=10000*apb111(1)/ev_joule
       if (prntflag.eq.1) then
        print*,'APB(111):',apb111(1)
       endif
c      vmdfile='ni3al-apb.lammpstrj'
c      call writevmd(stacks,vmdfile)
       stacks(3)=stacks(3)-1
       call loadfcc111(t1,t2,l1,latis,stacks)
c      call calcper()
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms9:',rmssum
        endif
         return
       endif
       call calctote(etot1,neigh,rdis)
       xshift = rv(1,1) - rv(1,stacks(2)+2)
       yshift = rv(2,1) - rv(2,stacks(2)+2)
       shifts(1)=2.0*abs(xshift)/3.0
       shifts(2)=2.0*abs(yshift)/3.0
       shifts(3)=0.0
        zplanes=4
       if (prntflag.eq.1) then
        print*,'Shifting',zplanes,'planes of FCC111 to x:',
     .shifts(1),'and y:',shifts(2)
        print*,'area:',perlen(1)*perlen(2)
       endif
       call shiftplanes(stacks,zplanes,shifts)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms10:',rmssum
        endif
         return
       endif
       call calctote(etot2,neigh,rdis)
       sisf111(1)=dabs(etot2-etot1)/dabs(perlen(1)*perlen(2))
       sisf111(1)=sisf111(1)/2.0
       sisf111(1)=10000*sisf111(1)/ev_joule
       if (prntflag.eq.1) then
        print*,'SISF(111):',sisf111(1)
       endif
       stacks(3)=stacks(3)+1
c      vmdfile='ni3al-sisf.lammpstrj'
c      call writevmd(stacks,vmdfile)
       endif

       endif
      

       if( ( ntypes .lt. 3 .and. ( at1 .eq. 1 .and. at2 .eq. 3 ) )
     $ .or. ( ntypes .gt. 2 ) ) then
c      ---------------------------------------------               
c      AlCo B2 BCC
       if (prntflag.eq.1) then
        write(*,*)'------- AlCo BCC B2 -------'
       endif
       if (prmid.gt.0) then
        write(*,*)'PROC:',procnum,' part:',prmid,' AlCo BCC B2'
       endif
       b2=1
       if( ntypes .eq. 2 ) then
        t1=1
        t2=2
       else
        t1=1
        t2=3
       endif
       call loadbcc100(t1,t2,b2,lata0,stacks)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms11:',rmssum
        endif
         return
       endif
       if (prntflag.eq.1) then
        write(*,*)'AlCo Test a_0:',lata0
       endif
       call calclatcon(minenergy,latis,lata0,neigh,rdis,iter,prntflag)
       b2alat(2)=latis
       if (prntflag.eq.1) then
        write(*,*)'AlCo B2 a_0:',latis
       endif
       call loadbcc100(t1,t2,b2,latis,stacks)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms12:',rmssum
        endif
         return
       endif
       call calctote(etot1,neigh,rdis)
       etot1=etot1/natoms
       b2ecoh(2)=etot1
       if (prntflag.eq.1) then
        print*,'AlCo B2 BCC Ecoh:',etot1
       endif
       tn=2
       typs(1)=t1
       typs(2)=t2
       call calctotetype(tn,typs,tnum,embed,etot1,neigh,rdis)
       etot1=etot1/natoms
       if (prntflag.eq.1) then
        print*,'AlCo B2 Ecoh:',etot1
        print*,'AlCo B2 Embed1:',embed(1),'num:',tnum(1)
        print*,'AlCo B2 Embed2:',embed(2),'num:',tnum(2)
       endif
c      vmdfile='AlCo-B2.lammpstrj'
c      call writevmd(stacks,vmdfile)
        
c   *****************************************************
c   * EAM fcn file ready to use in Sublimation Energy,  *
c   *     Vacancy Formation Energy, Bulk Modulus,       *
c   *     and Elastic Constant Calculations for FCC     *
c   *****************************************************
    
      do 1841 iat=1,natoms
        e(iat)=0.0
        fe(iat)=0.0
        fp(iat)=0.0
        fpp(iat)=0.0
        phie(iat)=0.0
        rho(iat)=0.0
        evacf(iat)=0.0
1841  continue
      
      rdrar = 1.0/drar 
      do 1871 iat=1,natoms
c
c     compute the contribution to rho(i) from particle j
c     Calculating the embedding energy F(Rho) in fe(i)
c      and Pair Energy contribution Phi(Rij) in phie(i)
c
ccdir$ novector
      do 1851 j=1,neigh(iat)
      rdist=rdis(iat,j,1)
      jat=rdis(iat,j,2)
      pij = rdist*rdrar + 1.0
      kij = pij
      kij = min0(kij,nrar-1)
      pij = pij - kij
c       make sure that p is less than 1.00
c       then if r is out of range, p = 1.00 and rho = last value of rhor
      pij = amin1(pij,1.0)
      ity=itype(iat)
      jty=itype(jat)
      rhojnei(iat,j) = ((rhorar3(kij,jty)*pij+
     $   rhorar2(kij,jty))*pij+
     $   rhorar1(kij,jty))*pij+
     $   rhorar(kij,jty)
      rhoinei(iat,j) = ((rhorar3(kij,ity)*pij+
     $   rhorar2(kij,ity))*pij+
     $   rhorar1(kij,ity))*pij+
     $   rhorar(kij,ity)
      rho(iat) = rho(iat) + rhojnei(iat,j)
      z2ij = ((z2rar3(kij,ity,jty)*pij+
     $          z2rar2(kij,ity,jty))*pij+
     $          z2rar1(kij,ity,jty))*pij+
     $          z2rar(kij,ity,jty)
      z2pij = (z2rar6(kij,ity,jty)*pij+z2rar5(kij,ity,jty))*pij+
     $         z2rar4(kij,ity,jty)
c     if(phitype) then
        phiij = z2ij
c     else
c       pij=1.0/rdist
c       phiij = z2ij * pij
c     endif
      phie(iat) = phie(iat) + 0.5 * phiij
1851  continue
ccdir$ vector
c        
1871  continue
c
c       store fsubi in e(i), the derivative of fsubi in fp
c       and the second derivative of fsubi in fpp
c
      rdrhoar = 1.0/drhoar
      do 2001 iat=1,natoms
      pij = rho(iat)*rdrhoar + 1.0
      kij = pij
      kij = max(1,min(kij,nrhoar-1))
      pij = pij - kij
      ity=itype(iat)
      fe(iat) = ((frhoar3(kij,ity)*pij+
     $        frhoar2(kij,ity))*pij+
     $        frhoar1(kij,ity))*pij+
     $        frhoar(kij,ity)
      fp(iat) = (frhoar6(kij,ity)*pij+
     $                  frhoar5(kij,ity))*pij+
     $                  frhoar4(kij,ity)
      fpp(iat) = (2.0*frhoar6(kij,ity)*pij+
     $                  frhoar5(kij,ity))/drhoar
2001  continue
c
c       store vacancy formation energy in evacf(i)
c
      rdrhoar = 1.0/drhoar
      do 3101 iat=1,natoms
        do 3111 j=1,neigh(iat)
        jat=rdis(iat,j,2)
        pij = (rho(jat)-rhoinei(iat,j))*rdrhoar + 1.0
        kij = pij
        kij = max(1,min(kij,nrhoar-1))
        pij = pij - kij
        ity=itype(iat)
        jty=itype(jat)
        evacf(iat) = evacf(iat) + ( fe(jat) - 
     $      (((frhoar3(kij,jty)*pij+
     $        frhoar2(kij,jty))*pij+
     $        frhoar1(kij,jty))*pij+
     $        frhoar(kij,jty))  )
3111    continue
      e(iat) = phie(iat) + fe(iat)
3101  continue

      do 3121 iat=1,natoms
      evacf(iat) = evacf(iat) + phie(iat)
3121  continue
      

      do 3153 iat=1,natoms
        ca(iat) = 0.0
        do 3152 cij=1,7
         econ(cij,iat) = 0.0
         bulkmod(cij,iat) = 0.0
         uu(cij) = 0.0
         ww(cij) = 0.0
         vv(cij) = 0.0
3152    continue
         do 3501 j=1,neigh(iat)
           jat=rdis(iat,j,2)
           do 3511 kc=1,3
3511       diski(kc) = rdisd(iat,j,kc)
c          call elascon(phitype,diski,jat,iat,uu,ww,vv,bulkmod) 
           call elastic(phitype,diski,jat,iat,uu,ww,vv,bulkmod) 
3501     continue
        econ(1,iat) = 0.5*uu(1) + fp(iat) *
     .                   ww(1) + fpp(iat) * vv(1)*vv(1)
        econ(2,iat) = 0.5*uu(2) + fp(iat) *
     .                   ww(2) + fpp(iat) * vv(1)*vv(2)
        econ(3,iat) = 0.5*uu(3) + fp(iat) *
     .                   ww(3) + fpp(iat) * vv(3)*vv(3)
        econ(4,iat) = 0.5*uu(4) + fp(iat) *
     .                   ww(4) + fpp(iat) * vv(1)*vv(4)
        econ(5,iat) = 0.5*uu(5) + fp(iat) *
     .                   ww(5) + fpp(iat) * vv(4)*vv(4)
        econ(6,iat) = 0.5*uu(6) + fp(iat) *
     .                   ww(6) + fpp(iat) * vv(5)*vv(6)
        econ(7,iat) = 0.5*uu(7) + fp(iat) *
     .                   ww(7) + fpp(iat) * vv(7)*vv(7)
        ca(iat)= fpp(iat)*(vv(1)*vv(2)-vv(7)*vv(7))
        ca(iat)= ca(iat)/ev_joule
        econ(1,iat) = econ(1,iat)/ev_joule
        econ(2,iat) = econ(2,iat)/ev_joule
        econ(3,iat) = econ(3,iat)/ev_joule
        econ(4,iat) = econ(4,iat)/ev_joule
        econ(5,iat) = econ(5,iat)/ev_joule
        econ(6,iat) = econ(6,iat)/ev_joule
        econ(7,iat) = econ(7,iat)/ev_joule
        bulkmod(1,iat) = 0.5*bulkmod(1,iat) + fp(iat) * 
     .    bulkmod(2,iat) + fpp(iat) * bulkmod(3,iat)**2
        bulkmod(4,iat) = 0.5*bulkmod(4,iat) + fp(iat) * 
     .    bulkmod(5,iat) + fpp(iat) * bulkmod(6,iat)**2
        bulkmod(1,iat) = bulkmod(1,iat)/ev_joule
        bulkmod(4,iat) = bulkmod(4,iat)/ev_joule
3153     continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         Calculation of Diatomic Strength and         c
c          Diatomic Length for Alloys                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if( ntypes .eq. 2 ) then
        t1=1
        t2=2
      else
        t1=1
        t2=3
      endif
      call calcdiatom(t1,t2,decal,recal)
      if (prntflag.eq.1) then
         write(*,*)'de:',decal,'re:',recal
      endif
      datome(2)=decal
      datomr(2)=recal
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     ------------------------------------------------------------
c       Generating bulk positions for FCC structure
c     ------------------------------------------------------------
c      l1fcc=0
c      t1=2
c      t2=3
c      call loadfcc100(t1,t2,l1fcc,lata0,stacks)
c      call loadneigh(rcutoff,neigh,rdis,rdisd)
c     ------------------------------------------------------------

c*****************************************************

        omega0fcc = lata0*lata0*lata0/4.d0
        omega0bcc = lata0*lata0*lata0/2.d0
        omega0hcp = 2.*sqrt3*hcpc*lata0*lata0/4.d0
        omega0 = latis*latis*latis/2.d0
        c11(2)=0.0
        c12(2)=0.0
        c13(2)=0.0
        c33(2)=0.0
        c44(2)=0.0
        c55(2)=0.0
        c66(2)=0.0
        c66ca(2)=0.0
        bbb(2)=0.0
        bbb2(2)=0.0
        bbb3(2)=0.0
        bbb4(2)=0.0
        embed1(2)=0.0
        embed2(2)=0.0
        vacform1(2)=0.0
        vacform2(2)=0.0
        vacnum1(2)=0
        vacnum2(2)=0
        do 3201 iat=1,natoms
        bbb2(2)=bbb2(2)+(econ(1,iat)+2*econ(2,iat))/(3*omega0)
        c11(2)=c11(2)+econ(1,iat)/omega0
        c12(2)=c12(2)+econ(2,iat)/omega0
        c44(2)=c44(2)+econ(3,iat)/omega0
        c13(2)=c13(2)+econ(4,iat)/omega0
        c33(2)=c33(2)+econ(5,iat)/omega0
        c55(2)=c55(2)+econ(6,iat)/omega0
        c66(2)=c66(2)+econ(7,iat)/omega0
        c66ca(2)=c66ca(2)+(econ(2,iat)-ca(iat))/omega0
        bbb(2)=bbb(2)+bulkmod(1,iat)/(9*omega0)
        bbb3(2)=bbb3(2)+bulkmod(4,iat)/(9.*omega0)
        bbb4(2)=bbb4(2)+(2.*econ(1,iat)+2.*econ(2,iat)+
     .             4.*econ(4,iat)+econ(5,iat))/(9.*omega0)
        if( itype(iat) .eq. 1 ) then
          vacform1(2)=vacform1(2)+evacf(iat)
          embed1(2)=embed1(2)+e(iat)
          vacnum1(2)=vacnum1(2)+1
        else
          vacform2(2)=vacform2(2)+evacf(iat)
          embed2(2)=embed2(2)+e(iat)
          vacnum2(2)=vacnum2(2)+1
        endif
c       write(77,*)iat
c       write(77,*)econ(1,iat)/omega0,econ(2,iat)/omega0,
c    .econ(3,iat)/omega0
c       write(77,*)bulkmod(1,iat)/(9*omega0),bbb2
3201     continue
        bbb(2)=bbb(2)/natoms
        bbb2(2)=bbb2(2)/natoms
        bbb3(2)=bbb3(2)/natoms
        bbb4(2)=bbb4(2)/natoms
        c11(2)=c11(2)/natoms
        c12(2)=c12(2)/natoms
        c44(2)=c44(2)/natoms
        c13(2)=c13(2)/natoms
        c33(2)=c33(2)/natoms
        c55(2)=c55(2)/natoms
        c66(2)=c66(2)/natoms
        c66ca(2)=c66ca(2)/natoms
        if( vacnum1(2) .gt. 0 ) then
           vacform1(2)=vacform1(2)/vacnum1(2)
           embed1(2)=embed1(2)/vacnum1(2)
        endif
        if( vacnum2(2) .gt. 0 ) then
           vacform2(2)=vacform2(2)/vacnum2(2)
           embed2(2)=embed2(2)/vacnum2(2)
        endif
c      -------------------------------------------------
       endif

       if ( ( ntypes .lt. 3 .and. ( at1 .eq. 1 .and. at2 .eq. 2 ))
     $ .or. ( ntypes .gt. 2 ) )then
c      -------------------------------------------------       
c      CoNi HEX
       if (prntflag.eq.1) then
        write(*,*)'------- CoNi HEX  -------'
       endif
        write(*,*)'PROC:',procnum,' part:',prmid,' CoNi HEX'
       l1=1
       if( ntypes .eq. 2 ) then
        t1=1
        t2=2
       else
        t1=1
        t2=2
       endif
       hcpc=caratio*lata0
c      if (prntflag.eq.1) then
        write(*,*)'CoNi Test c/a:',caratio
        write(*,*)'CoNi Test a_0:',lata0
c      endif
       call loadhcp111(t1,t2,l1,lata0,caratio,stacks)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms13:',rmssum
        endif
         return
       endif
c      if (prntflag.eq.1) then
c       write(*,*)'CoNi Test a_0:',lata0
c      endif
       rmssum=10000.0
       call calclatca(minenergy,latis,lata0,c,caratio,neigh,rdis,rdisd,
     &rmssum,prntflag)
       if (rmssum.gt.99999.0) then
          rmssum=100000.0
c         return
          write(*,*)'calclatca: no min latis,c/a'
          latis=2.0
          c=1.0
       endif
       hexalat(1)=latis
       hexca(1)=c
c      if (prntflag.eq.1) then
        write(*,*)'CoNi HEX rmssum:',rmssum
        write(*,*)'CoNi HEX a_0:',latis
        write(*,*)'CoNi HEX c/a:',c
c      endif
       hcpc=c*latis
       call loadhcp111(t1,t2,l1,latis,c,stacks)
       call loadneigh(rcutoff,neigh,rdis,rdisd,istop,rmssum,prntflag)
       if(istop.eq.1) then 
        if (prntflag.eq.1) then
         write(*,*)'rms14:',rmssum
        endif
         return
       endif
c      call calcminehex(etot1,latis,lata0,c,c0,neigh,rdis,rdisd)
       call calctote(etot1,neigh,rdis)
       etot1=etot1/natoms
       hexecoh(1)=etot1
       if (prntflag.eq.1) then
        print*,'CoNi HEX Ecoh:',etot1
       endif
       tn=2
       typs(1)=t1
       typs(2)=t2
       call calctotetype(tn,typs,tnum,embed,etot1,neigh,rdis)
       etot1=etot1/natoms
       if (prntflag.eq.1) then
        print*,'CoNi HEX Ecoh:',etot1
        print*,'CoNi HEX Embed1:',embed(1),'num:',tnum(1)
        print*,'CoNi HEX Embed2:',embed(2),'num:',tnum(2)
       endif
c      vmdfile='CoNi-HEX.lammpstrj'
c      call writevmd(stacks,vmdfile)
        
c   *****************************************************
c   * EAM fcn file ready to use in Sublimation Energy,  *
c   *     Vacancy Formation Energy, Bulk Modulus,       *
c   *     and Elastic Constant Calculations for FCC     *
c   *****************************************************
    
      do 1842 iat=1,natoms
        e(iat)=0.0
        fe(iat)=0.0
        fp(iat)=0.0
        fpp(iat)=0.0
        phie(iat)=0.0
        rho(iat)=0.0
        evacf(iat)=0.0
1842  continue
      
      rdrar = 1.0/drar 
      do 1872 iat=1,natoms
c
c     compute the contribution to rho(i) from particle j
c     Calculating the embedding energy F(Rho) in fe(i)
c      and Pair Energy contribution Phi(Rij) in phie(i)
c
ccdir$ novector
      do 1852 j=1,neigh(iat)
      rdist=rdis(iat,j,1)
      jat=rdis(iat,j,2)
      pij = rdist*rdrar + 1.0
      kij = pij
      kij = min0(kij,nrar-1)
      pij = pij - kij
c       make sure that p is less than 1.00
c       then if r is out of range, p = 1.00 and rho = last value of rhor
      pij = amin1(pij,1.0)
      ity=itype(iat)
      jty=itype(jat)
      rhojnei(iat,j) = ((rhorar3(kij,jty)*pij+
     $   rhorar2(kij,jty))*pij+
     $   rhorar1(kij,jty))*pij+
     $   rhorar(kij,jty)
      rhoinei(iat,j) = ((rhorar3(kij,ity)*pij+
     $   rhorar2(kij,ity))*pij+
     $   rhorar1(kij,ity))*pij+
     $   rhorar(kij,ity)
      rho(iat) = rho(iat) + rhojnei(iat,j)
      z2ij = ((z2rar3(kij,ity,jty)*pij+
     $          z2rar2(kij,ity,jty))*pij+
     $          z2rar1(kij,ity,jty))*pij+
     $          z2rar(kij,ity,jty)
      z2pij = (z2rar6(kij,ity,jty)*pij+z2rar5(kij,ity,jty))*pij+
     $         z2rar4(kij,ity,jty)
c     if(phitype) then
        phiij = z2ij
c     else
c       pij=1.0/rdist
c       phiij = z2ij * pij
c     endif
      phie(iat) = phie(iat) + 0.5 * phiij
1852  continue
ccdir$ vector
c        
1872  continue
c
c       store fsubi in e(i), the derivative of fsubi in fp
c       and the second derivative of fsubi in fpp
c
      rdrhoar = 1.0/drhoar
      do 2002 iat=1,natoms
      pij = rho(iat)*rdrhoar + 1.0
      kij = pij
      kij = max(1,min(kij,nrhoar-1))
      pij = pij - kij
      ity=itype(iat)
      fe(iat) = ((frhoar3(kij,ity)*pij+
     $        frhoar2(kij,ity))*pij+
     $        frhoar1(kij,ity))*pij+
     $        frhoar(kij,ity)
      fp(iat) = (frhoar6(kij,ity)*pij+
     $                  frhoar5(kij,ity))*pij+
     $                  frhoar4(kij,ity)
      fpp(iat) = (2.0*frhoar6(kij,ity)*pij+
     $                  frhoar5(kij,ity))/drhoar
2002  continue
c
c       store vacancy formation energy in evacf(i)
c
      rdrhoar = 1.0/drhoar
      do 3102 iat=1,natoms
        do 3112 j=1,neigh(iat)
        jat=rdis(iat,j,2)
        pij = (rho(jat)-rhoinei(iat,j))*rdrhoar + 1.0
        kij = pij
        kij = max(1,min(kij,nrhoar-1))
        pij = pij - kij
        ity=itype(iat)
        jty=itype(jat)
        evacf(iat) = evacf(iat) + ( fe(jat) - 
     $      (((frhoar3(kij,jty)*pij+
     $        frhoar2(kij,jty))*pij+
     $        frhoar1(kij,jty))*pij+
     $        frhoar(kij,jty))  )
3112    continue
      e(iat) = phie(iat) + fe(iat)
3102  continue

      do 3122 iat=1,natoms
      evacf(iat) = evacf(iat) + phie(iat)
3122  continue
      

      do 3154 iat=1,natoms
        ca(iat) = 0.0
        do 3155 cij=1,7
         econ(cij,iat) = 0.0
         bulkmod(cij,iat) = 0.0
         uu(cij) = 0.0
         ww(cij) = 0.0
         vv(cij) = 0.0
3155    continue
         do 3502 j=1,neigh(iat)
           jat=rdis(iat,j,2)
           do 3512 kc=1,3
3512       diski(kc) = rdisd(iat,j,kc)
           call elastic(phitype,diski,jat,iat,uu,ww,vv,bulkmod) 
3502     continue
        econ(1,iat) = 0.5*uu(1) + fp(iat) *
     .                   ww(1) + fpp(iat) * vv(1)*vv(1)
        econ(2,iat) = 0.5*uu(2) + fp(iat) *
     .                   ww(2) + fpp(iat) * vv(1)*vv(2)
        econ(3,iat) = 0.5*uu(3) + fp(iat) *
     .                   ww(3) + fpp(iat) * vv(3)*vv(3)
        econ(4,iat) = 0.5*uu(4) + fp(iat) *
     .                   ww(4) + fpp(iat) * vv(1)*vv(4)
        econ(5,iat) = 0.5*uu(5) + fp(iat) *
     .                   ww(5) + fpp(iat) * vv(4)*vv(4)
        econ(6,iat) = 0.5*uu(6) + fp(iat) *
     .                   ww(6) + fpp(iat) * vv(5)*vv(6)
        econ(7,iat) = 0.5*uu(7) + fp(iat) *
     .                   ww(7) + fpp(iat) * vv(7)*vv(7)
        ca(iat)= fpp(iat)*(vv(1)*vv(2)-vv(7)*vv(7))
        ca(iat)= ca(iat)/ev_joule
        econ(1,iat) = econ(1,iat)/ev_joule
        econ(2,iat) = econ(2,iat)/ev_joule
        econ(3,iat) = econ(3,iat)/ev_joule
        econ(4,iat) = econ(4,iat)/ev_joule
        econ(5,iat) = econ(5,iat)/ev_joule
        econ(6,iat) = econ(6,iat)/ev_joule
        econ(7,iat) = econ(7,iat)/ev_joule
        bulkmod(1,iat) = 0.5*bulkmod(1,iat) + fp(iat) * 
     .    bulkmod(2,iat) + fpp(iat) * bulkmod(3,iat)**2
        bulkmod(4,iat) = 0.5*bulkmod(4,iat) + fp(iat) * 
     .    bulkmod(5,iat) + fpp(iat) * bulkmod(6,iat)**2
        bulkmod(1,iat) = bulkmod(1,iat)/ev_joule
        bulkmod(4,iat) = bulkmod(4,iat)/ev_joule
3154     continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         Calculation of Diatomic Strength and         c
c          Diatomic Length for Alloys                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if( ntypes .eq. 2 ) then
        t1=1
        t2=2
       else
        t1=1
        t2=2
       endif
      call calcdiatom(t1,t2,decal,recal)
      if (prntflag.eq.1) then
       write(*,*)'de:',decal,'re:',recal
      endif
      datome(3)=decal
      datomr(3)=recal
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     ------------------------------------------------------------
c       Generating bulk positions for FCC structure
c     ------------------------------------------------------------
c      l1fcc=0
c      t1=2
c      t2=3
c      call loadfcc100(t1,t2,l1fcc,lata0,stacks)
c      call loadneigh(rcutoff,neigh,rdis,rdisd)
c     ------------------------------------------------------------

c*****************************************************

        omega0fcc = lata0*lata0*lata0/4.d0
        omega0bcc = lata0*lata0*lata0/2.d0
        omega0hcp = 2.*sqrt3*hcpc*lata0*lata0/4.d0
        omega0 = sqrt3*hcpc*latis*latis/4.d0
        c11(3)=0.0
        c12(3)=0.0
        c13(3)=0.0
        c33(3)=0.0
        c44(3)=0.0
        c55(3)=0.0
        c66(3)=0.0
        c66ca(3)=0.0
        bbb(3)=0.0
        bbb2(3)=0.0
        bbb3(3)=0.0
        bbb4(3)=0.0
        embed1(3)=0.0
        embed2(3)=0.0
        vacform1(3)=0.0
        vacform2(3)=0.0
        vacnum1(3)=0
        vacnum2(3)=0
        do 3202 iat=1,natoms
        bbb2(3)=bbb2(3)+(econ(1,iat)+2*econ(2,iat))/(3*omega0)
        c11(3)=c11(3)+econ(1,iat)/omega0
        c12(3)=c12(3)+econ(2,iat)/omega0
        c44(3)=c44(3)+econ(3,iat)/omega0
        c13(3)=c13(3)+econ(4,iat)/omega0
        c33(3)=c33(3)+econ(5,iat)/omega0
        c55(3)=c55(3)+econ(6,iat)/omega0
        c66(3)=c66(3)+econ(7,iat)/omega0
        c66ca(3)=c66ca(3)+(econ(2,iat)-ca(iat))/omega0
        bbb(3)=bbb(3)+bulkmod(1,iat)/(9*omega0)
        bbb3(3)=bbb3(3)+bulkmod(4,iat)/(9.*omega0)
        bbb4(3)=bbb4(3)+(2.*econ(1,iat)+2.*econ(2,iat)+
     .             4.*econ(4,iat)+econ(5,iat))/(9.*omega0)
        if( itype(iat) .eq. 1 ) then
          vacform1(3)=vacform1(3)+evacf(iat)
          embed1(3)=embed1(3)+e(iat)
          vacnum1(3)=vacnum1(3)+1
        else
          vacform2(3)=vacform2(3)+evacf(iat)
          embed2(3)=embed2(3)+e(iat)
          vacnum2(3)=vacnum2(3)+1
        endif
c       write(77,*)iat
c       write(77,*)econ(1,iat)/omega0,econ(2,iat)/omega0,
c    .econ(3,iat)/omega0
c       write(77,*)bulkmod(1,iat)/(9*omega0),bbb2
3202     continue
        bbb(3)=bbb(3)/natoms
        bbb2(3)=bbb2(3)/natoms
        bbb3(3)=bbb3(3)/natoms
        bbb4(3)=bbb4(3)/natoms
        c11(3)=c11(3)/natoms
        c12(3)=c12(3)/natoms
        c44(3)=c44(3)/natoms
        c13(3)=c13(3)/natoms
        c33(3)=c33(3)/natoms
        c55(3)=c55(3)/natoms
        c66(3)=c66(3)/natoms
        c66ca(3)=c66ca(3)/natoms
        if( vacnum1(3) .gt. 0 ) then
           vacform1(3)=vacform1(3)/vacnum1(3)
           embed1(3)=embed1(3)/vacnum1(3)
        endif
        if( vacnum2(3) .gt. 0 ) then
           vacform2(3)=vacform2(3)/vacnum2(3)
           embed2(3)=embed2(3)/vacnum2(3)
        endif

      do 3203 i=1,3
      if( b2alat(i) .lt. 0.5 ) then
          b2alat(i)=0.d0
      endif
      if( l12alat(i) .lt. 0.5 ) then
          l12alat(i)=0.d0
      endif
      if( hexalat(i) .lt. 0.5 ) then
          hexalat(i)=0.d0
      endif
3203  continue

c     ----------------------------------------------
      endif
c     ----------------------------------------------
      fitnums=0
      if( ntypes .lt. 3 .and. ( at1 .eq. 2 .and. at2 .eq. 4 ) ) then
c     CuNi NiCu
         if (weightdata(1).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(1)*weightdata(1)
         fitvals(fitnums)=l10alat(1)*weightdata(1)
         endif
         if (weightdata(2).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(2)*weightdata(2)
         fitvals(fitnums)=-l10ecoh(1)*weightdata(2)
         endif
         if (weightdata(3).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(3)*weightdata(3)
         fitvals(fitnums)=-embed1(3)*weightdata(3)
         endif
         if (weightdata(4).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(4)*weightdata(4)
         fitvals(fitnums)=-embed2(3)*weightdata(4)
         endif
         if (weightdata(5).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(5)*weightdata(5)
         fitvals(fitnums)=l12alat(1)*weightdata(5)
         endif
         if (weightdata(6).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(6)*weightdata(6)
         fitvals(fitnums)=-l12ecoh(1)*weightdata(6)
         endif
         if (weightdata(7).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(7)*weightdata(7)
         fitvals(fitnums)=-embed1(1)*weightdata(7)
         endif
         if (weightdata(8).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(8)*weightdata(8)
         fitvals(fitnums)=-embed2(1)*weightdata(8)
         endif
         if (weightdata(9).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(9)*weightdata(9)
         fitvals(fitnums)=bbb(1)*weightdata(9)
         endif
         if (weightdata(10).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(10)*weightdata(10)
         fitvals(fitnums)=c11(1)*weightdata(10)
         endif
         if (weightdata(11).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(11)*weightdata(11)
         fitvals(fitnums)=c12(1)*weightdata(11)
         endif
         if (weightdata(12).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(12)*weightdata(12)
         fitvals(fitnums)=c44(1)*weightdata(12)
         endif
         if (weightdata(13).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(13)*weightdata(13)
         fitvals(fitnums)=-vacform1(1)*weightdata(13)
         endif
         if (weightdata(14).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(14)*weightdata(14)
         fitvals(fitnums)=-vacform2(1)*weightdata(14)
         endif
         if (weightdata(15).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(15)*weightdata(15)
         fitvals(fitnums)=l13alat(1)*weightdata(15)
         endif
         if (weightdata(16).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(16)*weightdata(16)
         fitvals(fitnums)=-l13ecoh(1)*weightdata(16)
         endif
         if (weightdata(17).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(17)*weightdata(17)
         fitvals(fitnums)=-embed1(4)*weightdata(17)
         endif
         if (weightdata(18).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(18)*weightdata(18)
         fitvals(fitnums)=-embed2(4)*weightdata(18)
         endif
         if (weightdata(19).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(19)*weightdata(19)
         fitvals(fitnums)=l11alat(1)*weightdata(19)
         endif
         if (weightdata(20).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(20)*weightdata(20)
         fitvals(fitnums)=-l11ecoh(1)*weightdata(20)
         endif
         if (weightdata(21).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(21)*weightdata(21)
         fitvals(fitnums)=-embed1(2)*weightdata(21)
         endif
         if (weightdata(22).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(22)*weightdata(22)
         fitvals(fitnums)=-embed2(2)*weightdata(22)
         endif
         if (weightdata(23).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(23)*weightdata(23)
         fitvals(fitnums)=bbb(2)*weightdata(23)
         endif
         if (weightdata(24).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(24)*weightdata(24)
         fitvals(fitnums)=c11(2)*weightdata(24)
         endif
         if (weightdata(25).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(25)*weightdata(25)
         fitvals(fitnums)=c12(2)*weightdata(25)
         endif
         if (weightdata(26).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(26)*weightdata(26)
         fitvals(fitnums)=c44(2)*weightdata(26)
         endif
         if (weightdata(27).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(27)*weightdata(27)
         fitvals(fitnums)=-vacform1(2)*weightdata(27)
         endif
         if (weightdata(28).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(28)*weightdata(28)
         fitvals(fitnums)=-vacform2(2)*weightdata(28)
         endif
         if (weightdata(29).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(29)*weightdata(29)
         fitvals(fitnums)=-datome(1)*weightdata(29)
         endif
         if (weightdata(30).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(30)*weightdata(30)
         fitvals(fitnums)=datomr(1)*weightdata(30)
         endif
         if (weightdata(31).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(31)*weightdata(31)
         fitvals(fitnums)=-l12ecoh(2)*weightdata(31)
         endif
         if (weightdata(32).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(32)*weightdata(32)
         fitvals(fitnums)=-l11ecoh(2)*weightdata(32)
         endif
         if (weightdata(33).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(33)*weightdata(33)
         fitvals(fitnums)=sisf111(1)*weightdata(33)
         endif
         if (weightdata(34).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(34)*weightdata(34)
         fitvals(fitnums)=apb100(1)*weightdata(34)
         endif
         if (weightdata(35).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(35)*weightdata(35)
         fitvals(fitnums)=apb111(1)*weightdata(35)
         endif
      endif
      if( ( ntypes .lt. 3 .and. ( at1 .eq. 2 .and. at2 .eq. 3 ) ) 
     $ .or. ( ntypes .gt. 2 ) ) then
         if (weightdata(1).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(1)*weightdata(1)
         fitvals(fitnums)=b2alat(1)*weightdata(1)
         endif
         if (weightdata(2).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(2)*weightdata(2)
         fitvals(fitnums)=-b2ecoh(1)*weightdata(2)
         endif
         if (weightdata(5).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(5)*weightdata(5)
         fitvals(fitnums)=l12alat(1)**weightdata(5)
         endif
         if (weightdata(6).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(6)*weightdata(6)
         fitvals(fitnums)=-l12ecoh(1)*weightdata(6)
         endif
         if (weightdata(9).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(9)*weightdata(9)
         fitvals(fitnums)=bbb(1)*weightdata(9)
         endif
         if (weightdata(10).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(10)*weightdata(10)
         fitvals(fitnums)=c11(1)*weightdata(10)
         endif
         if (weightdata(11).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(11)*weightdata(11)
         fitvals(fitnums)=c12(1)*weightdata(11)
         endif
         if (weightdata(12).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(12)*weightdata(12)
         fitvals(fitnums)=c44(1)*weightdata(12)
         endif
         if (weightdata(13).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(13)*weightdata(13)
         fitvals(fitnums)=-vacform1(1)*weightdata(13)
         endif
         if (weightdata(14).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(14)*weightdata(14)
         fitvals(fitnums)=-vacform2(1)*weightdata(14)
         endif
         if (weightdata(15).gt.0) then 
c        fitnums=fitnums+1
c        expvals(fitnums)=expdata(15)*weightdata(15)
c        fitvals(fitnums)=sisf111*weightdata(15)
         endif
         if (weightdata(16).gt.0) then 
c        fitnums=fitnums+1
c        expvals(fitnums)=expdata(16)*weightdata(16)
c        fitvals(fitnums)=apb100*weightdata(16)
         endif
         if (weightdata(17).gt.0) then 
c        fitnums=fitnums+1
c        expvals(fitnums)=expdata(17)*weightdata(17)
c        fitvals(fitnums)=apb111*weightdata(17)
         endif
       endif
       if( ( ntypes .lt. 3 .and. ( at1 .eq. 1 .and. at2 .eq. 3 ) )
     $ .or. ( ntypes .gt. 2 ) ) then
         if (weightdata(18).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(18)*weightdata(18)
         fitvals(fitnums)=b2alat(2)*weightdata(18)
         endif
         if (weightdata(19).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(19)*weightdata(19)
         fitvals(fitnums)=-b2ecoh(2)*weightdata(19)
         endif
         if (weightdata(20).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(20)*weightdata(20)
         fitvals(fitnums)=-embed1(2)*weightdata(20)
         endif
         if (weightdata(21).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(21)*weightdata(21)
         fitvals(fitnums)=-embed2(2)*weightdata(21)
         endif
         if (weightdata(22).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(22)*weightdata(22)
         fitvals(fitnums)=bbb(2)*weightdata(22)
         endif
         if (weightdata(23).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(23)*weightdata(23)
         fitvals(fitnums)=c11(2)*weightdata(23)
         endif
         if (weightdata(24).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(24)*weightdata(24)
         fitvals(fitnums)=c12(2)*weightdata(24)
         endif
         if (weightdata(25).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(25)*weightdata(25)
         fitvals(fitnums)=c44(2)*weightdata(25)
         endif
         if (weightdata(26).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(26)*weightdata(26)
         fitvals(fitnums)=-vacform1(2)*weightdata(26)
         endif
         if (weightdata(27).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(27)*weightdata(27)
         fitvals(fitnums)=-vacform2(2)*weightdata(27)
         endif
       endif
       if ( ( ntypes .lt. 3 .and. ( at1 .eq. 1 .and. at2 .eq. 2 ))
     $ .or. ( ntypes .gt. 2 ) )then
         if (weightdata(28).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(28)*weightdata(28)
         fitvals(fitnums)=hexalat(1)*weightdata(28)
         endif
         if (weightdata(29).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(29)*weightdata(29)
         fitvals(fitnums)=hexca(1)*weightdata(29)
         endif
         if (weightdata(30).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(30)*weightdata(30)
         fitvals(fitnums)=-hexecoh(1)*weightdata(30)
         endif
         if (weightdata(31).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(31)*weightdata(31)
         fitvals(fitnums)=-embed1(3)*weightdata(31)
         endif
         if (weightdata(32).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(32)*weightdata(32)
         fitvals(fitnums)=-embed2(3)*weightdata(32)
         endif
         if (weightdata(33).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(33)*weightdata(33)
         fitvals(fitnums)=c11(3)*weightdata(33)
         endif
         if (weightdata(34).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(34)*weightdata(34)
         fitvals(fitnums)=c12(3)*weightdata(34)
         endif
         if (weightdata(35).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(35)*weightdata(35)
         fitvals(fitnums)=c33(3)*weightdata(35)
         endif
         if (weightdata(36).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(36)*weightdata(36)
         fitvals(fitnums)=c44(3)*weightdata(36)
         endif
         if (weightdata(37).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(37)*weightdata(37)
         fitvals(fitnums)=-vacform1(3)*weightdata(37)
         endif
         if (weightdata(38).gt.0) then 
         fitnums=fitnums+1
         expvals(fitnums)=expdata(38)*weightdata(38)
         fitvals(fitnums)=-vacform2(3)*weightdata(38)
         endif
       endif
c      if (fitvals(7).lt.fitvals(6)) then
c           fitvals(6)=fitvals(6)*10.0
c           fitvals(7)=fitvals(7)*10.0
c      endif
c     ----------------------------------------------
       rmssum=0.0
       if (prntflag.eq.1) then
        write(*,*)'fitnums:',fitnums
       endif
       do 4000 i=1,fitnums
        if (fitvals(i).eq.0.0) then
            fitvals(i)=1000
        endif
        if (fitvals(i).lt.0.01) then
            fitvals(i)=fitvals(i)*100000
c           fitvals(i)=fitvals(i)*100
        endif
        errors(i)=expvals(i)-fitvals(i)
        errors(i)=errors(i)*errors(i)
        rmssum=rmssum+errors(i)
c       write(*,*)expvals(i),fitvals(i),errors(i)
4000  continue
      if (-datome(1).gt.datomr(1)) then
        rmssum=rmssum+100
      endif
      if (-vacform1(1).lt.1.0 .or.
     .    -vacform1(1).gt.2.0) then
        rmssum=rmssum+100
      endif
      if (-vacform2(1).lt.1.0 .or.
     .    -vacform2(1).gt.2.0) then
        rmssum=rmssum+100
      endif
      if (-vacform1(2).lt.1.0 .or.
     .    -vacform1(2).gt.2.0) then
        rmssum=rmssum+100
      endif
      if (-vacform2(2).lt.1.0 .or.
     .    -vacform2(2).gt.2.0) then
        rmssum=rmssum+100
      endif
c     if (datomr(1).lt.2.1) then
c       rmssum=rmssum+10
c     endif
c     if (-datome(1).gt.2.1) then
c       rmssum=rmssum+10
c     endif
c     if (c11(2).gt.2.055 .or. 
c    .    c11(2).lt.2.035 ) then
c       rmssum=rmssum+100
c     endif
c     if (c12(2).gt.1.41 .or. 
c    .    c12(2).lt.1.39 ) then
c       rmssum=rmssum+100
c     endif
c     if (c44(2).gt.1.01 .or. 
c    .    c44(2).lt.0.98 ) then
c       rmssum=rmssum+100
c     endif
      if (l12alat(1).gt.3.60 .or. 
     .    l12alat(1).lt.3.58 ) then
        rmssum=rmssum+64
      endif
      if (l11alat(1).gt.3.577 .or.
     .    l11alat(1).lt.3.56 ) then
        rmssum=rmssum+64
      endif
      if (-l11ecoh(1).gt.3.976 .or. 
     .    -l11ecoh(1).lt.3.968 ) then
        rmssum=rmssum+200
      endif
      if (-l12ecoh(1).gt.3.754 .or.
     .    -l12ecoh(1).lt.3.7476 ) then
        rmssum=rmssum+200
      endif
      rmssum=sqrt(abs(rmssum))
c     ----------------------------------------------
      if (prntflag.eq.1) then
      open(unit=70,file=fitvalues)
      if(ntypes .lt. 3 .and. (at1 .eq. 2 .and. at2 .eq. 4))then
      write(*,*)'--------------------'
      write(70,7278)l10alat(1),-l10ecoh(1),-embed1(3),-embed2(3),
     .l12alat(1),-l12ecoh(1),-embed1(1),-embed2(1),
     .bbb(1),bbb3(1),c11(1),c12(1),c13(1),c33(1),c44(1),c55(1),c66(1),
     .-vacform1(1),-vacform2(1),l13alat(1),-l13ecoh(1),
     .-embed1(4),-embed2(4),l11alat(1),-l11ecoh(1),
     .-embed1(2),-embed2(2),bbb(2),bbb3(2),c11(2),c12(2),
     .c13(2),c33(2),c44(2),c55(2),c66(2),-vacform1(2),-vacform2(2)
     .-datome(1),datomr(1),-l12ecoh(2),-l11ecoh(2),
     .sisf111(1),apb100(1),apb111(1)
      write(*,7277)l10alat(1),-l10ecoh(1),-embed1(3),-embed2(3),
     .l12alat(1),-l12ecoh(1),-embed1(1),-embed2(1),
     .bbb(1),c11(1),c12(1),c44(1),c13(1),c33(1),c55(1),c66(1),
     .-vacform1(1),-vacform2(1),bbb3(1)
      write(*,5257)l13alat(1),-l13ecoh(1),-embed1(4),-embed2(4),
     .l11alat(1),-l11ecoh(1),-embed1(2),-embed2(2),
     .bbb(2),c11(2),c12(2),c44(2),c13(2),c33(2),c55(2),c66(2),
     .-vacform1(2),-vacform2(2),bbb3(2)
      write(*,7275)-datome(1),datomr(1),-l12ecoh(2),-l11ecoh(2),
     .sisf111(1),apb100(1),apb111(1),sisf111(2),apb100(2),apb111(2),
     .sisf111(3),apb100(3),apb111(3),sisf111(4),apb100(4),apb111(4)
      else
      write(*,*)'--------------------'
      write(70,7878)b2alat(1),-b2ecoh(1),-embed1(4),-embed2(4),
     .l12alat(1),-l12ecoh(1),-embed1(1),-embed2(1),bbb(1),bbb3(1),
     .c11(1),c12(1),c13(1),c33(1),c44(1),c66(1),-datome(1),datomr(1),
     .-vacform1(1),-vacform2(1),sisf111(1),apb100(1),apb111(1),
     .b2alat(2),
     .-b2ecoh(2),-embed1(2),-embed2(2),bbb(2),bbb3(2),c11(2),c12(2),
     .c13(2),c33(2),c44(2),c66(2),-datome(2),datomr(2),-vacform1(2),
     .-vacform2(2),hexalat(1),hexca(1),-hexecoh(1),-embed1(3),
     .-embed2(3),bbb(3),bbb3(3),c11(3),c12(3),c13(3),c33(3),c44(3),
     .c66(3),-datome(3),datomr(3),-vacform1(3),-vacform2(3)
      write(*,7777)l12alat(1),-l12ecoh(1),-embed1(1),-embed2(1),
     .bbb(1),c11(1),c12(1),c44(1),c13(1),c33(1),c55(1),c66(1),
     .-vacform1(1),-vacform2(1),bbb3(1),-datome(1),datomr(1)
      write(*,7575)sisf111(1),apb100(1),apb111(1)
      write(*,5757)b2alat(1),-b2ecoh(1),-embed1(4),-embed2(4)
      write(*,7775)b2alat(2),-b2ecoh(2),-embed1(2),-embed2(2),
     .bbb(2),c11(2),c12(2),c44(2),c13(2),c33(2),c55(2),c66(2),
     .-vacform1(2),-vacform2(2),bbb3(2),-datome(2),datomr(2)
      write(*,7757)hexalat(1),hexca(1),-hexecoh(1),-embed1(3),
     .-embed2(3),bbb(3),c11(3),c12(3),c13(3),c33(3),c44(3),c55(3),
     .c66(3),-vacform1(3),-vacform2(3),bbb3(3),-datome(3),datomr(3)
      endif
7777    format((' L1 Ni3Al',/,' a_0:',f9.5,/,' Ecoh:',f9.5,/,
     .' Ecoh1:',f9.5,/,' Ecoh2:',f9.5,/,' bmod:',f9.5,/,
     .' c11:',f9.5,/,' c12:',f9.5,/,' c44:',f9.5,/,
     .' c13:',f9.5,/,' c33:',f9.5,/,' c55:',f9.5,/,
     .' c66:',f9.5,/,' Evf1:',f9.5,/,' Evf2:',f9.5,/,' Bmod:',f9.5,/,
     .' De:',f9.5,/,' Re',f9.5,/))
7878    format(56(E25.16,1x))
7575    format((' SISF(111):',f10.5,/,' ABP(100):',f10.5,/,
     .' ABP(111):',f10.5,/))
5757    format((/,' B2 NiAl',/,' a_0:',f9.5,/,' Ecoh:',f9.5,/,
     .' Ecoh1:',f9.5,/,' Ecoh2:',f9.5,/))
7775    format((/,' B2 AlCo',/,' a_0:',f9.5,/,' Ecoh:',f9.5,/,
     .' Ecoh1:',f9.5,/,' Ecoh2:',f9.5,/,' bmod:',f9.5,/,
     .' c11:',f9.5,/,' c12:',f9.5,/,' c44:',f9.5,/,
     .' c13:',f9.5,/,' c33:',f9.5,/,' c55:',f9.5,/,
     .' c66:',f9.5,/,' Evf1:',f9.5,/,' Evf2:',f9.5,/,' Bmod:',f9.5,/,
     .' De:',f9.5,/,' Re',f9.5,/))
7757    format((/,' HEX CoNi',/,' a_0:',f9.5,/,'c/a:',f9.5,/,
     .' Ecoh:',f9.5,/,' Ecoh1:',f9.5,/,' Ecoh2:',f9.5,/,
     .' bmod:',f10.5,/,' c11:',f10.5,/,' c12:',f10.5,/,' c13:',f10.5,/,
     .' c33:',f10.5,/,' c44:',f10.5,/,' c55:',f10.5,/,' c66:',f10.5,/,
     .' Evf1:',f9.5,/,' Evf2:',f9.5,/,' Bmod:',f10.5,/,' De:',f9.5,/,
     .' Re',f9.5,/))
7277    format((' L10 CuNi',/,' a_0:',f10.5,/,' Ecoh:',f10.5,/,
     .' Embed1:',f10.5,/,' Embed2:',f10.5,/,
     .' L12 Cu3Ni',/,' a_0:',f10.5,/,' Ecoh1:',f10.5,/,
     .' Embed1:',f10.5,/,' Embed2:',f10.5,/,
     .' bmod:',f10.5,/,' c11:',f9.5,/,' c12:',f9.5,/,' c44:',f9.5,/,
     .' c13:',f9.5,/,' c33:',f9.5,/,' c55:',f9.5,/,
     .' c66:',f9.5,/,' Evf1:',f9.5,/,' Evf2:',f9.5,/,
     .' Bmod:',f10.5,/))
7278    format(45(E25.16,1x))
7275    format((/,' De:',f10.5,/,' Re:',f10.5,/,
     .' L10-L12:',f10.5,/,' L13-L11:',f10.5,/,
     .' Cu3Ni SISF(111):',f10.5,/,' Cu3Ni ABP(100):',f10.5,/,
     .' Cu3Ni ABP(111):',f10.5,/,
     .' CuNi SISF(111):',f10.5,/,' CuNi ABP(100):',f10.5,/,
     .' CuNi ABP(111):',f10.5,/,
     .' Cu SISF(111):',f10.5,/,' Cu ABP(100):',f10.5,/,
     .' Cu ABP(111):',f10.5,/,
     .' Ni SISF(111):',f10.5,/,' Ni ABP(100):',f10.5,/,
     .' Ni ABP(111):',f10.5,/))
5257    format((' L13 CuNi',/,' a_0:',f10.5,/,' Ecoh:',f10.5,/,
     .' Embed1:',f10.5,/,' Embed2:',f10.5,/,
     .' L11 Cu3Ni',/,' a_0:',f10.5,/,' Ecoh1:',f10.5,/,
     .' Embed1:',f10.5,/,' Embed2:',f10.5,/,
     .' bmod:',f10.5,/,' c11:',f9.5,/,' c12:',f9.5,/,' c44:',f9.5,/,
     .' c13:',f9.5,/,' c33:',f9.5,/,' c55:',f9.5,/,
     .' c66:',f9.5,/,' Evf1:',f9.5,/,' Evf2:',f9.5,/,
     .' Bmod:',f10.5,/))
c     write(*,*)'--------------------'
c     write(*,*)'RMS:',rmssum
        close(70)
      endif
c     ----------------------------------------------
c     ALLOY CALCULATION ENDS HERE
      endif
c     ----------------------------------------------
       if (rmssum.eq.0.0) then
           rmssum=100007.0
       endif
      return
      end


      subroutine splineinit(x,y,n,ypp)
      integer n,MAXN
      double precision yp1,ypn,x(n),y(n),ypp(n)
      PARAMETER (MAXN=10000)
c     INPUT: x,y,yp1,ypn
c     Given arrays x(1:n) and y(1:n) containing a tabulated function
c     yi = f(xi) with x1<x2<...<xn, and given values yp1 and ypn for
c     the first derivative of the interpolating function at points 1
c     and n
c     OUTPUT: yp2
c     this routine returns an array yp2(1:n) of length n which contains
c     the second derivatives of the interpolating function at the
c     tabulated points xi. If yp1 and/or ypn are equal to 1.dE30 or
c     larger, the routine is signaled to set the corresponding boundary
c     condiditon for "naturel spline", with zero second derivative of
c     that boundary.
c     Parameter: NMAX is the largest anticipated value of n
      integer i,k
      double precision p,qn,sig,un,u(MAXN)
c       if (yp1.gt..99e30) then
          ypp(1)=0.
          u(1)=0.
c       else
c         ypp(1)=-0.5
c         u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
c       endif
        do 11 i=2,n-1
          sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
          p=sig*ypp(i-1)+2.
          ypp(i)=(sig-1.)/p
          u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     .   / (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11      continue
c       if (ypn.gt..99e30) then
              qn=0.
              un=0.
c       else
c             qn=0.5
c             un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
c       endif
        ypp(n)=(un-qn*u(n-1))/(qn*ypp(n-1)+1.)
        do 12 k=n-1,1,-1
           ypp(k)=ypp(k)*ypp(k+1)+u(k)
12      continue
      return
      end

      subroutine splineinter(xa,ya,yppa,n,x,y)
      integer n
      double precision x,y,xa(n),yppa(n),ya(n)
c     INPUT: xa(n),ya(n),yppa(n),n,x
c     OUTPUT: y
c     Given the arrays xa(1:n) and ya(1:n) of length n, which tabulated
c     functions (with xa's order), and given the array yppa(1:n), which
c     is the output from splineinit subroutine, and given a value of x
c     this routine returns a cubic-spline interpolation value y.
      integer k,khi,klo
      double precision a,b,h
        klo=1
        khi=n
        do while ((khi-klo).gt.1)
              k=(khi+klo)/2
              if(xa(k).gt.x)then
                      khi=k
              else
                      klo=k
              endif
        enddo
        h=xa(khi)-xa(klo)
        if(h.eq.0.) then
              write(*,*)"Bad xa input in splineinter"
              pause
        endif
        a=(xa(khi)-x)/h
        b=(x-xa(klo))/h
        y=a*ya(klo)+b*ya(khi)+
     .     ((a**3-a)*yppa(klo)+(b**3-b)*yppa(khi))*(h**2)/6.
      return
      end

c*************************************************************
c                 Force Constant Calculation                 *
c*************************************************************
      subroutine forcon(neigh,rdis,rdisd,jj,iat,fp,fpp,dcon)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      integer neigh(nmax),jj,iat
      double precision rdis(nmax,nmax,2),rdisd(nmax,nmax,3)
      double precision fp(nmax),fpp(nmax),rhoip(neimax),rhoipp(neimax),
     . rhojp(neimax),rhojpp(neimax),diski(3),dcon(3,3),gvec(3),djk(3)
      common /grids/ frhoar(ngridar,nelmax),frhoar1(ngridar,nelmax),
     .  frhoar2(ngridar,nelmax),frhoar3(ngridar,nelmax),
     .  frhoar4(ngridar,nelmax),frhoar5(ngridar,nelmax),
     .  frhoar6(ngridar,nelmax),rhorar(ngridar,nelmax),
     .  rhorar1(ngridar,nelmax),rhorar2(ngridar,nelmax),
     .  rhorar3(ngridar,nelmax),rhorar4(ngridar,nelmax),
     .  rhorar5(ngridar,nelmax),rhorar6(ngridar,nelmax),
     .  z2rar(ngridar,nelmax,nelmax),z2rar1(ngridar,nelmax,nelmax),
     .  z2rar2(ngridar,nelmax,nelmax),z2rar3(ngridar,nelmax,nelmax),
     .  z2rar4(ngridar,nelmax,nelmax),z2rar5(ngridar,nelmax,nelmax),
     .  z2rar6(ngridar,nelmax,nelmax),drhoar,drar,nrhoar,nrar
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      real rij,pij,r(nmax),p(nmax),delta
      double precision rdist,rdisti
      double precision rdrar,roipij,rojpij,roippij,rojppij,
     .  z2ij,z2pij,z2ppij,phiij,phipij,phippij,dis(3,nmax)
      integer kij,phitype,mc,mr,m,indx(2),k(nmax)
        logical prtctl
       prtctl=.false.
c
      rdrar = 1.0/drar
      rdist=rdis(iat,jj,1)
      jat=rdis(iat,jj,2)
      rij=rdist
      pij = rij*rdrar + 1.0
      kij = pij
      kij = min0(kij,nrar-1)
      pij = pij - kij
      pij = amin1(pij,1.)
c
ccdir$ novector
      ity = itype(iat) 
      jty = itype(jat)
      z2ij = ((z2rar3(kij,ity,jty)*pij+
     $          z2rar2(kij,ity,jty))*pij+
     $          z2rar1(kij,ity,jty))*pij+
     $          z2rar(kij,ity,jty)
      z2pij = (z2rar6(kij,ity,jty)*pij+z2rar5(kij,ity,jty))*pij+
     $         z2rar4(kij,ity,jty)
      z2ppij=(2.0*z2rar6(kij,ity,jty)*pij+z2rar5(kij,ity,jty))/drar

      roipij = (rhorar6(kij,ity)*pij+rhorar5(kij,ity))*pij+
     $         rhorar4(kij,ity)
      rojpij = (rhorar6(kij,jty)*pij+rhorar5(kij,jty))*pij+
     $         rhorar4(kij,jty)

      roippij=(2.0*rhorar6(kij,ity)*pij+rhorar5(kij,ity))/drar
      rojppij=(2.0*rhorar6(kij,jty)*pij+rhorar5(kij,jty))/drar
ccdir$ vector
c
c     compute the contribution at
c      particle i due to particle j

      pij=1.0/rij
      diski(1) = rdisd(iat,jj,1)
      diski(2) = rdisd(iat,jj,2)
      diski(3) = rdisd(iat,jj,3)
      diski(1) = diski(1)*pij
      diski(2) = diski(2)*pij
      diski(3) = diski(3)*pij

      phiij = z2ij
      phipij = z2pij
      phippij = z2ppij

      phipij  = phipij  + fp(iat) * rojpij + fp(jat) * roipij 
      phippij = phippij + fp(iat) * rojppij+ fp(jat) * roippij 
        do mc=1,3
        do mr=1,3
        delta=0.0
        if(mc.eq.mr)delta=1.0
        dcon(mr,mc)=-( phippij*diski(mr)*diski(mc)
     .               + phipij*(delta-diski(mr)*diski(mc))/rij )
        enddo
        enddo
        indx(1)=iat
        indx(2)=jat
        
      do 3000 ij=1,2 
        i=indx(ij)

      do 3200 j=1,neigh(i)
      p(j) = rdis(i,j,1)*rdrar + 1.0
      k(j) = p(j)
      k(j) = min0(k(j),nrar-1)
      p(j) = p(j) - k(j)
      p(j) = amin1(p(j),1.)
3200     continue
c
      do 3500 j=1,neigh(i)
      jty = itype(rdis(i,j,2))
      rhojp(j) = (rhorar6(k(j),jty)*p(j)+rhorar5(k(j),jty))*p(j)+
     $         rhorar4(k(j),jty)
      rhojpp(j)=(2.0*rhorar6(k(j),jty)*p(j)+rhorar5(k(j),jty))/drar
3500  continue
c
c      particle i due to particle j
c
        do m=1,3
        gvec(m)=0.0
        enddo

      do 3700 j=1,neigh(i)
      p(j)=1.0/rdis(i,j,1)
      dis(1,j) = rdisd(i,j,1)*p(j)
      dis(2,j) = rdisd(i,j,2)*p(j)
      dis(3,j) = rdisd(i,j,3)*p(j)

        do m=1,3
        gvec(m) = gvec(m) + rhojp(j)* dis(m,j)
        enddo
3700  continue

        do mc=1,3
        do mr=1,3
        delta=0.0
        if(i.eq.iat)then
        dcon(mr,mc) = dcon(mr,mc) -
     .      fpp(iat)*rojpij * diski(mr)*gvec(mc)
        else
        dcon(mr,mc) = dcon(mr,mc) +
     .      fpp(jat)*roipij * gvec(mr)*diski(mc) 
        endif
        enddo
        enddo
3000  continue

        i=iat

      do 4100 j=1,neigh(i)

        do m=1,3
        djk(m)=rdisd(i,j,m)-diski(m)*rij
        enddo
        rjk=djk(1)**2+djk(2)**2+djk(3)**2
        if(rjk.le.1.0e-3.or.rjk.gt.rcutsq)goto 4100

      rik = rdis(i,j,1)
      pik = rik*rdrar + 1.0
      kik = pik
      kik = min0(kik,nrar-1)
      pik = pik - kik
      pik = amin1(pik,1.)
      jty = itype(iat)
      roipik = (rhorar6(kik,jty)*pik+rhorar5(kik,jty))*pik+
     $         rhorar4(kik,jty)

      rjk=sqrt(rjk)
      pik = rjk*rdrar + 1.0
      kik = pik
      kik = min0(kik,nrar-1)
      pik = pik - kik
      pik = amin1(pik,1.)
      jty = itype(jat)
      rojpjk = (rhorar6(kik,jty)*pik+rhorar5(kik,jty))*pik+
     $         rhorar4(kik,jty)
        do m=1,3
        djk(m)=djk(m)/rjk
        dis(m,j)=rdisd(i,j,m)/rik
        enddo

        do mc=1,3
        do mr=1,3

        dcon(mr,mc) = dcon(mr,mc) + fpp(rdis(i,j,2))*roipik*rojpjk
     .      * djk(mr)* dis(mc,j) 
        enddo
        enddo

4100  continue
        return
        end 



c*************************************************************
c     Elastic Constants for (HCP) C11, C12, C44, C13, C33, C66
c*************************************************************

      subroutine elastic(phitype,diski,jat,iat,uu,ww,vv,bulkmod)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      common /grids/ frhoar(ngridar,nelmax),frhoar1(ngridar,nelmax),
     .  frhoar2(ngridar,nelmax),frhoar3(ngridar,nelmax),
     .  frhoar4(ngridar,nelmax),frhoar5(ngridar,nelmax),
     .  frhoar6(ngridar,nelmax),rhorar(ngridar,nelmax),
     .  rhorar1(ngridar,nelmax),rhorar2(ngridar,nelmax),
     .  rhorar3(ngridar,nelmax),rhorar4(ngridar,nelmax),
     .  rhorar5(ngridar,nelmax),rhorar6(ngridar,nelmax),
     .  z2rar(ngridar,nelmax,nelmax),z2rar1(ngridar,nelmax,nelmax),
     .  z2rar2(ngridar,nelmax,nelmax),z2rar3(ngridar,nelmax,nelmax),
     .  z2rar4(ngridar,nelmax,nelmax),z2rar5(ngridar,nelmax,nelmax),
     .  z2rar6(ngridar,nelmax,nelmax),drhoar,drar,nrhoar,nrar
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      double precision diski(3),uu(7),ww(7),vv(7),bulkmod(7,nmax)
      real rij,pij
      double precision rdrar,roipij,rojpij,roippij,rojppij,
     .  z2ij,z2pij,z2ppij,phiij,phipij,phippij,dis(3)
      integer kij,phitype,lattype
        logical prtctl
       prtctl=.false.
c      prtctl=.true.
c
c  compute the reciprocal of the grid spacings
c
        if(prtctl)print *,'jat=',jat,'   iat=',iat
      rdrar = 1.0/drar 
      rij= sqrt(diski(1)**2 + diski(2)**2 + diski(3)**2)
      pij = rij*rdrar + 1.0
      kij = pij
      kij = min0(kij,nrar-1)
      pij = pij - kij
      pij = amin1(pij,1.)
c
ccdir$ novector
      ity = itype(iat) 
      jty = itype(jat)
      z2ij = ((z2rar3(kij,ity,jty)*pij+
     $          z2rar2(kij,ity,jty))*pij+
     $          z2rar1(kij,ity,jty))*pij+
     $          z2rar(kij,ity,jty)
      z2pij = (z2rar6(kij,ity,jty)*pij+z2rar5(kij,ity,jty))*pij+
     $         z2rar4(kij,ity,jty)
      z2ppij=(2.0*z2rar6(kij,ity,jty)*pij+z2rar5(kij,ity,jty))/drar

      roipij = (rhorar6(kij,ity)*pij+rhorar5(kij,ity))*pij+
     $         rhorar4(kij,ity)
      rojpij = (rhorar6(kij,jty)*pij+rhorar5(kij,jty))*pij+
     $         rhorar4(kij,jty)

      roippij=(2.0*rhorar6(kij,ity)*pij+rhorar5(kij,ity))/drar
      rojppij=(2.0*rhorar6(kij,jty)*pij+rhorar5(kij,jty))/drar
ccdir$ vector
c
      pij=1.0/rij
      dis(1) = diski(1)*pij
      dis(2) = diski(2)*pij
      dis(3) = diski(3)*pij

       phiij = z2ij
       phipij = z2pij
       phippij = z2ppij
c        c11  part of uu
         uu(1)= uu(1) + (phippij-(phipij*pij))*
     .       dis(1)*dis(1)*dis(1)*dis(1)*rij*rij
c        c12  part of uu
         uu(2)= uu(2) + (phippij-(phipij*pij))*
     .       dis(1)*dis(1)*dis(2)*dis(2)*rij*rij
c        c44  part of uu
         uu(3)= uu(3) + (phippij-(phipij*pij))*
     .       dis(2)*dis(3)*dis(2)*dis(3)*rij*rij
c        c13  part of uu
         uu(4)= uu(4) + (phippij-(phipij*pij))*
     .       dis(1)*dis(1)*dis(3)*dis(3)*rij*rij
c        c33  part of uu
         uu(5)= uu(5) + (phippij-(phipij*pij))*
     .       dis(3)*dis(3)*dis(3)*dis(3)*rij*rij
c        c55  part of uu
         uu(6)= uu(6) + (phippij-(phipij*pij))*
     .       dis(1)*dis(3)*dis(1)*dis(3)*rij*rij
c        c66  part of uu
         uu(7)= uu(7) + (phippij-(phipij*pij))*
     .       dis(1)*dis(2)*dis(1)*dis(2)*rij*rij
         uu(1)= uu(1)+phipij*dis(1)*dis(1)*rij
         uu(3)= uu(3)+phipij*dis(2)*dis(2)*rij
         uu(5)= uu(5)+phipij*dis(3)*dis(3)*rij
         uu(6)= uu(6)+phipij*dis(1)*dis(1)*rij
         uu(7)= uu(7)+phipij*dis(1)*dis(1)*rij

        bulkmod(1,iat) = bulkmod(1,iat) + phippij*(rij**2)
        bulkmod(4,iat) = bulkmod(4,iat) + phippij*(rij**2)
     .                                      - 2.*phipij*rij

c        c11  part of ww
         ww(1)= ww(1) + (rojppij-(rojpij*pij))*
     .       dis(1)*dis(1)*dis(1)*dis(1)*rij*rij
c        c12  part of ww
         ww(2)= ww(2) + (rojppij-(rojpij*pij))*
     .       dis(1)*dis(1)*dis(2)*dis(2)*rij*rij
c        c44  part of ww
         ww(3)= ww(3) + (rojppij-(rojpij*pij))*
     .       dis(2)*dis(3)*dis(2)*dis(3)*rij*rij
c        c13  part of ww
         ww(4)= ww(4) + (rojppij-(rojpij*pij))*
     .       dis(1)*dis(1)*dis(3)*dis(3)*rij*rij
c        c33  part of ww
         ww(5)= ww(5) + (rojppij-(rojpij*pij))*
     .       dis(3)*dis(3)*dis(3)*dis(3)*rij*rij
c        c55  part of ww
         ww(6)= ww(6) + (rojppij-(rojpij*pij))*
     .       dis(1)*dis(3)*dis(1)*dis(3)*rij*rij
c        c66  part of ww
         ww(7)= ww(7) + (rojppij-(rojpij*pij))*
     .       dis(1)*dis(2)*dis(1)*dis(2)*rij*rij
c        c11,c12,c44,c13,c33,c66  part of ww for hcp
c        c11,c44,c33,c66  part of ww for hcp
         ww(1)= ww(1)+rojpij*dis(1)*dis(1)*rij
         ww(3)= ww(3)+rojpij*dis(2)*dis(2)*rij
         ww(5)= ww(5)+rojpij*dis(3)*dis(3)*rij
         ww(6)= ww(6)+rojpij*dis(1)*dis(1)*rij
         ww(7)= ww(7)+rojpij*dis(1)*dis(1)*rij
c       
        bulkmod(2,iat) = bulkmod(2,iat) + rojppij*(rij**2)
        bulkmod(5,iat) = bulkmod(5,iat) + rojppij*(rij**2)
     .                                      - 2.*rojpij*rij
        
     
c        c11,c12,c44,c13,c33,c66  part of vv
         vv(1)= vv(1) + rojpij*dis(1)*dis(1)*rij
         vv(2)= vv(2) + rojpij*dis(2)*dis(2)*rij
         vv(3)= vv(3) + rojpij*dis(2)*dis(3)*rij
         vv(4)= vv(4) + rojpij*dis(3)*dis(3)*rij
         vv(5)= vv(5) + rojpij*dis(1)*dis(3)*rij
         vv(6)= vv(6) + rojpij*dis(1)*dis(3)*rij
         vv(7)= vv(7) + rojpij*dis(1)*dis(2)*rij
c       
        bulkmod(3,iat) = bulkmod(3,iat) + rojpij*rij
        bulkmod(6,iat) = bulkmod(6,iat) + rojpij*rij
       

        return
        end 


c*************************************************************
c       Elastic Constants C11, C12, C44 calculations
c*************************************************************

      subroutine elascon(phitype,diski,jat,iat,uu,ww,vv,bulkmod)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      common /grids/ frhoar(ngridar,nelmax),frhoar1(ngridar,nelmax),
     .  frhoar2(ngridar,nelmax),frhoar3(ngridar,nelmax),
     .  frhoar4(ngridar,nelmax),frhoar5(ngridar,nelmax),
     .  frhoar6(ngridar,nelmax),rhorar(ngridar,nelmax),
     .  rhorar1(ngridar,nelmax),rhorar2(ngridar,nelmax),
     .  rhorar3(ngridar,nelmax),rhorar4(ngridar,nelmax),
     .  rhorar5(ngridar,nelmax),rhorar6(ngridar,nelmax),
     .  z2rar(ngridar,nelmax,nelmax),z2rar1(ngridar,nelmax,nelmax),
     .  z2rar2(ngridar,nelmax,nelmax),z2rar3(ngridar,nelmax,nelmax),
     .  z2rar4(ngridar,nelmax,nelmax),z2rar5(ngridar,nelmax,nelmax),
     .  z2rar6(ngridar,nelmax,nelmax),drhoar,drar,nrhoar,nrar
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      double precision diski(3),uu(7),ww(7),vv(7),bulkmod(7,nmax)
      real rij,pij
      double precision rdrar,roipij,rojpij,roippij,rojppij,
     .  z2ij,z2pij,z2ppij,phiij,phipij,phippij
      integer kij,phitype
        logical prtctl
       prtctl=.false.
c      prtctl=.true.
c
        if(prtctl)print *,'jat=',jat,'   iat=',iat
      rdrar = 1.0/drar 
      rij= sqrt(diski(1)**2 + diski(2)**2 + diski(3)**2)
      pij = rij*rdrar + 1.0
      kij = pij
      kij = min0(kij,nrar-1)
      pij = pij - kij
      pij = amin1(pij,1.)
c
ccdir$ novector
      ity = itype(iat) 
      jty = itype(jat)
      z2ij = ((z2rar3(kij,ity,jty)*pij+
     $          z2rar2(kij,ity,jty))*pij+
     $          z2rar1(kij,ity,jty))*pij+
     $          z2rar(kij,ity,jty)
      z2pij = (z2rar6(kij,ity,jty)*pij+z2rar5(kij,ity,jty))*pij+
     $         z2rar4(kij,ity,jty)
      z2ppij=(2.0*z2rar6(kij,ity,jty)*pij+z2rar5(kij,ity,jty))/drar

      roipij = (rhorar6(kij,ity)*pij+rhorar5(kij,ity))*pij+
     $         rhorar4(kij,ity)
      rojpij = (rhorar6(kij,jty)*pij+rhorar5(kij,jty))*pij+
     $         rhorar4(kij,jty)

      roippij=(2.0*rhorar6(kij,ity)*pij+rhorar5(kij,ity))/drar
      rojppij=(2.0*rhorar6(kij,jty)*pij+rhorar5(kij,jty))/drar
ccdir$ vector
c
      pij=1.0/rij
      diski(1) = diski(1)*pij
      diski(2) = diski(2)*pij
      diski(3) = diski(3)*pij

c     if(phitype) then
       phiij = z2ij
       phipij = z2pij
       phippij = z2ppij
c     else
c      phiij = z2ij * pij
c      phipij = z2pij * pij - phiij * pij
c      phippij = z2ppij * pij - 2.0*phipij * pij
c     endif

        uu(1)= uu(1) + (phippij-(phipij*pij))*
     .       diski(1)*diski(1)*diski(1)*diski(1)*rij*rij
        uu(2)= uu(2) + (phippij-(phipij*pij))*
     .       diski(1)*diski(1)*diski(2)*diski(2)*rij*rij
        uu(3)= uu(3) + (phippij-(phipij*pij))*
     .       diski(2)*diski(3)*diski(2)*diski(3)*rij*rij
        bulkmod(1,iat) = bulkmod(1,iat) + phippij*(rij**2)


        ww(1)= ww(1) + (rojppij-(rojpij*pij))*
     .       diski(1)*diski(1)*diski(1)*diski(1)*rij*rij
        ww(2)= ww(2) + (rojppij-(rojpij*pij))*
     .       diski(1)*diski(1)*diski(2)*diski(2)*rij*rij
        ww(3)= ww(3) + (rojppij-(rojpij*pij))*
     .       diski(2)*diski(3)*diski(2)*diski(3)*rij*rij
        bulkmod(2,iat) = bulkmod(2,iat) + rojppij*(rij**2)
     
        
        vv(1)= vv(1) + rojpij*diski(1)*diski(1)*rij
        vv(2)= vv(2) + rojpij*diski(2)*diski(2)*rij
        vv(3)= vv(3) + rojpij*diski(2)*diski(3)*rij
        bulkmod(3,iat) = bulkmod(3,iat) + rojpij*rij
        
     
        return
        end 


      subroutine calctote(etot,neigh,rdis)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision etot,rdis(nmax,nmax,2)
      integer neigh(nmax)
      common /grids/ frhoar(ngridar,nelmax),frhoar1(ngridar,nelmax),
     .  frhoar2(ngridar,nelmax),frhoar3(ngridar,nelmax),
     .  frhoar4(ngridar,nelmax),frhoar5(ngridar,nelmax),
     .  frhoar6(ngridar,nelmax),rhorar(ngridar,nelmax),
     .  rhorar1(ngridar,nelmax),rhorar2(ngridar,nelmax),
     .  rhorar3(ngridar,nelmax),rhorar4(ngridar,nelmax),
     .  rhorar5(ngridar,nelmax),rhorar6(ngridar,nelmax),
     .  z2rar(ngridar,nelmax,nelmax),z2rar1(ngridar,nelmax,nelmax),
     .  z2rar2(ngridar,nelmax,nelmax),z2rar3(ngridar,nelmax,nelmax),
     .  z2rar4(ngridar,nelmax,nelmax),z2rar5(ngridar,nelmax,nelmax),
     .  z2rar6(ngridar,nelmax,nelmax),drhoar,drar,nrhoar,nrar
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      double precision rhoij,femb,phiij
      real rij,pij
      double precision rdrar,rdrhoar,z2ij,rdist
      integer kij,phitype,iat,j,jat,ity,jty
      
      etot=0.0
      rdrar = 1.0/drar 
      rdrhoar = 1.0/drhoar
c
c     compute the contribution to rho(i) from particle j
c     Calculating the embedding energy F(Rho) in fe(i)
c      and Pair Energy contribution Phi(Rij) in phie(i)
c
ccdir$ novector
      do 11 iat=1,natoms
       rhoij=0.0
       phiij=0.0
       ity=itype(iat)
      do 12 j=1,neigh(iat)
      rdist=rdis(iat,j,1)
      jat=rdis(iat,j,2)
      pij = rdist*rdrar + 1.0
      kij = pij
      kij = min0(kij,nrar-1)
      pij = pij - kij
      pij = amin1(pij,1.0)
      jty=itype(jat)
      rhoij = rhoij + ((rhorar3(kij,jty)*pij+
     $   rhorar2(kij,jty))*pij+
     $   rhorar1(kij,jty))*pij+
     $   rhorar(kij,jty)
      z2ij = ((z2rar3(kij,ity,jty)*pij+
     $          z2rar2(kij,ity,jty))*pij+
     $          z2rar1(kij,ity,jty))*pij+
     $          z2rar(kij,ity,jty)
c     if(phitype) then
        z2ij = z2ij
c     else
c       pij=1.0/rdist
c       z2ij = z2ij * pij
c     endif
      phiij = phiij + 0.5 * z2ij
12    continue
      pij = rhoij*rdrhoar + 1.0
      kij = pij
      kij = max(1,min(kij,nrhoar-1))
      pij = pij - kij
      femb = ((frhoar3(kij,ity)*pij+
     $        frhoar2(kij,ity))*pij+
     $        frhoar1(kij,ity))*pij+
     $        frhoar(kij,ity)
      etot = etot + phiij + femb
11    continue
ccdir$ vector
c        
      
      return
      end


      subroutine calctotetype(tn,typ,typn,embed,etot,neigh,rdis)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision etot,embed(nelmax),rdis(nmax,nmax,2)
      integer neigh(nmax),tn,typ(nelmax),typn(nelmax)
      common /grids/ frhoar(ngridar,nelmax),frhoar1(ngridar,nelmax),
     .  frhoar2(ngridar,nelmax),frhoar3(ngridar,nelmax),
     .  frhoar4(ngridar,nelmax),frhoar5(ngridar,nelmax),
     .  frhoar6(ngridar,nelmax),rhorar(ngridar,nelmax),
     .  rhorar1(ngridar,nelmax),rhorar2(ngridar,nelmax),
     .  rhorar3(ngridar,nelmax),rhorar4(ngridar,nelmax),
     .  rhorar5(ngridar,nelmax),rhorar6(ngridar,nelmax),
     .  z2rar(ngridar,nelmax,nelmax),z2rar1(ngridar,nelmax,nelmax),
     .  z2rar2(ngridar,nelmax,nelmax),z2rar3(ngridar,nelmax,nelmax),
     .  z2rar4(ngridar,nelmax,nelmax),z2rar5(ngridar,nelmax,nelmax),
     .  z2rar6(ngridar,nelmax,nelmax),drhoar,drar,nrhoar,nrar
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      double precision rhoij,femb,phiij
      real rij,pij
      double precision rdrar,rdrhoar,z2ij,rdist
      integer kij,phitype,iat,j,jat,ity,jty
      
      etot=0.0
      do 10 i=1,3
         embed(i)=0.0
         typn(i)=0
10    continue
      rdrar = 1.0/drar 
      rdrhoar = 1.0/drhoar
c
c     compute the contribution to rho(i) from particle j
c     Calculating the embedding energy F(Rho) in fe(i)
c      and Pair Energy contribution Phi(Rij) in phie(i)
c
ccdir$ novector
      do 11 iat=1,natoms
       rhoij=0.0
       phiij=0.0
       ity=itype(iat)
      do 12 j=1,neigh(iat)
      rdist=rdis(iat,j,1)
      jat=rdis(iat,j,2)
      pij = rdist*rdrar + 1.0
      kij = pij
      kij = min0(kij,nrar-1)
      pij = pij - kij
      pij = amin1(pij,1.0)
      jty=itype(jat)
      rhoij = rhoij + ((rhorar3(kij,jty)*pij+
     $   rhorar2(kij,jty))*pij+
     $   rhorar1(kij,jty))*pij+
     $   rhorar(kij,jty)
      z2ij = ((z2rar3(kij,ity,jty)*pij+
     $          z2rar2(kij,ity,jty))*pij+
     $          z2rar1(kij,ity,jty))*pij+
     $          z2rar(kij,ity,jty)
c     if(phitype) then
        z2ij = z2ij
c     else
c       pij=1.0/rdist
c       z2ij = z2ij * pij
c     endif
      phiij = phiij + 0.5 * z2ij
12    continue
      pij = rhoij*rdrhoar + 1.0
      kij = pij
      kij = max(1,min(kij,nrhoar-1))
      pij = pij - kij
      femb = ((frhoar3(kij,ity)*pij+
     $        frhoar2(kij,ity))*pij+
     $        frhoar1(kij,ity))*pij+
     $        frhoar(kij,ity)
      if(tn.eq.3) then
        if(ity.eq.typ(1)) then
          typn(1)=typn(1)+1
          embed(1)=embed(1)+phiij+femb
        else if(ity.eq.typ(2)) then
          typn(2)=typn(2)+1
          embed(2)=embed(2)+phiij+femb
        else if(ity.eq.typ(3)) then
          typn(3)=typn(3)+1
          embed(3)=embed(3)+phiij+femb
        endif
      else if(tn.eq.2) then
        if(ity.eq.typ(1)) then
          typn(1)=typn(1)+1
          embed(1)=embed(1)+phiij+femb
        else if(ity.eq.typ(2)) then
          typn(2)=typn(2)+1
          embed(2)=embed(2)+phiij+femb
        endif
      else
          typn(1)=typn(1)+1
          embed(1)=embed(1)+phiij+femb
      endif
      etot = etot + phiij + femb
11    continue
ccdir$ vector
c        
      if(tn.eq.3) then
          embed(1)=embed(1)/typn(1)
          embed(2)=embed(2)/typn(2)
          embed(3)=embed(3)/typn(3)
      else if(tn.eq.2) then
          embed(1)=embed(1)/typn(1)
          embed(2)=embed(2)/typn(2)
      else
          embed(1)=embed(1)/typn(1)
      endif
      
      return
      end



      subroutine calcmine(etot,lata,lata0,neigh,rdis)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision etot,lata,lata0,rdis(nmax,nmax,2)
      integer neigh(nmax)
      common /grids/ frhoar(ngridar,nelmax),frhoar1(ngridar,nelmax),
     .  frhoar2(ngridar,nelmax),frhoar3(ngridar,nelmax),
     .  frhoar4(ngridar,nelmax),frhoar5(ngridar,nelmax),
     .  frhoar6(ngridar,nelmax),rhorar(ngridar,nelmax),
     .  rhorar1(ngridar,nelmax),rhorar2(ngridar,nelmax),
     .  rhorar3(ngridar,nelmax),rhorar4(ngridar,nelmax),
     .  rhorar5(ngridar,nelmax),rhorar6(ngridar,nelmax),
     .  z2rar(ngridar,nelmax,nelmax),z2rar1(ngridar,nelmax,nelmax),
     .  z2rar2(ngridar,nelmax,nelmax),z2rar3(ngridar,nelmax,nelmax),
     .  z2rar4(ngridar,nelmax,nelmax),z2rar5(ngridar,nelmax,nelmax),
     .  z2rar6(ngridar,nelmax,nelmax),drhoar,drar,nrhoar,nrar
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      double precision rhoij,femb,phiij
      real rij,pij
      double precision rdrar,rdrhoar,z2ij,rdist
      integer kij,phitype,iat,j,jat,ity,jty
      
      etot=0.0
      rdrar = 1.0/drar 
      rdrhoar = 1.0/drhoar
c
c     compute the contribution to rho(i) from particle j
c     Calculating the embedding energy F(Rho) in fe(i)
c      and Pair Energy contribution Phi(Rij) in phie(i)
c
ccdir$ novector
      do 11 iat=1,natoms
       rhoij=0.0
       phiij=0.0
       ity=itype(iat)
      do 12 j=1,neigh(iat)
      rdist=(lata/lata0)*rdis(iat,j,1)
      jat=rdis(iat,j,2)
      pij = rdist*rdrar + 1.0
      kij = pij
      kij = min0(kij,nrar-1)
      pij = pij - kij
      pij = amin1(pij,1.0)
      jty=itype(jat)
      rhoij = rhoij + ((rhorar3(kij,jty)*pij+
     $   rhorar2(kij,jty))*pij+
     $   rhorar1(kij,jty))*pij+
     $   rhorar(kij,jty)
      z2ij = ((z2rar3(kij,ity,jty)*pij+
     $          z2rar2(kij,ity,jty))*pij+
     $          z2rar1(kij,ity,jty))*pij+
     $          z2rar(kij,ity,jty)
c     if(phitype) then
        z2ij = z2ij
c     else
c       pij=1.0/rdist
c       z2ij = z2ij * pij
c     endif
      phiij = phiij + 0.5 * z2ij
12    continue
      pij = rhoij*rdrhoar + 1.0
      kij = pij
      kij = max(1,min(kij,nrhoar-1))
      pij = pij - kij
      femb = ((frhoar3(kij,ity)*pij+
     $        frhoar2(kij,ity))*pij+
     $        frhoar1(kij,ity))*pij+
     $        frhoar(kij,ity)
      etot = etot + phiij + femb
11    continue
ccdir$ vector
c        
      
      return
      end


      subroutine calcminehex(etot,lata,lata0,ca,c0,neigh,rdis,rdisd,
     .rrms)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision etot,lata,lata0,rdis(nmax,nmax,2),rrms
      double precision ca,c0,rdisd(nmax,nmax,3)
      integer neigh(nmax)
      common /grids/ frhoar(ngridar,nelmax),frhoar1(ngridar,nelmax),
     .  frhoar2(ngridar,nelmax),frhoar3(ngridar,nelmax),
     .  frhoar4(ngridar,nelmax),frhoar5(ngridar,nelmax),
     .  frhoar6(ngridar,nelmax),rhorar(ngridar,nelmax),
     .  rhorar1(ngridar,nelmax),rhorar2(ngridar,nelmax),
     .  rhorar3(ngridar,nelmax),rhorar4(ngridar,nelmax),
     .  rhorar5(ngridar,nelmax),rhorar6(ngridar,nelmax),
     .  z2rar(ngridar,nelmax,nelmax),z2rar1(ngridar,nelmax,nelmax),
     .  z2rar2(ngridar,nelmax,nelmax),z2rar3(ngridar,nelmax,nelmax),
     .  z2rar4(ngridar,nelmax,nelmax),z2rar5(ngridar,nelmax,nelmax),
     .  z2rar6(ngridar,nelmax,nelmax),drhoar,drar,nrhoar,nrar
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      double precision rhoij,femb,phiij
      real rij,pij
      double precision rdrar,rdrhoar,z2ij,rdist
      integer kij,phitype,iat,j,jat,ity,jty,rmsflag
      
      rmsflag=0 
      etot=0.0
      rdrar = 1.0/drar 
      rdrhoar = 1.0/drhoar
c
c     compute the contribution to rho(i) from particle j
c     Calculating the embedding energy F(Rho) in fe(i)
c      and Pair Energy contribution Phi(Rij) in phie(i)
c
ccdir$ novector
      do 11 iat=1,natoms
       rhoij=0.0
       phiij=0.0
       ity=itype(iat)
      do 12 j=1,neigh(iat)
      rdist=((lata/lata0)*rdisd(iat,j,1))**2
      rdist=rdist+((lata/lata0)*rdisd(iat,j,2))**2
      rdist=rdist+((ca/c0)*rdisd(iat,j,3))**2
      if (rdist.lt.0.0) then
          rmsflag=1
          exit
      endif
      rdist=sqrt(rdist)
      if (rdist.gt.100000.0 .or. rdist.lt.0.0) then
          rmsflag=1
          exit
      endif
      jat=rdis(iat,j,2)
      pij = rdist*rdrar + 1.0
      kij = pij
      kij = min0(kij,nrar-1)
      pij = pij - kij
      pij = amin1(pij,1.0)
      jty=itype(jat)
      rhoij = rhoij + ((rhorar3(kij,jty)*pij+
     $   rhorar2(kij,jty))*pij+
     $   rhorar1(kij,jty))*pij+
     $   rhorar(kij,jty)
      z2ij = ((z2rar3(kij,ity,jty)*pij+
     $          z2rar2(kij,ity,jty))*pij+
     $          z2rar1(kij,ity,jty))*pij+
     $          z2rar(kij,ity,jty)
c     if(phitype) then
        z2ij = z2ij
c     else
c       pij=1.0/rdist
c       z2ij = z2ij * pij
c     endif
      phiij = phiij + 0.5 * z2ij
12    continue
      if (rmsflag.eq.1) then
          exit
      endif
      pij = rhoij*rdrhoar + 1.0
      kij = pij
      kij = max(1,min(kij,nrhoar-1))
      pij = pij - kij
      femb = ((frhoar3(kij,ity)*pij+
     $        frhoar2(kij,ity))*pij+
     $        frhoar1(kij,ity))*pij+
     $        frhoar(kij,ity)
      etot = etot + phiij + femb
11    continue
ccdir$ vector
c      
      if (rmsflag.eq.1) then
          rrms=100100.0
      endif
      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         Calculation of Diatomic Strength and         c
c          Diatomic Length for Type1 Element           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calcdiatom(t1,t2,decalc,recalc)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision demin1,demin2,demax,dedif
      double precision decalc,dephi,defe1,defe2,derho1,derho2,dez2ij
      double precision recalc,rei,ref,res
      integer ren,t1,t2
      common /grids/ frhoar(ngridar,nelmax),frhoar1(ngridar,nelmax),
     .  frhoar2(ngridar,nelmax),frhoar3(ngridar,nelmax),
     .  frhoar4(ngridar,nelmax),frhoar5(ngridar,nelmax),
     .  frhoar6(ngridar,nelmax),rhorar(ngridar,nelmax),
     .  rhorar1(ngridar,nelmax),rhorar2(ngridar,nelmax),
     .  rhorar3(ngridar,nelmax),rhorar4(ngridar,nelmax),
     .  rhorar5(ngridar,nelmax),rhorar6(ngridar,nelmax),
     .  z2rar(ngridar,nelmax,nelmax),z2rar1(ngridar,nelmax,nelmax),
     .  z2rar2(ngridar,nelmax,nelmax),z2rar3(ngridar,nelmax,nelmax),
     .  z2rar4(ngridar,nelmax,nelmax),z2rar5(ngridar,nelmax,nelmax),
     .  z2rar6(ngridar,nelmax,nelmax),drhoar,drar,nrhoar,nrar
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      double precision lrho,latfrho,lphi,latphi,totenergy
      double precision latai,lataf,latstep,lata
      real rij,pij,latno
      double precision rdrar,rdrhoar,z2ij,rdist
      integer kij,phitype,iat,j,jat,ity,jty
        rei=7.0
        ref=1.0
        ren=10001
        res=abs(rei-ref)/ren
        recalc=rei
        decalc=0.0 
        demax=0.0 
        demin1=0.0 
        demin2=0.0 
        dedif=0.0 
        rdrar = 1.0/drar 
      do 3700 i=1,ren
          pij = rei*rdrar + 1.0
          kij = pij
          kij = min0(kij,nrar-1)
          pij = pij - kij
          pij = amin1(pij,1.0)
      jty=t1
      derho1 = ((rhorar3(kij,jty)*pij+
     $   rhorar2(kij,jty))*pij+
     $   rhorar1(kij,jty))*pij+
     $   rhorar(kij,jty)
      jty=t2
      derho2 = ((rhorar3(kij,jty)*pij+
     $   rhorar2(kij,jty))*pij+
     $   rhorar1(kij,jty))*pij+
     $   rhorar(kij,jty)
      ity=t1
      jty=t2
      dephi = ((z2rar3(kij,ity,jty)*pij+
     $          z2rar2(kij,ity,jty))*pij+
     $          z2rar1(kij,ity,jty))*pij+
     $          z2rar(kij,ity,jty)
          pij = derho1*rdrhoar + 1.0
          kij = pij
          kij = max(1,min(kij,nrhoar-1))
          pij = pij - kij
      ity=t2
      defe2 = ((frhoar3(kij,ity)*pij+
     $        frhoar2(kij,ity))*pij+
     $        frhoar1(kij,ity))*pij+
     $        frhoar(kij,ity)
          pij = derho2*rdrhoar + 1.0
          kij = pij
          kij = max(1,min(kij,nrhoar-1))
          pij = pij - kij
      ity=t1
      defe1 = ((frhoar3(kij,ity)*pij+
     $        frhoar2(kij,ity))*pij+
     $        frhoar1(kij,ity))*pij+
     $        frhoar(kij,ity)
         if ((defe1+defe2+dephi).lt.demin) then
            demin = defe1 + defe2 + dephi
            decalc = defe1 + defe2 + dephi
            recalc = rei
         endif
         rei=rei-res
3700   continue

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccc


c****************************************************
c*   Calculation of Equilibrium Lattice Constant    *
c*                  for an Alloy                    *
c****************************************************
      subroutine calclatcon(minenergy,latis,lata0,neigh,rdis,iter,
     .prntflag)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision minenergy,latis,lata0,rdis(nmax,nmax,2)
      integer neigh(nmax),iter,prntflag
      common /grids/ frhoar(ngridar,nelmax),frhoar1(ngridar,nelmax),
     .  frhoar2(ngridar,nelmax),frhoar3(ngridar,nelmax),
     .  frhoar4(ngridar,nelmax),frhoar5(ngridar,nelmax),
     .  frhoar6(ngridar,nelmax),rhorar(ngridar,nelmax),
     .  rhorar1(ngridar,nelmax),rhorar2(ngridar,nelmax),
     .  rhorar3(ngridar,nelmax),rhorar4(ngridar,nelmax),
     .  rhorar5(ngridar,nelmax),rhorar6(ngridar,nelmax),
     .  z2rar(ngridar,nelmax,nelmax),z2rar1(ngridar,nelmax,nelmax),
     .  z2rar2(ngridar,nelmax,nelmax),z2rar3(ngridar,nelmax,nelmax),
     .  z2rar4(ngridar,nelmax,nelmax),z2rar5(ngridar,nelmax,nelmax),
     .  z2rar6(ngridar,nelmax,nelmax),drhoar,drar,nrhoar,nrar
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      double precision lrho,latfrho,lphi,latphi,totenergy
      double precision latai,lataf,latstep,mina
      double precision minai,minaf,minstep,lata
      double precision minua,minub,minda,mindb
      double precision minu,mind,latu,latd
      double precision latu1,latu2,latu3,latu4,latd1,latd2,latd3,latd4
      double precision totu1,totu2,totu3,totu4,totd1,totd2,totd3,totd4
      double precision lata1,lata2,lata3,lata4,tota1,tota2,tota3,tota4
      double precision emini,emint,latold,minold,tota
      double precision hnum,derivu,derivd,deriva
      double precision hitlata,hitca,hitmin,hitold,hitlatold
      real rij,pij,latno
      double precision rdrar,rdrhoar,z2ij,rdist
      integer kij,phitype,iat,j,jat,ity,jty,i,k
     
c
c     Sample some points to begin iteration
c
      latold=0.5
      latai=8.5
      lataf=0.5
      latno=17
      latstep=abs(latai-lataf)/(latno-1)
      emini=100000.0
      rdrar = 1.0/drar 
      rdrhoar = 1.0/drhoar
      do 47 i=1,latno
      mina=latai-(i-1)*latstep
      call calcmine(totenergy,mina,lata0,neigh,rdis)
        if(totenergy.lt.emini) then
                emini=totenergy
                lata=mina
        endif
47    continue
c
c     Searching with Bisection Method with First Derivatives
c
      iter=0
      latno=100
      minenergy=emini
      minold=0.0
      latold=0.0
      hnum=0.0001
      latstep=1.0
      do 5006 j=1,5
      latu=lata+latstep
      latd=lata-latstep
      if (prntflag.eq.1) then
       write(*,*)'   >>>Searching if minimum is inside region:'
       write(*,*)'      ',latd,'<lata<',latu
      endif
      latu1=latu+2*hnum
      latu2=latu+hnum
      latu3=latu-hnum
      latu4=latu-2*hnum
      call calcmine(totu1,latu1,lata0,neigh,rdis)
      call calcmine(totu2,latu2,lata0,neigh,rdis)
      call calcmine(totu3,latu3,lata0,neigh,rdis)
      call calcmine(totu4,latu4,lata0,neigh,rdis)
      latd1=latd+2*hnum
      latd2=latd+hnum
      latd3=latd-hnum
      latd4=latd-2*hnum
      call calcmine(totd1,latd1,lata0,neigh,rdis)
      call calcmine(totd2,latd2,lata0,neigh,rdis)
      call calcmine(totd3,latd3,lata0,neigh,rdis)
      call calcmine(totd4,latd4,lata0,neigh,rdis)
      derivu=(totu4-8.0*totu3+8.0*totu2-totu1)/(12.0*hnum)
      derivd=(totd4-8.0*totd3+8.0*totd2-totd1)/(12.0*hnum)
      if (prntflag.eq.1) then
       write(*,*)'      deriv up:',derivu
       write(*,*)'      deriv down:',derivd
      endif
      i=0
      if (derivu*derivd<0.and.derivd<0) then
        if (prntflag.eq.1) then
         write(*,*)'   >>>Minimum is inside region.'
        endif
        do 5007 i=1,latno
         lata=(latu+latd)/2.0
         lata1=lata+2*hnum
         lata2=lata+hnum
         lata3=lata-hnum
         lata4=lata-2*hnum
         call calcmine(tota1,lata1,lata0,neigh,rdis)
         call calcmine(tota2,lata2,lata0,neigh,rdis)
         call calcmine(tota3,lata3,lata0,neigh,rdis)
         call calcmine(tota4,lata4,lata0,neigh,rdis)
         call calcmine(minenergy,lata,lata0,neigh,rdis)
         deriva=(tota4-8.0*tota3+8.0*tota2-tota1)/(12.0*hnum)
c        write(*,*)'   >>>----ITERATING-----'
c        write(*,*)'      iter:',i
c        write(*,*)'      latis:',lata
c        write(*,*)'      minenergy:',minenergy
c        write(*,*)'      old energy:',minold
         latis=lata
         if (prntflag.eq.1) then
          if(abs(latold-latis).lt.1.0E-6) then
           if(minenergy.gt.minold) then
              write(*,*)'      Warning: Minimum energy is exceeded.'
              write(*,*)'      Old:',minold,'< New:',minenergy
              write(*,*)'      latis:',lata
           else
              write(*,*)'      Minimum energy found:',minenergy
              write(*,*)'      latis:',lata
           endif
           exit
          endif
         endif
         if(minenergy.lt.minold) then
           if (prntflag.eq.1) then
            write(*,*)'      Aprouching to minimum energy.'
            write(*,*)'      New:',minenergy,'< Old:',minold
            write(*,*)'      latis:',lata
           endif
c           hitold=minold
c           hitmin=minenergy
c           hitlata=lata
c           hitlatold=latold
c           latold=latis
c           minold=minenergy
         else
           if (prntflag.eq.1) then
            write(*,*)'      Minimum energy is exceeded:'
            write(*,*)'      Old:',minold,'< New:',minenergy
            write(*,*)'      latis:',lata
           endif
c           minold=minenergy
c           minenergy=hitmin
c           latold=lata
c           lata=hitlata
c           latis=hitlata
            hnum=0.000001
c           exit
         endif
         latold=latis
         minold=minenergy
         if (deriva>0) then
                 latu=lata
                 latd=latd
         else if (deriva<0) then
                 latu=latu
                 latd=lata
         endif
         if(minenergy.lt.emini) then
                emini=minenergy
                mina=lata
         endif
5007     continue
c        write(*,*)'   <<<--ITERATION END---'
      else
       if (prntflag.eq.1) then
        write(*,*)'   <<<Minimum is not inside region.'
        write(*,*)'   <<<Expanding region 0.1 out'
       endif
        latstep=latstep+0.5
c       exit 
      endif
      if (i.gt.0) then
        exit
      endif
5006  continue
        iter=i
       if (prntflag.eq.1) then
        write(*,*)'...iteration:',iter
        write(*,*)'Bulk search about min point:',mina
       endif
      latai=mina+0.05
      lataf=mina-0.05
      latno=201
      latstep=abs(latai-lataf)/(latno-1)
      do 57 i=1,latno
      mina=latai-(i-1)*latstep
      call calcmine(totenergy,mina,lata0,neigh,rdis)
        if(totenergy.lt.minenergy) then
                minenergy=totenergy
                latis=mina
        endif
57    continue
        if (prntflag.eq.1) then
         write(*,*)'----SEARCH RESULTS-----'
         write(*,*)'  latis:',latis
         write(*,*)'  minenergy:',minenergy
         write(*,*)'  First found minimum energy latis:',mina
         write(*,*)'  First found minimum energy:',emini
         write(*,*)'  Bisection Method latis:',latold
         write(*,*)'  Bisection Method energy:',minold
        endif
      return
      end
c*****************************************************


c****************************************************
c*   Calculation of Equilibrium Lattice Constant    *
c*                  for an Alloy                    *
c****************************************************
      subroutine calclatca(mine,latis,lata0,c,c0,neigh,rdis,rdisd,
     &rmss,prntflag)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision mine,latis,lata0,rdis(nmax,nmax,2),rmss
      double precision c,c0,rdisd(nmax,nmax,3)
      integer neigh(nmax),iter,prntflag
      common /grids/ frhoar(ngridar,nelmax),frhoar1(ngridar,nelmax),
     .  frhoar2(ngridar,nelmax),frhoar3(ngridar,nelmax),
     .  frhoar4(ngridar,nelmax),frhoar5(ngridar,nelmax),
     .  frhoar6(ngridar,nelmax),rhorar(ngridar,nelmax),
     .  rhorar1(ngridar,nelmax),rhorar2(ngridar,nelmax),
     .  rhorar3(ngridar,nelmax),rhorar4(ngridar,nelmax),
     .  rhorar5(ngridar,nelmax),rhorar6(ngridar,nelmax),
     .  z2rar(ngridar,nelmax,nelmax),z2rar1(ngridar,nelmax,nelmax),
     .  z2rar2(ngridar,nelmax,nelmax),z2rar3(ngridar,nelmax,nelmax),
     .  z2rar4(ngridar,nelmax,nelmax),z2rar5(ngridar,nelmax,nelmax),
     .  z2rar6(ngridar,nelmax,nelmax),drhoar,drar,nrhoar,nrar
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      double precision lrho,latfrho,lphi,latphi,totenergy
      double precision latai,lataf,latstep,mina
      double precision cai,caf,cstep,minc
      double precision minai,minaf,minstep,lata
      double precision minua,minub,minda,mindb
      double precision minu,mind,latu,latd,minenergy
      double precision latu1,latu2,latu3,latu4,latd1,latd2,latd3,latd4
      double precision totu1,totu2,totu3,totu4,totd1,totd2,totd3,totd4
      double precision lata1,lata2,lata3,lata4,tota1,tota2,tota3,tota4
      double precision emini,emina,eminc,latold,minold,tota
      double precision hnum,derivu,derivd,deriva
      double precision hitlata,hitca,hitmin,hitold,hitlatold
      real rij,pij,latno
      double precision rdrar,rdrhoar,z2ij,rdist,ctest
      double precision eps
      integer kij,phitype,iat,j,jat,ity,jty,i,k
     
c
c     Sample some points to begin iteration
c
      latold=0.5
      latai=8.5
      lataf=1.0
      latno=81
      cai=2.5
      caf=0.5
      cno=21
      latstep=abs(latai-lataf)/(latno-1)
      cstep=abs(cai-caf)/(cno-1)
      emini=100000.0
      rdrar = 1.0/drar 
      rdrhoar = 1.0/drhoar
      do 47 i=1,latno
      mina=latai-(i-1)*latstep
       do 57 j=1,cno
        minc=cai-(j-1)*cstep
        call calcminehex(totenergy,mina,lata0,minc,c0,neigh,rdis,
     &rdisd,rmss)
        if (rmss.gt.10000.0) then
            return
        endif
        if(totenergy.lt.emini) then
                emini=totenergy
                lata=mina
                ctest=minc
        endif
57     continue
47    continue
c
c     Searching with Steepest Descent Method with c/a values
c
      eps=10000
      iter=0
      latno=100
      minenergy=emini
      minold=0.0
      latold=0.0
      hnum=0.001
      latstep=0.01
      if (prntflag.eq.1) then
        write(*,*)'   >>>Minimum is near:'
        write(*,*)'      lata:',lata,'c:',ctest
      endif
      do 5006 j=1,latno
      latu=lata
      latd=ctest
      latu1=latu+2*hnum
      latu2=latu+hnum
      latu3=latu-hnum
      latu4=latu-2*hnum
      call calcminehex(totu1,latu1,lata0,latd,c0,neigh,rdis,rdisd,
     &rmss)
        if (rmss.gt.10000.0) then
            return
        endif
      call calcminehex(totu2,latu2,lata0,latd,c0,neigh,rdis,rdisd,
     &rmss)
        if (rmss.gt.10000.0) then
            return
        endif
      call calcminehex(totu3,latu3,lata0,latd,c0,neigh,rdis,rdisd,
     &rmss)
        if (rmss.gt.10000.0) then
            return
        endif
      call calcminehex(totu4,latu4,lata0,latd,c0,neigh,rdis,rdisd,
     &rmss)
        if (rmss.gt.10000.0) then
            return
        endif
      latd1=latd+2*hnum
      latd2=latd+hnum
      latd3=latd-hnum
      latd4=latd-2*hnum
      call calcminehex(totd1,latu,lata0,latd1,c0,neigh,rdis,rdisd,
     &rmss)
        if (rmss.gt.10000.0) then
            return
        endif
      call calcminehex(totd2,latu,lata0,latd2,c0,neigh,rdis,rdisd,
     &rmss)
        if (rmss.gt.10000.0) then
            return
        endif
      call calcminehex(totd3,latu,lata0,latd3,c0,neigh,rdis,rdisd,
     &rmss)
        if (rmss.gt.10000.0) then
            return
        endif
      call calcminehex(totd4,latu,lata0,latd4,c0,neigh,rdis,rdisd,
     &rmss)
        if (rmss.gt.10000.0) then
            return
        endif
      derivu=(totu4-8.0*totu3+8.0*totu2-totu1)/(12.0*hnum)
      derivd=(totd4-8.0*totd3+8.0*totd2-totd1)/(12.0*hnum)
c     write(*,*)'      deriv up:',derivu
c     write(*,*)'      deriv down:',derivd
      eps=abs(derivu)+abs(derivd)
      latu=latu-latstep*derivu
      latd=latd-latstep*derivd
      call calcminehex(minenergy,latu,lata0,latd,c0,neigh,rdis,rdisd,
     &rmss)
        if (rmss.gt.10000.0) then
            return
        endif
c        write(*,*)'   >>>----ITERATING-----'
c        write(*,*)'      iter:',i
c        write(*,*)'      latis:',lata
c        write(*,*)'      minenergy:',minenergy
c        write(*,*)'      old energy:',minold
         latis=latu
         if(minenergy.lt.minold) then
           if (prntflag.eq.1) then
            write(*,*)'      Aprouching to minimum energy.'
            write(*,*)'      New:',minenergy,'< Old:',minold
           endif
         else
           if (prntflag.eq.1) then
            write(*,*)'      Minimum energy is exceeded:'
            write(*,*)'      Old:',minold,'< New:',minenergy
           endif
            hnum=0.000001
            latstep=latstep*0.01
         endif
         if(minenergy.lt.emini) then
                emini=minenergy
                emina=latu
                eminc=latd
         endif
         if(eps.lt.1.0E-4) then
          if (prntflag.eq.1) then
           if(minenergy.gt.minold) then
              write(*,*)'      Warning: Minimum energy is exceeded.'
              write(*,*)'      Old:',minold,'< New:',minenergy
           else
              write(*,*)'      Minimum energy found:',minenergy
           endif
          endif
           latold=latis
           minold=minenergy
           exit
         endif
         latold=latis
         minold=minenergy
5006  continue
        iter=j
        if (prntflag.eq.1) then
         write(*,*)'...iteration:',iter
         write(*,*)'Bulk search about min point:',emina,'c',eminc
        endif
      latai=emina+0.001
      lataf=emina-0.001
      latno=11
      latstep=abs(latai-lataf)/(latno-1)
      cai=eminc+0.001
      caf=eminc-0.001
      cno=11
      cstep=abs(cai-caf)/(cno-1)
      do 67 i=1,latno
      mina=latai-(i-1)*latstep
       do 68 j=1,cno
        minc=cai-(j-1)*cstep
        call calcminehex(totenergy,mina,lata0,minc,c0,neigh,rdis,rdisd,
     &rmss)
        if (rmss.gt.10000.0) then
            return
        endif
        if(totenergy.lt.minenergy) then
                minenergy=totenergy
                latis=mina
                c=minc
        endif
68     continue
67    continue
      if (latis.lt.0.5) then
          latis=2.0
      endif
      if (c.lt.0.3) then
          c=0.5
      endif
        if (prntflag.eq.1) then
         write(*,*)'----SEARCH RESULTS-----'
         write(*,*)'  latis:',latis
         write(*,*)'  c/a:',c
         write(*,*)'  minenergy:',minenergy
         write(*,*)'  First found minimum energy latis:',emina
         write(*,*)'  First found minimum energy c:',eminc
         write(*,*)'  First found minimum energy:',emini
         write(*,*)'  Bisection Method latis:',latold
         write(*,*)'  Bisection Method energy:',minold
        endif

      return
      end
c*****************************************************


c****************************************************
c*   Calculation of Equilibrium Lattice Constant    *
c*                  for an Alloy                    *
c****************************************************
      subroutine calclathex(mine,latis,lata0,c,ccc,neigh,rdis,rdisd)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision mine,minenergy,latis,lata0,rdis(nmax,nmax,2)
      double precision lattest,c,ccc,rdisd(nmax,nmax,3)
      integer neigh(nmax)
      common /grids/ frhoar(ngridar,nelmax),frhoar1(ngridar,nelmax),
     .  frhoar2(ngridar,nelmax),frhoar3(ngridar,nelmax),
     .  frhoar4(ngridar,nelmax),frhoar5(ngridar,nelmax),
     .  frhoar6(ngridar,nelmax),rhorar(ngridar,nelmax),
     .  rhorar1(ngridar,nelmax),rhorar2(ngridar,nelmax),
     .  rhorar3(ngridar,nelmax),rhorar4(ngridar,nelmax),
     .  rhorar5(ngridar,nelmax),rhorar6(ngridar,nelmax),
     .  z2rar(ngridar,nelmax,nelmax),z2rar1(ngridar,nelmax,nelmax),
     .  z2rar2(ngridar,nelmax,nelmax),z2rar3(ngridar,nelmax,nelmax),
     .  z2rar4(ngridar,nelmax,nelmax),z2rar5(ngridar,nelmax,nelmax),
     .  z2rar6(ngridar,nelmax,nelmax),drhoar,drar,nrhoar,nrar
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      double precision lrho,latfrho,lphi,latphi,totenergy
      double precision latai,lataf,latstep,lata,cstep,ctest,cold
      real rij,pij,latno
      double precision rdrar,rdrhoar,z2ij,rdist
      integer kij,phitype,iat,j,jat,ity,jty
     
      cstep=0.1
      ctest=ccc
      mine=0.0
      minetest=0.0
      minenergy=100000.0
      latai=6.0
      lataf=1.0
      lattest=lataf
      latno=5001
      latstep=abs(latai-lataf)/(latno-1)
      rdrar = 1.0/drar 
      rdrhoar = 1.0/drhoar
      do 4007 i=1,latno
      lata=latai-(i-1)*latstep
      totenergy=0.0
      do 47 iat=1,natoms
      lrho=0.d0
      latphi=0.d0
      do 45 j=1,neigh(iat)
       rdist=((lata/lata0)*rdisd(iat,j,1))**2
       rdist=rdist+((lata/lata0)*rdisd(iat,j,2))**2
       rdist=rdist+((ctest/ccc)*rdisd(iat,j,3))**2
       rdist=sqrt(rdist)
       jat=rdis(iat,j,2)
c
        pij = rdist*rdrar + 1.0
        kij = pij
        kij = min0(kij,nrar-1)
        pij = pij - kij
        pij = amin1(pij,1.0)
        ity=itype(iat)
        jty=itype(jat)
        lrho = lrho +
     $  ((rhorar3(kij,jty)*pij+rhorar2(kij,jty))*pij+
     $    rhorar1(kij,jty))*pij+rhorar(kij,jty)
        lphi = ((z2rar3(kij,ity,jty)*pij+z2rar2(kij,ity,jty))*pij+
     $    z2rar1(kij,ity,jty))*pij+z2rar(kij,ity,jty)
          latphi=latphi+0.5*lphi
45      continue
          pij = lrho*rdrhoar + 1.0
          kij = pij
          kij = max(1,min(kij,nrhoar-1))
          pij = pij - kij
        latfrho = ((frhoar3(kij,ity)*pij+
     $        frhoar2(kij,ity))*pij+
     $        frhoar1(kij,ity))*pij+
     $        frhoar(kij,ity)
        totenergy = totenergy + latfrho + latphi
47      continue
        if(totenergy.lt.minenergy) then
                minenergy=totenergy
                lattest=lata
        endif
4007    continue
           latis=lattest
           c=ctest
      return
      end
c*****************************************************



      subroutine loadprimcell(pcell,ptyp,patm,pvec,stacks)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision pvec(3),pcell(3,32)
      integer natoms,l1,t1,t2,patm,ptyp(32)
      common /particle/ rv(6,nmax)
      common /unitcell/ ucell(nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      integer stacks(3),id,i,j,ii,jj,mx,my,mz
      double precision a,b,c,mins(3),xx(2),yy(2),zz(2)
c     ===========================================
c         
c       Input Parameters:
c       -----------------
c       pcell: Primitive cell atom coordinates in scaled unit cell length
c       patm : Number of atoms in primitive cell
c       pvec : Primitive Cell Vectors
c       zlc  : Lattice constant
c       stacks : Unit cell replica numbers in x,y,z directions
c
c       zlc:  lattice constant in angstrom (A)
c     '      metal  lattice constant (A)'
c     '      =====  ===================='
c     '      Ag    : 4.09               '
c     '      Au    : 4.08               '
c     '      Cu    : 3.615              '
c     '      Ni    : 3.52               '
c     '      Pd    : 3.89               '
c     '      Pt    : 3.92               '
c     '      Al    : 4.04               '
c
c       stacks(1): # of chains in x direction
c       stacks(2): # of atoms in y direction in each chain
c       stacks(3): # of layers in z direction
c
c       Outputs:
c       --------
c       natoms: Number of atoms in structure
c       Periodic boundary conditions:
c         perub(1-3): Upper bounds for 1=x,2=y,3=z
c         perlb(1-3): Lower bounds for 1=x,2=y,3=z
c       rv(6,natoms): x,y,z,vx,vy,vz of each atom
c       itype(natoms): atom type for each atom
c
c
c     ===========================================
       mx=stacks(1)
       my=stacks(2)
       mz=stacks(3)
       if (mx .lt. 1) mx=1
       if (my .lt. 1) my=1
       if (mz .lt. 1) mz=1
       
       mins(1)=0.
       mins(2)=0.
       mins(3)=0.
       a=pvec(1)
       b=pvec(2)
       c=pvec(3)
       do 10 ii=1,patm
         if (pcell(1,ii).lt.mins(1)) mins(1)=pcell(1,ii)
         if (pcell(2,ii).lt.mins(2)) mins(2)=pcell(2,ii)
         if (pcell(3,ii).lt.mins(3)) mins(3)=pcell(3,ii)
10     continue
       xx(1)=mins(1)
       xx(2)=a
       yy(1)=mins(2)
       yy(2)=b
       zz(1)=mins(3)
       zz(2)=c
       do 20 ii=1,patm
         if ( DABS(pcell(1,ii)-xx(1)) .lt. xx(2) .and.
     .                    pcell(1,ii) .gt. xx(1) ) then
              xx(2)=pcell(1,ii)-xx(1)
              xx(1)=pcell(1,ii)
         endif
         if ( DABS(pcell(2,ii)-yy(1)) .lt. yy(2) .and.
     .                    pcell(2,ii) .gt. yy(1) ) then
              yy(2)=pcell(2,ii)-yy(1)
              yy(1)=pcell(2,ii)
         endif
         if ( DABS(pcell(3,ii)-zz(1)) .lt. zz(2) .and.
     .                    pcell(3,ii) .gt. zz(1) ) then
              zz(2)=pcell(3,ii)-zz(1)
              zz(1)=pcell(3,ii)
         endif
20     continue
      perlb(1)=-xx(2)*0.25;
      perub(1)=a*(mx-1)+(a-xx(2)*0.25);
      perlb(2)=-yy(2)*0.25;
      perub(2)=b*(my-1)+(b-yy(2)*0.25);
      perlb(3)=-zz(2)*0.25;
      perub(3)=c*(mz-1)+(c-zz(2)*0.25);
      id=1
       do 30 k=1,mz
       do 40 i=1,mx
       do 50 j=1,my
       do 60 ii=1,patm
          itype(id)=ptyp(ii)
          ucell(id)=ptyp(ii)
          rv(1,id)=(pcell(1,ii)+i-1)*a
          rv(2,id)=(pcell(2,ii)+j-1)*b
          rv(3,id)=(pcell(3,ii)+k-1)*c
          rv(4,id)=0.
          rv(5,id)=0.
          rv(6,id)=0.
          id=id+1
60     continue
50     continue
40     continue
30     continue
       natoms=id-1
      return
      end


      subroutine loadfcc100(t1,t2,l1,zlc,stacks)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision zlc
      integer natoms,l1,t1,t2
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      integer stacks(3),nchain,nyatom,nlayer,izcnt
      integer ilayer,id,i,j,jj
      double precision znnd,height,dz,ychain
c     ===========================================
c         
c       Input Parameters:
c       -----------------
c       l1:  0 = Non L1 structure single type. All types=t1
c            1 = L1_0 structure for alloys (start with corner) 
c                          types=t1(corner+face),t2(facelayer)
c            2 = L1_1 structure for alloys
c                          types=t1(corner111),t2(other111)
c            3 = L1_2 structure for alloys 
c                          types=t1(corner),t2(face-center)
c            4 = L1_3 structure for alloys 
c                          types=t1(corner111),t2(other111)
c            5+ = Not used. Behave as 0; all types=t1 
c
c       zlc:  lattice constant in angstrom (A)
c     '      metal  lattice constant (A)'
c     '      =====  ===================='
c     '      Ag    : 4.09               '
c     '      Au    : 4.08               '
c     '      Cu    : 3.615              '
c     '      Ni    : 3.52               '
c     '      Pd    : 3.89               '
c     '      Pt    : 3.92               '
c     '      Al    : 4.04               '
c
c       stacks(1): # of chains in x direction
c       stacks(2): # of atoms in y direction in each chain
c       stacks(3): # of layers in z direction
c
c       Outputs:
c       --------
c       natoms: Number of atoms in structure
c       Periodic boundary conditions:
c         perub(1-3): Upper bounds for 1=x,2=y,3=z
c         perlb(1-3): Lower bounds for 1=x,2=y,3=z
c       rv(6,natoms): x,y,z,vx,vy,vz of each atom
c       itype(natoms): atom type for each atom
c
c
c     ===========================================
       nchain=stacks(1)
       nyatom=stacks(2)
       nlayer=stacks(3)
       if (l1.gt.4) then
           l1=0
       endif
c
c      FCC structure: x displacement lat/2
c                     y displacement lat
c                     z displacement lat/2
       znnd=zlc/2.
       height=zlc
       dz=zlc/2.
       id=0
       perlb(1)=-(1./4.)*znnd
       perub(1)=(3./4.)*znnd+(nchain-1)*znnd
       perlb(2)=-(1./3.)*height
       perub(2)=(2./3.)*height+(nyatom-1)*height
       perlb(3)=-(2./3.)*dz-(nlayer-1)*dz
       perub(3)=(1./3.)*dz
       izcnt=0
       do 10 ilayer=1,nlayer
c      construct the type 1 layers
          do 20 i=1,nchain
          ychain=MOD((ilayer-1)+i+1,2)*height/2.
             do 30 j=1,nyatom
                  id=id+1
                  rv(1,id)=(i-1)*znnd
                  rv(2,id)=(j-1)*height+ychain
                  rv(3,id)=-izcnt*dz
                  rv(4,id)=0.0
                  rv(5,id)=0.0
                  rv(6,id)=0.0
                  if(l1.gt.0) then
c                 L1 structures                          
                    if(l1.lt.2) then
c                   L1_0 structure                          
                      if(MOD(ilayer,2).eq.1) then
                        itype(id)=t1
                      else
                        itype(id)=t2
                      endif
                    else if (l1.lt.3) then
c                   L1_1 structure                          
                     if(MOD(ilayer,4).eq.1.or.MOD(ilayer,4).eq.3) then
                        if(MOD(ilayer,4).eq.1) then
                          jj=j
                        else
                          jj=j+1
                        endif
                        if(MOD(i,4).eq.1.or.MOD(i,4).eq.0) then
                          jj=jj
                        else
                          jj=jj+1
                        endif
                     else
                        if(MOD(ilayer,4).eq.2) then
                          jj=j
                        else
                          jj=j+1
                        endif
                        if(MOD(i,4).eq.1.or.MOD(i,4).eq.2) then
                          jj=jj
                        else
                          jj=jj+1
                        endif
                     endif
                     if(MOD(jj,2).eq.1) then
                        itype(id)=t1
                     else
                        itype(id)=t2
                     endif
                    else if (l1.lt.4) then
c                   L1_2 structure                          
                      if(MOD(ilayer,2).eq.1) then
                        if(MOD(i+1,2).eq.1) then
                           itype(id)=t2
                        else
                           itype(id)=t1
                        endif
                      else
                        itype(id)=t2
                      endif
                    else if (l1.lt.5) then
c                   L1_3 structure                          
                     if(MOD(ilayer,4).eq.1.or.MOD(ilayer,4).eq.3) then
                        if(MOD(ilayer,4).eq.1) then
                          jj=j
                        else
                          jj=j+1
                        endif
                        if(MOD(i,4).eq.1.or.MOD(i,4).eq.0) then
                          jj=jj
                        else
                          jj=jj+1
                        endif
                        if(MOD(jj,2).eq.1) then
                           itype(id)=t2
                        else
                           itype(id)=t1
                        endif
                     else
                        itype(id)=t2
                     endif
                    endif
                  else
c                 Single type FCC structure                          
                    itype(id)=t1
                  endif
30            continue
20        continue
          izcnt=izcnt+1
10      continue
c       count the number of atom
        natoms=id
      return
      end


      subroutine loadfcc100nnd(t1,t2,l1,zlc,stacks)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision zlc
      integer natoms,l1,t1,t2
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      integer stacks(3),nchain,nyatom,nlayer,izcnt
      integer ilayer,id,i,j
      double precision znnd,height,dz,ychain,sqrt2
c     ===========================================
c         
c       Input Parameters:
c       -----------------
c       l1:  0 = L1 structure for alloys (Ni3Al) types=1(corner),2(face)
c            1 = FCC structure for pure metals all types=1
c            2 = FCC structure for pure metals all types=2
c            any non zero = FCC structure all types=non zero number
c
c       zlc:  lattice constant in angstrom (A)
c     '      metal  lattice constant (A)'
c     '      =====  ===================='
c     '      Ag    : 4.09               '
c     '      Au    : 4.08               '
c     '      Cu    : 3.615              '
c     '      Ni    : 3.52               '
c     '      Pd    : 3.89               '
c     '      Pt    : 3.92               '
c     '      Al    : 4.04               '
c
c       stacks(1): # of chains in x direction
c       stacks(2): # of atoms in y direction in each chain
c       stacks(3): # of layers in z direction
c
c       Outputs:
c       --------
c       natoms: Number of atoms in structure
c       Periodic boundary conditions:
c         perub(1-3): Upper bounds for 1=x,2=y,3=z
c         perlb(1-3): Lower bounds for 1=x,2=y,3=z
c       rv(6,natoms): x,y,z,vx,vy,vz of each atom
c       itype(natoms): atom type for each atom
c
c
c     ===========================================
       nchain=stacks(1)
       nyatom=stacks(2)
       nlayer=stacks(3)
c
c      FCC structure: x displacement lat/2
c                     y displacement lat
c                     z displacement lat/2
       sqrt2=1.4142135623730950488016887242096980785696718753769480731767
       znnd=zlc/sqrt2
       height=zlc/sqrt2
       dz=zlc/2.
       id=0
       perlb(1)=-(1./4.)*znnd
       perub(1)=(3./4.)*znnd+(nchain-1)*znnd
       perlb(2)=-(1./3.)*height
       perub(2)=(2./3.)*height+(nyatom-1)*height
       perlb(3)=-(2./3.)*dz-(nlayer-1)*dz
       perub(3)=(1./3.)*dz
       izcnt=0
       do 10 ilayer=1,nlayer
c      construct the type 1 layers
          xchain=MOD(ilayer+1,2)*znnd/2.
          ychain=MOD(ilayer+1,2)*height/2.
          do 20 i=1,nchain
             do 30 j=1,nyatom
                  id=id+1
                  rv(1,id)=(i-1)*znnd+xchain
                  rv(2,id)=(j-1)*height+ychain
                  rv(3,id)=-izcnt*dz
                  rv(4,id)=0.0
                  rv(5,id)=0.0
                  rv(6,id)=0.0
                  if(l1.lt.1) then
                    if(MOD(ilayer,2)) then
                      if(MOD(i,2)) then
                        if(2-MOD(j+1,2).eq.1) then
                           itype(id)=t2
                        else
                           itype(id)=t1
                        endif
                      else
                        if(2-MOD(j,2).eq.1) then
                           itype(id)=t2
                        else
                           itype(id)=t1
                        endif
                      endif
                    else
                      itype(id)=t2
                    endif
                  else
                    itype(id)=l1
                  endif
30            continue
20        continue
          izcnt=izcnt+1
10      continue
c       count the number of atom
        natoms=id
      return
      end


      subroutine loadfcc111(t1,t2,l1,zlc,stacks)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision zlc
      integer natoms,l1,t1,t2
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      integer stacks(3),nchain,nyatom,nlayer,izcnt
      integer ilayer,id,i,j,ii,jj
      double precision znnd,height,dz,ychain,xlayer,sqrt3,sqrt4
c     ===========================================
c         
c       Input Parameters:
c       -----------------
c       L1:  0 = Pure metals all types=t1
c            1 = L1_0 alloy structure types=t1(corner),t2(face)
c            2 = L1_1 alloy structure types=t1(layer1),t2(layer2)
c            3 = L1_2 alloy structure types=t1(corner),t2(face)
c            4 = L1_3 alloy structure types=t1(cornerlayer),t2(alloylayer)
c            5+ = Pure metals all types=t1 (same as L1=0)
c
c       zlc:  lattice constant in angstrom (A)
c     '      metal  lattice constant (A)'
c     '      =====  ===================='
c     '      Ag    : 4.09               '
c     '      Au    : 4.08               '
c     '      Cu    : 3.615              '
c     '      Ni    : 3.52               '
c     '      Pd    : 3.89               '
c     '      Pt    : 3.92               '
c     '      Al    : 4.04               '
c
c       stacks(1): # of chains in x direction
c       stacks(2): # of atoms in y direction in each chain
c       stacks(3): # of layers in z direction
c
c       Outputs:
c       --------
c       natoms: Number of atoms in structure
c       Periodic boundary conditions:
c         perub(1-3): Upper bounds for 1=x,2=y,3=z
c         perlb(1-3): Lower bounds for 1=x,2=y,3=z
c       rv(6,natoms): x,y,z,vx,vy,vz of each atom
c       itype(natoms): atom type for each atom
c
c
c     ===========================================
       nchain=stacks(1)
       nyatom=stacks(2)
       nlayer=stacks(3)
       if (l1.gt.4) then
           l1=0
       endif
c
c      FCC 111 stacking structure: x displacement lat*sqrt(3/8)
c                                  y displacement lat/sqrt(2)
c                                  z displacement lat/sqrt(3)
       sqrt2=1.4142135623730950488016887242096980785696718753769480731767
       sqrt3=1.7320508075688772935274463415058723669428052538103806280558
       dz=zlc/sqrt3
       height=zlc/sqrt2
       znnd=zlc*sqrt3/(2.*sqrt2)
       id=0
       perlb(1)=-(1./4.)*znnd
       perub(1)=(3./4.)*znnd+(nchain-1)*znnd
       perlb(2)=-(1./3.)*height
       perub(2)=(2./3.)*height+(nyatom-1)*height
       perlb(3)=-(2./3.)*dz-(nlayer-1)*dz
       perub(3)=(1./3.)*dz
       izcnt=0
       do 10 ilayer=1,nlayer
       if(MOD(ilayer,3).eq.1) then
c      construct the type 1 layers
          izcnt=izcnt+1
          do 20 i=1,nchain
             ychain=MOD(i+1,2)*height/2.
             do 30 j=1,nyatom
                  id=id+1
                  rv(1,id)=(i-1)*znnd
                  rv(2,id)=(j-1)*height+ychain
                  rv(3,id)=-(izcnt-1)*dz
                  rv(4,id)=0.0
                  rv(5,id)=0.0
                  rv(6,id)=0.0
                  if(l1.gt.0) then
                   if(l1.lt.2) then
c                  L10 structure
                    if(MOD(i,2)) then
                      itype(id)=t1
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.3) then
c                  L11 structure
                    if(MOD(ilayer,2).eq.1) then
                      itype(id)=t1
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.4) then
c                  L12 structure
                    if(MOD(i,2)) then
                      if(MOD(i,4).eq.1) then
                         jj=2-MOD(j+1,2)
                      else if(MOD(i,4).eq.3) then
                         jj=2-MOD(j,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.5) then
c                  L13 structure
                    if(MOD(ilayer,6).eq.4) then
                       ii=i+2
                       jj=j
                    else
                       ii=i
                       jj=j
                    endif
                    if(MOD(ii,2)) then
                      if(MOD(ii,4).eq.1) then
                         jj=2-MOD(jj+1,2)
                      else if(MOD(ii,4).eq.3) then
                         jj=2-MOD(jj,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    else
                      itype(id)=t2
                    endif
                   else
                    itype(id)=t1
                   endif
                  else
                    itype(id)=t1
                  endif
30            continue
20        continue
        elseif(MOD(ilayer,3).eq.2) then
c       construct the type 2 layers which fill every other threefold hollow
c       of the type 1 layers.
          izcnt=izcnt+1
          do 21 i=1,nchain
             ychain=MOD(i,2)*height/2.
             xlayer=znnd/3.
             do 31 j=1,nyatom
                  id=id+1
                  rv(1,id)=(i-1)*znnd+xlayer
                  rv(2,id)=(j-1)*height+ychain
                  rv(3,id)=-(izcnt-1)*dz
                  rv(4,id)=0.0
                  rv(5,id)=0.0
                  rv(6,id)=0.0
                  if(l1.gt.0) then
                   if(l1.lt.2) then
c                  L10 structure
                    if(MOD(i,2)) then
                      itype(id)=t2
                    else
                      itype(id)=t1
                    endif
                   else if(l1.lt.3) then
c                  L11 structure
                    if(MOD(ilayer,2).eq.1) then
                      itype(id)=t1
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.4) then
c                  L12 structure
                    if(MOD(i,2)) then
                      itype(id)=t2
                    else
                      if(MOD(i,4).eq.2) then
                         jj=2-MOD(j+1,2)
                      else if(MOD(i,4).eq.0) then
                         jj=2-MOD(j,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    endif
                   else if(l1.lt.5) then
c                  L13 structure
                    if(MOD(ilayer,6).eq.5) then
                       ii=i+1
                       jj=j+1
                    else
                       ii=i+1
                       jj=j
                    endif
                    if(MOD(ii,2)) then
                      if(MOD(ii,4).eq.1) then
                         jj=2-MOD(jj+1,2)
                      else if(MOD(ii,4).eq.3) then
                         jj=2-MOD(jj,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    else
                      itype(id)=t2
                    endif
                   else
                    itype(id)=t1
                   endif
                  else
                    itype(id)=t1
                  endif
31            continue
21        continue
        else
c       construct the type 3 layers which fill every other threefold hollow
c       of type 2 layers.
          izcnt=izcnt+1
          do 22 i=1,nchain
             ychain=MOD(i+1,2)*height/2.
             xlayer=2*znnd/3.
             do 32 j=1,nyatom
                  id=id+1
                  rv(1,id)=(i-1)*znnd+xlayer
                  rv(2,id)=(j-1)*height+ychain
                  rv(3,id)=-(izcnt-1)*dz
                  rv(4,id)=0.0
                  rv(5,id)=0.0
                  rv(6,id)=0.0
                  if(l1.gt.0) then
                   if(l1.lt.2) then
c                  L10 structure
                    if(MOD(i,2)) then
                      itype(id)=t1
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.3) then
c                  L11 structure
                    if(MOD(ilayer,2).eq.1) then
                      itype(id)=t1
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.4) then
c                  L12 structure
                    if(MOD(i,2)) then
                      if(MOD(i,4).eq.1) then
                         jj=2-MOD(j,2)
                      else if(MOD(i,4).eq.3) then
                         jj=2-MOD(j+1,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.5) then
c                  L13 structure
                    if(MOD(ilayer,6).eq.0) then
                       ii=i
                       jj=j
                    else
                       ii=i
                       jj=j+1
                    endif
                    if(MOD(ii,2)) then
                      if(MOD(ii,4).eq.1) then
                         jj=2-MOD(jj+1,2)
                      else if(MOD(ii,4).eq.3) then
                         jj=2-MOD(jj,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    else
                      itype(id)=t2
                    endif
                   else
                    itype(id)=t1
                   endif
                  else
                    itype(id)=t1
                  endif
32            continue
22        continue
        endif
10      continue
c       count the number of atom
        natoms=id
      return
      end


      subroutine loadsfault111(t1,t2,l1,zlc,stacks,zplns)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision zlc
      integer natoms,l1,t1,t2
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      integer stacks(3),nchain,nyatom,nlayer,izcnt
      integer ilayer,id,i,j,ii,jj
      integer iilayer,zplns(2)
      double precision znnd,height,dz,ychain,xlayer,sqrt3,sqrt4
c     ===========================================
c         
c       Input Parameters:
c       -----------------
c       L1:  0 = Pure metals all types=t1
c            1 = L1_0 alloy structure types=t1(corner),t2(face)
c            2 = L1_1 alloy structure types=t1(layer1),t2(layer2)
c            3 = L1_2 alloy structure types=t1(corner),t2(face)
c            4 = L1_3 alloy structure types=t1(cornerlayer),t2(alloylayer)
c            5+ = Pure metals all types=t1 (same as L1=0)
c
c       zlc:  lattice constant in angstrom (A)
c     '      metal  lattice constant (A)'
c     '      =====  ===================='
c     '      Ag    : 4.09               '
c     '      Au    : 4.08               '
c     '      Cu    : 3.615              '
c     '      Ni    : 3.52               '
c     '      Pd    : 3.89               '
c     '      Pt    : 3.92               '
c     '      Al    : 4.04               '
c
c       stacks(1): # of chains in x direction
c       stacks(2): # of atoms in y direction in each chain
c       stacks(3): # of layers in z direction
c
c       Outputs:
c       --------
c       natoms: Number of atoms in structure
c       Periodic boundary conditions:
c         perub(1-3): Upper bounds for 1=x,2=y,3=z
c         perlb(1-3): Lower bounds for 1=x,2=y,3=z
c       rv(6,natoms): x,y,z,vx,vy,vz of each atom
c       itype(natoms): atom type for each atom
c
c
c     ===========================================
       nchain=stacks(1)
       nyatom=stacks(2)
       nlayer=stacks(3)
       if (l1.gt.4) then
           l1=0
       endif
c
c      FCC 111 stacking structure: x displacement lat*sqrt(3/8)
c                                  y displacement lat/sqrt(2)
c                                  z displacement lat/sqrt(3)
       sqrt2=1.4142135623730950488016887242096980785696718753769480731767
       sqrt3=1.7320508075688772935274463415058723669428052538103806280558
       dz=zlc/sqrt3
       height=zlc/sqrt2
       znnd=zlc*sqrt3/(2.*sqrt2)
       id=0
       perlb(1)=-(1./4.)*znnd
       perub(1)=(3./4.)*znnd+(nchain-1)*znnd
       perlb(2)=-(1./3.)*height
       perub(2)=(2./3.)*height+(nyatom-1)*height
       perlb(3)=-(2./3.)*dz-(nlayer-1)*dz
       perub(3)=(1./3.)*dz
       izcnt=0
       do 10 iilayer=1,nlayer
        if(iilayer.gt.zplns(2)) then
         ilayer=iilayer-1
        else if(iilayer.ge.zplns(1)) then
         ilayer=iilayer+1
        else
         ilayer=iilayer
        endif
       if(MOD(ilayer,3).eq.1) then
c      construct the type 1 layers
          izcnt=izcnt+1
          do 20 i=1,nchain
             ychain=MOD(i+1,2)*height/2.
             do 30 j=1,nyatom
                  id=id+1
                  rv(1,id)=(i-1)*znnd
                  rv(2,id)=(j-1)*height+ychain
                  rv(3,id)=-(izcnt-1)*dz
                  rv(4,id)=0.0
                  rv(5,id)=0.0
                  rv(6,id)=0.0
                  if(l1.gt.0) then
                   if(l1.lt.2) then
c                  L10 structure
                    if(MOD(i,2)) then
                      itype(id)=t1
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.3) then
c                  L11 structure
                    if(MOD(ilayer,2).eq.1) then
                      itype(id)=t1
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.4) then
c                  L12 structure
                    if(MOD(i,2)) then
                      if(MOD(i,4).eq.1) then
                         jj=2-MOD(j+1,2)
                      else if(MOD(i,4).eq.3) then
                         jj=2-MOD(j,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.5) then
c                  L13 structure
                    if(MOD(ilayer,6).eq.4) then
                       ii=i+2
                       jj=j
                    else
                       ii=i
                       jj=j
                    endif
                    if(MOD(ii,2)) then
                      if(MOD(ii,4).eq.1) then
                         jj=2-MOD(jj+1,2)
                      else if(MOD(ii,4).eq.3) then
                         jj=2-MOD(jj,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    else
                      itype(id)=t2
                    endif
                   else
                    itype(id)=t1
                   endif
                  else
                    itype(id)=t1
                  endif
30            continue
20        continue
        elseif(MOD(ilayer,3).eq.2) then
c       construct the type 2 layers which fill every other threefold hollow
c       of the type 1 layers.
          izcnt=izcnt+1
          do 21 i=1,nchain
             ychain=MOD(i,2)*height/2.
             xlayer=znnd/3.
             do 31 j=1,nyatom
                  id=id+1
                  rv(1,id)=(i-1)*znnd+xlayer
                  rv(2,id)=(j-1)*height+ychain
                  rv(3,id)=-(izcnt-1)*dz
                  rv(4,id)=0.0
                  rv(5,id)=0.0
                  rv(6,id)=0.0
                  if(l1.gt.0) then
                   if(l1.lt.2) then
c                  L10 structure
                    if(MOD(i,2)) then
                      itype(id)=t2
                    else
                      itype(id)=t1
                    endif
                   else if(l1.lt.3) then
c                  L11 structure
                    if(MOD(ilayer,2).eq.1) then
                      itype(id)=t1
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.4) then
c                  L12 structure
                    if(MOD(i,2)) then
                      itype(id)=t2
                    else
                      if(MOD(i,4).eq.2) then
                         jj=2-MOD(j+1,2)
                      else if(MOD(i,4).eq.0) then
                         jj=2-MOD(j,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    endif
                   else if(l1.lt.5) then
c                  L13 structure
                    if(MOD(ilayer,6).eq.5) then
                       ii=i+1
                       jj=j+1
                    else
                       ii=i+1
                       jj=j
                    endif
                    if(MOD(ii,2)) then
                      if(MOD(ii,4).eq.1) then
                         jj=2-MOD(jj+1,2)
                      else if(MOD(ii,4).eq.3) then
                         jj=2-MOD(jj,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    else
                      itype(id)=t2
                    endif
                   else
                    itype(id)=t1
                   endif
                  else
                    itype(id)=t1
                  endif
31            continue
21        continue
        else
c       construct the type 3 layers which fill every other threefold hollow
c       of type 2 layers.
          izcnt=izcnt+1
          do 22 i=1,nchain
             ychain=MOD(i+1,2)*height/2.
             xlayer=2*znnd/3.
             do 32 j=1,nyatom
                  id=id+1
                  rv(1,id)=(i-1)*znnd+xlayer
                  rv(2,id)=(j-1)*height+ychain
                  rv(3,id)=-(izcnt-1)*dz
                  rv(4,id)=0.0
                  rv(5,id)=0.0
                  rv(6,id)=0.0
                  if(l1.gt.0) then
                   if(l1.lt.2) then
c                  L10 structure
                    if(MOD(i,2)) then
                      itype(id)=t1
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.3) then
c                  L11 structure
                    if(MOD(ilayer,2).eq.1) then
                      itype(id)=t1
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.4) then
c                  L12 structure
                    if(MOD(i,2)) then
                      if(MOD(i,4).eq.1) then
                         jj=2-MOD(j,2)
                      else if(MOD(i,4).eq.3) then
                         jj=2-MOD(j+1,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.5) then
c                  L13 structure
                    if(MOD(ilayer,6).eq.0) then
                       ii=i
                       jj=j
                    else
                       ii=i
                       jj=j+1
                    endif
                    if(MOD(ii,2)) then
                      if(MOD(ii,4).eq.1) then
                         jj=2-MOD(jj+1,2)
                      else if(MOD(ii,4).eq.3) then
                         jj=2-MOD(jj,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    else
                      itype(id)=t2
                    endif
                   else
                    itype(id)=t1
                   endif
                  else
                    itype(id)=t1
                  endif
32            continue
22        continue
        endif
10      continue
c       count the number of atom
        natoms=id
      return
      end

      

      subroutine loadtwin111(t1,t2,l1,zlc,stacks,zplns)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision zlc
      integer natoms,l1,t1,t2
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      integer stacks(3),nchain,nyatom,nlayer,izcnt
      integer ilayer,id,i,j,ii,jj
      integer iilayer,zplns(2)
      double precision znnd,height,dz,ychain,xlayer,sqrt3,sqrt4
c     ===========================================
c         
c       Input Parameters:
c       -----------------
c       L1:  0 = Pure metals all types=t1
c            1 = L1_0 alloy structure types=t1(corner),t2(face)
c            2 = L1_1 alloy structure types=t1(layer1),t2(layer2)
c            3 = L1_2 alloy structure types=t1(corner),t2(face)
c            4 = L1_3 alloy structure types=t1(cornerlayer),t2(alloylayer)
c            5+ = Pure metals all types=t1 (same as L1=0)
c
c       zlc:  lattice constant in angstrom (A)
c     '      metal  lattice constant (A)'
c     '      =====  ===================='
c     '      Ag    : 4.09               '
c     '      Au    : 4.08               '
c     '      Cu    : 3.615              '
c     '      Ni    : 3.52               '
c     '      Pd    : 3.89               '
c     '      Pt    : 3.92               '
c     '      Al    : 4.04               '
c
c       stacks(1): # of chains in x direction
c       stacks(2): # of atoms in y direction in each chain
c       stacks(3): # of layers in z direction
c
c       Outputs:
c       --------
c       natoms: Number of atoms in structure
c       Periodic boundary conditions:
c         perub(1-3): Upper bounds for 1=x,2=y,3=z
c         perlb(1-3): Lower bounds for 1=x,2=y,3=z
c       rv(6,natoms): x,y,z,vx,vy,vz of each atom
c       itype(natoms): atom type for each atom
c
c
c     ===========================================
       nchain=stacks(1)
       nyatom=stacks(2)
       nlayer=stacks(3)
       if (l1.gt.4) then
           l1=0
       endif
c
c      FCC 111 stacking structure: x displacement lat*sqrt(3/8)
c                                  y displacement lat/sqrt(2)
c                                  z displacement lat/sqrt(3)
       sqrt2=1.4142135623730950488016887242096980785696718753769480731767
       sqrt3=1.7320508075688772935274463415058723669428052538103806280558
       dz=zlc/sqrt3
       height=zlc/sqrt2
       znnd=zlc*sqrt3/(2.*sqrt2)
       id=0
       perlb(1)=-(1./4.)*znnd
       perub(1)=(3./4.)*znnd+(nchain-1)*znnd
       perlb(2)=-(1./3.)*height
       perub(2)=(2./3.)*height+(nyatom-1)*height
       perlb(3)=-(2./3.)*dz-(nlayer-1)*dz
       perub(3)=(1./3.)*dz
       izcnt=0
       do 10 iilayer=1,nlayer
       if(iilayer.ge.zplns(1) .and. 
     .    iilayer.le.zplns(2)) then
         ilayer=nlayer-iilayer
       else
         ilayer=iilayer
       endif
       if(MOD(ilayer,3).eq.1) then
c      construct the type 1 layers
          izcnt=izcnt+1
          do 20 i=1,nchain
             ychain=MOD(i+1,2)*height/2.
             do 30 j=1,nyatom
                  id=id+1
                  rv(1,id)=(i-1)*znnd
                  rv(2,id)=(j-1)*height+ychain
                  rv(3,id)=-(izcnt-1)*dz
                  rv(4,id)=0.0
                  rv(5,id)=0.0
                  rv(6,id)=0.0
                  if(l1.gt.0) then
                   if(l1.lt.2) then
c                  L10 structure
                    if(MOD(i,2)) then
                      itype(id)=t1
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.3) then
c                  L11 structure
                    if(MOD(ilayer,2).eq.1) then
                      itype(id)=t1
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.4) then
c                  L12 structure
                    if(MOD(i,2)) then
                      if(MOD(i,4).eq.1) then
                         jj=2-MOD(j+1,2)
                      else if(MOD(i,4).eq.3) then
                         jj=2-MOD(j,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.5) then
c                  L13 structure
                    if(MOD(ilayer,6).eq.4) then
                       ii=i+2
                       jj=j
                    else
                       ii=i
                       jj=j
                    endif
                    if(MOD(ii,2)) then
                      if(MOD(ii,4).eq.1) then
                         jj=2-MOD(jj+1,2)
                      else if(MOD(ii,4).eq.3) then
                         jj=2-MOD(jj,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    else
                      itype(id)=t2
                    endif
                   else
                    itype(id)=t1
                   endif
                  else
                    itype(id)=t1
                  endif
30            continue
20        continue
        elseif(MOD(ilayer,3).eq.2) then
c       construct the type 2 layers which fill every other threefold hollow
c       of the type 1 layers.
          izcnt=izcnt+1
          do 21 i=1,nchain
             ychain=MOD(i,2)*height/2.
             xlayer=znnd/3.
             do 31 j=1,nyatom
                  id=id+1
                  rv(1,id)=(i-1)*znnd+xlayer
                  rv(2,id)=(j-1)*height+ychain
                  rv(3,id)=-(izcnt-1)*dz
                  rv(4,id)=0.0
                  rv(5,id)=0.0
                  rv(6,id)=0.0
                  if(l1.gt.0) then
                   if(l1.lt.2) then
c                  L10 structure
                    if(MOD(i,2)) then
                      itype(id)=t2
                    else
                      itype(id)=t1
                    endif
                   else if(l1.lt.3) then
c                  L11 structure
                    if(MOD(ilayer,2).eq.1) then
                      itype(id)=t1
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.4) then
c                  L12 structure
                    if(MOD(i,2)) then
                      itype(id)=t2
                    else
                      if(MOD(i,4).eq.2) then
                         jj=2-MOD(j+1,2)
                      else if(MOD(i,4).eq.0) then
                         jj=2-MOD(j,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    endif
                   else if(l1.lt.5) then
c                  L13 structure
                    if(MOD(ilayer,6).eq.5) then
                       ii=i+1
                       jj=j+1
                    else
                       ii=i+1
                       jj=j
                    endif
                    if(MOD(ii,2)) then
                      if(MOD(ii,4).eq.1) then
                         jj=2-MOD(jj+1,2)
                      else if(MOD(ii,4).eq.3) then
                         jj=2-MOD(jj,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    else
                      itype(id)=t2
                    endif
                   else
                    itype(id)=t1
                   endif
                  else
                    itype(id)=t1
                  endif
31            continue
21        continue
        else
c       construct the type 3 layers which fill every other threefold hollow
c       of type 2 layers.
          izcnt=izcnt+1
          do 22 i=1,nchain
             ychain=MOD(i+1,2)*height/2.
             xlayer=2*znnd/3.
             do 32 j=1,nyatom
                  id=id+1
                  rv(1,id)=(i-1)*znnd+xlayer
                  rv(2,id)=(j-1)*height+ychain
                  rv(3,id)=-(izcnt-1)*dz
                  rv(4,id)=0.0
                  rv(5,id)=0.0
                  rv(6,id)=0.0
                  if(l1.gt.0) then
                   if(l1.lt.2) then
c                  L10 structure
                    if(MOD(i,2)) then
                      itype(id)=t1
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.3) then
c                  L11 structure
                    if(MOD(ilayer,2).eq.1) then
                      itype(id)=t1
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.4) then
c                  L12 structure
                    if(MOD(i,2)) then
                      if(MOD(i,4).eq.1) then
                         jj=2-MOD(j,2)
                      else if(MOD(i,4).eq.3) then
                         jj=2-MOD(j+1,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    else
                      itype(id)=t2
                    endif
                   else if(l1.lt.5) then
c                  L13 structure
                    if(MOD(ilayer,6).eq.0) then
                       ii=i
                       jj=j
                    else
                       ii=i
                       jj=j+1
                    endif
                    if(MOD(ii,2)) then
                      if(MOD(ii,4).eq.1) then
                         jj=2-MOD(jj+1,2)
                      else if(MOD(ii,4).eq.3) then
                         jj=2-MOD(jj,2)
                      endif
                      if(MOD(jj,2).eq.1) then
                         itype(id)=t2
                      else
                         itype(id)=t1
                      endif
                    else
                      itype(id)=t2
                    endif
                   else
                    itype(id)=t1
                   endif
                  else
                    itype(id)=t1
                  endif
32            continue
22        continue
        endif
10      continue
c       count the number of atom
        natoms=id
      return
      end


      subroutine loadbcc100(t1,t2,b2,zlc,stacks)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision zlc
      integer natoms,b2,t1,t2
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      integer stacks(3),nchain,nyatom,nlayer,izcnt
      integer ilayer,id,i,j
      double precision znnd,height,dz,ychain,xlayer
c     ===========================================
c         
c       Input Parameters:
c       -----------------
c       B2:  0 = BCC structure all types=t1
c            0<= B2 structure for alloys (NiAl) 
c                     types=1(corner),2(center)
c
c       zlc:  lattice constant in angstrom (A)
c     '      metal  lattice constant (A)'
c     '      =====  ===================='
c     '      Ag    : 4.09               '
c     '      Au    : 4.08               '
c     '      Cu    : 3.615              '
c     '      Ni    : 3.52               '
c     '      Pd    : 3.89               '
c     '      Pt    : 3.92               '
c     '      Al    : 4.04               '
c
c       stacks(1): # of chains in x direction
c       stacks(2): # of atoms in y direction in each chain
c       stacks(3): # of layers in z direction
c
c       Outputs:
c       --------
c       natoms: Number of atoms in structure
c       Periodic boundary conditions:
c         perub(1-3): Upper bounds for 1=x,2=y,3=z
c         perlb(1-3): Lower bounds for 1=x,2=y,3=z
c       rv(6,natoms): x,y,z,vx,vy,vz of each atom
c       itype(natoms): atom type for each atom
c
c
c     ===========================================
       nchain=stacks(1)
       nyatom=stacks(2)
       nlayer=stacks(3)
c
c      BCC structure: x displacement lat
c                     y displacement lat
c                     z displacement lat/2
       znnd=zlc
       height=zlc
       dz=zlc/2.
       id=0
       perlb(1)=-(1./4.)*znnd
       perub(1)=(3./4.)*znnd+(nchain-1)*znnd
       perlb(2)=-(1./3.)*height
       perub(2)=(2./3.)*height+(nyatom-1)*height
       perlb(3)=-(2./3.)*dz-(nlayer-1)*dz
       perub(3)=(1./3.)*dz
       izcnt=0
       do 10 ilayer=1,nlayer
c      construct the type 1 layers
          xlayer=MOD(ilayer+1,2)*znnd/2.
          ychain=MOD(ilayer+1,2)*height/2.
          izcnt=izcnt+1
          do 20 i=1,nchain
             do 30 j=1,nyatom
                  id=id+1
                  rv(1,id)=(i-1)*znnd+xlayer
                  rv(2,id)=(j-1)*height+ychain
                  rv(3,id)=-(izcnt-1)*dz
                  rv(4,id)=0.0
                  rv(5,id)=0.0
                  rv(6,id)=0.0
                  if(b2.gt.0) then
                    if(MOD(ilayer,2).eq.1) then
                      itype(id)=t1
                    else
                      itype(id)=t2
                    endif
                  else
                    itype(id)=t1
                  endif
30            continue
20        continue
10      continue
c       count the number of atom
        natoms=id
      return
      end
      
      subroutine loadhcp111(t1,t2,l1,zlc,c,stacks)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision zlc,c
      integer natoms,l1,t1,t2
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      character*3 lat(nelmax),latty
      integer stacks(3),nchain,nyatom,nlayer,izcnt
      integer ilayer,id,i,j
      double precision znnd,height,dz,ychain,xlayer,sqrt3
c     ===========================================
c         
c       Input Parameters:
c       -----------------
c       L1:  0 = L1 structure for hcp alloy types=1(base1),2(base2)
c            1 = HCP structure for pure metals 2 bravais lattice types=1,2 for each A,B layers
c            2 = HCP structure for pure metals all types=t1
c            2<= HCP structure for pure metals all types=t1
c
c       zlc:  lattice constant in angstrom (A)
c     '      metal  lattice constant (A)'
c     '      =====  ===================='
c     '      Ag    : 4.09               '
c     '      Au    : 4.08               '
c     '      Cu    : 3.615              '
c     '      Ni    : 3.52               '
c     '      Pd    : 3.89               '
c     '      Pt    : 3.92               '
c     '      Al    : 4.04               '
c
c       c/a:  c/a constant in angstrom (A)
c     '      metal  c/lattice constant (A)'
c     '      =====  ===================='
c     '      Co    : 1.623              '
c
c
c       stacks(1): # of chains in x direction
c       stacks(2): # of atoms in y direction in each chain
c       stacks(3): # of layers in z direction
c
c       Outputs:
c       --------
c       natoms: Number of atoms in structure
c       Periodic boundary conditions:
c         perub(1-3): Upper bounds for 1=x,2=y,3=z
c         perlb(1-3): Lower bounds for 1=x,2=y,3=z
c       rv(6,natoms): x,y,z,vx,vy,vz of each atom
c       itype(natoms): atom type for each atom
c
c
c     ===========================================
       nchain=stacks(1)
       nyatom=stacks(2)
       nlayer=stacks(3)
c
c      HCP stacking structure: x displacement lat*sqrt(3/4)    lat*sqrt(3/8)
c                              y displacement lat              lat/sqrt(2)
c                              z displacement c                lat/sqrt(3)
       dz=c*zlc/2.
       height=zlc
       sqrt3=1.7320508075688772935274463415058723669428052538103806280558
c      znnd=zlc*sqrt(3.)/2.
       znnd=zlc*sqrt3/2.
       id=0
       perlb(1)=-(1./4.)*znnd
       perub(1)=(3./4.)*znnd+(nchain-1)*znnd
       perlb(2)=-(1./3.)*height
       perub(2)=(2./3.)*height+(nyatom-1)*height
       perlb(3)=-(2./3.)*dz-(nlayer-1)*dz
       perub(3)=(1./3.)*dz
       izcnt=0
       do 10 ilayer=1,nlayer
       if(MOD(ilayer,2).eq.1) then
c      construct the type 1 layers
          izcnt=izcnt+1
          do 20 i=1,nchain
             ychain=MOD(i+1,2)*height/2.
             do 30 j=1,nyatom
                  id=id+1
                  rv(1,id)=(i-1)*znnd
                  rv(2,id)=(j-1)*height+ychain
                  rv(3,id)=-(izcnt-1)*dz
                  rv(4,id)=0.0
                  rv(5,id)=0.0
                  rv(6,id)=0.0
                  if(l1.lt.1) then
                    if(MOD(i,2).eq.1) then
                      if(MOD(i,4).eq.1) then
                        if(2-MOD(j+1,2).eq.1) then
                           itype(id)=t1
                        else
                           itype(id)=t2
                        endif
                      else if(MOD(i,4).eq.3) then
                        if(2-MOD(j,2).eq.1) then
                           itype(id)=t1
                        else
                           itype(id)=t2
                        endif
                      endif
                    else
                      itype(id)=t1
                    endif
                  else
                    if(l1.lt.2) then
                      itype(id)=t1
                    else
                      itype(id)=t1
                    endif
                  endif
30            continue
20        continue
        else
c       construct the type 2 layers which fill every other threefold hollow
c       of type 1 layers.
          izcnt=izcnt+1
          do 22 i=1,nchain
             ychain=MOD(i+1,2)*height/2.
             xlayer=2*znnd/3.
             do 32 j=1,nyatom
                  id=id+1
                  rv(1,id)=(i-1)*znnd+xlayer
                  rv(2,id)=(j-1)*height+ychain
                  rv(3,id)=-(izcnt-1)*dz
                  rv(4,id)=0.0
                  rv(5,id)=0.0
                  rv(6,id)=0.0
                  if(l1.lt.1) then
                    if(MOD(i,2).eq.1) then
                      if(MOD(i,4).eq.1) then
                        if(2-MOD(j,2).eq.1) then
                           itype(id)=t1
                        else
                           itype(id)=t2
                        endif
                      else if(MOD(i,4).eq.3) then
                        if(2-MOD(j+1,2).eq.1) then
                           itype(id)=t1
                        else
                           itype(id)=t2
                        endif
                      endif
                    else
                      itype(id)=t1
                    endif
                  else
                    if(l1.lt.2) then
                      itype(id)=t2
                    else
                      itype(id)=t1
                    endif
                  endif
32            continue
22        continue
        endif
10      continue
c       count the number of atom
        natoms=id
      return
      end

      
      subroutine printtyp(stacks)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      integer natoms
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      integer stacks(3),nchain,nyatom,nlayer
      integer i,j,k,id,ii
      nchain=stacks(1)
      nyatom=stacks(2)
      nlayer=stacks(3)
      id=1
      do 10 k= 1,nlayer
        write(*,80)'-------Layer ',k,' -------'
        do 20 j = 1,nchain
            if(nyatom.eq.1) then
               write(*,91)itype(id)
            else if(nyatom.eq.2) then
               write(*,92)(itype(id+ii),ii=0,nyatom-1)
            else if(nyatom.eq.3) then
               write(*,93)(itype(id+ii),ii=0,nyatom-1)
            else if(nyatom.eq.4) then
               write(*,94)(itype(id+ii),ii=0,nyatom-1)
            else if(nyatom.eq.5) then
               write(*,95)(itype(id+ii),ii=0,nyatom-1)
            else if(nyatom.eq.6) then
               write(*,96)(itype(id+ii),ii=0,nyatom-1)
            else if(nyatom.eq.7) then
               write(*,97)(itype(id+ii),ii=0,nyatom-1)
            else if(nyatom.eq.8) then
               write(*,98)(itype(id+ii),ii=0,nyatom-1)
            else if(nyatom.eq.9) then
               write(*,99)(itype(id+ii),ii=0,nyatom-1)
            else
               write(*,90)(itype(id+ii),ii=0,nyatom-1)
            endif
            id=id+nyatom
20    continue
10    continue
91    format(1(i1,1x))
92    format(2(i1,1x))
93    format(3(i1,1x))
94    format(4(i1,1x))
95    format(5(i1,1x))
96    format(6(i1,1x))
97    format(7(i1,1x))
98    format(8(i1,1x))
99    format(9(i1,1x))
90    format(10(i1,1x))
80    format(a13,i2,a8)
      return
      end


      subroutine writevmd(filename)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      integer natoms
      character*80 filename
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      integer stacks(3),nchain,nyatom,nlayer
      integer i,j,k,id,ii
      filename=trim(filename)
      open(unit=55,file=filename)
      write(*,*)'Writing atom positions in LAMMPS dump format to ',
     .filename,'file.'
      write(55,7770)'ITEM: TIMESTEP'
      write(55,7771)'1'
      write(55,7772)'ITEM: NUMBER OF ATOMS'
      write(55,*)natoms
      write(55,7774)'ITEM: BOX BOUNDS'
      write(55,7770)'-10.000 10.000'
      write(55,7770)'-10.000 10.000'
      write(55,7770)'-10.000 10.000'
      write(55,7776)'ITEM: ATOMS'
      id=0
      do 10 i= 1,natoms
            id=id+1
            if(id.lt.10) then
              write(55,8881)id,itype(id),(rv(ii,id),ii=1,3)
            else if(id.lt.100) then
              write(55,8882)id,itype(id),(rv(ii,id),ii=1,3)
            else if(id.lt.1000) then
              write(55,8883)id,itype(id),(rv(ii,id),ii=1,3)
            else if(id.lt.10000) then
              write(55,8884)id,itype(id),(rv(ii,id),ii=1,3)
            else if(id.lt.100000) then
              write(55,8885)id,itype(id),(rv(ii,id),ii=1,3)
            else
              write(55,8886)id,itype(id),(rv(ii,id),ii=1,3)
            endif
10    continue
8881    format(i1,1x,i1,3g15.5)
8882    format(i2,1x,i1,3g15.5)
8883    format(i3,1x,i1,3g15.5)
8884    format(i4,1x,i1,3g15.5)
8885    format(i5,1x,i1,3g15.5)
8886    format(i6,1x,i1,3g15.5)
8887    format(a10,a22,a6)
7770    format(a14)
7771    format(a1)
7772    format(a21)
7773    format(a5)
7774    format(a16)
7776    format(a11)
      return
      end

      subroutine writedyn86(filename)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      integer natoms
      character*80 filename
      common /particle/ rv(6,nmax)
      common /unitcell/ ucell(nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      data conmas/1.0365e-4/
      integer stacks(3),nchain,nyatom,nlayer
      integer i,j,k,id,ii
      filename=trim(filename)
      open(unit=55,file=filename)
      write(*,*)'Writing DYN86 PHONON format to ',
     .filename,'file.'
      write(55,*)'DYN86 PHONON format to file.'
      write(55,9502)natoms,ntypes
9502  format(2i10,e15.8)
      write(55,9503) (perub(i),i=1,3),(perlb(i),i=1,3)
9503  format(3e25.16)
      write(55,9504) (amass(i)*conmas,ielement(i),i=1,ntypes)
9504  format(e25.16,i10)
      write(55,9505) ((rv(i,j),i=1,6),itype(j),
     .                IDINT(ucell(j)),j=1,natoms)
9505  format(3e25.16/3e25.16/2i10)
      return
      end

      subroutine removeatoms(atomnos)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      integer natoms
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      integer i,j,zplanes,atomnos(2),ii
      ii=0
      do 10 i=1,natoms
       if(i.lt.atomnos(2) .or.
     .    i.gt.atomnos(1) ) then
          ii=ii+1
          rv(1,ii) = rv(1,i)
          rv(2,ii) = rv(2,i)
          rv(3,ii) = rv(3,i)
       endif
10    continue
       natoms=ii
      return
      end

      subroutine shiftatoms(atomnos,shiftys)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      integer natoms
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      integer i,atomnos(2)
      double precision shiftys(3)
      do 10 i=1,natoms
       if(i.le.atomnos(2) .and.
     .    i.ge.atomnos(1) ) then
          rv(1,i) = rv(1,i) + shiftys(1)
          rv(2,i) = rv(2,i) + shiftys(2)
          rv(3,i) = rv(3,i) + shiftys(3)
       endif
10    continue
      return
      end

      subroutine shiftplanes(stacks,zplanes,shifts)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      integer natoms
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      integer stacks(3),nchain,nyatom,nlayer
      integer i,zplanes
      double precision newx,newy,newz,shifts(3)
      nchain=stacks(1)
      nyatom=stacks(2)
      nlayer=stacks(3)
      do 10 i= 1,natoms
       if(i.lt.(nchain*nyatom*zplanes+1)) then
        rv(1,i) = rv(1,i) + shifts(1)
        rv(2,i) = rv(2,i) + shifts(2)
        rv(3,i) = rv(3,i) + shifts(3)
       endif
10    continue
      call putbacktobox()
      return
      end


      subroutine shiftplanest(stacks,zplanes,shifts)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      integer natoms
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      integer stacks(3),nchain,nyatom,nlayer
      integer i,zplanes(2)
      double precision newx,newy,newz,shifts(3)
      nchain=stacks(1)
      nyatom=stacks(2)
      nlayer=stacks(3)
      do 10 i= 1,natoms
       if(i.lt.(nchain*nyatom*zplanes(2)+1) .and. 
     .    i.gt.(nchain*nyatom*zplanes(1))) then
        rv(1,i) = rv(1,i) + shifts(1)
        rv(2,i) = rv(2,i) + shifts(2)
        rv(3,i) = rv(3,i) + shifts(3)
       endif
10    continue
      call putbacktobox()
      return
      end



      subroutine shiftplanesx(stacks,xplanes,shifts)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      integer natoms
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      integer stacks(3),nchain,nyatom,nlayer
      integer i,j,k,id,xplanes
      double precision newx,newy,newz,shifts(3)
      nchain=stacks(1)
      nyatom=stacks(2)
      nlayer=stacks(3)
      id=0
      do 10 i= 1,nlayer
       do 11 j= 1,nchain
        do 12 k= 1,nyatom
         id=id+1
         if(j.le.xplanes) then
            rv(1,id) = rv(1,id) + shifts(1)
            rv(2,id) = rv(2,id) + shifts(2)
            rv(3,id) = rv(3,id) + shifts(3)
         endif
12      continue
11     continue
10    continue
      call putbacktobox()
      return
      end


      subroutine shiftplanesy(stacks,yplanes,shifts)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      integer natoms
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      common /types/ amass(nelmax),ielement(nelmax),netype(nelmax),
     .               ntypes,itype(nmax)
      integer stacks(3),nchain,nyatom,nlayer
      integer i,j,k,id,yplanes
      double precision newx,newy,newz,shifts(3)
      nchain=stacks(1)
      nyatom=stacks(2)
      nlayer=stacks(3)
      id=0
      do 10 i= 1,nlayer
       do 11 j= 1,nchain
        do 12 k= 1,nyatom
         id=id+1
         if(k.le.yplanes) then
            rv(1,id) = rv(1,id) + shifts(1)
            rv(2,id) = rv(2,id) + shifts(2)
            rv(3,id) = rv(3,id) + shifts(3)
         endif
12      continue
11     continue
10    continue
      call putbacktobox()
      return
      end


      subroutine putbacktobox()
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      integer natoms
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      integer i,j
c  this subroutine puts all the particles back into the basic
c  periodic box or cell.
      do 1100 j = 1,3
      do 1100 i = 1,natoms
      if (rv(j,i).gt.perub(j)) rv(j,i) = rv(j,i) - perlen(j)
      if (rv(j,i).lt.perlb(j)) rv(j,i) = rv(j,i) + perlen(j)
1100  continue
      return
      end

      subroutine calcper()
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      integer i
c
c  compute perlen
c
      do 100 i = 1,3
100   perlen(i) = perub(i) - perlb(i)
c
      return
      end

      subroutine checkper(rcutoff,istop,rms,prntflag)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision rcutoff,rms
      integer prntflag
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      double precision permin,rcutoffsq
      integer i,j,k,istop,ii
      logical onlyoi,twoi,threei,fouri,notinlist

      onlyoi=.true.
      twoi=.false.
      threei=.false.
      fouri=.false.
c
c  compute perlen
c
      do 100 i = 1,3
100   perlen(i) = perub(i) - perlb(i)
c
      rcutoffsq=sqrt(rcutoff)
      permin = 2.*rcutoffsq
      istop = 0
      do 200 ii=1,2
      if(perlen(ii).ge.permin)go to 200
      if(prntflag.eq.1) then
         write(6,9230)ii
         write(6,9240)permin,perlen(ii)
      endif
 9230 format('   periodicity is too short in the ',i2,'  direction')
 9240 format('   permin = ',e15.5,'  periodicity = ',e15.5)
      rms=110000.0+rcutoff
      istop = 1
 200   continue
      if(istop.eq.1) return
c
c       z direction:
c       four levels
c
      if(perlen(3).lt.2.*rcutoffsq)then
          onlyoi=.false.
          twoi=.true.
      endif
      if(onlyoi)go to 55 
      if(perlen(3).lt.rcutoffsq)threei=.true.
      if(perlen(3).lt.2.*rcutoffsq/3.)fouri=.true.
      if(perlen(3).lt.0.5*rcutoffsq)then
         ii=3
         permin = 0.5*rcutoffsq
         if(prntflag.eq.1) then
           write(6,9230)ii
           write(6,9240)permin,perlen(ii)
         endif
         rms=120000.0+rcutoff
         istop=1
         return
c        stop
      endif
 55   continue
      return
      end


      subroutine loadneigh(rcutoff,neigh,rdis,rdisd,istop,rms,prntflag)
      parameter (nmax=4000, neimax=3000, nelmax=3, 
     . ngrid=20000, ngridar=20000)
      double precision rcutoff,rdis(nmax,nmax,2),rdisd(nmax,nmax,3),rms
      integer neigh(nmax),prntflag
      double precision dis(3,nmax),r(nmax),rneigh(nmax),
     . dneigh(3,nmax)
      integer nneighs(nmax),jneighs(nmax*nmax/2),neiind(0:nmax)
      common /particle/ rv(6,nmax)
      common /lattice/ perub(3),perlb(3),perlen(3),alat,
     . xbound(2),ybound(2),zbound(2),natoms,latty
      double precision permin,rcutoffsq,disz,rim
      double precision dneighs(3,nmax*nmax/2)
      integer i,j,k,nneigh,nc,jneigh(nmax)
      integer istop,kc,kcoord,neitot,nneips
      integer iat,jat,jj,ii
      logical onlyoi,twoi,threei,fouri,notinlist

      call checkper(rcutoff,istop,rms,prntflag)
      if(istop.eq.1) then
        if (prntflag.eq.1) then
         write(*,*)'rms1:',rms
        endif
         return
      endif
     
      if(rcutoff.gt.0.0) then
      rcutoffsq=rcutoff**2
      neitot = 0
      neiind(0) = neitot
      do 1000 i= 1,natoms
c
      do 1100 j = 1,i-1
c
c     compute the square of the distance to the closest periodic image
c
      dis(1,j) = rv(1,i) - rv(1,j)
      dis(1,j) = dis(1,j) - perlen(1)*nint(dis(1,j)/perlen(1))
      dis(2,j) = rv(2,i) - rv(2,j)
      dis(2,j) = dis(2,j) - perlen(2)*nint(dis(2,j)/perlen(2))
      dis(3,j) = rv(3,i) - rv(3,j)
      dis(3,j) = dis(3,j) - perlen(3)*nint(dis(3,j)/perlen(3))
      r(j) = dis(1,j)**2 + dis(2,j)**2 + dis(3,j)**2
1100  continue
    
c     ----------------------------------      
c
c     determine which pairs are separated by less than rcut
c     and store the needed information about these pairs
c
      nneigh = 0
      do 1200 j = 1,i-1
c
c     if nearest periodic image is out of range, then all images
c     will be
c
      if (r(j).gt.rcutoffsq) go to 1200
      nneigh = nneigh + 1
      rneigh(nneigh) = r(j)
      jneigh(nneigh) = j
      do 1250 kcoord = 1,3
1250  dneigh(kcoord,nneigh) = dis(kcoord,j)
c     check periodic images in z direction
      if(onlyoi)go to 1200
      disz = -sign(perlen(3),dis(3,j))
      rim = r(j) + 2.*dis(3,j)*disz + disz**2
c     if next nearest image is out of range, subsequent ones will also be
      if (rim.gt.rcutoffsq) go to 1200
      nneigh = nneigh + 1
      rneigh(nneigh) = rim
      jneigh(nneigh) = j
      dneigh(1,nneigh) = dis(1,j)
      dneigh(2,nneigh) = dis(2,j)
      dneigh(3,nneigh) = dis(3,j) + disz
      if(.not.threei)go to 1200
      disz = -disz
      rim = r(j) + 2.*dis(3,j)*disz + disz**2
      if (rim.gt.rcutoffsq) go to 1200
      nneigh = nneigh + 1
      rneigh(nneigh) = rim
      jneigh(nneigh) = j
      dneigh(1,nneigh) = dis(1,j)
      dneigh(2,nneigh) = dis(2,j)
      dneigh(3,nneigh) = dis(3,j) + disz
      if(.not.fouri)go to 1200
      disz = -2.*sign(perlen(3),dis(3,j))
      rim = r(j) + 2.*dis(3,j)*disz + disz**2
      if (rim.gt.rcutoffsq) go to 1200
      nneigh = nneigh + 1
      rneigh(nneigh) = rim
      jneigh(nneigh) = j
      dneigh(1,nneigh) = dis(1,j)
      dneigh(2,nneigh) = dis(2,j)
      dneigh(3,nneigh) = dis(3,j) + disz
1200  continue
c
c       now do diagonal (i=j) term for three or four images of self
c       both cases produce two images
c
      nneips = nneigh
      if(threei)then
c       first image of self
         nneips = nneips + 1
         disz = perlen(3)
         rim = disz**2
c     don't need to check against range here
         rneigh(nneips) = rim
         jneigh(nneips) = i
         dneigh(1,nneips) = 0.0
         dneigh(2,nneips) = 0.0
         dneigh(3,nneips) = disz
c       second image of self
         nneips = nneips + 1
         disz = -disz
         rneigh(nneips) = rim
         jneigh(nneips) = i
         dneigh(1,nneips) = 0.0
         dneigh(2,nneips) = 0.0
         dneigh(3,nneips) = disz
      endif
c
c     end of inline gneigh(i)
c
      nneighs(i) = nneigh
      do 1300 j=1,nneips
      neitot = neitot + 1
      jneighs(neitot) = jneigh(j)
      dneighs(1,neitot) = dneigh(1,j)
      dneighs(2,neitot) = dneigh(2,j)
      dneighs(3,neitot) = dneigh(3,j)
1300  continue
      neiind(i) = neitot
1000  continue
c     ----------------------------------      
      else
        write(6,9013)rcutoff
9013    format(' Cut-off radius',e25.16,' is not bigger than 0')
        stop
      endif
     
      do 1400 iat=1,natoms
c
c  obtain the information about the neighbors of the given atom
c
      nneips = neiind(iat) - neiind(iat-1)

       nneigh=nneighs(iat)
      do 1500 j=1,nneips
      jj = neiind(iat-1) + j
      jneigh(j) = jneighs(jj)
      dneigh(1,j) = dneighs(1,jj)
      dneigh(2,j) = dneighs(2,jj)
      dneigh(3,j) = dneighs(3,jj)
1500     continue
      nc=nneips
      do 3000 i=iat+1,natoms
c
c  obtain the information about the neighbors of the given atom
c
      nneips = neiind(i) - neiind(i-1)

      do 3100 j=1,nneips
      jj = neiind(i-1) + j
      if(jneighs(jj).eq.iat)then
         notinlist=.true.
         do 3200 k=1,nc
           if(jneigh(k).eq.i)then
              notinlist=.false.
           endif
3200     continue
         if(notinlist)then
         nc=nc+1
         jneigh(nc) = i
         dneigh(1,nc) = -dneighs(1,jj)
         dneigh(2,nc) = -dneighs(2,jj)
         dneigh(3,nc) = -dneighs(3,jj)
          endif
        endif
3100     continue
3000     continue
c************************************************************************

c   All neighbors set       
       neigh(iat)=nc
       do 1600 j=1,nc
       jat = jneigh(j)
       rdis(iat,j,1)=sqrt(dneigh(1,j)**2+dneigh(2,j)**2+dneigh(3,j)**2)
       rdis(iat,j,2)=jat
       rdisd(iat,j,1)=dneigh(1,j)
       rdisd(iat,j,2)=dneigh(2,j)
       rdisd(iat,j,3)=dneigh(3,j)
1600   continue
1400   continue
       return
       end


