c     Predicting phi/psi
c	Eshel Faraggi -- Dec. 2007
c	21 window 
c	g77 -ffixed-line-length-none -Wall -W -Wsurprising  -O2 -o phipsi_phipsi0.e phipsi_phipsi0.f
c	ifort -132 -watch all -warn all -O2 -o phipsi_phipsi0.e phipsi_phipsi0.f
c	ifort -132 -O2 -o phipsi_phipsi0.e phipsi_phipsi0.f
c	ifort -132 -g -o phipsi_phipsi0.e phipsi_phipsi0.f
      program TP_r
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=31,nw1=101,nw2=101,nrxo=1,numav=5)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2),wo2(nw2+1,nrxo,2)
      common/filin/ protin(nrdx,nip+3), pk1, pk2
      common/filout/ phia(nrdx),icphi(nrdx),psia(nrdx), phiav(nrdx,numav)
      common/inout/ resv(nrx+1,nip)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(nrxo,2)
      common/edges/ ibegia(nrdx),iendia(nrdx),ibegoa(nrdx),iendoa(nrdx), nisf, nosf
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2)
      common/peaks/ peks(nrdx)
      common/misc/ alpha, itsd
      character seqres(nrdx)
      common seqres
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf
      character*150 spec, spxdir, workdir, listfl, flout
      call getarg(1,spec)
      call getarg(2,spxdir)
      call getarg(3,workdir)
      call getarg(4,listfl)
c      read(spec,'(I20)') inum
c* Definition of the general input/output files
      call getdefs(listfl)
      call netinit()
      call setio(spec,workdir)	! initialize input/output
CCCCC  PSI:
      pk1 = -0.5d0
      pk2 = 0.5d0
C      "-------------------PEAKS-------------------"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do 287 iav=1,numav
        call initphi(phia,icphi)
        call readwei("psi",0,iav,spxdir)
        do ix=1,nta
            call readwin(ix,0,0)
            call calcnet(ix,0,2.92d0)
        enddo
        call normphi(phia,phiav,iav)
287   continue
      call getpks()
C      "-------------------ANGLES-------------------"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do iav=1,numav
        call initphi(phia,icphi)
        call readwei("psi",10,iav,spxdir)
        do ix=1,nta
            call readwin(ix,0,1)
            call calcnet(ix,1,2.92d0)
        enddo
        call normphiz(iav)
      enddo
      call getangs()
      do ix=1,nta
          psia(ix) = phia(ix)
      enddo
CCCCC  PHI:
      pk1 = -0.34722d0
      pk2 = 0.31944d0
C      "-------------------PEAKS-------------------"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do iav=1,numav
        call initphi(phia,icphi)
        call readwei("phi",0,iav,spxdir)
        do ix=1,nta
            call readwin(ix,0,0)
            call calcnet(ix,0,2.92d0)
        enddo
        call normphi(phia,phiav,iav)
      enddo
      call getpks()
C      "-------------------ANGLES-------------------"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do iav=1,numav
        call initphi(phia,icphi)
        call readwei("phi",10,iav,spxdir)
        do ix=1,nta
            call readwin(ix,0,1)
            call calcnet(ix,1,2.92d0)
        enddo
        call normphiz(iav)
      enddo
      call getangs()
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC--- OUTPUT ANGLE FILES  ---CCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      flout = trim(workdir)//"/out_phipsi0"
      open(unit=99,file=flout,status='unknown')
      do i=1,nta
        tphia = 180.d0 * phia(i)
        if (psia(i).lt.0.72222) then
          tpsia = psia(i) + 0.27778
        else
          tpsia = psia(i) - 1.72222
        endif
        tpsia = 180.d0 * tpsia
        write(99,499) i,seqres(i),tphia,tpsia,phia(i),psia(i)
      enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
499   format(I9,3x,A1,1x,F7.1,1x,F7.1,1x,F9.5,1x,F9.5)
      stop
      end
c_______________________________________________________________________
      
c------------- subroutine getdefs-------------------------------
c--
c* Get initial definitions
c--
      subroutine getdefs(listfl)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=31,nw1=101,nw2=101,nrxo=1,numav=5)
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2)
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf
      common/misc/ alpha, itsd
      character listfl*150
c* Definitions
c*********** alpha -- parameter for transfer function
      alpha = 0.2d0
c*********** dfac -- parameter for distance decay matrix
      dfac = 2.d0
      open(unit=3,file=listfl,status='old')
c* read names of files containing training proteins
      i = 1
      do 245 while (i.gt.0)
        read(3,'(A)',end=246)rtrainf(i)
	i = i + 1
245   continue
246   continue
      close(3)
      ntf = i-1
      end
c___________________________________________________________________
c------------- subroutine netinit -------------------------------
c--
c* Initialize the neural network (the weights are randomized)
c--
      subroutine netinit()
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=31,nw1=101,nw2=101,nrxo=1,numav=5)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2),wo2(nw2+1,nrxo,2)
      common/inout/ resv(nrx+1,nip)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(nrxo,2)
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2)
      common/misc/ alpha, itsd
      
c***   setting up the biases
      resv(nrx+1,1) = 1.d0
      hnu2(nw2+1) = 1.d0
      hnu1(nw1+1) = 1.d0
        delhi = (nrx-1.d0)/(nw1-1.d0)
c***    initialize distance matrix
	do i=1,nrx
          do k=1,nw1
	    dij1(i,k) = dfac / dsqrt(1.d0+((k-1)*delhi-i+1)*
     *                          ((k-1)*delhi-i+1))
	  enddo
	enddo
        jc = (nw1+1)/2
        jc2 = (nw2+1)/2
        mc = (nrxo+1)/2
	do i=1,nw1
          do j=1,nw2
	    dijh(i,j) = dfac / dsqrt(1.d0+delhi*delhi*(i-jc-j+jc2)*
     *                          (i-jc-j+jc2))
          enddo
	enddo
	do i=1,nw2
          do j=1,nrxo
	    dij2(i,j) = dfac / dsqrt(1.d0+((i-jc)*delhi-j+mc)*
     *                          ((i-jc)*delhi-j+mc))
          enddo
	enddo
      end
c___________________________________________________________________
c------------- subroutine setio -------------------------------
c--
c* Setup the input/output arrays
c--
      subroutine setio(ftfl,workdir)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=31,nw1=101,nw2=101,nrxo=1,numav=5)
      common/filin/ protin(nrdx,nip+3), pk1, pk2
      common/edges/ ibegia(nrdx),iendia(nrdx),ibegoa(nrdx),iendoa(nrdx), nisf, nosf
      character seqres(nrdx)
      common seqres
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf
      
      character ftfl*150, inflnm*150, workdir*150, tmpchr*1
      double precision tprot(49)
c* midpoint for input/output windows
      nisf = (nrx+1)/2
      nosf = (nrxo+1)/2
      inflnm = trim(workdir)//"/out_ss0"
      open(unit=44,file=inflnm,status='old')
      inflnm = trim(workdir)//"/out_rsa0"
      open(unit=45,file=inflnm,status='old')
      nta = 0
        i = 1
        do 265 while (i.gt.0)
          nta = nta + 1
          read(44,*,end=266) ntmp, tmpchr, tmpchr, protin(nta,nip+1), protin(nta,nip+2), protin(nta,nip+3), (tprot(j),j=1,27),
     *      ibegia(nta),iendia(nta),ibegoa(nta),iendoa(nta)
          read(45,*,end=266) ntmp, tmpchr, atmp, protin(nta,nip)
          seqres(nta) = tmpchr
          do j=1,27
            protin(nta,j) = tprot(j)
          enddo
265     continue
266     continue
      close(44)
      close(45)
      nta = nta - 1
      end
c___________________________________________________________________
c------------- subroutine initphi -------------------------------
c--
c* Initialize output angle array
c--
      subroutine initphi(phit,icphit)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=31,nw1=101,nw2=101,nrxo=1,numav=5)
      common/filin/ protin(nrdx,nip+3), pk1, pk2
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf
      double precision phit(nrdx)
      integer icphit(nrdx)
      do i=1,nta
        phit(i) = 0.d0
        icphit(i) = 0
      enddo
      end
c___________________________________________________________________
c------------- subroutine readwei -------------------------------
c--
c* Read optimized weights
c--
      subroutine readwei(angid,ishift,iav2,spxdir)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=31,nw1=101,nw2=101,nrxo=1,numav=5)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2),wo2(nw2+1,nrxo,2)
      
      character angid*3,weifl*150,weishft*2,avo,weidir*150, spxdir*150
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCC--- INPUT WEIGHTS  ---CCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      itmp = ishift + 77
      write(weishft,'(I2)') itmp
      write(avo,'(I1)') iav2
      
      weidir = trim(spxdir)//"/weights/phipsi0/"
      weifl = trim(weidir)//angid//"/av"//avo//"/fort."//weishft
      open(unit=77,file=weifl,status='unknown')
      do i=1,nrx+1
	do j=1,nip
	  do k=1,nw1
	    read(77,*) wi(i,j,k)
	  enddo
	enddo
      enddo
      close(77)
      itmp = itmp + 1
      write(weishft,'(I2)') itmp
      weifl = trim(weidir)//angid//"/av"//avo//"/fort."//weishft
      open(unit=77,file=weifl,status='unknown')
      do i=1,nw1+1
        do j=1,nw2
	  read(77,*) wh(i,j)
        enddo
      enddo
      close(77)
      itmp = itmp + 1
      write(weishft,'(I2)') itmp
      weifl = trim(weidir)//angid//"/av"//avo//"/fort."//weishft
      open(unit=77,file=weifl,status='unknown')
      do i=1,nw2+1
        do j=1,nrxo
	  read(77,*) (wo2(i,j,k),k=1,2)
        enddo
      enddo
      close(77)
      end
c___________________________________________________________________
c------------- subroutine readwin -------------------------------
c--
c* Read the current protein window
c--
      subroutine readwin(iwin,itp,itp2)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=31,nw1=101,nw2=101,nrxo=1,numav=5)
      common/filin/ protin(nrdx,nip+3), pk1, pk2
      common/inout/ resv(nrx+1,nip)
      common/edges/ ibegia(nrdx),iendia(nrdx),ibegoa(nrdx),iendoa(nrdx), nisf, nosf
     
      integer iwin,itp,itp2
      do i=ibegia(iwin),iendia(iwin)
        irw = iwin-nisf+i
        resv(i,1) = protin(irw,nip+1)
        resv(i,2) = protin(irw,nip+2)
        resv(i,3) = protin(irw,nip+3)
        do j=1,20
          resv(i,j+3) = protin(irw,j)
        enddo
        do j=1,7
          resv(i,j+23) = protin(irw,j+20)
        enddo
        resv(i,31) = protin(irw,nip)
      enddo
      
      end
c___________________________________________________________________
c------------- subroutine calcnet -------------------------------
c--
c* Calculate the neural network's output angle given input residue
c* sequence parameters in resv matrix
c--
      subroutine calcnet(iwin,itp,aafac)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=31,nw1=101,nw2=101,nrxo=1,numav=5)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2),wo2(nw2+1,nrxo,2)
      common/filout/ phia(nrdx),icphi(nrdx),psia(nrdx), phiav(nrdx,numav)
      common/inout/ resv(nrx+1,nip)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(nrxo,2)
      common/edges/ ibegia(nrdx),iendia(nrdx),ibegoa(nrdx),iendoa(nrdx), nisf, nosf
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2)
      common/peaks/ peks(nrdx)
      integer iwin,itp
      double precision aafac
      double precision tfunc !transfer function
c** Calculate nueron values of the first hidden layer
      do i=1,nw1
        ztemp=0.d0
	do j=ibegia(iwin),iendia(iwin)
	  do k=1,nip
	    ztemp = ztemp + resv(j,k) * wi(j,k,i) * dij1(j,i)
	  enddo
	enddo
        ztemp = ztemp + resv(nrx+1,1) * wi(nrx+1,1,i)
	hnu1(i) = tfunc(ztemp)
      enddo
      do i=1,nw2
        ztemp=0.d0
	do j=1,nw1
          ztemp = ztemp + hnu1(j) * wh(j,i) * dijh(j,i)
	enddo
        ztemp = ztemp + hnu1(nw1+1) * wh(nw1+1,i)
	hnu2(i) = tfunc(ztemp)
      enddo
c** Calculate output nueron values
      do i=ibegoa(iwin),iendoa(iwin)
        ztemp=0.d0
        do j=1,nw2
          ztemp = ztemp + hnu2(j) * wo2(j,i,1) * dij2(j,i)
        enddo
        ztemp = ztemp + hnu2(nw2+1) * wo2(nw2+1,i,1)
        onu2(i,1) = tfunc(ztemp)
        ztemp=0.d0
        do j=1,nw2
          ztemp = ztemp + hnu2(j) * wo2(j,i,2) * dij2(j,i)
        enddo
        ztemp = ztemp + hnu2(nw2+1) * wo2(nw2+1,i,2)
        onu2(i,2) = tfunc(ztemp)
CCCCC! PEAK (itp=0)
        if(itp.le.0) then
          if(onu2(i,1).ge.onu2(i,2)) then
            i2phin = -1
          else
            i2phin = 1
          endif
          phia(iwin+i-nosf) = phia(iwin+i-nosf) + i2phin
CCCCC! dPHI (itp=1)
       else
          phia(iwin+i-nosf) = phia(iwin+i-nosf) + peks(iwin+i-nosf) + (onu2(i,1)+onu2(i,2))/aafac
          icphi(iwin+i-nosf) = icphi(iwin+i-nosf) + 1
        endif
      enddo
      end
c___________________________________________________________________
c------------- subroutine normphi -------------------------------
c--
c* Normalize peak assignement (averaging)
c--
      subroutine normphi(phit,phitv,iav2)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=31,nw1=101,nw2=101,nrxo=1,numav=5)
      common/filin/ protin(nrdx,nip+3), pk1, pk2
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf
      double precision phit(numav), phitv(nrdx,numav)
      do i=1,nta
        if (phit(i).lt.0) then
          phitv(i,iav2) = -1.d0
        else
          phitv(i,iav2) = 1.d0
        endif
      enddo
      end
c___________________________________________________________________
c------------- subroutine normphiz -------------------------------
c--
c* Normalize output angle (averaging)
c--
      subroutine normphiz(iav2)
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=31,nw1=101,nw2=101,nrxo=1,numav=5)
      common/filin/ protin(nrdx,nip+3), pk1, pk2
      common/filout/ phia(nrdx),icphi(nrdx),psia(nrdx), phiav(nrdx,numav)
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf
      do i=1,nta
        phiav(i,iav2) = phia(i) / icphi(i)
      enddo
      end
c___________________________________________________________________
c------------- subroutine getpks -------------------------------
c--
c* Get the peak prediction
c--
      subroutine getpks()
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=31,nw1=101,nw2=101,nrxo=1,numav=5)
      common/filin/ protin(nrdx,nip+3), pk1, pk2
      common/filout/ phia(nrdx),icphi(nrdx),psia(nrdx), phiav(nrdx,numav)
      common/peaks/ peks(nrdx)
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf
      
c** Calculate peaks
      do i=1,nta
        phia(i) = phiav(i,1)
        do j=2,numav
          phia(i) = phia(i) + phiav(i,j)
        enddo
        if (phia(i).le.0.d0) then
          peks(i) = pk1
        else
          peks(i) = pk2
        endif
      enddo
      end
c___________________________________________________________________
c------------- subroutine getangs -------------------------------
c--
c* Get the peak prediction
c--
      subroutine getangs()
      implicit real*8(a-h,o-z)
      PARAMETER(nrpx=8000,nrdx=2000000,nrx=21,nip=31,nw1=101,nw2=101,nrxo=1,numav=5)
      common/filin/ protin(nrdx,nip+3), pk1, pk2
      common/filout/ phia(nrdx),icphi(nrdx),psia(nrdx), phiav(nrdx,numav)
      common/peaks/ peks(nrdx)
      character rtrainf(nrpx)*150
      common/sinout/ ntf,nta,rtrainf
      
c** Calculate peaks
      do i=1,nta
        phia(i) = phiav(i,1)
        do j=2,numav
          phia(i) = phia(i) + phiav(i,j)
        enddo
        phia(i) = phia(i) / numav
      enddo
      end
c___________________________________________________________________
c------------- function tfunc(xoutp) -------------------------------
c--
c* This is the transfer function f(x) = 1/(1+exp(alpha*x))
c* Note: the derivative f'(x) = alpha*f(x)*(1-f(x)). This fact is used
c* in the computation
c--
      real*8 function tfunc(xoutp)
      implicit real*8(a-h,o-z)
      common/misc/ alpha, itsd
      tfunc = tanh(alpha * xoutp)
      return
      end function tfunc
c___________________________________________________________________
c------------- function randnm() -------------------------------
c--
c* LCG psuedorandom 
c--
      real*8 function randnm()
      implicit real*8(a-h,o-z)
      common/misc/ alpha, itsd
      
      dvs = 2147483647.d0
      dnm = 2147483711.d0
      dml = 16807.d0
      if(itsd .lt. 1) itsd = itsd + dvs
      ditsd = itsd * dml
      ditsd = mod(ditsd, dvs)
      randnm = itsd / dnm
      itsd = ditsd
      return
      end function randnm
c___________________________________________________________________
