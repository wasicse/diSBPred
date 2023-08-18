c     General Artificial Neural Network
c	g77 -ffixed-line-length-none -Wall -W -Wsurprising  -O2 -o genn.e genn.f
c	ifort -132 -O2 -o genn.e genn.f
c       EF 2009
c
      program GENN
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1900,nrpx=2640,naax=60000,nrx=15,nip=34,nw1=51,nw2=51,nrxo=1)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2)
      common/oldwei/ oldwi(nrx+1,nip,nw1),oldwh(nw1+1,nw2)
      common/pwei/ wo2(nw2+1,nrxo,2),oldwo2(nw2+1,nrxo,2)
      common/filin/ prt1(naax,100),prt2(naax,101:200),prt3(naax,201:300)
      common/filin2/ prt4(naax,301:400),prt5(naax,401:500),prt6(naax,501:600),phir(naax,nrxo)
      common/filout/ phia(naax,nrxo)
      common/inout/ resv(nrx+1,nip),phi2(nrxo,2)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(nrxo,2)
      common/learn/ alpha, ratel, rmoment
      common/edges/ ibegia(naax),iendia(naax),ibegoa(naax),iendoa(naax), nisf, nosf
      common/misc/ rmae(nrxo), ccoef(nrxo), cvfrac, bsint, iqscale, itsd
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2)
      character rtrainf(nrpx)*120
      common/sinout/ ntotf,ntf,numtra,nta,idxprt(nrpx),rtrainf
      common/ipin/ ipind(naax,2)
      character fylevld*120
c* Definition of the general input/output files
      call getdefs(fylevld,ibinv,nte,npe,icv)
      call netinit()    ! initialize weights
      call setio()	! initialize input/output
      if (iqscale.eq.1) call scalin()	! scale input parameters
      call readwei(10)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC--- OUTPUT CROSS VALIDATION  ---CCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        do ix=1,nta
          call readwin(ix)
          call calcnet(ix)
        enddo
        call calcor(phir,phia,1,nta)
        call calcerr(phir,phia,1,nta,tserr,idupl)
        write(*,*) "MAE: ", (rmae(io),io=1,nrxo)
        write(*,*) "CORC: ", (ccoef(io),io=1,nrxo)
        write(*,*) "ID_max,duplications: ", tserr,idupl
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC--- OUTPUT prediction FILES  ---CCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      fylevld = "pred_test.out"
      open(unit=97,file=fylevld)
      do i=1,nta
        write(97,499) ipind(i,1),ipind(i,2),(phir(i,j),j=1,nrxo),(phia(i,j),j=1,nrxo)
      enddo
      close(97)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
499   format(I6,2x,I6,2x,100(D13.4))
      stop
      end
c_______________________________________________________________________
      
c------------- subroutine getdefs-------------------------------
c--
c* Get initial definitions
c--
      subroutine getdefs(fylevld,ibinv,nte,npe,icv)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1900,nrpx=2640,naax=60000,nrx=15,nip=34,nw1=51,nw2=51,nrxo=1)
      common/learn/ alpha, ratel, rmoment
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2)
      common/misc/ rmae(nrxo), ccoef(nrxo), cvfrac, bsint, iqscale, itsd
      character rtrainf(nrpx)*120
      common/sinout/ ntotf,ntf,numtra,nta,idxprt(nrpx),rtrainf
      character fyle*6, fylein*120, fylevld*120,tmp*120
      character rtrain*120, rtraing(nrpx)*120
c* Definition of the general input/output files
      fyle='spiner'
      fylein=fyle//'.in'
      fylevld=fyle//'.valdav'
      open(unit=2,file=fylein)
c****  ibinv is the number of epochs to skip before bin output
      ibinv = 50
c****  The following parameters are read from the input file progname.in
c*********** nte -- number of training epochs
      read(2,*)nte
c*********** npe -- size of protein set
      read(2,*)npe
c*********** ipe -- iterations per exemplar
      read(2,*)ipe
c*********** ratel -- learning rate
      read(2,*)ratel
c*********** rmoment -- initial momentum
      read(2,*)rmoment
c*********** alpha -- parameter for transfer function
      read(2,*)alpha
c*********** dfac -- parameter for distance decay matrix
      read(2,*)dfac
c*********** icv -- set for 10 fold cross validation
      read(2,*)icv
c*********** cvfrac -- fraction for cross validation
      read(2,*)cvfrac
c*********** iqscale -- Question: scale inputs? (0-no / 1-yea)
      read(2,*)iqscale
c*********** bsint -- boundary of symmetric interval to scale inputs into
      read(2,*)bsint
c*********** rtrain -- file containing names of training proteins
      read(2,*)rtrain,tmp
      close(2)
      open(unit=3,file=rtrain)
c* read names of files containing training proteins
      i = 1
      do 245 while (i.gt.0)
        read(3,'(A)',end=246)rtraing(i)
	i = i + 1
245   continue
246   continue
      close(3)
      ntotp = i-1
      if (npe.le.ntotp) then
        ntotf = npe
      else
        write(*,*) "Protein set too small for number of proteins selected, using full set."
        ntotf = ntotp
      endif
      ntf = ntotf
      do i=1,ntf
        rtrainf(i) = rtraing(i)
      enddo
      end
c___________________________________________________________________
c------------- subroutine netinit -------------------------------
c--
c* Initialize the neural network (the weights are randomized)
c--
      subroutine netinit()
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1900,nrpx=2640,naax=60000,nrx=15,nip=34,nw1=51,nw2=51,nrxo=1)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2)
      common/oldwei/ oldwi(nrx+1,nip,nw1),oldwh(nw1+1,nw2)
      common/pwei/ wo2(nw2+1,nrxo,2),oldwo2(nw2+1,nrxo,2)
      common/inout/ resv(nrx+1,nip),phi2(nrxo,2)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(nrxo,2)
      common/misc/ rmae(nrxo), ccoef(nrxo), cvfrac, bsint, iqscale, itsd
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2)
      
      integer*4 timeArray(3)    ! Holds the hour, minute, and second
      double precision randnm !random LCG
      call itime(timeArray)     ! Get the current time for seed
      ta1 = timeArray(1) / 1.74
      ta2 = timeArray(2) / 4.35
      ta3 = timeArray(3) / 2.8
      jta3 = timeArray(3)
      ita3 = exp(ta3)
      ta3 = (mod(ita3,61)/2.8)
      ta4 = exp(ta1) + exp(ta2) + exp(ta3)
      ita4log = log10(ta4)
      ita4 = ta4
      dta4 = 10**ita4log * (ta4 - ita4)
      itsd = abs(dta4 - ta4) + ta3
c      itsd = 44027
      write(*,*) "random seed = ",itsd
c***   setting up the biases
      resv(nrx+1,1) = 1.d0
      hnu2(nw2+1) = 1.d0
      hnu1(nw1+1) = 1.d0
c***    randomizing weights
	do i=1,nrx+1
	  do j=1,nip
	    do k=1,nw1
	      atmp = randnm() - 0.5
	      wi(i,j,k) = atmp
	      oldwi(i,j,k) = atmp / 10.d0
	    enddo
	  enddo
	enddo
	do i=1,nw1+1
          do j=1,nw2
	    atmp = randnm() - 0.5
	    wh(i,j) = atmp
	    oldwh(i,j) = atmp / 10.d0
          enddo
	enddo
	do i=1,nw2+1
          do j=1,nrxo
	    atmp = randnm() - 0.5
	    wo2(i,j,1) = atmp
	    oldwo2(i,j,1) = atmp / 10.d0
	    atmp = randnm() - 0.5
	    wo2(i,j,2) = atmp
	    oldwo2(i,j,2) = atmp / 10.d0
          enddo
	enddo
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
c------------- subroutine outwei -------------------------------
c--
c* Output weights
c--
      subroutine outwei(ishift)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1900,nrpx=2640,naax=60000,nrx=15,nip=34,nw1=51,nw2=51,nrxo=1)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2)
      common/pwei/ wo2(nw2+1,nrxo,2),oldwo2(nw2+1,nrxo,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCC--- OUTPUT WEIGHTS  ---CCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do i=1,nrx+1
	do j=1,nip
	  do k=1,nw1
	    write(ishift+37,*) wi(i,j,k)
	  enddo
	enddo
      enddo
      close(ishift+37)
      do i=1,nw1+1
        do j=1,nw2
	  write(ishift+38,*) wh(i,j)
        enddo
      enddo
      close(ishift+38)
      do i=1,nw2+1
        do j=1,nrxo
	  write(ishift+39,*) (wo2(i,j,k),k=1,2)
        enddo
      enddo
      close(ishift+39)
      end
c___________________________________________________________________
c------------- subroutine readwei -------------------------------
c--
c* Read optimized weights
c--
      subroutine readwei(ishift)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1900,nrpx=2640,naax=60000,nrx=15,nip=34,nw1=51,nw2=51,nrxo=1)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2)
      common/pwei/ wo2(nw2+1,nrxo,2),oldwo2(nw2+1,nrxo,2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCC--- INPUT WEIGHTS  ---CCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do i=1,nrx+1
	do j=1,nip
	  do k=1,nw1
	    read(ishift+37,*) wi(i,j,k)
	  enddo
	enddo
      enddo
      close(ishift+37)
      do i=1,nw1+1
        do j=1,nw2
	  read(ishift+38,*) wh(i,j)
        enddo
      enddo
      close(ishift+38)
      do i=1,nw2+1
        do j=1,nrxo
	  read(ishift+39,*) (wo2(i,j,k),k=1,2)
        enddo
      enddo
      close(ishift+39)
      end
c___________________________________________________________________
c------------- subroutine setio -------------------------------
c--
c* Setup the input/output arrays
c--
      subroutine setio()
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1900,nrpx=2640,naax=60000,nrx=15,nip=34,nw1=51,nw2=51,nrxo=1)
      common/filin/ prt1(naax,100),prt2(naax,101:200),prt3(naax,201:300)
      common/filin2/ prt4(naax,301:400),prt5(naax,401:500),prt6(naax,501:600),phir(naax,nrxo)
      common/edges/ ibegia(naax),iendia(naax),ibegoa(naax),iendoa(naax), nisf, nosf
      character rtrainf(nrpx)*120
      common/sinout/ ntotf,ntf,numtra,nta,idxprt(nrpx),rtrainf
      common/ipin/ ipind(naax,2)
      
      double precision tprot(700)
      integer iopos(nrxo)
      data iopos /35/
c* midpoint for input/output windows
      nisf = (nrx+1)/2
      nosf = (nrxo+1)/2
      nta = 0
      do ix=1,ntf
        ntao = nta
        open(unit=4,file=rtrainf(ix))
        i = 1
        do 265 while (i.gt.0)
          read(4,*,end=266) (tprot(j),j=1,36)
          nta = nta + 1
          ipind(nta,1) = idxprt(ix)
          ipind(nta,2) = i
          do j=1,nip
            if (j.le.100) then
              prt1(nta,j) = tprot(j)
            elseif (j.le.200) then
              prt2(nta,j) = tprot(j)
            elseif (j.le.300) then
              prt3(nta,j) = tprot(j)
            elseif (j.le.400) then
              prt4(nta,j) = tprot(j)
            elseif (j.le.500) then
              prt5(nta,j) = tprot(j)
            elseif (j.le.600) then
              prt6(nta,j) = tprot(j)
            else
              write(*,*) "Too many input features (limit: 600), aborting"
              stop
            endif
          enddo
          do j=1,nrxo
            phir(nta,j) = tprot(iopos(j))
          enddo
          i = i + 1
265     continue
266     continue
        close(4)
        numres = i-1  ! number of residues in current protein
c* Setup the edges for each amino acid
        do i=1,numres
          inta = ntao + i
          ib = 1 + nisf - i
          if(ib.lt.1) ib = 1
          ibegia(inta) = ib
          ib = numres + nisf - i
          if(ib.gt.nrx) ib = nrx
          iendia(inta) = ib
          ib = 1 + nosf - i
          if(ib.lt.1) ib = 1
          ibegoa(inta) = ib
          ib = numres + nosf - i
          if(ib.gt.nrxo) ib = nrxo
          iendoa(inta) = ib
        enddo
      enddo
      end
c___________________________________________________________________
c------------- subroutine scalin -------------------------------
c--
c* scale the input parameters
c--
      subroutine scalin()
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1900,nrpx=2640,naax=60000,nrx=15,nip=34,nw1=51,nw2=51,nrxo=1)
      common/filin/ prt1(naax,100),prt2(naax,101:200),prt3(naax,201:300)
      common/filin2/ prt4(naax,301:400),prt5(naax,401:500),prt6(naax,501:600),phir(naax,nrxo)
      common/misc/ rmae(nrxo), ccoef(nrxo), cvfrac, bsint, iqscale, itsd
      character rtrainf(nrpx)*120
      common/sinout/ ntotf,ntf,numtra,nta,idxprt(nrpx),rtrainf
      
      double precision tav(nip),tstd(nip)
      do j=1,nip
        tav(j) = 0.d0
        tstd(j) = 0.d0
      enddo
      do j=1,nip
        if (j.le.100) then
          do i=1,nta
            tav(j) = tav(j) + prt1(i,j)
          enddo
          tav(j) = tav(j) / nta
          do i=1,nta
            tstd(j) = tstd(j) + (prt1(i,j)-tav(j))*(prt1(i,j)-tav(j))
          enddo
        elseif (j.le.200) then
          do i=1,nta
            tav(j) = tav(j) + prt2(i,j)
          enddo
          tav(j) = tav(j) / nta
          do i=1,nta
            tstd(j) = tstd(j) + (prt2(i,j)-tav(j))*(prt2(i,j)-tav(j))
          enddo
        elseif (j.le.300) then
          do i=1,nta
            tav(j) = tav(j) + prt3(i,j)
          enddo
          tav(j) = tav(j) / nta
          do i=1,nta
            tstd(j) = tstd(j) + (prt3(i,j)-tav(j))*(prt3(i,j)-tav(j))
          enddo
        elseif (j.le.400) then
          do i=1,nta
            tav(j) = tav(j) + prt4(i,j)
          enddo
          tav(j) = tav(j) / nta
          do i=1,nta
            tstd(j) = tstd(j) + (prt4(i,j)-tav(j))*(prt4(i,j)-tav(j))
          enddo
        elseif (j.le.500) then
          do i=1,nta
            tav(j) = tav(j) + prt5(i,j)
          enddo
          tav(j) = tav(j) / nta
          do i=1,nta
            tstd(j) = tstd(j) + (prt5(i,j)-tav(j))*(prt5(i,j)-tav(j))
          enddo
        else
          do i=1,nta
            tav(j) = tav(j) + prt6(i,j)
          enddo
          tav(j) = tav(j) / nta
          do i=1,nta
            tstd(j) = tstd(j) + (prt6(i,j)-tav(j))*(prt6(i,j)-tav(j))
          enddo
        endif
        tstd(j) = tstd(j) / nta
        tstd(j) = dsqrt(tstd(j))
      enddo
      write(*,*) "Average of parameters:"
      write(*,'(700(D13.4))') (tav(j),j=1,nip)
      write(*,*) "STD of parameters:"
      write(*,'(700(D13.4))') (tstd(j),j=1,nip)
      do i=1,nta
        do j=1,nip
          if (j.le.100) then
            prt1(i,j) = bsint * (prt1(i,j)-tav(j))/tstd(j)
          elseif (j.le.200) then
            prt2(i,j) = bsint * (prt2(i,j)-tav(j))/tstd(j)
          elseif (j.le.300) then
            prt3(i,j) = bsint * (prt3(i,j)-tav(j))/tstd(j)
          elseif (j.le.400) then
            prt4(i,j) = bsint * (prt4(i,j)-tav(j))/tstd(j)
          elseif (j.le.500) then
            prt5(i,j) = bsint * (prt5(i,j)-tav(j))/tstd(j)
          else
            prt6(i,j) = bsint * (prt6(i,j)-tav(j))/tstd(j)
          endif
        enddo
      enddo
      end
c___________________________________________________________________
c------------- subroutine readwin -------------------------------
c--
c* Read the current protein window
c--
      subroutine readwin(iwin)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1900,nrpx=2640,naax=60000,nrx=15,nip=34,nw1=51,nw2=51,nrxo=1)
      common/filin/ prt1(naax,100),prt2(naax,101:200),prt3(naax,201:300)
      common/filin2/ prt4(naax,301:400),prt5(naax,401:500),prt6(naax,501:600),phir(naax,nrxo)
      common/inout/ resv(nrx+1,nip),phi2(nrxo,2)
      common/edges/ ibegia(naax),iendia(naax),ibegoa(naax),iendoa(naax), nisf, nosf
      integer iwin
      do i=ibegia(iwin),iendia(iwin)
        irw = iwin-nisf+i
        do j=1,nip
          if (j.le.100) then
            resv(i,j) = prt1(irw,j)
          elseif (j.le.200) then
            resv(i,j) = prt2(irw,j)
          elseif (j.le.300) then
            resv(i,j) = prt3(irw,j)
          elseif (j.le.400) then
            resv(i,j) = prt4(irw,j)
          elseif (j.le.500) then
            resv(i,j) = prt5(irw,j)
          else
            resv(i,j) = prt6(irw,j)
          endif
        enddo
      enddo
      do i=1,nrxo
        phi2(i,1) = phir(iwin,i)
        phi2(i,2) = phir(iwin,i)
      enddo
      end
c___________________________________________________________________
c------------- subroutine calcnet -------------------------------
c--
c* Calculate the neural network's output angle given input residue
c* sequence parameters in resv matrix
c--
      subroutine calcnet(iwin)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1900,nrpx=2640,naax=60000,nrx=15,nip=34,nw1=51,nw2=51,nrxo=1)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2)
      common/pwei/ wo2(nw2+1,nrxo,2),oldwo2(nw2+1,nrxo,2)
      common/filout/ phia(naax,nrxo)
      common/inout/ resv(nrx+1,nip),phi2(nrxo,2)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(nrxo,2)
      common/edges/ ibegia(naax),iendia(naax),ibegoa(naax),iendoa(naax), nisf, nosf
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2)
      integer iwin
      
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
      do i=1,nrxo
        do i2=1,2
          ztemp=0.d0
          do j=1,nw2
            ztemp = ztemp + hnu2(j) * wo2(j,i,1) * dij2(j,i)
          enddo
          ztemp = ztemp + hnu2(nw2+1) * wo2(nw2+1,i,1)
          onu2(i,i2) = tfunc(ztemp)
        enddo
        phia(iwin,i) = (onu2(i,1)+onu2(i,2))/2.d0
      enddo
      end
c___________________________________________________________________
c------------- subroutine backprop -------------------------------
c--
c* Calculate the neural network's output error given training output
c* phi. Backpropagate this error to the hidden neurons
c* and update the weights.
c* Subroutine also caluculates the RMSE as rmse
c--
      subroutine backprop(iwin)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1900,nrpx=2640,naax=60000,nrx=15,nip=34,nw1=51,nw2=51,nrxo=1)
      common/weights/ wi(nrx+1,nip,nw1),wh(nw1+1,nw2)
      common/oldwei/ oldwi(nrx+1,nip,nw1),oldwh(nw1+1,nw2)
      common/pwei/ wo2(nw2+1,nrxo,2),oldwo2(nw2+1,nrxo,2)
      common/inout/ resv(nrx+1,nip),phi2(nrxo,2)
      common/neurons/ hnu1(nw1+1),hnu2(nw2+1),onu2(nrxo,2)
      common/learn/ alpha, ratel, rmoment
      common/edges/ ibegia(naax),iendia(naax),ibegoa(naax),iendoa(naax), nisf, nosf
      common/dcays/ dij1(nrx,nw1),dij2(nw2,nrxo),dfac,dijh(nw1,nw2)
      integer iwin
      
      real*8 errh1(nrx,nw1),birh1(nw1),errh2(nw1,nw2),birh2(nw2),
     *biro2(nrxo,2),erro2(nw2,nrxo,2)
      
      ibegi = ibegia(iwin)
      iendi = iendia(iwin)
      ibego = ibegoa(iwin)
      iendo = iendoa(iwin)
c** Calculate errors at the output layer and the RMSE
      do i=1,nrxo
        errotmp = phi2(i,1) - onu2(i,1)
        biro2(i,1) = errotmp*alpha*(1.d0-onu2(i,1)*onu2(i,1))
        errotmp = phi2(i,2) - onu2(i,2)
        biro2(i,2) = errotmp*alpha*(1.d0-onu2(i,2)*onu2(i,2))
        do j=1,nw2
          erro2(j,i,1) = biro2(i,1) * dij2(j,i)
          erro2(j,i,2) = biro2(i,2) * dij2(j,i)
        enddo
      enddo
     
c** Calculate the errors at the second hidden layer
      do k=1,nw2
        ztemp=0.d0
        do j=1,nrxo
	  ztemp = ztemp + wo2(k,j,1)*erro2(k,j,1) + 
     *                    wo2(k,j,2)*erro2(k,j,2) 
        enddo
        birh2(k) = ztemp * alpha  * (1 - hnu2(k) * hnu2(k))
        do j=1,nw1
	  errh2(j,k) = birh2(k) * dijh(j,k)
        enddo
      enddo
c** Calculate the errors at the first hidden layer
      do k=1,nw1
        ztemp=0.d0
        do j=1,nw2
	  ztemp = ztemp + wh(k,j)*errh2(k,j)
        enddo
        birh1(k) = ztemp * alpha  * (1 - hnu1(k) * hnu1(k))
        do j=ibegi,iendi
	  errh1(j,k) = birh1(k) * dij1(j,k)
        enddo
c** Update weights:
c** Update input/hidden layer weights (wi)
        do i=ibegi,iendi
          do j=1,nip
            ztemp = ratel*errh1(i,k)*resv(i,j) + rmoment*oldwi(i,j,k)
            wi(i,j,k) = wi(i,j,k) + ztemp
            oldwi(i,j,k) = ztemp
          enddo
        enddo
c** Update input bias weights (wi(nrx+1,1,1-nw1))
        ztemp = ratel*birh1(k) + rmoment*oldwi(nrx+1,1,k)
        wi(nrx+1,1,k) = wi(nrx+1,1,k) + ztemp
        oldwi(nrx+1,1,k) = ztemp
      enddo
c** Update hidden layer weights (wh)
      do k=1,nw2
        do j=1,nw1
          ztemp = ratel*errh2(j,k)*hnu1(j) + rmoment*oldwh(j,k)
          wh(j,k) = wh(j,k) + ztemp 
          oldwh(j,k) = ztemp
        enddo
c** Update hidden bias weights (wh(nw1+1,1-nw2))
        ztemp = ratel*birh2(k) + rmoment*oldwh(nw1+1,k)
        wh(nw1+1,k) = wh(nw1+1,k) + ztemp
        oldwh(nw1+1,k) = ztemp
c** Update hidden/output layer weights (wo)
        do i=1,nrxo
	  ztemp = ratel*erro2(k,i,1)*hnu2(k) + rmoment*oldwo2(k,i,1)
	  wo2(k,i,1) = wo2(k,i,1) + ztemp
	  oldwo2(k,i,1) = ztemp
	  ztemp = ratel*erro2(k,i,2)*hnu2(k) + rmoment*oldwo2(k,i,2)
	  wo2(k,i,2) = wo2(k,i,2) + ztemp
	  oldwo2(k,i,2) = ztemp
        enddo
      enddo
c** Update hidden/output layer bias weights (wo)
        k=nw2+1
        do i=1,nrxo
	  ztemp = ratel*biro2(i,1) + rmoment*oldwo2(k,i,1)
	  wo2(k,i,1) = wo2(k,i,1) + ztemp
	  oldwo2(k,i,1) = ztemp
	  ztemp = ratel*biro2(i,2) + rmoment*oldwo2(k,i,2)
	  wo2(k,i,2) = wo2(k,i,2) + ztemp
	  oldwo2(k,i,2) = ztemp
        enddo
      end
c___________________________________________________________________
c------------- subroutine calcerr -------------------------------
c--
c* Calculate the neural network's error
c--
      subroutine calcerr(phin,phip,istart,ifinish,cerr,idup)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1900,nrpx=2640,naax=60000,nrx=15,nip=34,nw1=51,nw2=51,nrxo=1)
      
      double precision phin(naax,nrxo),phip(naax,nrxo),atmpa(nrxo)
c** Calculate error
      cerr = 0.d0
      idup = 0
      do i=istart,ifinish
        do j=1,nrxo
          atmpa(j) = phin(i,j)
        enddo
        call argmax(atmpa,nrxo,iamax,rmax,ierr)
        idup = idup + ierr
        niamax = iamax
        do j=1,nrxo
          atmpa(j) = phip(i,j)
        enddo
        call argmax(atmpa,nrxo,iamax,rmax,ierr)
        idup = idup + ierr
        if (niamax.eq.iamax) then
          cerr = cerr + 1.d0
        endif
      enddo
      cerr = cerr / (ifinish - istart + 1)
      end
c___________________________________________________________________
c------------- subroutine calcor -------------------------------
c--
c* Calculate the neural network's correlation and Q10 errors
c--
      subroutine calcor(phin2,phip2,istart,ifinish)
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1900,nrpx=2640,naax=60000,nrx=15,nip=34,nw1=51,nw2=51,nrxo=1)
      common/misc/ rmae(nrxo), ccoef(nrxo), cvfrac, bsint, iqscale, itsd
      
      double precision phin2(naax,nrxo),phip2(naax,nrxo),phin(naax),phip(naax)
      double precision tphi(naax),tphip(naax)
c** Calculate error, correlation and Q10
      do io=1,nrxo
        rmae(io) = 0.d0
        avphi = 0.d0
        avphip = 0.d0
        do i=istart,ifinish
          phin(i) = phin2(i,io)
          phip(i) = phip2(i,io)
          erro1 = dabs(phin(i) - phip(i))
          rmae(io) = rmae(io) + erro1
          avphi = avphi + phin(i)
          avphip = avphip + phip(i)
        enddo
        rmae(io) = rmae(io) / (ifinish - istart + 1)
        avphi = avphi / (ifinish - istart + 1)
        avphip = avphip / (ifinish - istart + 1)
        do i=istart,ifinish
          tphi(i) = phin(i) - avphi
          tphip(i) = phip(i) - avphip
        enddo
        crhn = 0.d0
        cnor1hn = 0.d0
        cnor2hn = 0.d0
        do i=istart,ifinish
          crhn = crhn + tphi(i)*tphip(i)
          cnor1hn = cnor1hn + tphi(i)*tphi(i)
          cnor2hn = cnor2hn + tphip(i)*tphip(i)
        enddo
        ccoef(io) = crhn / (sqrt(cnor1hn)*sqrt(cnor2hn))
      enddo
      end
c___________________________________________________________________
c------------- function tfunc(xoutp) -------------------------------
c--
c* This is the transfer function. Its derivative f'(x) is used
c* in the computation
c--
      real*8 function tfunc(xoutp)
      implicit real*8(a-h,o-z)
      common/learn/ alpha, ratel, rmoment
      tfunc = tanh(alpha * xoutp)
      return
      end function tfunc
c___________________________________________________________________
c------------- subroutine argmax -------------------------------
c--
c* Return argument iax for maximal element of array arx
c-- arx - aray to work on; narx - length of array; iax - argmax; axmax - maximum, ier - index for error output
      subroutine argmax(arx,narx,iax,axmax,ier)
      implicit real*8(a-h,o-z)
      
      dimension arx(*)
      ier = 0
      iax = 1
      axmax = arx(1)
      if (narx.ge.2) then
        do i=2,narx
          if (arx(i).ge.axmax) then
            if (arx(i).gt.axmax) then
              iax = i
              axmax = arx(i)
              ier = 0
            else
              ier = ier + 1
            endif
          endif
        enddo
      endif
      end
c___________________________________________________________________
c------------- function randnm() -------------------------------
c--
c* LCG psuedorandom 
c--
      real*8 function randnm()
      implicit real*8(a-h,o-z)
      PARAMETER(nrdx=1900,nrpx=2640,naax=60000,nrx=15,nip=34,nw1=51,nw2=51,nrxo=1)
      common/misc/ rmae(nrxo), ccoef(nrxo), cvfrac, bsint, iqscale, itsd
      
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
