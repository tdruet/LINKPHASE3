! PROGRAM : LINKPHASE3
! Author  : Tom DRUET
! Copyright (C) 2014-2020

!This program is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!any later version.

!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.

!You should have received a copy of the GNU General Public License
!along with this program.  If not, see <http://www.gnu.org/licenses/>.

program linkphase3
implicit none
! version 23/02/2020


integer ::parent1,parent2,nmarq,nhap,i,j,k,l,pos(2),round,ori,num,useold,verif,maxinfo,maxk
integer ::all1,all2,sall1,sall2,dall1,dall2,pall1,pall2,io,maxid,n1,n2,ntemplates,autosome,parstart,parend,sex
integer*1, allocatable ::hap(:,:,:),typ(:,:),hapin(:),prephaseinfo(:,:)
real*8 ::val,pmax,lik,lik2,bestlik,conv
integer ::nround,ani,nani,k1,k2,print_val,freqprint,pere,mere,add1,add2,parent,round2,nround2,nseed,nroundextra,seed,map_option
integer, allocatable ::sire(:),dam(:),sexes(:),newid(:),oldid(:),byparent(:,:)
real*8 ::position,position2,ph1,ph2,epsi,gerr
real*8,allocatable ::posi(:,:)
logical*1, allocatable ::genotyped(:),haplotyped(:)
character*80 ::genofile,pedfile,oldfile,markfile,step
character*1000000 ::genoline1
integer ::ori0,ori1,ori2,mk1,mk2,n_switched,n_erased,n_phased
real*8, allocatable ::alpha(:,:),scaling(:),bjk(:,:,:),beta(:,:),gamma(:,:),L11(:),L21(:),L12(:),L22(:)
real*8 ::pjump(0:1),pi,prob_emission
real*8, allocatable ::emap(:,:),probrecs(:,:),numrecs(:,:)
real*8 ::edenom1,edenom2,econv1,econv2,totlik1,totlik2
integer ::noffspring,firsto,lasto,oldparent,mc1,parsex,offspring,p12,p11,neighbour,r12,r11,firstm,lastm,c1,c2,sextemplate
integer ::ninfo,nprogeny(4),maxoffspring,numoff,inphase,outphase,m1,m2,nheter,nhomoz,checkrec,switch,reading
integer ::expanded,limit1,ninfol,ninfor,inphasel,inphaser,outphasel,outphaser,lastknown,lastori,famsize,nclust,nclust2,clust
integer ::lasttophase,stopt,ntested,refoff,phase_option,mate,nprinted,gender,norigins1,norigins2,nemission
real*8    ::rrate,endpos,startpos,dist,rratel,rrater,ngam(0:2),nh1,nh2,entro1,entro2,nall11,nall12,nall21,nall22
real*8, allocatable ::gammas(:,:),expco(:)
real*8 ::b11,b12,b21,b22,PTOT,p_error,p_correct
real*8, allocatable ::par_error(:),npar(:),noff(:),off_error(:),entropy(:),nentropy(:)
integer, allocatable ::phase1(:),phase2(:),info_offspring(:),haplotypes(:,:),hap2(:,:)
integer, allocatable ::moutphase(:,:),minphase(:,:),info_marker(:,:),nrec(:),origins(:),nall(:,:),flanking(:,:,:),phased(:),tab(:,:)
integer, allocatable ::tophase(:)
logical ::prephased
integer ::mapiter,niter_map,sexmap


 call read_data
 print*,'Finished reading data'
 call homozygotes
 print*,'Finished homozygotes'
 call mendelian
 print*,'Finished mendelian'
 call halfsib_phasing
 if(phase_option>1)print*,'Finished phasing families'


 call write_data
 if(phase_option>2 .and. nprinted>0)then
   do gender=0,2
     call requiem(gender)
   enddo
 endif

contains

! **************************************************************
! ********* Reading data (observations) and parameters *********
! **************************************************************

subroutine read_data
implicit none
integer, allocatable ::ids(:),genoin(:)
character*20 ::paramline
character*1 ::ispar

phase_option=1;reading=0;autosome=1;niter_map=1;sexmap=0;gerr=1.0d-3
open(100,file='linkin.txt')
read(100,*,iostat=io)paramline
if(io/=0)stop 'Parameter file too short !'
if(paramline/='#PEDIGREE_FILE')stop 'Error in parameter file, expected :: #PEDIGREE_FILE'
read(100,*,iostat=io)pedfile
if(io/=0)stop 'Parameter file too short !'
read(100,*,iostat=io)paramline
if(io/=0)stop 'Parameter file too short !'
if(paramline/='#GENOTYPE_FILE')stop 'Error in parameter file, expected :: #GENOTYPE_FILE'
read(100,*,iostat=io)genofile
if(io/=0)stop 'Parameter file too short !'
read(100,*,iostat=io)paramline
if(io/=0)stop 'Parameter file too short !'
if(paramline/='#MARKER_FILE')stop 'Error in parameter file, expected :: #MARKER_FILE'
read(100,*,iostat=io)markfile
if(io/=0)stop 'Parameter file too short !'
read(100,*,iostat=io)paramline
if(io/=0)stop 'Parameter file too short !'
if(paramline/='#HALFSIB_PHASING')stop 'Error in parameter file, expected :: #HALFSIB_PHASING'
read(100,*,iostat=io)paramline
if(io/=0)stop 'Parameter file too short !'
if(paramline=='yes')phase_option=2
read(100,*,iostat=io)paramline
if(io/=0)stop 'Parameter file too short !'
if(paramline/='#HMM_PHASING')stop 'Error in parameter file, expected :: #HMM_PHASING'
read(100,*,iostat=io)paramline
if(io/=0)stop 'Parameter file too short !'
if(paramline=='yes' .and. phase_option==2)phase_option=3
if(paramline=='yes' .and. phase_option==1)stop 'Error in parameter file, HMM_PHASING must be preceded by HALFSIB_PHASING'
read(100,*,iostat=io)paramline
if(io/=0)stop 'Parameter file too short !'
if(paramline/='#N_TEMPLATES')stop 'Error in parameter file, expected :: #N_TEMPLATES'
map_option=0
read(100,*,iostat=io)ntemplates
if(io/=0)stop 'Parameter file too short !'
read(100,*,iostat=io)paramline
if(io/=0)stop 'Parameter file too short !'
if(paramline/='#CHECK_PREPHASING')stop 'Error in parameter file, expected :: #CHECK_PREPHASING'
read(100,*,iostat=io)paramline
if(io/=0)stop 'Parameter file too short !'
if(paramline=='yes')map_option=1
do 
 read(100,*,iostat=io)paramline
 if(io/=0)exit
 if(paramline=='#COLUMNS')reading=1
 if(paramline=='#SEXCHROM')autosome=0
 if(paramline=='#SEXMAP')sexmap=1
 if(paramline=='#ITERATIONS')then
  read(100,*,iostat=io)niter_map
  if(io/=0)stop 'Number of iterations expected after #ITERATIONS option!'
 endif
 if(paramline=='#GENO_ERROR')then
  read(100,*,iostat=io)gerr
  if(io/=0)stop 'Value of genotyping error expected after #GENO_ERROR option!'
 endif
enddo

open(50,file=markfile)

nmarq=0
do
read(50,*,iostat=io)num
if(io/=0)exit
nmarq=nmarq+1
enddo

print*,'Number of markers ::',nmarq
rewind(50)

if(nmarq>=10000000)then
 print*,'More than 10,000,000 markers!'
 print*,'The output format must be update.'
 print*,'The program will now stop.'
 stop
endif

allocate(posi(nmarq,2))
posi=0.0

num=0;parend=0;parstart=nmarq;ispar='A'
do
if(autosome==1 .and. sexmap==0)read(50,*,iostat=io)k,markfile,position
if(autosome==0 .and. sexmap==0)read(50,*,iostat=io)k,markfile,position,ispar
if(autosome==1 .and. sexmap==1)read(50,*,iostat=io)k,markfile,position,position2
if(autosome==0 .and. sexmap==1)read(50,*,iostat=io)k,markfile,position,position2,ispar
if(io/=0)exit
num=num+1
if(num<=nmarq)posi(num,1)=position/100.0
if(sexmap==0)posi(num,2)=posi(num,1)
if(sexmap==1)posi(num,2)=position2/100.0
if(autosome==1 .or. ispar .ne. 'X')then
 parstart=min(parstart,num)
 parend=max(parend,num)
endif
if(num==1)cycle
if(posi(num,1) < posi(num-1,1) .or. posi(num,2) < posi(num-1,2))then
  print*,'Problem in the map file: marker positions are not consecutive !!!'
  stop
else if(posi(num,1)==posi(num-1,1) .or. posi(num,2)==posi(num-1,2))then
  print*,'Warning: some markers with identical positions ::',num,(num-1)
endif
enddo
print*,parstart,parend

if(sexmap==0)then
 print*,'Length of marker map (in cM or Mb) ::',(posi(nmarq,1)-posi(1,1))*100.0
else
 print*,'Length of first marker map (in cM or Mb) ::',(posi(nmarq,1)-posi(1,1))*100.0
 print*,'Length of second marker map (in cM or Mb) ::',(posi(nmarq,2)-posi(1,2))*100.0
endif

verif=0

! **************************************
! ***** Reading data (observations) ****
! **************************************
nhap=0;maxid=0;nani=0

open(9,file=pedfile)

do
if(autosome==1)read(9,*,iostat=io)ani,pere,mere
if(autosome==0)read(9,*,iostat=io)ani,pere,mere,sex
if(io/=0)exit
 maxid=max(maxid,ani)
 maxid=max(maxid,pere)
 maxid=max(maxid,mere)
enddo

print*,'Maximum number of animal red in pedigree file ::',maxid

if(maxid>=10000000)then
 print*,'Maxid ID larger than 10,000,000!'
 print*,'The output format must be update.'
 print*,'The program will now stop.'
 stop
endif

rewind(9)

allocate(newid(0:maxid),oldid(0:maxid))
newid=0;oldid=0

print*,'Name of genotypes file ::',genofile
open(13,file=genofile)

select case(reading)
case(0) ! one individual per line
 do 
!read(13,'(i6)',iostat=io)num
 read(13,*,iostat=io)num
 if(io/=0)exit
  if(num>maxid)print*,'Error: genotyped animal not in pedigree file ::',num
  if(newid(num)==0)then
    nani=nani+1
    newid(num)=nani
  endif
 enddo
 rewind(13)
case(1) ! one individual per column
 read(13,'(a100000000)',iostat=io)genoline1
 genoline1=adjustl(genoline1)
 k=len(trim(genoline1))
 nani=1
 do i=2,k
  if(genoline1(i:i)==" " .and. genoline1((i-1):(i-1))/=" ")nani=nani+1
 enddo
 rewind(13)
 allocate(ids(nani),genoin(nani))
 read(13,*,iostat=io)(ids(i),i=1,nani)
 do i=1,nani
  num=ids(i)
  if(num>maxid)print*,'Error: genotyped animal not in pedigree file ::',num
  if(newid(num)==0)then
    newid(num)=i
  endif
 enddo
end select

print*,'Number of genotyped animals ::',nani

! old animals with smallest number
! number must not be function of order in input file
! order must be the same as in initial pedigree file

nani=0
do i=1,maxid
 if(newid(i)/=0)then ! individual has genotype or haplotype
   nani=nani+1
   newid(i)=nani
   oldid(nani)=i
 endif
enddo

print*,'Confirmation - Number of genotyped or phased animals ::',nani

nhap=2*nani

allocate(hap(nani,2,nmarq),typ(nani,2*nmarq),genotyped(nani))
allocate(sire(nani),dam(nani),haplotyped(nani),hapin(nmarq),byparent(2*nani,2),prephaseinfo(nani,nmarq))
allocate(sexes(nani))

sire=0;dam=0;hap=0;byparent=0;prephaseinfo=0;sexes=0

! saved pedigree contains only genotyped animals

do
if(autosome==1)read(9,*,iostat=io)ani,pere,mere
if(autosome==0)read(9,*,iostat=io)ani,pere,mere,sex
if(io/=0)exit
 if(pere>ani .or. mere>ani)then
  print*,'Error: pedigree is not sorted - parent IDs must be lower than offspring IDs !!!'
  stop
 endif
 ani=newid(ani);pere=newid(pere);mere=newid(mere)
 if(ani==0)cycle
 sire(ani)=pere
 dam(ani)=mere
 sexes(ani)=sex
enddo

typ=0;genotyped=.false.;haplotyped=.false.

select case(reading)
case(0)
 do
 !read(13,'(i6,<nmarq>(i2,i2))',iostat=io)ani,(typ(newid(ani),j),j=1,2*nmarq)
 read(13,*,iostat=io)ani,(typ(newid(ani),j),j=1,2*nmarq)
 if(io/=0)exit
 ani=newid(ani)
 if(autosome==0 .and. sexes(ani)==0)cycle
 do i=1,nmarq
  if(typ(ani,2*i-1)/=0)then
   genotyped(ani)=.true.
  endif
 enddo
 enddo
case(1)
 k=0
 do
  read(13,*,iostat=io)(genoin(i),i=1,nani)
  if(io/=0)exit
  k=k+1
  do i=1,nani
   ani=newid(ids(i))
   if(autosome==0 .and. sexes(ani)==0)cycle
   if(genoin(i)/=-1)genotyped(ani)=.true.
   select case(genoin(i))
    case(-1)
     typ(ani,2*k-1)=0;typ(ani,2*k)=0
    case(0)
     typ(ani,2*k-1)=1;typ(ani,2*k)=1
    case(1)
     typ(ani,2*k-1)=1;typ(ani,2*k)=2
    case(2)
     typ(ani,2*k-1)=2;typ(ani,2*k)=2
   end select
  enddo
 enddo
end select

print*,'End of reading genotyped and pedigree files '

num=0
do i=1,nani
 if(genotyped(i))verif=verif+1
 do j=i+1,nani
   if(sire(j)==i .or. dam(j)==i)then
     num=num+1
     byparent(num,1)=i
     byparent(num,2)=j
   endif
 enddo
enddo

print*,'Number of animals with genotype ::',verif


end subroutine

!******************************************************************************************
!********************* inference tanks to homozygous markers ******************************
!******************************************************************************************

subroutine homozygotes
implicit none


do i=1,nani
 do k=1,nmarq
  if(typ(i,2*k-1)/=0)then
   if(typ(i,2*k-1)==typ(i,2*k))then
     hap(i,1,k)=typ(i,2*k-1)
     hap(i,2,k)=typ(i,2*k-1)
     haplotyped(i)=.true.
     prephaseinfo(i,k)=1
     if((k < parstart .or. k > parend) .and. sexes(i)==1)hap(i,1,k)=9
   endif
  endif
 enddo
enddo

end subroutine

!******************************************************************************************
!********************* inference tanks to mendelian rules *********************************
!******************************************************************************************


subroutine mendelian
implicit none

do i=1,nani
 if(.not.genotyped(i))cycle
 pere=sire(i)
 mere=dam(i)
 if(pere==0 .and. mere==0)cycle


 if(pere/=0 .and. mere==0)then
    if(.not.genotyped(pere))cycle
 else if(pere==0 .and. mere/=0)then
    if(.not.genotyped(mere))cycle
 else if(pere/=0 .and. mere/=0)then
    if(.not.genotyped(pere) .and. .not.genotyped(mere))cycle
 endif

 do k=1,nmarq
  if(typ(i,2*k-1)==0)cycle
  if(hap(i,1,k)/=0)cycle

  sall1=0;sall2=0;dall1=0;dall2=0  
  all1=typ(i,2*k-1);all2=typ(i,2*k)
  if(pere/=0)then
    sall1=typ(pere,2*k-1);sall2=typ(pere,2*k)
  endif  
  if(mere/=0)then
    dall1=typ(mere,2*k-1);dall2=typ(mere,2*k)
  endif  
  
  if(sall1/=0 .and. dall1/=0)then

    if(sall1==sall2 .and. dall1==dall2)then     ! both parent homozygous
      if(all1==sall1 .and. all2==dall1)then
         hap(i,1,k)=all1;hap(i,2,k)=all2
         haplotyped(i)=.true.
         prephaseinfo(i,k)=4 ! both informative
         cycle
      else if(all2==sall1 .and. all1==dall1)then
         hap(i,1,k)=all2;hap(i,2,k)=all1
         haplotyped(i)=.true.
         prephaseinfo(i,k)=4 ! both informative
         cycle
      endif
      cycle ! incompatibility - don't prephase
    endif


    if(sall1==sall2 .and. dall1/=dall2)then                         ! sire homozygous
      if(all1==sall1)then
         hap(i,1,k)=all1;hap(i,2,k)=all2
         haplotyped(i)=.true.
         prephaseinfo(i,k)=2
         cycle
      else if(all2==sall1)then
         hap(i,1,k)=all2;hap(i,2,k)=all1
         haplotyped(i)=.true.
         prephaseinfo(i,k)=2
         cycle
      endif
      cycle ! incompatibility - but offspring should be prephased as homozygous (in case of SNPs)
    endif

    if(dall1==dall2)then                           ! dam homozygous  
      if(all2==dall1)then
         hap(i,1,k)=all1;hap(i,2,k)=all2
         haplotyped(i)=.true.
         prephaseinfo(i,k)=3
         cycle
      else if(all1==dall1)then
         hap(i,1,k)=all2;hap(i,2,k)=all1
         haplotyped(i)=.true.
         prephaseinfo(i,k)=3
         cycle
      endif
      cycle
    endif

    ! parents both heterozygous (with SNPs then stop)
    ! test identical genotype for both parents
   if(sall1==dall1 .and. sall2==dall2)cycle
   if(sall1==dall2 .and. sall2==dall1)cycle
   
   if(all1==sall1 .and. all2==dall1)then
         hap(i,1,k)=all1;hap(i,2,k)=all2
         haplotyped(i)=.true.
         cycle
   else if(all1==sall1 .and. all2==dall2)then
         hap(i,1,k)=all1;hap(i,2,k)=all2
         haplotyped(i)=.true.
         cycle
   else if(all1==sall2 .and. all2==dall1)then
         hap(i,1,k)=all1;hap(i,2,k)=all2
         haplotyped(i)=.true.
         cycle
   else if(all1==sall2 .and. all2==dall2)then
         hap(i,1,k)=all1;hap(i,2,k)=all2
         haplotyped(i)=.true.
         cycle
   else if(all2==sall1 .and. all1==dall1)then
         hap(i,1,k)=all2;hap(i,2,k)=all1
         haplotyped(i)=.true.
         cycle
   else if(all2==sall1 .and. all1==dall2)then
         hap(i,1,k)=all2;hap(i,2,k)=all1
         haplotyped(i)=.true.
         cycle
   else if(all2==sall2 .and. all1==dall1)then
         hap(i,1,k)=all2;hap(i,2,k)=all1
         haplotyped(i)=.true.
         cycle
   else if(all2==sall2 .and. all1==dall2)then
         hap(i,1,k)=all2;hap(i,2,k)=all1
         haplotyped(i)=.true.
         cycle
   endif

  else if(sall1/=0 .and. dall1==0)then

    if(sall1==sall2)then                           ! sire homozygous
      prephaseinfo(i,k)=2
      if(all1==sall1)then
         hap(i,1,k)=all1;hap(i,2,k)=all2
         haplotyped(i)=.true.
         cycle
      else
         hap(i,1,k)=all2;hap(i,2,k)=all1
         haplotyped(i)=.true.
         cycle
      endif
    endif

    ! sire heterozygous
    ! test identical genotype for sire and progeny
   if(sall1==all1 .and. sall2==all2)cycle
   if(sall1==all2 .and. sall2==all1)cycle
    
   if(all1==sall1)then
         hap(i,1,k)=all1;hap(i,2,k)=all2
         haplotyped(i)=.true.
         cycle
   else if(all1==sall2)then
         hap(i,1,k)=all1;hap(i,2,k)=all2
         haplotyped(i)=.true.
         cycle
   else if(all2==sall1)then
         hap(i,1,k)=all2;hap(i,2,k)=all1
         haplotyped(i)=.true.
         cycle
   else if(all2==sall2)then
         hap(i,1,k)=all2;hap(i,2,k)=all1
         haplotyped(i)=.true.
         cycle
   endif


  else if(sall1==0 .and. dall1/=0)then

    if(dall1==dall2)then                           ! dam homozygous
      prephaseinfo(i,k)=3
      if(all2==dall1)then
         hap(i,1,k)=all1;hap(i,2,k)=all2
         haplotyped(i)=.true.
         cycle
      else
         hap(i,1,k)=all2;hap(i,2,k)=all1
         haplotyped(i)=.true.
         cycle
      endif
    endif

    ! dam heterozygous
    ! test identical genotype for dam and progeny
   if(dall1==all1 .and. dall2==all2)cycle
   if(dall1==all2 .and. dall2==all1)cycle
    
   if(all2==dall1)then
         hap(i,1,k)=all1;hap(i,2,k)=all2
         haplotyped(i)=.true.
         cycle
   else if(all2==dall2)then
         hap(i,1,k)=all1;hap(i,2,k)=all2
         haplotyped(i)=.true.
         cycle
   else if(all1==dall1)then
         hap(i,1,k)=all2;hap(i,2,k)=all1
         haplotyped(i)=.true.
         cycle
   else if(all1==dall2)then
         hap(i,1,k)=all2;hap(i,2,k)=all1
         haplotyped(i)=.true.
         cycle
   endif


  endif
 enddo
enddo
   

end subroutine


!******************************************************************************************
!*************************** phase in halfsib families ************************************
!******************************************************************************************

subroutine halfsib_phasing
implicit none
! p11: phase1 has allele 1 for next marker / p12: phase1 has allele 2 for next marker

if(phase_option==3)then
 open(102,file='genotyping_errors.txt')
 open(104,file='convergence_per_parent.txt')

 open(200,file='origins_hmm.txt')
 open(201,file='recombinations_hmm')
 open(202,file='nrec_hmm.txt')
 open(203,file='emission_parents.txt')
endif

open(300,file='origins.txt')
open(301,file='recombinations')
open(302,file='nrec.txt')
if(niter_map > 1)open(103,file='emap.txt')

allocate(phase1(nmarq),phase2(nmarq),info_marker(nmarq,3),par_error(nmarq),entropy(nmarq),nentropy(nmarq))
allocate(phased(nmarq),npar(nmarq),noff(nmarq),off_error(nmarq),tophase(nmarq))
allocate(emap(nmarq,2),probrecs(nmarq,2),numrecs(nmarq,2))

emap=0.d0
do k=1,nmarq-1
 emap(k,1)=probrec(k,k+1,1)
 emap(k,2)=probrec(k,k+1,2)
enddo
do mapiter=1,niter_map
probrecs=0.d0;numrecs=0.d0;totlik1=0.d0;totlik2=0.d0

num=1;noffspring=0;firsto=1;oldparent=byparent(1,1);maxoffspring=0; info_marker=0;ngam=0.0;par_error=0.d0
npar=0.d0;noff=0.d0;off_error=0.d0;tophase=1;nprinted=0;nemission=0;norigins1=0;norigins2=0
entropy=0.00;nentropy=0.00

!*****************************************
!******** loop: for each parent **********
!*****************************************

do while(byparent(num,1)/=0)
 do while(byparent(num,1)==oldparent)
  lasto=num
  num=num+1
  noffspring=noffspring+1
 enddo

 prephased=.FALSE.
 if(sire(oldparent)/=0)then
   if(genotyped(sire(oldparent)))prephased=.TRUE.
 endif
 if(dam(oldparent)/=0)then
   if(genotyped(dam(oldparent)))prephased=.TRUE.
 endif

 maxoffspring=max(noffspring,maxoffspring)
 phase1=0;phase2=0;mc1=0;phased=0
 if(sire(byparent(firsto,2))==oldparent)parsex=1
 if(dam(byparent(firsto,2))==oldparent)parsex=2

 
 allocate(nrec(noffspring))
 nrec=0
 
 if(noffspring<11)famsize=noffspring
 if(noffspring>10 .and. noffspring<21)famsize=11
 if(noffspring>20 .and. noffspring<51)famsize=12
 if(noffspring>50 .and. noffspring<101)famsize=13
 if(noffspring>100)famsize=14

 call listtophase
! run first time HMM to remove prephasing errors / done only in normal phasing (not when focus in on detecting map errors)
 step='prephasing'
 if(phase_option==3 .and. prephased .and. map_option==1 .and. mapiter==1)call phase_hmm(1,gerr,0) 

 ! impute SNPs to increase informativity
   ! store imputed haplotypes in hap2
   allocate(hap2(noffspring,nmarq))
   do offspring=firsto,lasto
     l=offspring-firsto+1
     hap2(l,:)=hap(byparent(offspring,2),parsex,:)
   enddo
   n_phased=0;l=0
   do offspring=firsto,lasto
     do k=1,nmarq
       if(hap2(offspring-firsto+1,k)/=0)n_phased=n_phased+1
     enddo
   enddo

 if((noffspring>2 .or. prephased) .and. phase_option>1)then

 ! initialize hmm with homozygotes markers and offspring
  ntested=0
  do refoff=firsto,lasto
   ntested=ntested+1
   if(ntested>ntemplates)exit
   phase1=0;phase2=0
   sextemplate=sexes(byparent(refoff,2))
   do k=1,nmarq
     if(prephaseinfo(byparent(refoff,2),k)==0 .or. prephaseinfo(byparent(refoff,2),k)>4)cycle ! use only mendelian markers
     if(parsex==1 .and. prephaseinfo(byparent(refoff,2),k)==2)cycle
     if(parsex==2 .and. prephaseinfo(byparent(refoff,2),k)==3)cycle
     if(typ(oldparent,2*k-1)/=typ(oldparent,2*k) .and. hap(byparent(refoff,2),parsex,k)/=0)then
       phase1(k)=hap(byparent(refoff,2),parsex,k) ! as offspring
       phase2(k)=3-hap(byparent(refoff,2),parsex,k)
     endif
     if(hap(oldparent,1,k)==9 .and. hap(byparent(refoff,2),parsex,k)/=0)then
       phase1(k)=hap(byparent(refoff,2),parsex,k) ! as offspring
       phase2(k)=3-hap(byparent(refoff,2),parsex,k)
     endif
   enddo
   call phase_hmm(0,0.0d0,1)
!   call phase_hmm(0,1.0d-3,1) 
   n_phased=0
   do offspring=firsto,lasto
     do k=1,nmarq
       if(hap2(offspring-firsto+1,k)/=0)n_phased=n_phased+1
     enddo
   enddo
  enddo
! delete imputed markers which are inconsistent (different imputed alleles with different references)
! these cases are flagged with a value 3
  do offspring=firsto,lasto
   do k=1,nmarq
    if(hap2(offspring-firsto+1,k)==3)hap2(offspring-firsto+1,k)=0
   enddo
  enddo
 
 phase1=0;phase2=0;mc1=0;phased=0;sextemplate=0

   call phase_byone(7)
   call phase_byone(7)
   do i=1,4
    call phase_byone(4)
   enddo
   call phase_byone(5)
!   call phase_byone(6)
!   call phase_hmm(0,0.d0,0) ! include ibd groups
 endif ! offspring > 2 or prephased

   if(phase_option==2)then ! store phased markers
     do k=1,nmarq 
       if(hap(oldparent,1,k)==0 .and. phase1(k)/=0 .and. phased(k)==1)then
        hap(oldparent,1,k)=phase1(k)
        hap(oldparent,2,k)=phase2(k)
       endif
     enddo
   endif
   step='hmmphasing'
   if(phase_option==3)call phase_hmm(2,gerr,0) 
   if((noffspring>2 .or. prephased) .and. mapiter == niter_map)call printrec

 deallocate(nrec)
 deallocate(hap2)
 noffspring=0;firsto=num;oldparent=byparent(num,1);info_marker=0
enddo

if(phase_option==3 .and. niter_map > 1)then ! estimate and print only if more than one iteration required
 econv1=0.d0;econv2=0.d0
 edenom1=dot_product(emap(:,1),emap(:,1))
 edenom2=dot_product(emap(:,2),emap(:,2))
 do k=1,nmarq-1
  econv1=econv1+(emap(k,1)-probrecs(k,1)/numrecs(k,1))**2
  emap(k,1)=probrecs(k,1)/numrecs(k,1)
  econv2=econv2+(emap(k,2)-probrecs(k,2)/numrecs(k,2))**2
  emap(k,2)=probrecs(k,2)/numrecs(k,2)
  if(mapiter==niter_map)write(103,'(i6,2(1x,f11.6),2(1x,f11.8))')k,posi(k,1),posi(k+1,1),emap(k,1),emap(k,2)
 enddo
 if(edenom1 > 0.00)econv1=econv1/edenom1
 if(edenom2 > 0.00)econv2=econv2/edenom2
 print'(a76,1x,i4,2(1x,ES11.3),2(1x,f11.6),2(1x,f16.2))', &
    ' Iteration for map estimation, convergence, map estimates and likelihoods ::', &
    mapiter,econv1,econv2,sum(emap(:,1)),sum(emap(:,2)),totlik1,totlik2
endif

enddo ! end iteration on map


deallocate(phase1,phase2,phased,info_marker,tophase)

if(phase_option>2)then
 close(102);close(104)
 close(200);close(201);close(202);close(203)
 !close(204)
endif
close(300);close(301);close(302)

end subroutine

subroutine write_data
implicit none
character*50 ::hfile

open(101,file='phases')

if(reading==0)then
 do i=1,nani
  if(.not.genotyped(i) .and. .not.haplotyped(i))cycle
  do j=1,2
!  write(101,'(i6,i2,1x,<nmarq>i2)')oldid(i),j,(hap(i,j,k),k=1,nmarq)
    write(101,'(i7,i2,1x)',advance='no')oldid(i),j
    do k=1,nmarq
      if(k<nmarq)write(101,'(i2)',advance='no')hap(i,j,k)
      if(k==nmarq)write(101,'(i2)',advance='yes')hap(i,j,k)
    enddo
  enddo
 enddo
 close(101)
else if(reading==1)then
 do i=1,nani
  if(.not.genotyped(i) .and. .not.haplotyped(i))cycle
  write(101,'(i7,1x,i7)',advance='no')oldid(i),oldid(i)
 enddo
 write(101,*)
 do i=1,nani
  if(.not.genotyped(i) .and. .not.haplotyped(i))cycle
  j=1
  write(101,'(i2)',advance='no')j
  j=2
  write(101,'(i2)',advance='no')j
 enddo
 write(101,*)
 do k=1,nmarq
  do i=1,nani
  if(.not.genotyped(i) .and. .not.haplotyped(i))cycle
  write(101,'(2i2)',advance='no')hap(i,1,k),hap(i,2,k)
  enddo
  write(101,*)
 enddo
 close(101)
 deallocate(hap,typ,genotyped,hapin,prephaseinfo)
 hfile='origins.txt'
 if(phase_option>1)call transpose(hfile,norigins1,1)
 hfile='origins_hmm.txt'
 if(phase_option>2)call transpose(hfile,norigins2,2)
 hfile='emission_parents.txt'
 if(phase_option>2)call transpose(hfile,nemission,2)
endif



end subroutine

!#### coding = 1 for integers and coding = 2 for reals
subroutine transpose(hfilename,nlines,coding)
implicit none
integer ::coding,io,i,j,k,l,nmark,nlines,nline,num1,num2
integer, allocatable :: table1(:,:),table2i(:,:),vali(:)
real*8, allocatable :: table2r(:,:),valr(:)
character*50 ::hfilename

allocate(table1(2,nlines))
if(coding==1)allocate(table2i(nmarq,nlines),vali(nmarq))
if(coding==2)allocate(table2r(nmarq,nlines),valr(nmarq))

open(101,file=hfilename)
nline=1

do
if(coding==1)read(101,*,iostat=io)num1,num2,(vali(i),i=1,nmarq)
if(coding==2)read(101,*,iostat=io)num1,num2,(valr(i),i=1,nmarq)
if(io/=0)exit
table1(1,nline)=num1
table1(2,nline)=num2
if(coding==1)table2i(:,nline)=vali
if(coding==2)table2r(:,nline)=valr
nline=nline+1
enddo

close(101)

open(101,file=hfilename)
do j=1,nlines
 if(j<nlines)write(101,'(i7,1x)',advance='no')table1(1,j)
 if(j==nlines)write(101,'(i7)',advance='yes')table1(1,j)
enddo
do j=1,nlines
 if(j<nlines)write(101,'(i7,1x)',advance='no')table1(2,j)
 if(j==nlines)write(101,'(i7)',advance='yes')table1(2,j)
enddo
do i=1,nmarq
 do j=1,nlines
 if(coding==1 .and. j<nlines)write(101,'(i2)',advance='no')table2i(i,j)
 if(coding==1 .and. j==nlines)write(101,'(i2)',advance='yes')table2i(i,j)
 if(coding==2 .and. j<nlines)write(101,'(1x,f8.5)',advance='no')table2r(i,j)
 if(coding==2 .and. j==nlines)write(101,'(1x,f8.5)',advance='yes')table2r(i,j)
 enddo
enddo

close(101)

deallocate(table1)
if(coding==1)deallocate(table2i,vali)
if(coding==2)deallocate(table2r,valr)

end subroutine

function probrec(m1,m2,mapi)
implicit none
integer ::m1,m2,mapi
real*8 ::probrec

probrec=0.5d0*(1.d0-dexp(-2.d0*abs(posi(m1,mapi)-posi(m2,mapi))))

end function

!****************************************************************************************************
!***************************** create list of markers to phase **************************************
!****************************************************************************************************

! this subroutine count number of heterozygotes, homozygotes
! counts number of informative offspring per marker
! identify markers to phase (1), markers already phased by Mendelian rules (2)
! excludes markers no to phase (0): homozygotes or only one informative offspring (and no Mendelian information)

subroutine listtophase
implicit none

 tophase=1;info_marker=0;nhomoz=0;nheter=0;lasttophase=0
 do k=1,nmarq
  if(parsex==1 .and. (k<parstart .or. k>parend))then
    tophase(k)=0
    if(typ(oldparent,2*k)/=0)then
      tophase(k)=2
      lasttophase=k
      phased(k)=1;phase1(k)=hap(oldparent,1,k);phase2(k)=hap(oldparent,2,k)
    endif
    cycle
  endif
  do offspring=firsto,lasto
   if(hap(byparent(offspring,2),parsex,k)/=0)info_marker(k,1)=info_marker(k,1)+1
   if(hap(byparent(offspring,2),parsex,k)==1)info_marker(k,2)=info_marker(k,2)+1
   if(hap(byparent(offspring,2),parsex,k)==2)info_marker(k,3)=info_marker(k,3)+1
  enddo
  if(typ(oldparent,2*k-1)==typ(oldparent,2*k))then
    tophase(k)=0
    if(typ(oldparent,2*k)/=0)nhomoz=nhomoz+1
  endif
  if(typ(oldparent,2*k-1)/=typ(oldparent,2*k))then
    nheter=nheter+1 
    if(hap(oldparent,1,k)/=0)then
      tophase(k)=2
      lasttophase=k
      phased(k)=1;phase1(k)=hap(oldparent,1,k);phase2(k)=hap(oldparent,2,k)
      cycle ! marker already phased by Mendelian rules, keep it
    endif
    if(info_marker(k,1)<2)tophase(k)=0 ! remove markers informative in only one offspring -- ideally, even if two informative must be connected to other markers
  endif
  if(tophase(k)>0)lasttophase=k
 enddo

end subroutine


subroutine phase_hmm(printout,epsil,impute)
implicit none
integer ::printout,impute
real*8 ::epsil

epsi=epsil

allocate(alpha(2,nmarq),beta(2,nmarq),gamma(2,nmarq),bjk(nmarq,2,2),scaling(nmarq))
allocate(L11(nmarq),L12(nmarq),L21(nmarq),L22(nmarq))
allocate(gammas(noffspring,nmarq))
allocate(expco(noffspring))

if(nhomoz>0 .and. nheter==0)then ! parent is 'inbred'

 ! output gammas (0.5) / estimation of recombination rate but number of heterozygote markers indicated in output
 ! remove conflicting opposite homozygotes (if not done in GWASH)

 do k=1,nmarq
  if(hap(oldparent,1,k)==0)cycle
  do offspring=firsto,lasto
     if(hap(byparent(offspring,2),parsex,k)==0)cycle
     if(hap(byparent(offspring,2),parsex,k)/=hap(oldparent,1,k))then
       hap(byparent(offspring,2),1,k)=0
       hap(byparent(offspring,2),2,k)=0
     endif   
  enddo
 enddo

endif
! else ! perform also hmm for inbred

 call initialize_hmm(impute)

 do while(conv>1e-16)
  round=round+1
  call reset_lik

  numoff=0
   do offspring=firsto,lasto
    call forward(impute,printout)
    call backward(impute)
    if(impute==0)call add_contribution
   enddo ! offspring
  if(impute==0)call update_bjk
  if(impute==1)call impute_hmm
  if(round==1000 .and. impute==0)call erase_errors
  if(impute==1)conv=1e-17
  !conv=1e-17

  if(conv<1e-16 .or. round==2000)exit
 
 enddo ! new parameters
 if(parsex==1 .and. printout==2 .and. (noffspring>2 .or. prephased))then
  if(parstart < parend)totlik1=totlik1+lik
 endif
 if(parsex==2 .and. printout==2 .and. (noffspring>2 .or. prephased))totlik2=totlik2+lik


  if(impute == 0)then
  numoff=0
   do offspring=firsto,lasto
    call forward(impute,printout)
    call backward(impute)
    if(printout==2 .and. (noffspring>2 .or. prephased))call ejump !### new subroutine with expected jumps prob 
   enddo ! offspring
  endif


 if(printout>0 .and. mapiter==niter_map)call output_hmm(printout)

! endif ! inbred parent or not


deallocate(alpha,beta,gamma,bjk,scaling,L11,L12,L21,L22)
deallocate(gammas,expco)

end subroutine

subroutine initialize_hmm(impute)
implicit none
integer ::impute

!************************************************************************************
!***************************** INITIALIZE *******************************************
!************************************************************************************
! either prephased, choose one offspring or random?

 if(impute==0)then
 do k=1,nmarq ! unphase non-phased 
  if(phased(k)==0)then
   phase1(k)=0
   phase2(k)=0
  endif
 enddo

 do k=1,nmarq ! phase homozygotes (not used up to know)
  if(hap(oldparent,1,k)==hap(oldparent,2,k))then
    if(hap(oldparent,1,k)==1)then
      phase1(k)=1;phase2(k)=1
    endif
    if(hap(oldparent,1,k)==2)then
      phase1(k)=2;phase2(k)=2
    endif
  endif
 enddo
endif

do k=1,nmarq
 if(phase1(k)==1)then
   bjk(k,1,1)=1.d0-epsi;bjk(k,1,2)=epsi
 else if(phase1(k)==2)then
   bjk(k,1,1)=epsi;bjk(k,1,2)=1.d0-epsi
 else if(phase1(k)==0)then
   bjk(k,1,1)=0.500;bjk(k,1,2)=0.500
 else if(phase1(k)==9)then ! sire on X (bjk are not used in emission)
   bjk(k,1,1)=0.500;bjk(k,1,2)=0.500
 endif

 if(phase2(k)==1)then
   bjk(k,2,1)=1.d0-epsi;bjk(k,2,2)=epsi
 else if(phase2(k)==2)then
   bjk(k,2,1)=epsi;bjk(k,2,2)=1.d0-epsi
 else if(phase2(k)==0)then
   bjk(k,2,1)=0.500;bjk(k,2,2)=0.500
 endif

enddo

conv=1.d0;round=0

end subroutine

subroutine reset_lik
implicit none

lik=0
L11=epsi/3;L12=epsi/3;L21=epsi/3;L22=epsi/3
do k=1,nmarq
 if(hap(oldparent,1,k)==1 .and. hap(oldparent,2,k)==1)then
  L11(k)=1.d0-epsi
 else if (hap(oldparent,1,k)==1 .and. hap(oldparent,2,k)==2)then
  L12(k)=1.d0-epsi
 else if (hap(oldparent,1,k)==2 .and. hap(oldparent,2,k)==1)then
  L21(k)=1.d0-epsi
 else if (hap(oldparent,1,k)==2 .and. hap(oldparent,2,k)==2)then
  L22(k)=1.d0-epsi
 else if (hap(oldparent,1,k)==9 .and. hap(oldparent,2,k)==1)then ! sires on X 
  L11(k)=1.d0-epsi
 else if (hap(oldparent,1,k)==9 .and. hap(oldparent,2,k)==2)then ! sires on X
  L22(k)=1.d0-epsi
 else if (hap(oldparent,1,k)==0 .and. hap(oldparent,2,k)==0 .and. typ(oldparent,2*k-1)/=0)then ! should be heterozygous
  L12(k)=0.5000-epsi/2.d0;L21(k)=0.5000-epsi/2.d0;L11(k)=epsi/2.d0;L22(k)=epsi/2.d0
 else if (hap(oldparent,1,k)==0 .and. hap(oldparent,2,k)==0 .and. typ(oldparent,2*k-1)==0)then 
  L12(k)=0.25;L21(k)=0.25;L11(k)=0.25;L22(k)=0.25
 endif  
enddo

end subroutine

!************************************************************************************
!***************************** FORWARD **********************************************
!************************************************************************************

subroutine forward(ibd,printlik)
implicit none
integer ::ibd,printlik

 alpha=0.0;scaling=0.0;pi=0.5

 ! initialisation

 do i=1,2
  prob_emission=emission(byparent(offspring,2),i,1,parsex)
  alpha(i,1)=pi*prob_emission
  scaling(1)=scaling(1)+alpha(i,1)
 enddo

 scaling(1)=1.0/scaling(1)
 alpha(:,1)=alpha(:,1)*scaling(1)

! induction: to get to i, two ways:
! 1) transition + jump into cluster i
! 2) no transition and in i at previous position
 
do k=2,nmarq
 pjump=jump(k-1,ibd)
 do i=1,2
    prob_emission=emission(byparent(offspring,2),i,k,parsex)

    alpha(i,k)=alpha(i,k-1)*pjump(0)   
    alpha(i,k)=alpha(i,k)+alpha(3-i,k-1)*pjump(1)

    alpha(i,k)=alpha(i,k)*prob_emission 
    scaling(k)=scaling(k)+alpha(i,k)

 enddo ! end i
scaling(k)=1.0/scaling(k)
alpha(:,k)=alpha(:,k)*scaling(k)
enddo ! end marker 

! termination

val=-sum(log(scaling(:)))
lik=lik+val

end subroutine

!************************************************************************************
!***************************** BACKWARD *********************************************
!************************************************************************************

subroutine backward(ibd)
implicit none
integer ::ibd
real*8 :: trans,notrans

beta=0.0;gamma=0.0

! initialisation

do i=1,2
 gamma(i,nmarq)=alpha(i,nmarq)*1.0  ! beta(i,nmarq)=1.0
 beta(i,nmarq)=1.0*scaling(nmarq)
enddo

! induction
! to arrive in k: with or without transition

do k=nmarq-1,1,-1
 pjump=jump(k,ibd)
 do i=1,2
    prob_emission=emission(byparent(offspring,2),i,k+1,parsex)
 
    beta(i,k)=beta(i,k)+pjump(0)*prob_emission*beta(i,k+1)
    beta(3-i,k)=beta(3-i,k)+pjump(1)*prob_emission*beta(i,k+1)
 enddo

 do i=1,2
  beta(i,k)=beta(i,k)*scaling(k)
  gamma(i,k)=alpha(i,k)*beta(i,k)/scaling(k)
 enddo
enddo

numoff=numoff+1
gammas(numoff,:)=gamma(1,:)

end subroutine


!************************************************************************************
!************************  Expected transition  *************************************
!************************************************************************************

subroutine ejump
implicit none
integer ::l0
real*8 ::pt11,pt12,pt21,pt22

! estimate RECS and NON-RECS

l0=offspring-firsto+1
expco(l0)=0.000
do k=1,nmarq-1
 pjump=jump(k,0) !# recombination prob
 prob_emission=emission(byparent(offspring,2),2,k+1,parsex)
 pt12=alpha(1,k)*pjump(1)*prob_emission*beta(2,k+1)
 pt22=alpha(2,k)*pjump(0)*prob_emission*beta(2,k+1)
 prob_emission=emission(byparent(offspring,2),1,k+1,parsex)
 pt21=alpha(2,k)*pjump(1)*prob_emission*beta(1,k+1)
 pt11=alpha(1,k)*pjump(0)*prob_emission*beta(1,k+1)
 probrecs(k,parsex)=probrecs(k,parsex)+pt12+pt21
 numrecs(k,parsex)=numrecs(k,parsex)+pt11+pt12+pt21+pt22
 expco(l0)=expco(l0)+(pt12+pt21)
enddo

end subroutine


! ******************************** estimation of bjk ************************************ 
! **********  = expected counting in i and allele 1 / expected counting in i ************
! ***************************************************************************************

subroutine add_contribution
implicit none

do k=1,nmarq
  if(parsex==1 .and. sexes(byparent(offspring,2))==1 .and. (k<parstart .or. k>parend))cycle
  if(prephaseinfo(byparent(offspring,2),k)==5)cycle ! don't use prephased based on Linkage
  if(parsex==1 .and. prephaseinfo(byparent(offspring,2),k)==2)cycle ! if parent sire, don't use prephased based on sire-mendelian
  if(parsex==2 .and. prephaseinfo(byparent(offspring,2),k)==3)cycle ! if parent dam, don't use prephased based on dam-mendelian
  if(hap(byparent(offspring,2),parsex,k)==1)then
    L11(k)=L11(k)*(1.d0-epsi)
    L22(k)=L22(k)*epsi
    L12(k)=L12(k)*(gamma(1,k)*(1.d0-epsi)+gamma(2,k)*epsi)
    L21(k)=L21(k)*(gamma(1,k)*epsi+gamma(2,k)*(1.d0-epsi))
! rescale to 1.000
    PTOT=L11(k)+L12(k)+L21(k)+L22(k)
    L11(k)=L11(k)/PTOT;L12(k)=L12(k)/PTOT;L21(k)=L21(k)/PTOT;L22(k)=L22(k)/PTOT
  else if(hap(byparent(offspring,2),parsex,k)==2)then
    L11(k)=L11(k)*epsi 
    L22(k)=L22(k)*(1.d0-epsi) 
    L12(k)=L12(k)*(gamma(1,k)*epsi+gamma(2,k)*(1.d0-epsi))
    L21(k)=L21(k)*(gamma(1,k)*(1.d0-epsi)+gamma(2,k)*epsi)
    PTOT=L11(k)+L12(k)+L21(k)+L22(k)
    L11(k)=L11(k)/PTOT;L12(k)=L12(k)/PTOT;L21(k)=L21(k)/PTOT;L22(k)=L22(k)/PTOT
  endif
enddo

end subroutine

! ******************************** estimation of bjk ************************************ 
! **********  = expected counting in i and allele 1 / expected counting in i ************
! ***************************************************************************************

subroutine update_bjk
implicit none

conv=0.d0
do k=1,nmarq
 if(phase1(k)==0)cycle !!! work only on prephased markers / HMM doesn't phase nicely in non-informative regions => inflated errors
 b11=L11(k)+L12(k)
 b12=L22(k)+L21(k)
 b21=L11(k)+L21(k)
 b22=L22(k)+L12(k)
 if(0.5*(bjk(k,1,1)+b11)>1e-16)conv=max(conv,((bjk(k,1,1)-b11)**2)/((0.5*bjk(k,1,1)+0.5*b11)**2)) !!! take 1e-16 to be sure bjk**2 is still computable
 if(0.5*(bjk(k,2,1)+b21)>1e-16)conv=max(conv,((bjk(k,2,1)-b21)**2)/((0.5*bjk(k,2,1)+0.5*b21)**2))
 if(round==2000)print'(2i6,5(1x,f16.9))',oldid(oldparent),k,b11,b21,bjk(k,1,1),bjk(k,2,1),conv
 bjk(k,1,1)=b11
 bjk(k,1,2)=b12
 bjk(k,2,1)=b21
 bjk(k,2,2)=b22
enddo

end subroutine

! erase possible errors to improve convergence

subroutine erase_errors
! remove potential errors in offspring to improve convergence after 1000 iterations
! parent genotypes are not changed (because the program estimates the emission probabilities)

do k=1,nmarq
 if(parsex==1 .and. (k < parstart .or. k > parend))cycle
 if(phase1(k)==0)cycle !!! work only on prephased markers / HMM doesn't phase nicely in non-informative regions => inflated errors
 do l=1,noffspring
   if(hap(byparent(firsto+l-1,2),parsex,k)/=0)then
     if(hap(byparent(firsto+l-1,2),parsex,k)==1)p_error=gammas(l,k)*bjk(k,1,2)+(1-gammas(l,k))*bjk(k,2,2)
     if(hap(byparent(firsto+l-1,2),parsex,k)==2)p_error=gammas(l,k)*bjk(k,1,1)+(1-gammas(l,k))*bjk(k,2,1)
     if(L11(K)>0.95 .or. (L12(K)+L21(k))>0.95 .or. L22(k)<0.95)then ! counts errors only when likelihood for one configuration is high / avoid uninformative positions
      if(famsize>9)noff(k)=noff(k)+1.00
      if(famsize>9)off_error(k)=off_error(k)+p_error
     endif
     if(p_error>0.95)then
      if(typ(byparent(firsto+l-1,2),2*k)/=0 .and. typ(byparent(firsto+l-1,2),2*k)==typ(byparent(firsto+l-1,2),2*k-1))then ! prephasing and genotyping errors
       write(102,'(i7,1x,i7,1x,a10,1x,a7,1x,2i1,5(1x,f7.4),1x,i5)')oldid(byparent(firsto+l-1,2)),k,step,'progeny',&
          typ(byparent(firsto+l-1,2),2*k),typ(byparent(firsto+l-1,2),2*k-1),p_error,L11(k),L12(k),L21(k),L22(k),noffspring 
       typ(byparent(firsto+l-1,2),2*k-1)=0
       typ(byparent(firsto+l-1,2),2*k)=0
      endif
      hap(byparent(firsto+l-1,2),1,k)=0
      hap(byparent(firsto+l-1,2),2,k)=0
      prephaseinfo(byparent(firsto+l-1,2),k)=0
     endif
   endif
 enddo
enddo

end subroutine

!************************************************************************************
!***************************** output all results ***********************************
!************************************************************************************


subroutine output_hmm(printout)
implicit none
integer ::printout,ninfor
! printout = 1 / print only errors (during prephasing checking)
! printout = 2 / print also convergence, emission probabilities, recombinations, origins in HMM

 if(printout==2)write(104,*)oldid(oldparent),noffspring,prephased,round,conv
 do k=1,nmarq
 if(parsex==1 .and. (k<parstart .or. k>parend))cycle
  if(typ(oldparent,2*k)/=0 .and. famsize>9 .and. phase1(k)/=0)then
    npar(k)=npar(k)+1.d00
    if((typ(oldparent,2*k-1)+typ(oldparent,2*k))==2)par_error(k)=par_error(k)+1.d0-L11(k)
    if((typ(oldparent,2*k-1)+typ(oldparent,2*k))==3)par_error(k)=par_error(k)+1.d0-L12(k)-L21(k)
    if((typ(oldparent,2*k-1)+typ(oldparent,2*k))==4)par_error(k)=par_error(k)+1.d0-L22(k)
!    if((typ(oldparent,2*k-1)+typ(oldparent,2*k))==2)par_error(k)=par_error(k)+1.d0-bjk(k,1,1)*bjk(k,2,1)
!    if((typ(oldparent,2*k-1)+typ(oldparent,2*k))==3)par_error(k)=par_error(k)+1.d0-(bjk(k,1,1)*bjk(k,2,2)+bjk(k,1,2)*bjk(k,2,1))
!    if((typ(oldparent,2*k-1)+typ(oldparent,2*k))==4)par_error(k)=par_error(k)+1.d0-bjk(k,1,2)*bjk(k,2,2)
  endif
  nh1=0.00;nh2=0.00;nall11=0.00;nall12=0.00;nall21=0.00;nall22=0.00
  do l=1,noffspring
!   if(k==1 .and. famsize>9)ngam=ngam+1.00
   if(hap(byparent(firsto+l-1,2),parsex,k)/=0 .and. phase1(k)/=0)then
    if(gammas(l,k)>0.95)then
      if(hap(byparent(firsto+l-1,2),parsex,k)==1)nall11=nall11+1.00
      if(hap(byparent(firsto+l-1,2),parsex,k)==2)nall12=nall12+1.00
      nh1=nh1+1.00
    else if(gammas(l,k)<0.05)then
      if(hap(byparent(firsto+l-1,2),parsex,k)==1)nall21=nall21+1.00
      if(hap(byparent(firsto+l-1,2),parsex,k)==2)nall22=nall22+1.00
      nh2=nh2+1.00
    endif
    if(hap(byparent(firsto+l-1,2),parsex,k)==1)p_error=gammas(l,k)*bjk(k,1,2)+(1-gammas(l,k))*bjk(k,2,2)
    if(hap(byparent(firsto+l-1,2),parsex,k)==2)p_error=gammas(l,k)*bjk(k,1,1)+(1-gammas(l,k))*bjk(k,2,1)
     if(famsize>9 .and. round<1000)then ! avoid counting twice
       off_error(k)=off_error(k)+p_error 
       noff(k)=noff(k)+1.00
    endif
    if(p_error>0.95)then
       if(typ(byparent(firsto+l-1,2),2*k)/=0 .and. typ(byparent(firsto+l-1,2),2*k)==typ(byparent(firsto+l-1,2),2*k-1))then ! prephasing and genotyping errors
        write(102,'(i7,1x,i7,1x,a10,1x,a7,1x,2i1,5(1x,f7.4),1x,i5)')oldid(byparent(firsto+l-1,2)),k,step,'progeny',&
          typ(byparent(firsto+l-1,2),2*k),typ(byparent(firsto+l-1,2),2*k-1),p_error,L11(k),L12(k),L21(k),L22(k),noffspring
        typ(byparent(firsto+l-1,2),2*k-1)=0
        typ(byparent(firsto+l-1,2),2*k)=0
       endif 
       hap(byparent(firsto+l-1,2),1,k)=0
       hap(byparent(firsto+l-1,2),2,k)=0
       prephaseinfo(byparent(firsto+l-1,2),k)=0
    endif
   endif
  enddo
  entro1=0.00;entro2=0.00
  if(nh1>0.00 .and. nall11>0.00)entro1=entro1+nall11/nh1*log(nall11/nh1)
  if(nh1>0.00 .and. nall12>0.00)entro1=entro1+nall12/nh1*log(nall12/nh1)
  if(nh2>0.00 .and. nall21>0.00)entro2=entro2+nall21/nh2*log(nall21/nh2)
  if(nh2>0.00 .and. nall22>0.00)entro2=entro2+nall22/nh2*log(nall22/nh2)
  entropy(k)=entropy(k)+nh1*entro1+nh2*entro2
  nentropy(k)=nentropy(k)+nh1+nh2
 enddo

!**** output prob_origins.txt, recombinations and number of recombinations / offspring from HMM

 if(printout==2 .and. (noffspring>2 .or. prephased))then
 do l=1,noffspring
   mate=0
   if(parsex==1 .and. dam(byparent(firsto+l-1,2))/=0)then
     if(genotyped(dam(byparent(firsto+l-1,2))))mate=1
   endif
   if(parsex==2 .and. sire(byparent(firsto+l-1,2))/=0)then
     if(genotyped(sire(byparent(firsto+l-1,2))))mate=1
   endif

   lastknown=0;lastori=0;ninfor=0
   write(200,'(i7,1x,i7)',advance='no')oldid(byparent(firsto+l-1,2)),oldid(oldparent)
   norigins2=norigins2+1
   do k=1,nmarq
     if(k<nmarq)write(200,'(1x,f7.4)',advance='no')gammas(l,k)
     if(k==nmarq)write(200,'(1x,f7.4)',advance='yes')gammas(l,k)
      if(gammas(l,k)>0.999 .and. (L12(k)>0.95 .or. L21(k)>0.95 .or. (parsex==1 .and. (k<parstart .or. k>parend))) &
          .and. hap(byparent(firsto+l-1,2),parsex,k)/=0)then
     if(lastknown/=0)then
       if(lastori==2)then
         write(201,'(4(1x,i7))')oldid(byparent(firsto+l-1,2)),oldid(oldparent),lastknown,k
         nrec(l)=nrec(l)+1
       endif
      endif
      lastknown=k;lastori=1
      ninfor=ninfor+1
     endif
      if(gammas(l,k)<0.001 .and. (L12(k)>0.95 .or. L21(k)>0.95 .or. (parsex==1 .and. (k<parstart .or. k>parend))) &
        .and. hap(byparent(firsto+l-1,2),parsex,k)/=0)then
      if(lastknown/=0)then
       if(lastori==1)then
         write(201,'(4(1x,i7))')oldid(byparent(firsto+l-1,2)),oldid(oldparent),lastknown,k
         nrec(l)=nrec(l)+1
       endif
      endif
      lastknown=k;lastori=2
      ninfor=ninfor+1
     endif
   enddo
   k=0
   if(prephased)k=1
   if(prephased .or. noffspring>2)write(202,'(2(i6,1x),i1,1x,i4,2(1x,i1),3(1x,i6),1x,i4,1x,f8.3)')&
     oldid(byparent(firsto+l-1,2)),oldid(oldparent),parsex,noffspring,k,mate,nheter,nhomoz,ninfor,nrec(l),expco(l)
enddo
endif


! keep as phased markers with high emission probability / don't change markers which are non prephased (model works poorly in less informative regions)
! work with both alleles simultaneously because imputation and phasing are not symmetric (both alleles are not always informative)
 do k=1,nmarq 
 if(phase1(k)==0)cycle !!! work only on prephased markers / HMM doesn't phase nicely in non-informative regions => inflated errors
 hap(oldparent,:,k)=0 ! 
! if(L11(k)>=(1.d0-epsi))hap(oldparent,1,k)=1
! if(L11(k)>=(1.d0-epsi))hap(oldparent,2,k)=1
! if(L12(k)>=(1.d0-epsi))hap(oldparent,1,k)=1
! if(L12(k)>=(1.d0-epsi))hap(oldparent,2,k)=2
! if(L21(k)>=(1.d0-epsi))hap(oldparent,1,k)=2
! if(L21(k)>=(1.d0-epsi))hap(oldparent,2,k)=1
! if(L22(k)>=(1.d0-epsi))hap(oldparent,1,k)=2
! if(L22(k)>=(1.d0-epsi))hap(oldparent,2,k)=2
 if(L11(k)>0.95)hap(oldparent,1,k)=1
 if(L11(k)>0.95)hap(oldparent,2,k)=1
 if(L12(k)>0.95)hap(oldparent,1,k)=1
 if(L12(k)>0.95)hap(oldparent,2,k)=2
 if(L21(k)>0.95)hap(oldparent,1,k)=2
 if(L21(k)>0.95)hap(oldparent,2,k)=1
 if(L22(k)>0.95)hap(oldparent,1,k)=2
 if(L22(k)>0.95)hap(oldparent,2,k)=2
 if(parsex==1 .and. (k<parstart .or. k>parend) .and. hap(oldparent,2,k)/=0)hap(oldparent,1,k)=9
 if((typ(oldparent,2*k-1)+typ(oldparent,2*k))==2 .and. (1.d0-L11(k))>0.95)&
      write(102,'(i7,1x,i7,1x,a10,1x,a7,1x,2i1,5(1x,f7.4),1x,i5)')oldid(oldparent),k,step,'parents',typ(oldparent,2*k-1),&
      typ(oldparent,2*k),(1.d0-L11(k)),L11(k),L12(k),L21(k),L22(k),noffspring
 if((typ(oldparent,2*k-1)+typ(oldparent,2*k))==2 .and. (1.d0-L11(k))>0.95)typ(oldparent,(2*k-1):(2*k))=0
 if((typ(oldparent,2*k-1)+typ(oldparent,2*k))==3 .and. (1.d0-L12(k)-L21(k))>0.95)&
   write(102,'(i7,1x,i7,1x,a10,1x,a7,1x,2i1,5(1x,f7.4),1x,i5)')oldid(oldparent),k,step,'parents',typ(oldparent,2*k-1),&
    typ(oldparent,2*k),(1.d0-L12(k)-L21(k)),L11(k),L12(k),L21(k),L22(k),noffspring
 if((typ(oldparent,2*k-1)+typ(oldparent,2*k))==3 .and. (1.d0-L12(k)-L21(k))>0.95)typ(oldparent,(2*k-1):(2*k))=0
 if((typ(oldparent,2*k-1)+typ(oldparent,2*k))==4 .and. (1.d0-L22(k))>0.95)&
    write(102,'(i7,1x,i7,1x,a10,1x,a7,1x,2i1,5(1x,f7.4),1x,i5)')oldid(oldparent),k,step,'parents',typ(oldparent,2*k-1),&
    typ(oldparent,2*k),(1.d0-L22(k)),L11(k),L12(k),L21(k),L22(k),noffspring
 if((typ(oldparent,2*k-1)+typ(oldparent,2*k))==4 .and. (1.d0-L22(k))>0.95)typ(oldparent,(2*k-1):(2*k))=0

 if(typ(oldparent,2*k-1)==0)then ! remove prephasing in progeny
  do l=1,noffspring
   if(parsex==1 .and. prephaseinfo(byparent(firsto+l-1,2),k)==2)then
      hap(byparent(firsto+l-1,2),1,k)=0
      hap(byparent(firsto+l-1,2),2,k)=0
      prephaseinfo(byparent(firsto+l-1,2),k)=0
   endif
   if(parsex==2 .and. prephaseinfo(byparent(firsto+l-1,2),k)==3)then
      hap(byparent(firsto+l-1,2),1,k)=0
      hap(byparent(firsto+l-1,2),2,k)=0
      prephaseinfo(byparent(firsto+l-1,2),k)=0
   endif
  enddo
 endif

 enddo

! output haplotypes of parents as emission probabilities

 if(printout==2)then
 ori=1
 write(203,'(i7,1x,i1)',advance='no')oldid(oldparent),ori
 nemission=nemission+1
 do k=1,nmarq
  if(parsex==1 .and. (k<parstart .or. k>parend))bjk(k,1,1)=9.d0
  if(k<nmarq)write(203,'(1x,f8.5)',advance='no')bjk(k,1,1)
  if(k==nmarq)write(203,'(1x,f8.5)',advance='yes')bjk(k,1,1)
 enddo

 ori=2
 write(203,'(i7,1x,i1)',advance='no')oldid(oldparent),ori
 nemission=nemission+1
 do k=1,nmarq
  if(k<nmarq)write(203,'(1x,f8.5)',advance='no')bjk(k,2,1)
  if(k==nmarq)write(203,'(1x,f8.5)',advance='yes')bjk(k,2,1)
 enddo
 endif

end subroutine


! *************************************************************************************** 
! ******************  impute based on IBDs with a give offspring ************************
! ***************************************************************************************

subroutine impute_hmm
implicit none
integer :: ii,jj,allibd1,allibd2

 do k=1,nmarq
! first check if all 'ibd' haplotypes carry the same allele
   if(parsex==1 .and. (k<parstart .or. k>parend))cycle
   do offspring=firsto,lasto
    allibd1=1;allibd2=1
    l=offspring-firsto+1
    if(gammas(l,k)>0.999999 .and. hap(byparent(offspring,2),parsex,k)/=phase1(k))allibd1=0
    if(gammas(l,k)<0.000001 .and. hap(byparent(offspring,2),parsex,k)/=phase2(k))allibd2=0
!    if(gammas(l,k)>0.999 .and. hap(byparent(offspring,2),parsex,k)/=phase1(k))allibd1=0
!    if(gammas(l,k)<0.001 .and. hap(byparent(offspring,2),parsex,k)/=phase2(k))allibd2=0
   enddo

   do offspring=firsto,lasto
    l=offspring-firsto+1
    if(hap2(l,k)>2 .or. tophase(k)==0 .or. hap(byparent(offspring,2),parsex,k)/=0 .or. phase1(k)==0)cycle ! if phase2(k)=0 then phase1(k)=0
    if(hap2(l,k)>0)then ! allready imputed once
     if(gammas(l,k)>0.999999 .and. hap2(l,k)/=phase1(k))hap2(l,k)=3 ! 3 => as been imputed with two different alleles; reset to 0 latter
     if(gammas(l,k)<0.000001 .and. hap2(l,k)/=phase2(k))hap2(l,k)=3
!     if(gammas(l,k)>0.999 .and. hap2(l,k)/=phase1(k))hap2(l,k)=3 ! 3 => as been imputed with two different alleles; reset to 0 latter
!     if(gammas(l,k)<0.001 .and. hap2(l,k)/=phase2(k))hap2(l,k)=3
    else ! not already imputed
     if(gammas(l,k)>0.999999 .and. allibd1==1)hap2(l,k)=phase1(k)
     if(gammas(l,k)<0.000001 .and. allibd2==1)hap2(l,k)=phase2(k)
!     if(gammas(l,k)>0.999 .and. allibd1==1)hap2(l,k)=phase1(k)
!     if(gammas(l,k)<0.001 .and. allibd2==1)hap2(l,k)=phase2(k)
    endif
   enddo
 enddo

end subroutine



function jump(marker,ibd)
implicit none
integer ::marker,ibd
real*8 ::jump(0:1)

! first indice 0 or 1 tells if the first haplotype jumped or not
! second indice for the second haplotype

if(mapiter==1)then
 if(ibd==0)then
  jump(0)=1.0-probrec(marker,marker+1,parsex)
  jump(1)=probrec(marker,marker+1,parsex)
 else if(ibd==1)then
  jump(0)=(1.0-probrec(marker,marker+1,parsex))**2+probrec(marker,marker+1,parsex)**2
  jump(1)=2.0*probrec(marker,marker+1,parsex)*(1.0-probrec(marker,marker+1,parsex))
 endif
else
 if(ibd==0)then
  jump(0)=1.0-emap(marker,parsex)
  jump(1)=emap(marker,parsex)
 else if(ibd==1)then
  jump(0)=(1.0-emap(marker,parsex))**2+emap(marker,parsex)**2
  jump(1)=2.0*emap(marker,parsex)*(1.0-emap(marker,parsex))
 endif
endif

end function

function emission(animal,cluster1,marker,parentalorigin)
implicit none
integer ::animal,cluster1,marker,genocompa,parentalorigin
real*8 ::emission

 if(parsex==1 .and. (marker < parstart .or. marker > parend))then
  emission=1.d0
  if(sextemplate==0)then ! for hmm_phasing
   if(sexes(animal)==1 .and. cluster1==2)emission=0.d0
   if(sexes(animal)==2 .and. cluster1==1)emission=0.d0
  else ! in case of 'imputation' HMM on X 
   if(sexes(animal)/=sextemplate)emission=0.d0
  endif
  return
 endif

! if not phased based on homozygous or mendelian, non-informative
! use only those prephased based on homozygous or mendelian from second parent
! prephasing from same parent introduces a bias (if genotyping error in that parent)
 emission=0.d0
 if(hap(animal,parentalorigin,marker)==0 .or. prephaseinfo(animal,marker)==5)emission=1.d0
 if(parentalorigin==1 .and. prephaseinfo(animal,marker)==2)emission=1.d0
 if(parentalorigin==2 .and. prephaseinfo(animal,marker)==3)emission=1.d0
 if(emission==1.d0)then
  return
 endif
 emission=1.d0

! else, emission is simply bjk
 if(hap(animal,parentalorigin,marker)==1)then
  emission=(1.d0-epsi)*bjk(marker,cluster1,1)+epsi*bjk(marker,cluster1,2)
 else if (hap(animal,parentalorigin,marker)==2)then
  emission=(1.d0-epsi)*bjk(marker,cluster1,2)+epsi*bjk(marker,cluster1,1)
 endif
 return
 
end function

subroutine requiem(gender)
implicit none
integer ::i,j,k,l,io,nmark,mark,gender
integer ::ani,parent,m1,m2,nrec,iter
real ::rec,scaling,control,composite,posi1,posi2
real,allocatable ::initrec(:),postrec(:),recounts(:)
integer, allocatable ::rectable(:,:)

open(512,file='recombinations')
if(gender==0)open(513,file='MCS.txt') ! map confidence score for males and females
if(gender==1)open(513,file='RRm.txt') ! map confidence score for males
if(gender==2)open(513,file='RRf.txt') ! map confidence score for females

! this is not the number of markers but intervals = nmark-1
nmark=nmarq

allocate(initrec(nmark),postrec(nmark),recounts(nmark))
initrec=0.00

do k=1,nmark-1
 if(gender==0)initrec(k)=0.5*(probrec(k,k+1,1)+probrec(k,k+1,2))
 if(gender==1)initrec(k)=probrec(k,k+1,1)
 if(gender==2)initrec(k)=probrec(k,k+1,2)
 if(initrec(k)<0.000001)initrec(k)=0.000001
enddo

nrec=0
do
read(512,*,iostat=io)ani,parent,m1,m2
if(io/=0)exit
if(gender==2 .and. newid(parent)==sire(newid(ani)))cycle
if(gender==1 .and. newid(parent)==dam(newid(ani)))cycle
nrec=nrec+1
enddo

rewind(512)

allocate(rectable(nrec,3))
rectable=0

nrec=0
do
read(512,*,iostat=io)ani,parent,m1,m2
if(io/=0)exit
if(gender==2 .and. newid(parent)==sire(newid(ani)))cycle
if(gender==1 .and. newid(parent)==dam(newid(ani)))cycle
nrec=nrec+1
rectable(nrec,1)=parent
rectable(nrec,2)=m1
rectable(nrec,3)=m2
enddo

postrec=initrec

do iter=1,100
recounts=0.00
do i=1,nrec
 m1=rectable(i,2);m2=rectable(i,3)
 scaling=sum(postrec(m1:(m2-1)))
 control=0.00
 do k=m1,m2-1
  control=control+postrec(k)/scaling
  recounts(k)=recounts(k)+postrec(k)/scaling
 enddo
enddo
postrec=recounts/ngam(gender)
enddo

 
entropy=-entropy/log(2.d0)
do k=1,nmark
 if(npar(k)/=0.00)par_error(k)=par_error(k)/npar(k)
! if(noff(k)/=0.00)off_error(k)=off_error(k)/noff(k)
 if(nentropy(k)/=0.00)entropy(k)=entropy(k)/nentropy(k)
 composite=(1.d00-postrec(k))*(1.d00-entropy(k))*(1.d00-par_error(k))
 posi1=posi(k,1)
 if(k<nmarq)posi2=posi(k+1,1)
 if(k==nmarq)posi2=posi(k,1)+1.00
 if(gender==0)then
  write(513,'(i7,2(1x,f11.6),1x,f16.13,2(1x,f8.1,1x,f8.5),1x,f8.5)')k,posi1,posi2,postrec(k), &
    npar(k),par_error(k),nentropy(k),entropy(k),composite
 else
  write(513,'(i7,2(1x,f11.6),1x,f16.13)')k,posi1,posi2,postrec(k)
 endif
enddo

if(gender==0)then
 print*,''
 print*,'Gender = 0 for males and females'
 print*,'Gender = 1 for males'
 print*,'Gender = 2 for females'
 print*,''
endif
print'(a12,i1,a66,i10,1x,f11.1,1x,f10.4)',' For gender ',gender, &
 ', number of recombinations, gametes and length of the map (in cM) ::',nrec, &
 ngam(gender),100.0*sum(postrec(:))

if(gender==2)then
 print*,''
 print*,'The genetic lengths obtained with the counts are approximative.'
 print*,'They are not corrected for informativeness.'
 print*,'Better genetic maps can be estimated with the #ITERATIONS option.'
 print*,''
endif

deallocate(initrec,postrec,recounts,rectable)
close(511);close(512);close(513)

end subroutine

subroutine phase_byone(rules)
implicit none
integer ::rules,mininfo,maxdiff,i,info1
real*8 ::delta,maxdist,distortion,distor1,distor2

 limit1=ceiling(0.50*maxval(info_marker)) ! accept also markers with 90% of informative progeny 
 if(limit1 > 50)limit1=50
 if(limit1 < 2)limit1=2 ! keep unphased regions with only one informative offspring


! rules 1 : phase only reliable markers based on reliable markers - no differences
! rules 2 : phase only reliable markers based on reliable markerss - 0.999
! rules 3 : phase all based on reliable markers - 0.999
! rules 4 : phase all based on all markers - 0.999
! rules 5 : phase all based on all markers - 0.99
! rules 6 : phase all based on all markers - 0.99 (min-informative reduced)
! rules 7 : rules 6 but markers must be "balanced"
  
  distortion=0.00
  select case(rules)
    case(1)
      delta=6.906755;mininfo=2;maxdiff=0;maxdist=1000.d0
    case(2)
      delta=6.906755;mininfo=2;maxdiff=100000;maxdist=1000.d0
    case(3)
      delta=6.906755;mininfo=2;maxdiff=100000;maxdist=1000.d0
    case(4)
      delta=6.906755;mininfo=2;maxdiff=100000;maxdist=1000.d0
    case(5)
      delta=4.59511985;mininfo=2;maxdiff=100000;maxdist=1000.d0 ! P1=0.99 et P2=0.01
    case(6)
      delta=4.59511985;mininfo=1;maxdiff=100000;maxdist=1000.d0 ! P1=0.99 et P2=0.01
    case(7)
      delta=6.906755;mininfo=2;maxdiff=100000;maxdist=1000.d0;distortion=0.25
  end select

! indicate markers as phased when mendelian information available

 do k=1,nmarq
  if(tophase(k)==2)then
    phased(k)=1
    phase1(k)=hap(byparent(firsto,1),1,k)
    phase2(k)=hap(byparent(firsto,1),2,k)
   endif
 enddo


! fix a first markers in cases no phased at start

 if(sum(phased)==0)then
   do k=1,nmarq
     if(info_marker(k,2)>0 .or. info_marker(k,3)>0)distor1=1.d0*info_marker(k,2)/(1.d0*info_marker(k,2)+1.d0*info_marker(k,3))
     distor2=1.d0-distor1
     if(info_marker(k,1)>=limit1 .and. tophase(k)==1 .and. distor1 >= distortion .and. distor2 >=distortion)then
      phased(k)=1
      phase1(k)=1;phase2(k)=2
      exit
     endif
   enddo
 endif


do k=1,nmarq
 if(phased(k)==1)cycle
 if(rules<3 .and. info_marker(k,1)<limit1)cycle
 if(typ(oldparent,2*k-1)==typ(oldparent,2*k))cycle ! homozygote
 if(hap(oldparent,1,k)/=0)cycle ! already phased
 if(info_marker(k,2)>0 .or. info_marker(k,3)>0)distor1=1.d0*info_marker(k,2)/(1.d0*info_marker(k,2)+1.d0*info_marker(k,3))
 distor2=1.d0-distor1
 if(distor1 < distortion .or. distor2 < distortion)cycle

  phase1(k)=1;phase2(k)=2 ! start with random values

  numoff=0;ninfo=0;outphase=0;inphase=0;ph1=0.d0;ph2=0.d0

  do offspring=firsto,lasto ! cylce around offspring
   numoff=numoff+1
!   if(hap(byparent(offspring,2),parsex,k)==0)cycle
   if(hap2(numoff,k)==0)cycle


! search flanking phased markers
   m2=0;info1=0
   do l=k+1,nmarq
     if(phased(l)==0)cycle ! should skip homozygotes
     if(rules<4 .and. info_marker(l,1)<limit1)cycle
!     if(hap(byparent(offspring,2),parsex,l)==0)cycle
     if(hap2(numoff,l)==0)cycle
     m2=l
     exit
   enddo

   if(m2/=0)then
    m1=k 
    info1=1;ninfo=ninfo+1
!    if(hap(byparent(offspring,2),parsex,k)==1 .and. hap(byparent(offspring,2),parsex,m2)==phase1(m2))then
    if(hap2(numoff,k)==1 .and. hap2(numoff,m2)==phase1(m2))then
       inphase=inphase+1
       ph1=ph1+log(1.d0-probrec(m1,m2,parsex))
       ph2=ph2+log(probrec(m1,m2,parsex))
!    else if(hap(byparent(offspring,2),parsex,k)==1 .and. hap(byparent(offspring,2),parsex,m2)==phase2(m2))then
    else if(hap2(numoff,k)==1 .and. hap2(numoff,m2)==phase2(m2))then
       outphase=outphase+1
       ph2=ph2+log(1.d0-probrec(m1,m2,parsex))
       ph1=ph1+log(probrec(m1,m2,parsex))
!    else if(hap(byparent(offspring,2),parsex,k)==2 .and. hap(byparent(offspring,2),parsex,m2)==phase1(m2))then
    else if(hap2(numoff,k)==2 .and. hap2(numoff,m2)==phase1(m2))then
       outphase=outphase+1
       ph2=ph2+log(1.d0-probrec(m1,m2,parsex))
       ph1=ph1+log(probrec(m1,m2,parsex))
!    else if(hap(byparent(offspring,2),parsex,k)==2 .and. hap(byparent(offspring,2),parsex,m2)==phase2(m2))then
    else if(hap2(numoff,k)==2 .and. hap2(numoff,m2)==phase2(m2))then
       inphase=inphase+1
       ph1=ph1+log(1.d0-probrec(m1,m2,parsex))
       ph2=ph2+log(probrec(m1,m2,parsex))
    endif
   endif

! search flanking phased markers
   m2=0
   do l=k-1,1,-1
     if(phased(l)==0)cycle
     if(rules<4 .and. info_marker(l,1)<limit1)cycle
!     if(hap(byparent(offspring,2),parsex,l)==0)cycle
     if(hap2(numoff,l)==0)cycle
     m2=l
     exit
   enddo

   if(m2/=0)then
    m1=k 
    if(info1==0)ninfo=ninfo+1
!    if(hap(byparent(offspring,2),parsex,k)==1 .and. hap(byparent(offspring,2),parsex,m2)==phase1(m2))then
    if(hap2(numoff,k)==1 .and. hap2(numoff,m2)==phase1(m2))then
       inphase=inphase+1
       ph1=ph1+log(1.d0-probrec(m2,m1,parsex))
       ph2=ph2+log(probrec(m2,m1,parsex))
!    else if(hap(byparent(offspring,2),parsex,k)==1 .and. hap(byparent(offspring,2),parsex,m2)==phase2(m2))then
    else if(hap2(numoff,k)==1 .and. hap2(numoff,m2)==phase2(m2))then
       outphase=outphase+1
       ph2=ph2+log(1.d0-probrec(m2,m1,parsex))
       ph1=ph1+log(probrec(m2,m1,parsex))
!    else if(hap(byparent(offspring,2),parsex,k)==2 .and. hap(byparent(offspring,2),parsex,m2)==phase1(m2))then
    else if(hap2(numoff,k)==2 .and. hap2(numoff,m2)==phase1(m2))then
       outphase=outphase+1
       ph2=ph2+log(1.d0-probrec(m2,m1,parsex))
       ph1=ph1+log(probrec(m2,m1,parsex))
!    else if(hap(byparent(offspring,2),parsex,k)==2 .and. hap(byparent(offspring,2),parsex,m2)==phase2(m2))then
    else if(hap2(numoff,k)==2 .and. hap2(numoff,m2)==phase2(m2))then
       inphase=inphase+1
       ph1=ph1+log(1.d0-probrec(m2,m1,parsex))
       ph2=ph2+log(probrec(m2,m1,parsex))
    endif
   endif

  enddo ! offspring


!  if(oldid(oldparent)==821)print*,k,inphase,outphase

  if((ph1-ph2)>delta .and. ninfo>=mininfo .and. inphase>=mininfo .and. outphase<=maxdiff)then ! P1=0.999 et P2=0.001
      phased(k)=1 ! don't switch phase1
      n_phased=n_phased+1
      phase1(k)=1;phase2(k)=2
  else if((ph2-ph1)>delta .and. ninfo>=mininfo .and. outphase>=mininfo .and. inphase<=maxdiff)then
      phased(k)=1 ! switch phase1
      n_phased=n_phased+1
      phase1(k)=2;phase2(k)=1
  endif

enddo ! marker
 
end subroutine
subroutine printrec
implicit none
integer ::ninfor

 do l=1,noffspring

   mate=0
   if(parsex==1 .and. dam(byparent(firsto+l-1,2))/=0)then
     if(genotyped(dam(byparent(firsto+l-1,2))))mate=1
   endif
   if(parsex==2 .and. sire(byparent(firsto+l-1,2))/=0)then
     if(genotyped(sire(byparent(firsto+l-1,2))))mate=1
   endif

   ngam(0)=ngam(0)+1.00;ngam(parsex)=ngam(parsex)+1.00
   ninfor=0
   lastknown=0;lastori=0;nrec=0
   write(300,'(i7,1x,i7)',advance='no')oldid(byparent(firsto+l-1,2)),oldid(oldparent)
   norigins1=norigins1+1
   do k=1,nmarq
     ori1=0
     if(hap(oldparent,1,k)/=hap(oldparent,2,k) .and. hap(oldparent,1,k)/=0 .and. hap(oldparent,2,k)/=0)then
       if(hap(oldparent,1,k)==hap(byparent(firsto+l-1,2),parsex,k))ori1=1
       if(hap(oldparent,2,k)==hap(byparent(firsto+l-1,2),parsex,k))ori1=2
       if(ori1>0)ninfor=ninfor+1
     endif
     if(k<nmarq)write(300,'(1x,i1)',advance='no')ori1
     if(k==nmarq)write(300,'(1x,i1)',advance='yes')ori1
     if(ori1/=0)then
      if(lastknown/=0)then
       if(ori1/=lastori)then
         write(301,'(4(1x,i7))')oldid(byparent(firsto+l-1,2)),oldid(oldparent),lastknown,k
         nrec(l)=nrec(l)+1
         nprinted=nprinted+1
       endif
      endif
      lastknown=k;lastori=ori1
     endif
   enddo
   k=0
   if(prephased)k=1
   write(302,'(2(i7,1x),i1,1x,i4,2(1x,i1),3(1x,i7),1x,i4)')&
     oldid(byparent(firsto+l-1,2)),oldid(oldparent),parsex,noffspring,k,mate,nheter,nhomoz,ninfor,nrec(l)
 enddo

end subroutine


end program

