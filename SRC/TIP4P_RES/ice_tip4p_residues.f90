program ice_tip4p_residues

implicit none

integer :: nargc, iargc, nat, nmol, no, i, j, k, l, mind_IDX, c_o, c_v, c_h1, c_h2
integer,allocatable :: resnum(:), idx(:), rn(:,:)
integer, parameter :: cart=3, zeroi=0, one=1, two=2
double precision :: xdf, ydf, zdf, rsqdf, mind, icell(cart*cart)
double precision, allocatable :: pos(:,:)
character*2, allocatable :: sym(:)
character*5, allocatable :: resname(:), aname(:)
character*256 :: buffer, infile, line, boxfile, command


nargc=iargc()

if (nargc.lt.1) then
   write(*,*) "Usage:"
   write(*,*) "ice_tip4p_residues.f90 [.gro file]"
   stop
else
   call getarg(1,buffer)
   read(buffer,*) infile
endif

! open .gro file
nmol=0
open(unit=101, file=trim(infile), status='old')
read(101,*)
read(101,*) nat
allocate(pos(nat,cart),sym(nat),resnum(nat),resname(nat),aname(nat),idx(nat))
do i=1,nat
   read(101,'(i5,2a5,i6,3f8.3)') resnum(i), resname(i), aname(i), idx(i), (pos(i,j), j=1,3)
   if (trim(adjustl(aname(i))).eq.'O') then
      nmol=nmol+1
   endif
enddo
icell(:)=0.0d0
read(101,'(3f10.5)') icell(1), icell(5), icell(9)
close(101)

allocate(rn(nmol,4))

open(unit=102, file='tip4p.gro', status='unknown')
write(102,*) "Carefully crafted by hand"
write(102,*) nmol*4

! 
c_o=0
do i=1,nat
   if (trim(adjustl(aname(i))).eq.'O') then
      c_o=c_o+1
      rn(c_o,1)=i

      mind=1.0d30
      do j=1,nat
         if (trim(adjustl(aname(j))).eq.'X') then
            xdf=pos(j,1)-pos(i,1)
            ydf=pos(j,2)-pos(i,2)
            zdf=pos(j,3)-pos(i,3)
            call images(cart,0,1,1,icell,xdf,ydf,zdf)
            rsqdf=dsqrt(xdf**2.0d0+ydf**2.0d0+zdf**2.0d0)
            if (rsqdf.lt.mind) then
               mind=rsqdf
               mind_IDX=j
            endif
         endif
      enddo
      !write(*,*) mind, i, mind_IDX
      rn(c_o,2)=mind_IDX

      mind=1.0d30
      do j=1,nat
         if (trim(adjustl(aname(j))).eq.'H') then
            xdf=pos(j,1)-pos(i,1)
            ydf=pos(j,2)-pos(i,2)
            zdf=pos(j,3)-pos(i,3)
            call images(cart,0,1,1,icell,xdf,ydf,zdf)
            rsqdf=dsqrt(xdf**2.0d0+ydf**2.0d0+zdf**2.0d0)
            if (rsqdf.lt.mind) then
               mind=rsqdf
               mind_IDX=j
            endif
         endif
      enddo
      !write(*,*) mind, i, mind_IDX
      rn(c_o,3)=mind_IDX

      mind=1.0d30
      do j=1,nat
         if (trim(adjustl(aname(j))).eq.'H'.and.j.ne.rn(c_o,3)) then
            xdf=pos(j,1)-pos(i,1)
            ydf=pos(j,2)-pos(i,2)
            zdf=pos(j,3)-pos(i,3)
            call images(cart,0,1,1,icell,xdf,ydf,zdf)
            rsqdf=dsqrt(xdf**2.0d0+ydf**2.0d0+zdf**2.0d0)
            if (rsqdf.lt.mind) then
               mind=rsqdf
               mind_IDX=j
            endif
         endif
      enddo
      !write(*,*) mind, i, mind_IDX
      rn(c_o,4)=mind_IDX

      write(102,'(i5,2a5,i5,3f8.3)')   c_o, "SOL", "OW" , (c_o*4)-3, pos(rn(c_o,1),:) 
      write(102,'(i5,2a5,i5,3f8.3)')   c_o, "SOL", "HW1", (c_o*4)-1, pos(rn(c_o,3),:)
      write(102,'(i5,2a5,i5,3f8.3)')   c_o, "SOL", "HW2" , c_o*4,    pos(rn(c_o,4),:)
      write(102,'(i5,2a5,i5,3f8.3)')   c_o, "SOL", "MW",  (c_o*4)-2, pos(rn(c_o,2),:) 

   endif
enddo

write(102,'(3f10.5)') icell(1), icell(5), icell(9)

end program ice_tip4p_residues
