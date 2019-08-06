subroutine images (imcon,idnode,mxnode,natm,cell,xxx,yyy,zzz)

implicit real*8 (a-h,o-z)

dimension xxx(*),yyy(*),zzz(*)
dimension cell(9),rcell(9)

   call invert(cell,rcell,det)
   if(abs(det).lt.1.d-6) stop "zero determinant cell matrix"
   
   do i=idnode+1,natm,mxnode
      
      ssx=(rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i))
      ssy=(rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i))
      ssz=(rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i))
      
      xss=ssx-nint(ssx)
      yss=ssy-nint(ssy)
      zss=ssz-nint(ssz)
      
      xxx(i)=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
      yyy(i)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
      zzz(i)=(cell(3)*xss+cell(6)*yss+cell(9)*zss)
      
   enddo
   
return
end

! invert subroutine
subroutine invert(a,b,d)

real*8 a,b,d,r

dimension a(9),b(9)

b(1)=a(5)*a(9)-a(6)*a(8)
b(2)=a(3)*a(8)-a(2)*a(9)
b(3)=a(2)*a(6)-a(3)*a(5)
b(4)=a(6)*a(7)-a(4)*a(9)
b(5)=a(1)*a(9)-a(3)*a(7)
b(6)=a(3)*a(4)-a(1)*a(6)
b(7)=a(4)*a(8)-a(5)*a(7)
b(8)=a(2)*a(7)-a(1)*a(8)
b(9)=a(1)*a(5)-a(2)*a(4)

d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
r=0.d0
if(abs(d).gt.0.d0)r=1.d0/d

b(1)=r*b(1)
b(2)=r*b(2)
b(3)=r*b(3)
b(4)=r*b(4)
b(5)=r*b(5)
b(6)=r*b(6)
b(7)=r*b(7)
b(8)=r*b(8)
b(9)=r*b(9)
return
end
