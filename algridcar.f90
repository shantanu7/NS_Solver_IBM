program algridcar

 real, dimension(:), allocatable :: x, y, z
 integer mx, my, mz
 integer ndim
  
 integer ncoord, sti, edi, cdir, opt, cont
 real stx, edx, dx, dsmin

 
 cont = 1
 open(11,file='gscript.dat')
 read(11,*) ndim, mx, my, mz
 
 write(*,*) ndim
 write(*,*) mx, my, mz
 
 allocate(x(mx),y(my),z(mz))
 
 x = 0.0
 y = 0.0
 z = 0.0
 
 
 do while(cont .eq. 1)
 read(11,*) ncoord, stx, edx, dx, sti, edi, cdir, opt, cont
 
 select case ( ncoord )
 
 case (1)

 if(opt .eq. 2) then
  
  dx = (edx-stx)/(edi-sti)
  do i=sti,edi
  x(i) = stx + (i-sti)*dx
  enddo

 else 
 
  if(opt .eq. 0) then
   dsmin = dx
  elseif( opt .eq. 1 ) then
   dsmin = abs(x(edi)-x(edi+1))
  elseif( opt .eq. -1 ) then
   dsmin = abs(x(sti)-x(sti-1))
  endif
  
  call ali(stx,edx,dsmin,cdir,sti,edi,x,mx)
 
 endif
 
 case (2)

 if(opt .eq. 2) then
  
  dx = (edx-stx)/(edi-sti)
  do i=sti,edi
  y(i) = stx + (i-sti)*dx
  enddo

 else 
 
  if(opt .eq. 0) then
   dsmin = dx
  elseif( opt .eq. 1 ) then
   dsmin = abs(y(edi)-y(edi+1))
  elseif( opt .eq. -1 ) then
   dsmin = abs(y(sti)-y(sti-1))
  endif
  
  call ali(stx,edx,dsmin,cdir,sti,edi,y,my)
 
 endif

 case (3)

 if(opt .eq. 2) then
  
  dx = (edx-stx)/(edi-sti)
  do i=sti,edi
  z(i) = stx + (i-sti)*dx
  enddo

 else 
 
  if(opt .eq. 0) then
   dsmin = dx
  elseif( opt .eq. 1 ) then
   dsmin = abs(z(edi)-z(edi+1))
  elseif( opt .eq. -1 ) then
   dsmin = abs(z(sti)-z(sti-1))
  endif
  
  call ali(stx,edx,dsmin,cdir,sti,edi,z,mz)
 
 endif 


 end select
 
 enddo


 write(*,*) x(1), y(1), z(1)

 
 select case ( ndim )
 
 case (2)
 
  open (unit=8, file='grid.dat')
  write (8,*) 'variables = x, y'
  write (8,*) 'zone i =',mx,' j =',my,' f = block'
  write (8,'(10(E13.6,2X))') ((x(i),i=1,mx),j=1,my)
  write (8,'(10(E13.6,2X))') ((y(j),i=1,mx),j=1,my)
  close (8)

  open (unit=9, file='Agrid.dat')
  write (9,*) '1'
  write (9,*) mx, my, mz
  write (9,*) (x(i),i=1,mx)
  write (9,*) (y(j),j=1,my)
  close (9)
  
  open(unit=11, file='xgrid.dat')
  do i=1,mx
  write(11,*) i, x(i)
  enddo
  close(11)
  
  open(unit=12, file='ygrid.dat')
  do j=1,my
  write(12,*) j, y(j)
  enddo
  close(12)
 
 case (3)
 
  open (unit=8, file='grid_xy.dat')
  write (8,*) 'variables = x, y'
  write (8,*) 'zone i =',mx,' j =',my,' f = block'
  write (8,*) ((x(i),i=1,mx),j=1,my)
  write (8,*) ((y(j),i=1,mx),j=1,my)
  close (8)

  open (unit=8, file='grid_xz.dat')
  write (8,*) 'variables = x, z'
  write (8,*) 'zone i =',mx,' j =',mz,' f = block'
  write (8,*) ((x(i),i=1,mx),k=1,mz)
  write (8,*) ((z(k),i=1,mx),k=1,mz)
  close (8)
  
  
  open (unit=9, file='Agrid.dat')
  write (9,*) '1'
  write (9,*) mx, my, mz
  write (9,*) (x(i),i=1,mx)
  write (9,*) (y(j),j=1,my)
  write (9,*) (z(k),k=1,mz)
  close (9)
  
  open(unit=11, file='xgrid.dat')
  do i=1,mx
  write(11,*) i, x(i)
  enddo
  close(11)
  
  open(unit=12, file='ygrid.dat')
  do j=1,my
  write(12,*) j, y(j)
  enddo
  close(12)

  open(unit=12, file='zgrid.dat')
  do k=1,mz
  write(12,*) k, z(k)
  enddo
  close(12)
  
end select

 deallocate(x,y,z)

end program
  

  subroutine ali(stx,edx,dsmin,cdir,sti,edi,x,mx)
 
  real stx,edx,dsmin
  integer cdir,sti,edi,i,n,mx
  real dx,xx,a,filter,root,length
  real x(mx)

  n=(edi-sti)+1
  dx=1.0/(n-1)
  
  length = abs(edx-stx)
  ddsmin = dsmin/length
  
  a=root(ddsmin,dx)
  do i=1,n
  if (cdir.eq.1) then
  xx=dx*(n-i)
  filter=1-E(xx,a)
  else
  xx=dx*(i-1)
  filter=E(xx,a)
  endif
  x(sti+i-1)=filter*(edx-stx)+stx
  enddo
  return
  end


  real function root(dsmin,dx)
  real dsmin,dx
  real Es, a, err 
  real a0, df, f
  Es=0.0001
  a=2.
  err=1.
  do while(err.gt.Es)
  a0=a;
  f=tanh(a)*(dsmin-1.)+tanh(a*(1.-dx))
  df=(dsmin-1.)/(cosh(a)**2)+(1.-dx)/(cosh(a*(1.-dx))**2)
  a=a-f/df
  err=abs(a-a0)
  enddo
  root=a
  return
  end

  real function E(x,a)
  real x,a
  E=1.-tanh(a*(1.-x))/tanh(a)
  return
  end



