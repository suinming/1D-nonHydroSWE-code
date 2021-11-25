    program case2
      implicit none
      integer i,j,m,step,step_n,total_steps,total_nodes,n_near,n_correct,i_correct
      real*8, allocatable :: x(:),h(:,:),u(:,:),w(:,:),h_x(:,:),u_x(:,:),uh_x(:,:),x_c(:),q(:),q_x(:)
      real*8, allocatable :: weighting(:,:),lambda(:,:,:),zb(:),u_bc(:),h_bc(:),fig1(:,:)
      integer, allocatable :: n_loc(:),loc_index(:,:),kin(:)
      real*8 dt,g,t,ub,hb
      character*22 filename
      parameter(g=9.81d0)
      m=3
      n_near=11
	  n_correct=5
      open(1,file='input.prn')
      read(1,*) total_nodes,total_steps,dt
      close(1)
      allocate(kin(1:total_nodes))
      allocate(x(1:total_nodes))
      allocate(x_c(1:total_nodes))
      allocate(h(1:total_nodes,0:2))
      allocate(u(1:total_nodes,0:2))
      allocate(w(1:total_nodes,0:2))
      allocate(h_x(1:total_nodes,0:2))
      allocate(u_x(1:total_nodes,0:2))
      allocate(uh_x(1:total_nodes,0:2))
      allocate(zb(1:total_nodes))
	  allocate(q(1:total_nodes))
	  allocate(q_x(1:total_nodes))
      allocate(weighting(1:total_nodes,1:n_near*5))
      allocate(lambda(1:total_nodes,1:m,1:n_near*5))
      allocate(n_loc(1:total_nodes))
      allocate(loc_index(1:total_nodes,1:n_near*5))
      allocate(u_bc(0:total_steps))
      allocate(h_bc(0:total_steps))
      allocate(fig1(1:total_steps, 1:2))
      open(1,file='bathymetry.prn')
      do j=1,total_nodes
        read(1,*) kin(j)  
      end do
      rewind(1)
      do j=1,total_nodes
        read(1,*) kin(j),x(j),zb(j)
      end do
      close(1)

      open(1,file='wmbc.prn')
      do step=0,total_steps
        read(1,*) u_bc(step)  
      end do
      close(1)

      open(1,file='IC.prn')
      step=0  
      do j=1,total_nodes
        read(1,*) h(j,step),u(j,step)
        h(j,step)=h(j,step)-zb(j)
        if(h(j,step) .lt. 0) then
          h(j,step)=0
        end if
      end do
      close(1)

      call pre_MLS(total_nodes,m,n_near,x,weighting,lambda,n_loc,loc_index,g,h,dt,kin) 
      call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),u(:,step),w(:,step) &
              ,h_x(:,step),u_x(:,step),uh_x(:,step))    

!      write(filename,'(i5.5,a4)') step,'.dat'  
!      Open(2,file=filename)
!	  do j=1,total_nodes
!        write(2,'(i6,11f16.5)') kin(j),x(j),h(j,step)+zb(j),u(j,step),w(j,step)
!	  end do
!      close(2)

      step=1  
      dt=dt/2
      hb=(h_bc(0)+h_bc(1))/2
      ub=(u_bc(0)+u_bc(1))/2
      call predict(total_nodes,step,dt,g,h,u,h_x,u_x,uh_x)
	  call NF_BC(total_nodes,step,u,kin)
	  call WM_BC(total_nodes,step,h,u,kin,hb,ub)
      call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),u(:,step),w(:,step) &
              ,h_x(:,step),u_x(:,step),uh_x(:,step))     
	  do i_correct=1,n_correct
	    call calculate_q(total_nodes,step,dt,h,w,q)
	    call MLS_q(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,q,q_x)     
        call correct(total_nodes,step,dt,g,h,u,h_x,u_x,uh_x,q,q_x,n_correct,i_correct)
	    call NF_BC(total_nodes,step,u,kin)
	    call WM_BC(total_nodes,step,h,u,kin,hb,ub)
        call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),u(:,step),w(:,step) &
                ,h_x(:,step),u_x(:,step),uh_x(:,step))     
	  end do
      step=2  
      hb=h_bc(1)
	  ub=u_bc(1)
      call predict(total_nodes,step,dt,g,h,u,h_x,u_x,uh_x)
	  call NF_BC(total_nodes,step,u,kin)
	  call WM_BC(total_nodes,step,h,u,kin,hb,ub)
      call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),u(:,step),w(:,step) &
              ,h_x(:,step),u_x(:,step),uh_x(:,step))     
	  do i_correct=1,n_correct
	    call calculate_q(total_nodes,step,dt,h,w,q)
	    call MLS_q(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,q,q_x)   
        call correct(total_nodes,step,dt,g,h,u,h_x,u_x,uh_x,q,q_x,n_correct,i_correct)
	    call NF_BC(total_nodes,step,u,kin)
	    call WM_BC(total_nodes,step,h,u,kin,hb,ub)
        call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),u(:,step),w(:,step) &
                ,h_x(:,step),u_x(:,step),uh_x(:,step))  
	  end do
      do j=1,total_nodes
        h(j,step-1)=h(j,step)
        u(j,step-1)=u(j,step)
		w(j,step-1)=w(j,step)
        h_x(j,step-1)=h_x(j,step)
        u_x(j,step-1)=u_x(j,step)
        uh_x(j,step-1)=uh_x(j,step)
      end do
      step=1  
!      write(filename,'(i5.5,a4)') step,'.dat'  
!      Open(2,file=filename)
!	  do j=1,total_nodes
!        write(2,'(i6,11f16.5)') kin(j),x(j),h(j,step)+zb(j)
!	  end do
!      close(2)
      fig1(1, 1) = h(50,step) + zb(50)
      fig1(1, 2) = h(88,step) + zb(88)
      dt=dt*2
      step=2 
      do step_n=2,total_steps  
        hb=h_bc(step_n)
	    ub=u_bc(step_n)
!        write(*,'(a5,i5)') 'step=',step_n
        call predict(total_nodes,step,dt,g,h,u,h_x,u_x,uh_x)
	    call NF_BC(total_nodes,step,u,kin)
	    call WM_BC(total_nodes,step,h,u,kin,hb,ub)
        call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),u(:,step),w(:,step) &
                ,h_x(:,step),u_x(:,step),uh_x(:,step))     
		do i_correct=1,n_correct
	      call calculate_q(total_nodes,step,dt,h,w,q)
	      call MLS_q(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,q,q_x)   
          call correct(total_nodes,step,dt,g,h,u,h_x,u_x,uh_x,q,q_x,n_correct,i_correct)
	      call NF_BC(total_nodes,step,u,kin)
	      call WM_BC(total_nodes,step,h,u,kin,hb,ub)
          call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),u(:,step),w(:,step) &
                  ,h_x(:,step),u_x(:,step),uh_x(:,step))    
		end do

2       do j=1,total_nodes
          h(j,0)=h(j,1)
          u(j,0)=u(j,1)
          w(j,0)=w(j,1)
          h_x(j,0)=h_x(j,1)
          u_x(j,0)=u_x(j,1)
          uh_x(j,0)=uh_x(j,1)
          h(j,1)=h(j,2)
          u(j,1)=u(j,2)
          w(j,1)=w(j,2)
          h_x(j,1)=h_x(j,2)
          u_x(j,1)=u_x(j,2)
          uh_x(j,1)=uh_x(j,2)
        end do
      end do

    end



    subroutine pre_MLS(N,m,n_near,x,wt,lambda,n_loc,loc_index,g,h,dt,kin)
      implicit none
      integer i,j,k,l,m,N,n_near,n_loc(1:N),loc_index(1:N,1:n_near*5),kin(1:N)
      real*8 x(1:N),wt(1:N,1:n_near*5),lambda(1:N,1:m,1:n_near*5),d_near(1:n_near),p(1:m)
      real*8 r,x_loc,ro,g,rn,dt,epsilon,h(1:N)
      real*8 a(1:n_near*5,1:m),at(1:m,1:n_near*5),ata(1:m,1:m),ata_inv(1:m,1:m),aa(1:10),xx(1:10)
	  aa(1)=-13.287599973968
	  aa(2)=257.139674858618
	  aa(3)=45.8730830833093
	  aa(4)=-1152.47047375008
	  aa(5)=-641.980091013938
	  aa(6)=-90.7921629936896
	  aa(7)=2234.95009387434
	  aa(8)=1416.73353942545
	  aa(9)=759.11459081992
	  aa(10)=40.4045499548725
  
      do j=1,N
  
        call KNN(N,x,d_near,n_near,j)
        ro=d_near(n_near)*1.01
		if(kin(j) .eq. 0) then
		  rn=(d_near(2)+d_near(3))/2
		else
		  rn=d_near(2)
		end if
		xx(1)=1
		xx(2)=rn/h(j)
		xx(3)=(g*h(j))**.5/(rn/dt)
		xx(4)=xx(2)**2
		xx(5)=xx(2)*xx(3)
		xx(6)=xx(3)**2
		xx(7)=xx(2)**3
		xx(8)=xx(2)**2*xx(3)
		xx(9)=xx(2)*xx(3)**2
		xx(10)=xx(3)**3
		epsilon=.0d0
		do i=1,10
		  epsilon=epsilon+aa(i)*xx(i)
		end do
		epsilon=max(epsilon,0.d0)
  
        n_loc(j)=0
        do l=1,N 
          x_loc=x(l)-x(j) 
          r=abs(x_loc)  
          if(r .ge. ro) goto 1  
          n_loc(j)=n_loc(j)+1
          loc_index(j,n_loc(j))=l 
          wt(j,n_loc(j))=((1-r/ro)**epsilon)**.5 
1       end do
      end do

      do j=1,N
 
        do k=1,n_loc(j)
          x_loc=x(loc_index(j,k))-x(j)  
          call local_poly(x_loc,p,m) Local Polynomial
          do i=1,m
            a(k,i)=wt(j,k)*p(i)
          end do
        end do
  ! A inverse
        do k=1,n_loc(j)
          do i=1,m
            at(i,k)=a(k,i)
          end do
        end do
  ! Psuudo matrix
        do i=1,m
          do l=1,m
            ata(i,l)=0
            do k=1,n_loc(j)
              ata(i,l)=ata(i,l)+at(i,k)*a(k,l)    !atxa=at*a
            end do
          end do
        end do
  ! Psuudo matrix inverse
        call matrix_inv(m,ata,ata_inv)
  ! lambda matrix
        do i=1,m
          do k=1,n_loc(j)
            lambda(j,i,k)=0
            do l=1,m
              lambda(j,i,k)=lambda(j,i,k)+ata_inv(i,l)*at(l,k)
            end do
          end do
        end do
      end do
    end




    subroutine KNN(N,x,d_near,n_near,star) ! K Nearest Neighbors
      implicit none
      integer N,i,j,k,n_near,star
      real*8 x(1:N),r,d_near(1:n_near)
  
      do i=1,n_near
        d_near(i)=1.d99
      end do
  
      do j=1,N
        r=abs(x(j)-x(star))
        do i=1,n_near
          if(r .lt. d_near(i)) then
            do k=n_near-1,i,-1
              d_near(k+1)=d_near(k)
            end do
            d_near(i)=r
            goto 1
          end if
        end do
1     end do
    end




    subroutine local_poly(x,p,m)
      implicit none
      integer m
      real*8 x,p(1:m)
      p(1)=1.0d0
      p(2)=x
      p(3)=x**2
!      p(4)=x**3
    end




    subroutine matrix_inv(m,a,a_inv)
      implicit none
      integer i,j,k,m
      real*8 a(1:m,1:m),a_inv(1:m,1:m),z
      do i=1,m
        do j=1,m
          if(i .eq. j) then
            a_inv(i,j)=1.d0
          else
            a_inv(i,j)=0
          end if
        end do
      end do
      do i=1,m-1
        z=a(i,i)
        do j=i,m
          a(i,j)=a(i,j)/z
        end do
        do j=1,m
          a_inv(i,j)=a_inv(i,j)/z
        end do
        do k=i+1,m
          z=a(k,i)
          do j=1,m
            a(k,j)=a(k,j)-z*a(i,j)
            a_inv(k,j)=a_inv(k,j)-z*a_inv(i,j)
          end do
        end do
      end do
      do j=1,m
        a_inv(m,j)=a_inv(m,j)/a(m,m)
      end do
      do i=m-1,1,-1
        do j=m,i+1,-1
          do k=1,m
            a_inv(i,k)=a_inv(i,k)-a_inv(j,k)*a(i,j)
          end do
        end do
        do k=1,m
          a_inv(i,k)=a_inv(i,k)/a(i,i)
        end do
      end do
    end




    subroutine MLS(N,x,m,n_near,wt,lambda,n_loc,loc_index,h,u,w,h_x,u_x,uh_x)
      implicit none
      integer i,j,k,l,m,N,n_near,n_loc(1:N),loc_index(1:N,1:n_near*5)
      real*8 x(1:N),h(1:N),u(1:N),w(1:N),new_h(1:N),new_u(1:N),smooth_w(1:N),aa &
            ,h_x(1:N),u_x(1:N),uh_x(1:N),lambda(1:N,1:m,1:n_near*5),wt(1:N,1:n_near*5)
      real*8 alpha(1:m),beta(1:n_near*5)
  
    ! h
      do j=1,N
      ! beta_h
        do k=1,n_loc(j)
          beta(k)=wt(j,k)*h(loc_index(j,k))
        end do
      ! alpha_h
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      ! ouput h,h_x
        new_h(j)=alpha(1)
        h_x(j)=alpha(2)
      end do
    ! u
      do j=1,N
      ! beta_u
        do k=1,n_loc(j)
          beta(k)=wt(j,k)*u(loc_index(j,k))
        end do
      ! alpha_u
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      ! u,u_x
        new_u(j)=alpha(1)
        u_x(j)=alpha(2)
    ! uh
      ! beta_uh
        do k=1,n_loc(j)
          beta(k)=wt(j,k)*u(loc_index(j,k))*h(loc_index(j,k))
        end do
      ! alpha_uh
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      ! uh_x
        uh_x(j)=alpha(2)
      end do
  ! substitude
      do j=1,N
        h(j)=new_h(j)
        u(j)=new_u(j)
		w(j)=.5d0*(u(j)*h_x(j)-uh_x(j))
      end do
    end





    subroutine predict(total_nodes,step,dt,g,h,u,h_x,u_x,uh_x)
      implicit none
      integer j,total_nodes,step
      real*8 h(1:total_nodes,0:2),u(1:total_nodes,0:2) &
            ,h_x(1:total_nodes,0:2),u_x(1:total_nodes,0:2),uh_x(1:total_nodes,0:2) &
            ,dt,g
      do j=1,total_nodes
          h(j,step)=h(j,step-1)-dt*uh_x(j,step-1)
          u(j,step)=u(j,step-1)-dt*(g*h_x(j,step-1)+u(j,step-1)*u_x(j,step-1))
      end do
    end




    subroutine correct(total_nodes,step,dt,g,h,u,h_x,u_x,uh_x,q,q_x,n_correct,i_correct)
      implicit none
      integer j,total_nodes,step,n_correct,i_correct
      real*8 h(1:total_nodes,0:2),u(1:total_nodes,0:2),q(1:total_nodes),q_x(1:total_nodes) &
            ,h_x(1:total_nodes,0:2),u_x(1:total_nodes,0:2),uh_x(1:total_nodes,0:2) &
            ,dt,g,s1,s2,s3,s4,rr
	  rr=1.d0*i_correct/n_correct
      do j=1,total_nodes
        h(j,step)=h(j,step-1)-dt*(uh_x(j,step-1)+uh_x(j,step))/2
		s1=g*(h_x(j,step-1)+h_x(j,step))/2
		s2=(u(j,step-1)*u_x(j,step-1)+u(j,step)*u_x(j,step))/2
        s3=.5d0*q_x(j)*(1-abs(cos(3.14159265358979*rr/2))**2.d0)
		s4=.5d0*q(j)*(h_x(j,step-1)+h_x(j,step))/(h(j,step-1)+h(j,step))*(1-abs(cos(3.14159265358979*rr/2))**2.d0)
        u(j,step)=u(j,step-1)-dt*(s1+s2+s3+s4)
      end do
    end




    subroutine NF_BC(total_nodes,step,u,kin)
      implicit none
      integer j,total_nodes,step,kin(1:total_nodes)
      real*8 h(1:total_nodes,0:2),u(1:total_nodes,0:2)
      do j=1,total_nodes
        if(kin(j) .eq. 1) then
          u(j,step)=0
        end if
      end do
    end




    subroutine calculate_q(total_nodes,step,dt,h,w,q)
      implicit none
      integer j,total_nodes,step
      real*8 h(1:total_nodes,0:2),w(1:total_nodes,0:2),q(1:total_nodes),dt,w_t
      do j=1,total_nodes
	    w_t=(w(j,step)-w(j,step-1))/dt
        q(j)=(h(j,step)+h(j,step-1))/2*w_t
      end do
    end




    subroutine MLS_q(N,x,m,n_near,wt,lambda,n_loc,loc_index,q,q_x)
      implicit none
      integer i,j,k,l,m,N,n_near,n_loc(1:N),loc_index(1:N,1:n_near*5),n_smooth
      real*8 x(1:N),q(1:N),q_x(1:N),lambda(1:N,1:m,1:n_near*5),wt(1:N,1:n_near*5)
      real*8 alpha(1:m),beta(1:n_near*5),smooth_q(1:N),aa
	  n_smooth=1
  ! smooth the q
    do i=1,n_smooth
      do j=1,N
	    aa=0
		smooth_q(j)=0
        do k=1,n_loc(j)
		  aa=aa+wt(j,k)**2
          smooth_q(j)=smooth_q(j)+wt(j,k)**2*q(loc_index(j,k))
        end do
        smooth_q(j)=smooth_q(j)/aa
	  end do
      do j=1,N
        q(j)=smooth_q(j)
	  end do
	end do
  
    ! q
      do j=1,N
      ! beta_h
        do k=1,n_loc(j)
          beta(k)=wt(j,k)*q(loc_index(j,k))
        end do
      ! alpha_h
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      ! q_x
        q_x(j)=alpha(2)
      end do
    end




    subroutine WM_BC(total_nodes,step,h,u,kin,hb,ub)
      implicit none
      integer j,total_nodes,step,kin(1:total_nodes)
      real*8 h(1:total_nodes,0:2),u(1:total_nodes,0:2),hb,ub
      do j=1,total_nodes
        if(kin(j) .eq. 2) then
!		  h(j,step)=hb
          u(j,step)=ub
        end if
      end do
    end


