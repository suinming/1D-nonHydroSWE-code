    program SWE1D_wetDry_nonhydro
      implicit none
      integer i,j,m,step,step_n,total_steps,total_nodes,n_near,n_flux_BC,n_eta_BC,points,n_step_save
      integer count
      real*8 x_point(1:6)
      real*8, allocatable :: x(:),h(:,:),uh(:,:),h_x(:,:),uh_x(:,:),x_c(:) &
                            ,u2h_x(:,:),n_x(:),n_y(:),zb(:),zb_x(:),eta(:,:),eta1(:) &
                            ,flux_x(:,:),flux_y(:,:),flux1_x(:),flux1_y(:),coriolis(:) &
                            ,uh_xx(:,:),ra(:),wavemaker_uhArr(:),wavemaker_hArr(:)
      integer, allocatable :: kin(:),cor(:),region(:),nearest_b(:),nearest_i(:)
      real*8, allocatable :: weighting(:,:),lambda(:,:,:)
      integer, allocatable :: n_loc(:),loc_index(:,:)
      real*8 dt,dx,g,Cz,h_min,t_ramp,t,u_max,volume,h0,slope
      real*8 amplitude,x0,c,wavemaker_eta,wavemaker_u,term1,term2,wavemaker_uh_bc,wavemaker_h_bc
      real*8, allocatable :: compareWithexperiment(:),animation(:)
      real*8, allocatable :: w_velocity(:,:),q(:),q_x(:)
      integer n_correct,i_correct
      integer gaugePosition
      character*22 filename
      parameter(g=9.81d0)


	  ! parameter control -------------------------

      m=3
      n_near=10
	  total_nodes = 3000     ! 1376 2902
	  total_steps = 1000
	  n_correct = 1
	  n_eta_BC = 0
	  n_flux_BC = 0
	  t_ramp = 0
	  dt = 0.0005
	  dx = 0.002
	  h0 = 0.218
	  h_min = 0.001
	  slope = 2.75
	  !Cz = 55.56
	  Cz = 0
	  u_max = 25
	  points = 0
	  n_step_save = 250
	  gaugePosition = 1635

	  ! wave maker condition
	  amplitude = h0 * 0.12
	  x0 = -1.5
	  c = (g * (h0 + amplitude))**0.5

	  !--------------------------------------------

      allocate(kin(1:total_nodes))
      allocate(cor(1:total_nodes))
      allocate(region(1:total_nodes))
      allocate(x(1:total_nodes))
      allocate(x_c(1:total_nodes))
      allocate(h(1:total_nodes,0:2))
      allocate(uh(1:total_nodes,0:2))
      allocate(h_x(1:total_nodes,0:2))
      allocate(uh_x(1:total_nodes,0:2))
      allocate(u2h_x(1:total_nodes,0:2))
      allocate(n_x(1:total_nodes))
      allocate(zb(1:total_nodes))
      allocate(zb_x(1:total_nodes))

      ! non-hydrostatic
      allocate(w_velocity(1:total_nodes,0:2))
      allocate(q(1:total_nodes))
      allocate(q_x(1:total_nodes))

      allocate(eta(1001:1000+n_eta_BC,0:total_steps))
      allocate(eta1(1001:1000+n_eta_BC))
      allocate(flux_x(5001:5000+n_flux_BC,0:total_steps))
      allocate(flux1_x(5001:5000+n_flux_BC))
      allocate(weighting(1:total_nodes,1:n_near*5))
      allocate(lambda(1:total_nodes,1:m,1:n_near*5))
      allocate(n_loc(1:total_nodes))
      allocate(loc_index(1:total_nodes,1:n_near*5))
      allocate(coriolis(1:total_nodes))
      allocate(uh_xx(1:total_nodes,0:2))
      allocate(ra(1:total_nodes))
      allocate(nearest_b(1:total_nodes))
      allocate(nearest_i(1:total_nodes))
      allocate(wavemaker_uhArr(0:total_steps))
      allocate(wavemaker_hArr(0:total_steps))
      allocate(compareWithexperiment( 1:(total_steps) ) )
      allocate(animation(1:(total_steps * total_nodes)))

      ! wave maker on the boundary (where x = 0)

      do j = 0,total_steps

        term1 = ( (3 * amplitude) / (4 * (h0**3) ) )**0.5
        term2 =  x0 + (c * j * dt)
        wavemaker_eta = amplitude * ( 1 / cosh(term1 * term2) )**2
        wavemaker_u = (c * wavemaker_eta) / (wavemaker_eta + h0)
        wavemaker_uhArr(j) = wavemaker_u * (wavemaker_eta + h0)
        wavemaker_hArr(j) = wavemaker_eta + h0

      end do

	  ! initial data

	  do j=1,total_nodes
	      cor(j) = 0
		  region(j) = 0
		  coriolis(j) = 0
		  n_x(j) = 0

		  ! x
		  if(j .eq. 1)then
            x(j) = 0
          else
            x(j) = x(j - 1) + dx
          end if

		  ! kin
		  if(j .eq. 1)then
		    kin(j) = 2
		  else if(j .eq. total_nodes-1) then
		    kin(j) = 1
		  else if(j .eq. total_nodes) then
		    kin(j) = -1
		  else
		    kin(j) = 0
		  end if

		  ! zb and h
		  if(j .le. ceiling(22.87 * h0 / dx) + 1)then  ! 8.87 22.87

            h(j,0) = h0
            zb(j) = -h(j,0)

          else if(j .gt. ceiling(22.87 * h0 / dx) + 1 .and. j .le. ceiling(25.62 * h0 / dx) + 1)then  ! 11.62 25.62

            h(j,0) = h(j - 1,0) - (dx / slope)
            !zb(j) = zb(j - 1) + (dx / slope)
            zb(j) = -h(j,0)

          else

            h(j,0) = 0
            zb(j) = zb(j - 1) + (dx / slope)

          end if

		  uh(j,0) = 0

	  end do

	  ! start calculation -----------
      call pre_MLS(total_nodes,m,n_near,x,weighting,lambda,n_loc,loc_index,kin,cor,x_c,region,ra,nearest_b,nearest_i)
      call MLS_zb(total_nodes,m,n_near,weighting,lambda,n_loc,loc_index,zb,zb_x,kin,nearest_i)

      step=0
      call ghost_BC(total_nodes,kin,h(:,step),uh(:,step),nearest_b,nearest_i,x)
      call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),uh(:,step) &
              ,h_x(:,step),uh_x(:,step),u2h_x(:,step),zb,zb_x,h_min,uh_xx(:,step),u_max,ra,kin,w_velocity(:,step))

      step=1  !�b�Ģ��ɶ��B�M�Ģ��ɶ��B�������J step=1/2�A�W�[�Ҧ�í�w
!      dt=dt/2
!      wavemaker_uh_bc = (wavemaker_uhArr(0) + wavemaker_uhArr(1)) / 2
!      wavemaker_h_bc = (wavemaker_hArr(0) + wavemaker_hArr(1)) / 2

      call predict(total_nodes,step,dt,g,zb_x,h,uh,h_x,uh_x,u2h_x,Cz,h_min,coriolis,uh_xx,u_max,kin)
!print*, 'after predict'
!print*, h(2494,step),h_x(2494,step),zb(2792)
      call NF_bc(total_nodes,h(:,step),uh(:,step),kin,wavemaker_uh_bc,wavemaker_h_bc)
      call ghost_BC(total_nodes,kin,h(:,step),uh(:,step),nearest_b,nearest_i,x)
      call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),uh(:,step) &
              ,h_x(:,step),uh_x(:,step),u2h_x(:,step),zb,zb_x,h_min,uh_xx(:,step),u_max,ra,kin,w_velocity(:,step))
!print*, 'after MLS'
!print*, h(2494,step),h_x(2494,step)
      do i_correct = 1,n_correct

          call calculate_q(total_nodes,step,dt,h,w_velocity,q)
          call MLS_q(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,q,q_x,h(:,step),h_min,zb,kin)
          call correct(i_correct,n_correct,total_nodes,step,dt,g,zb_x,h,uh,h_x,uh_x,u2h_x,Cz&
                        ,h_min,coriolis,uh_xx,u_max,kin,q,q_x)
!print*, 'after correct'
!print*, h(2494,step),h_x(2494,step)

          call NF_bc(total_nodes,h(:,step),uh(:,step),kin,wavemaker_uh_bc,wavemaker_h_bc)
          call ghost_BC(total_nodes,kin,h(:,step),uh(:,step),nearest_b,nearest_i,x)
          call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),uh(:,step)&
                  ,h_x(:,step),uh_x(:,step),u2h_x(:,step),zb,zb_x,h_min,uh_xx(:,step),u_max,ra,kin,w_velocity(:,step))


      end do


!      step=2  !���O�Ģ��ɶ��B�Astep=1
!      wavemaker_uh_bc = wavemaker_uhArr(1)
!      wavemaker_h_bc = wavemaker_hArr(1)
!      !write(*,'(a5,i5)') 'step=',step/2
!      call predict(total_nodes,step,dt,g,zb_x,h,uh,h_x,uh_x,u2h_x,Cz,h_min,coriolis,uh_xx,u_max,kin)
!!print*, uh(2487,step),uh_x(2487,step),h(2487,step),h_x(2487,step)
!      call NF_bc(total_nodes,h(:,step),uh(:,step),kin,wavemaker_uh_bc,wavemaker_h_bc)
!      call ghost_BC(total_nodes,kin,h(:,step),uh(:,step),nearest_b,nearest_i,x)
!      call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),uh(:,step) &
!              ,h_x(:,step),uh_x(:,step),u2h_x(:,step),zb,zb_x,h_min,uh_xx(:,step),u_max,ra,kin,w_velocity(:,step))
!!print*, uh(2487,step),uh_x(2487,step),h(2487,step),h_x(2487,step)
!      do i_correct = 1,n_correct
!
!          call calculate_q(total_nodes,step,dt,h,w_velocity,q)
!          call MLS_q(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,q,q_x,h(:,step),h_min,zb,kin)
!!print*, uh(2487,step),uh_x(2487,step),h(2487,step),h_x(2487,step)
!          call correct(i_correct,n_correct,total_nodes,step,dt,g,zb_x,h,uh,h_x,uh_x,u2h_x,Cz &
!                        ,h_min,coriolis,uh_xx,u_max,kin,q,q_x)
!!print*, uh(2487,step),uh_x(2487,step),h(2487,step),h_x(2487,step)
!          call NF_bc(total_nodes,h(:,step),uh(:,step),kin,wavemaker_uh_bc,wavemaker_h_bc)
!          call ghost_BC(total_nodes,kin,h(:,step),uh(:,step),nearest_b,nearest_i,x)
!          call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),uh(:,step) &
!                  ,h_x(:,step),uh_x(:,step),u2h_x(:,step),zb,zb_x,h_min,uh_xx(:,step),u_max,ra,kin,w_velocity(:,step))
!
!      end do

      ! for compare with experiment
      compareWithexperiment(1) = ( h(gaugePosition,step) + zb(gaugePosition) ) / h0

      do j=1,total_nodes

        h(j,step-1)=h(j,step)
        uh(j,step-1)=uh(j,step)
        h_x(j,step-1)=h_x(j,step)
        uh_x(j,step-1)=uh_x(j,step)
        u2h_x(j,step-1)=u2h_x(j,step)
        uh_xx(j,step-1)=uh_xx(j,step)
        w_velocity(j,step-1) = w_velocity(j,step)

        ! data for animation
        animation(j) = h(j,step) + zb(j)

      end do


      !dt=dt*2
      step=2 !�p��Ģ��ɶ��B�H��A�Ystep>=2

      do step_n=2,total_steps  !�s�[�@�ӭp��index, step_n

        wavemaker_uh_bc = wavemaker_uhArr(step_n)
        wavemaker_h_bc = wavemaker_hArr(step_n)

        write(*,'(a5,i5)') 'step=',step_n   ! print the step
        call predict(total_nodes,step,dt,g,zb_x,h,uh,h_x,uh_x,u2h_x,Cz,h_min,coriolis,uh_xx,u_max,kin)
!print*, uh(2487,step),uh_x(2487,step),h(2487,step),h_x(2487,step)
        call NF_bc(total_nodes,h(:,step),uh(:,step),kin,wavemaker_uh_bc,wavemaker_h_bc)
        call ghost_BC(total_nodes,kin,h(:,step),uh(:,step),nearest_b,nearest_i,x)
        call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),uh(:,step)&
                ,h_x(:,step),uh_x(:,step),u2h_x(:,step),zb,zb_x,h_min,uh_xx(:,step),u_max,ra,kin,w_velocity(:,step))
!print*, uh(2487,step),uh_x(2487,step),h(2487,step),h_x(2487,step)
        do i_correct =1,n_correct

            call calculate_q(total_nodes,step,dt,h,w_velocity,q)
            call MLS_q(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,q,q_x,h(:,step),h_min,zb,kin)
!print*, uh(2487,step),uh_x(2487,step),h(2487,step),h_x(2487,step)
            call correct(i_correct,n_correct,total_nodes,step,dt,g,zb_x,h,uh,h_x,uh_x,u2h_x,Cz &
                            ,h_min,coriolis,uh_xx,u_max,kin,q,q_x)
!print*, uh(2487,step),uh_x(2487,step),h(2487,step),h_x(2487,step)
            call NF_bc(total_nodes,h(:,step),uh(:,step),kin,wavemaker_uh_bc,wavemaker_h_bc)
            call ghost_BC(total_nodes,kin,h(:,step),uh(:,step),nearest_b,nearest_i,x)
            call MLS(total_nodes,x,m,n_near,weighting,lambda,n_loc,loc_index,h(:,step),uh(:,step) &
                    ,h_x(:,step),uh_x(:,step),u2h_x(:,step),zb,zb_x,h_min,uh_xx(:,step),u_max,ra,kin,w_velocity(:,step))

        end do

        ! for compare with experiment
        compareWithexperiment(step_n) = ( h(gaugePosition,step) + zb(gaugePosition) ) / h0


2       do j=1,total_nodes
          h(j,0)=h(j,1)
          uh(j,0)=uh(j,1)
          h_x(j,0)=h_x(j,1)
          uh_x(j,0)=uh_x(j,1)
          u2h_x(j,0)=u2h_x(j,1)
          uh_xx(j,0)=uh_xx(j,1)
          w_velocity(j,0) = w_velocity(j,1)

          h(j,1)=h(j,2)
          uh(j,1)=uh(j,2)
          h_x(j,1)=h_x(j,2)
          uh_x(j,1)=uh_x(j,2)
          u2h_x(j,1)=u2h_x(j,2)
          uh_xx(j,1)=uh_xx(j,2)
          w_velocity(j,1) = w_velocity(j,2)

          ! data for animation
          animation((step_n - 1) * total_nodes + j) = h(j,step) + zb(j)    ! 1376 2902

        end do



        ! step data output

        if(step_n.eq.total_steps)then

          open(3,file='C:\Users\laserlab04\Desktop\zeta.txt')

          do j = 1,total_nodes
              write(3,'(11f16.5)')x(j) , h(j,step)+zb(j)
          end do

          close(3)


        end if


     end do

    ! output experiment data
!    open(3,file='C:\Users\laserlab04\Desktop\compareWithexperiment.txt')
!    do j = 1,total_steps
!        write(3,'(11f16.5)')compareWithexperiment(j)
!    end do
!    close(3)

    ! data for animation
!    open(3,file='C:\Users\laserlab04\Desktop\animation.txt')
!    do j = 1,total_nodes * total_steps
!
!        write(3,'(11f16.5)')animation(j)
!
!    end do
!    close(3)

    end

	! function ------------------------

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

    subroutine local_poly(x,p,m)
      implicit none
      integer m
      real*8 x,p(1:m)
      p(1)=1.0d0
      p(2)=x
      p(3)=x**2
    end

    subroutine KNN(N,x,d_near,n_near,cor,x_c,region,star)
      implicit none
      integer N,i,j,k,n_near,cor(1:N),region(1:N),star
      real*8 x(1:N),r,d_near(1:n_near),x_c(1:N)
  !  �@�}�l�̪�Z�����]1.d99
      do i=1,n_near
        d_near(i)=1.d99
      end do
  !  �}�l�j�M
      do j=1,N
        if(cor(star)*cor(j) .lt. 0) then
          if((x_c(star) .eq. x_c(j))) then
            if(((x_c(star)-x(star))-(x(j)-x(star)))*cor(star) .lt. 0) then
              goto 1  !���Q�B���쪺�I�n�ư�
            end if
          end if
        end if
        if(abs(region(star)-region(j)) .gt. 1) then
          goto 1  !���b�P�@�Ϫ��I�]�n�ư�
        end if
        if(abs(x(j)-x(star)) .gt. 1.d5) goto 1
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

    subroutine MLS(N,x,m,n_near,w,lambda,n_loc,loc_index,h,uh,h_x,uh_x,u2h_x,zb,zb_x,h_min &
        ,uh_xx,u_max,ra,kin,w_velocity)
      implicit none
      integer i,j,k,l,m,N,n_near,n_loc(1:N),loc_index(1:N,1:n_near*5),wet_neighbors(1:N),wet(1:N),kin(1:N)&
             ,shore(1:N),n_nearest_wet
      real*8 x(1:N),y(1:N),h(1:N),uh(1:N),new_h(1:N),new_uh(1:N)&
            ,h_x(1:N),uh_x(1:N),u2h_x(1:N)&
            ,lambda(1:N,1:m,1:n_near*5),w(1:N,1:n_near*5),zb(1:N),zb_x(1:N),h_min &
            ,uh_xx(1:N),theta(1:N),eta(1:N),ww &
            ,smooth_eta(1:N),smooth_uh(1:N),err_eta(1:N),err_uh(1:N),err_vh(1:N)
      real*8 alpha(1:m),beta(1:n_near*5),p(1:m) &
            ,a(1:n_near*5,1:m),at(1:m,1:n_near*5),ata(1:m,1:m),ata_inv(1:m,1:m)
      real*8 x_loc,y_loc,r,u,v,u_max,abs_u,ra(1:N),r_nearest_wet,lambda_e(1:m,1:n_near*5)&
             ,w_e(1:n_near*5),r_nearest_dry(1:N)
      real*8 w_velocity(1:N),new_w_velocity(1:N),aa

      do j=1,N
        new_h(j)=h(j)
        new_uh(j)=uh(j)
        if(h(j) .le. h_min) then   !�ʲ��P�_����
          wet(j)=0
        else
          wet(j)=1
        end if
      end do
  !��X��F���I�����I�A�íp��Z�̪��I���h��
      do j=1,N
        shore(j)=0  !�����]��0
        r_nearest_dry(j)=1.d99 !�����]���ܤj
        do k=1,n_loc(j)
          x_loc=x(loc_index(j,k))-x(j)
          r=(x_loc**2)**.5d0
          if(r .lt. ra(j)*4.d0/3.d0) then                                      !��F���I�����A
            if(wet(j) .eq. 0 .and. wet(loc_index(j,k)) .eq. 1) then            !�����I�A�B�o�@�I�����O���I
              if(zb(j)+h(j) .lt. zb(loc_index(j,k))+h(loc_index(j,k))) then    !�ӥ������I�����{���F���I�������٧C
                shore(j)=1                                                     !�h���O��F���I�����I�A�n�����I�B�z
              end if
            end if
          end if
          if(r .lt. r_nearest_dry(j) .and. wet(loc_index(j,k)) .eq. 0) then    !��F���I�����A�O���I�B�Z�����̪��I�٤p
            r_nearest_dry(j)=r                                                 !�h��s��̪��I���Z��
          end if
        end do
      end do
  !���s�w�q���I
      do j=1,N
        wet(j)=wet(j)+shore(j) !���O���I�A�h����1�F�Y�O���I�A���n�����I�B�z�A�h0+1=1�F���I�Z���I�ܻ��A�h��0
      end do
  !�p��ƶq�Τ�v
      do j=1,N
        wet_neighbors(j)=0
        do k=1,n_loc(j)
          if(wet(loc_index(j,k)) .eq. 1) then
            wet_neighbors(j)=wet_neighbors(j)+1
          end if
        end do
        theta(j)=1.d0*wet_neighbors(j)/n_loc(j) ! theta �N��P�����I����v�A�Y��1�h��ܩP�򳣬O���I�A�Y��0�h��ܩP�򳣬O���I

      end do

  !�����I�w�qeta�ȡC���I��eta�ä��s�b�A�γ̪����I��eta�ȷ���eta�C�Y�̪����I����ӥH�W�A�h�������C�Y�d�򤺧����L���I�A�h�Ω��ɰ��{��eta�ȡC
      do j=1,N
        if(theta(j) .eq. 0) then
          eta(j)=1.d99  !�Y�d�򤺧����L���I�A�h�]eta=1.99�A�b�p��eta���ɼƮɡA�o�����I�n�ư��A�Y�n���s�غc�����x�}
          ! change this part to eta(j) = zb(j)
          !eta(j) = zb(j)
          goto 11
        end if
        if(wet(j) .eq. 1) then
          eta(j)=zb(j)+h(j)  !�Y�����O���I�A�h��eta=zb+h
          goto 11
        end if
        eta(j)=0   !�ѤU���N�O���I
        ww=0
        do k=1,n_loc(j)
          if(wet(loc_index(j,k)) .eq. 1) then
            ww=ww+w(j,k)**2
            eta(j)=eta(j)+(h(loc_index(j,k))+zb(loc_index(j,k)))*w(j,k)**2
          end if
        end do
        eta(j)=eta(j)/ww
11    end do
  !�}�l�p��
    ! h������
      do j=1,N
        if(kin(j) .eq. -1) goto 21            !ghost point�Τ���U���z�q�����ɼơC
        if(wet(j) .eq. 0 .or. theta(j) .lt. .15d0) then !�Y�o�@�I�O���I�A�νd�����I�ƶq�ܤ֡A�h�]�����y�ʡC
          h_x(j)=-zb_x(j)
          goto 21
        end if
      !�p��beta_eta
        do k=1,n_loc(j)
          if(eta(loc_index(j,k)) .gt. 1.d98) then
            goto 22 !�Y�d�򤺦�eta=1.d99���I�A�h�n���s�غc�����x�}
          end if
          !beta(k)=w(j,k)*h(loc_index(j,k))
          beta(k)=w(j,k)*eta(loc_index(j,k))
        end do
      !�p��alpha_eta
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      !��X h,h_x,h_y
        !new_h(j)=alpha(1)  ! h=eta-zb
        !new_h(j)=h(j)
        !h_x(j)=alpha(2)
         new_h(j)=eta(j) - zb(j)
         h_x(j)=alpha(2) - zb_x(j)

        goto 21
      !���s�غc�����x�}
      !�p��A�x�}�����e
22      do k=1,n_loc(j)
          if(eta(loc_index(j,k)) .gt. 1.d98) then
            w_e(k)=0 !�Yeta=1.d99�A�h�ư��A�]�v����0
          else
            w_e(k)=w(j,k) !���M�N�ӭ쥻���v��
          end if
          x_loc=x(loc_index(j,k))-x(j)
          call local_poly(x_loc,p,m) !�p��Local Polynomial
          do i=1,m
            a(k,i)=w_e(k)*p(i)  ! ��eta�ɡA�u��eta������1.d99���I�Ӻ�C
          end do
        end do
      !����A�x�}����m�x�}
        do k=1,n_loc(j)
          do i=1,m
            at(i,k)=a(k,i)
          end do
        end do
      !����Psuudo�x�}
        do i=1,m
          do l=1,m
            ata(i,l)=0
            do k=1,n_loc(j)
              ata(i,l)=ata(i,l)+at(i,k)*a(k,l)    !atxa=at*a
            end do
          end do
        end do
        call matrix_inv(m,ata,ata_inv)
      !����lambda�x�}
        do i=1,m
          do l=1,n_loc(j)
            lambda_e(i,l)=0
            do k=1,m
              lambda_e(i,l)=lambda_e(i,l)+ata_inv(i,k)*at(k,l)
            end do
          end do
        end do
      !�p��beta_eta
        do k=1,n_loc(j)
          beta(k)=w_e(k)*h(loc_index(j,k))
        end do
      !�p��alpha_eta
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda_e(i,k)*beta(k)
          end do
        end do
      !��X h,h_x,h_y
        !new_h(j)=alpha(1) ! h=eta-zb
        new_h(j)=h(j)
        h_x(j)=alpha(2)

21    end do


    ! uh,u2h
      do j=1,N
        if(kin(j) .eq. -1) goto 31            !ghost point�Τ���U���z�q�����ɼơC
        if(wet(j) .eq. 0 .or. theta(j) .lt. .15d0) then !�Y�o�@�I�O���I�A�νd�����I�ƶq�ܤ֡A�h�]�����y�ʡC
          new_uh(j)=0
          uh_x(j)=0
          uh_xx(j)=0
          u2h_x(j)=0
          goto 31
        end if
      !�p��beta_uh
        do k=1,n_loc(j)
          beta(k)=w(j,k)*uh(loc_index(j,k))
        end do
      !�p��alpha_uh
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      !��X uh,uh_x,uh_y
        !new_uh(j)=alpha(1)
        new_uh(j)=uh(j)
        uh_x(j)=alpha(2)
        uh_xx(j)=2*alpha(4)
      !�p��beta_u2h
        do k=1,n_loc(j)
          if(h(loc_index(j,k)) .le. h_min) then
            beta(k)=0
          else
            beta(k)=w(j,k)*uh(loc_index(j,k))**2/h(loc_index(j,k))
          end if
        end do
      !�p��alpha_u2h
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      !��X u2h_x
        u2h_x(j)=alpha(2)
31    end do

!-----------------------------------------------------------------------------
    !�ˬd�O�_�����`���t�ȩάy�t�L�j
      ww=3.5d0
      do j=1,N
        if(kin(j) .eq. -1) goto 32

        if(r_nearest_dry(j) .gt. ww*ra(j)) then !��̪��I�Z���W�Lww����ra
          h(j)=new_h(j)
          uh(j)=new_uh(j)
        else
          h(j)=new_h(j)*(r_nearest_dry(j)/(ww*ra(j)))+h(j)*(1-r_nearest_dry(j)/(ww*ra(j)))
          uh(j)=new_uh(j)
        end if

        if(h(j) .lt. 0) then
          h(j)=0
        end if

        if(h(j) .le. h_min) then
          uh(j)=0
          goto 32
        end if

        u=uh(j)/h(j)
        abs_u=(u**2)**.5d0

!        if(abs_u .gt. u_max) then
!          uh(j)=uh(j)*u_max/abs_u
!        end if

!------------------------------------------------------------------

32    end do

!      print*, 'MLS '
!      print*, h(2792),h_x(2792),zb(2792),zb_x(2792)
      !print*, 'MLS uh uh_x'
      !print*, uh(2792),uh_x(2792)

        do j = 1,N

            if(wet(j) .eq. 0 .or. theta(j) .lt. .15d0 ) then

               w_velocity(j) = 0

            else

               w_velocity(j) = 0.5 * ( -uh_x(j) + ( uh(j) / h(j) ) * (h_x(j) + 2 * zb_x(j) ) )

            end if

        end do

!      ! smooth w
!           do i=1,8
!            do j=1,N
!              aa=0
!              new_w_velocity(j)=0
!              do k=1,n_loc(j)
!                aa=aa+w(j,k)**2
!                new_w_velocity(j)=new_w_velocity(j)+w(j,k)**2*w_velocity(loc_index(j,k))
!              end do
!              new_w_velocity(j)=new_w_velocity(j)/aa
!            end do
!            do j=1,N
!              w_velocity(j)=new_w_velocity(j)
!            end do
!          end do

!      !�}�l�p��
!          do j=1,N
!          !�p��beta_w
!            do k=1,n_loc(j)
!              beta(k)=w(j,k)*w_velocity(loc_index(j,k))
!            end do
!          !�p��alpha_w
!            do i=1,m
!              alpha(i)=0.d0
!              do k=1,n_loc(j)
!                alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
!              end do
!            end do
!          !��X smooth_w
!            new_w_velocity(j)=alpha(1)
!            new_w_velocity(j)=w_velocity(j)
!          end do
!      !�N��
!          do j=1,N
!            w_velocity(j)=new_w_velocity(j)
!          end do

    end


    subroutine predict(total_nodes,step,dt,g,zb_x,h,uh,h_x,uh_x,u2h_x,Cz0,h_min,coriolis,uh_xx,u_max,kin)
      implicit none
      integer j,total_nodes,step,kin(1:total_nodes)
      real*8 h(1:total_nodes,0:2),uh(1:total_nodes,0:2),h_x(1:total_nodes,0:2),uh_x(1:total_nodes,0:2) &
            ,u2h_x(1:total_nodes,0:2),dt,g,zb_x(1:total_nodes),Cz,Cz0,u,abs_u,h_min,coriolis(1:total_nodes) &
            ,Eddy,uh_xx(1:total_nodes,0:2),u_max

      Eddy=0.d0

      do j=1,total_nodes
        if(kin(j) .eq. -1) goto 1

        if(step .eq. 1) then

          ! h------
          h(j,step)=h(j,step-1)-dt * uh_x(j,step-1)

          if(h(j,step) .lt. 0) then
            h(j,step)=0
          end if

          ! uh------
          if(h(j,step) .le. h_min) then
            uh(j,step)=0
            goto 1
          end if

          ! step - 1
          if(h(j,step-1) .le. h_min) then
            u=0
          else
            u=uh(j,step-1)/h(j,step-1)
          end if

          uh(j,step)=uh(j,step-1)-dt * u2h_x(j,step-1)                             !convective terms
          uh(j,step)=uh(j,step)-dt*g*h(j,step-1)*(zb_x(j)+h_x(j,step-1))           !variable bottom elevation and gravity terms

        else

          ! h------
          h(j,step)=h(j,step-1)-dt * (1.5d0 * uh_x(j,step-1) - .5d0 * uh_x(j,step-2))

          if(h(j,step) .lt. 0) then
            h(j,step)=0
          end if

          ! uh-----
          if(h(j,step) .le. h_min) then
            uh(j,step)=0
            goto 1
          end if

          ! step - 1
          if(h(j,step-1) .le. h_min) then
            u=0
          else
            u=1.5d0*uh(j,step-1)/h(j,step-1)
          end if

          ! step - 2
          if(h(j,step-2) .le. h_min) then
            u=u
          else
            u=u-.5d0*uh(j,step-2)/h(j,step-2)
          end if


          uh(j,step)=uh(j,step-1)-dt*( 1.5d0 * u2h_x(j,step-1) - .5d0 * u2h_x(j,step-2) )
          uh(j,step)=uh(j,step)-dt*g*(1.5d0 * h(j,step-1) * (zb_x(j)+h_x(j,step-1)) &
                                    - .5d0 * h(j,step-2) * (zb_x(j)+h_x(j,step-2)))

        end if

        u=uh(j,step)/h(j,step)


1     end do

!    print*, 'in predict'
!    print*, zb_x(2792),h_x(2792,step-1)
!    print*, zb_x(2792),h_x(2792,step-2)

    end


    subroutine correct(i_correct,n_correct,total_nodes,step,dt,g,zb_x,h,uh,h_x,uh_x,u2h_x,Cz0,h_min&
        ,coriolis,uh_xx,u_max,kin,q,q_x)
      implicit none
      integer j,total_nodes,step,kin(1:total_nodes),i_correct,n_correct
      real*8 h(1:total_nodes,0:2),uh(1:total_nodes,0:2),h_x(1:total_nodes,0:2),uh_x(1:total_nodes,0:2)&
            ,u2h_x(1:total_nodes,0:2),dt,g,zb_x(1:total_nodes),Cz,Cz0,u,v,abs_u,h_min,coriolis(1:total_nodes) &
            ,Eddy,uh_xx(1:total_nodes,0:2),u_max
      real*8 q(1:total_nodes),q_x(1:total_nodes),rr

      Eddy=0.d0
      rr=1.d0*i_correct/n_correct

      do j=1,total_nodes
        if(kin(j) .eq. -1) goto 1
        ! h------
        h(j,step)=h(j,step-1)-dt * 0.5 * (uh_x(j,step)+uh_x(j,step-1))

        if(h(j,step) .lt. 0) then
          h(j,step)=0
        end if

        ! uh------
        if(h(j,step) .le. h_min) then
          uh(j,step)=0
          goto 1
        end if
        u=uh(j,step)/h(j,step)

        if(h(j,step-1) .gt. h_min) then
          u=u+uh(j,step-1)/h(j,step-1)
        end if

        u=u/2

        uh(j,step)=uh(j,step-1)- 0.5 * dt*(u2h_x(j,step)+u2h_x(j,step-1))                       !convective terms

        uh(j,step)=uh(j,step)- 0.5 * dt*g*(h(j,step)*(zb_x(j)+h_x(j,step)) &
                                   +h(j,step-1)*(zb_x(j)+h_x(j,step-1)))                 !variable bottom elevation and gravity terms

        ! non-hydrostatic terms


        uh(j,step) = uh(j,step) - 0.5 * dt * 0.5 * (h(j,step-1) + h(j,step)) * q_x(j)

        uh(j,step) = uh(j,step) - 0.5 * dt * q(j) * 0.5 * ( h_x(j,step-1) + h_x(j,step) + 2*zb_x(j) )


        u=uh(j,step)/h(j,step)


1     end do


!        print*, 'in correct'
!        print*, zb_x(2792),h_x(2792,step-1)
!        print*, zb_x(2792),h_x(2792,step)


    end


    subroutine NF_bc(N,h,uh,kin,wavemaker_uh_bc,wavemaker_h_bc)
      implicit none
      integer j,N,kin(1:N)
      real*8 h(1:N),uh(1:N),wavemaker_uh_bc,wavemaker_h_bc

      do j=1,N
        ! no flow bc

        if(kin(j) .eq. 1) then
            uh(j) = 0
            h(j) = 0
        end if

        ! wave maker bc

        if(kin(j) .eq. 2) then
            uh(j) = 0
            h(j) = 0.218
            !uh(j) = wavemaker_uh_bc
            !h(j) = wavemaker_h_bc
        end if

      end do
    end


    subroutine NS_bc(N,uh,vh,kin)  ! not using this function
      implicit none
      integer j,N,kin(N)
      real*8 uh(1:N),vh(1:N)
      do j=1,N
        if(kin(j) .eq. 2) then
          uh(j)=0
          vh(j)=0
        end if
      end do
    end


!    subroutine FF_bc(N,n_x,h,uh,kin,h_min)  ! not using this function
!      implicit none
!      integer j,N,kin(1:N)
!      real*8 n_x(1:N),h(1:N),uh(1:N),u,v,abs_u,cr_u,h_min
!      do j=1,N
!        if(kin(j) .eq. 3) then
!          if(h(j) .gt. h_min) then
!            u=uh(j)/h(j)
!            v=vh(j)/h(j)  !----v can be deleted----
!            abs_u=(u**2+v**2)**.5  !----abs_u = abs(u)----
!            cr_u=(9.81*h(j))**.5
!            if(abs_u .lt. cr_u) then
!              uh(j)=n_x(j)*cr_u*h(j)
!              vh(j)=n_y(j)*cr_u*h(j)
!            end if
!          end if
!        end if
!      end do
!
!      do j = 1,N
!
!        if(kin(j) .eq. 2)then
!
!        end if
!
!      end do
!
!    end


    subroutine eta_bc(N,h,zb,kin,eta,n_eta_BC)
      implicit none
      integer j,N,kin(N),n_eta_BC
      real*8 h(1:N),zb(1:N),eta(1001:1000+n_eta_BC)
      do j=1,N
        if((kin(j) .gt. 1.d3) .and. (kin(j) .le. 5.d3)) then
          h(j)=eta(kin(j))-zb(j)
        end if
      end do
    end


    subroutine flux_bc(N,uh,vh,kin,flux_x,flux_y,n_flux_BC)  ! not using this function
      implicit none
      integer j,N,kin(N),n_flux_BC
      real*8 uh(1:N),vh(1:N),flux_x(5001:5000+n_flux_BC),flux_y(5001:5000+n_flux_BC)
      do j=1,N
        if(kin(j) .gt. 5.d3) then
          uh(j)=flux_x(kin(j))
          vh(j)=flux_y(kin(j)) !----vh----
        end if
      end do
    end


    subroutine MLS_zb(N,m,n_near,w,lambda,n_loc,loc_index,zb,zb_x,kin,nearest_i)
      implicit none
      integer i,j,k,m,N,n_near,n_loc(1:N),loc_index(1:N,1:n_near*5),kin(1:N),nearest_i(1:N)
      real*8 zb(1:N),zb_x(1:N),new_zb(1:N),lambda(1:N,1:m,1:n_near*5),w(1:N,1:n_near*5)
      real*8 alpha(1:m),beta(1:n_near*5)
! ����ghost point��zb�������������I��zb
      do j=1,N
        if(kin(j) .eq. -1) then
          zb(j)=zb(nearest_i(j))
        end if
      end do
!�}�l�p��
      do j=1,N
  !�p��beta_zb
        do k=1,n_loc(j)
          beta(k)=w(j,k)*zb(loc_index(j,k))
        end do
  !�p��alpha_zb
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
  !��X zb,zb_x,zb_y
        !new_zb(j)=alpha(1)
        new_zb(j)=zb(j)

        if(j .eq. 3000)then
          zb_x(j)=0
        else
          zb_x(j)=alpha(2)
        end if

      end do
      do j=1,N
        zb(j)=new_zb(j)
      end do
    end


    subroutine pre_MLS(N,m,n_near,x,w,lambda,n_loc,loc_index,kin,cor,x_c,region,ra,nearest_b,nearest_i)
      implicit none
      integer i,j,k,l,m,N,n_near,kin(1:N),n_loc(1:N),loc_index(1:N,1:n_near*5),cor(1:N),region(1:N) &
             ,nearest_b(1:N),nearest_i(1:N)
      real*8 x(1:N),w(1:N,1:n_near*5),lambda(1:N,1:m,1:n_near*5),d_near(1:n_near),p(1:m),x_c(1:N)
      real*8 r,x_loc,ro,epsilon,ra(1:N),r_nearest_b
      real*8 a(1:n_near*5,1:m),at(1:m,1:n_near*5),ata(1:m,1:m),ata_inv(1:m,1:m),x_i
      epsilon=-10d0
  !�}�l��Collocation
      do j=1,N
  !���n_near�̪��I���Z��
        call KNN(N,x,d_near,n_near,cor,x_c,region,j)
        ra(j)=(d_near(2)+d_near(3)+d_near(4)+d_near(5))/4.d0 ! ra �j���Olocal nodal spacing
        ro=max(ra(j)*4.5d0,d_near(n_near))*1.001
  !�j�I�íp���v��
        n_loc(j)=0
        do l=1,N !�v�@�j�I
          if(cor(l)*cor(j) .lt. 0) then
            if((x_c(l) .eq. x_c(j))) then
              if(((x_c(j)-x(j))-(x(l)-x(j)))*cor(j) .lt. 0) then   ! might cause problems
                goto 1 !���Q�B���쪺�I�n�ư�
              end if
            end if
          end if
          if(abs(region(l)-region(j)) .gt. 1) then
            goto 1 !���b�P�@�Ϫ��I�]�n�ư�
          end if
          x_loc=x(l)-x(j) !�p��۹理��
          r=abs(x_loc) !�p���focused node�������Z��  !----r = x_loc----
          if(r .ge. ro) goto 1  !�Y�Zfocused node�ӻ��A�]�n�ư�
          n_loc(j)=n_loc(j)+1
          loc_index(j,n_loc(j))=l !�o�̪�local_index�N�O�פ�̪�k�C�Y��l��node�b��j�ӧ����d�򤺡A��s����k
          w(j,n_loc(j))=(1-r/ro)**3.65d0 !abs((exp(epsilon*(r/ro)**2)-exp(epsilon))/(1.d0-exp(epsilon)))**.5d0 !�p���v���Y��
1       end do
      end do
! �j�Mghost point�ҹ���������I
      do j=1,N
        if(kin(j) .ge. 0) goto 2
        r_nearest_b=1.d99
        do k=1,n_loc(j)
          x_loc=x(loc_index(j,k))-x(j) !�p��۹理��
          r=abs(x_loc)  !�p���focused node�������Z��  !----r = x_loc----
          if(kin(loc_index(j,k)) .gt. 0) then !������k�I�O����I�A�����Ikin�O0�Aghost point��kin�O-1�A���O�����]���Oghost�N�O����I
            if(r .lt. r_nearest_b) then !�Z����쥻��r_nearest_b�٤p
              r_nearest_b=r
              nearest_b(j)=loc_index(j,k) !�o��̪�����I�������s��
            end if
          end if
        end do
2     end do
! �j�Mghost point�ҹ����������I
      do j=1,N
        if(kin(j) .ge. 0) goto 3
        x_i=2*x(nearest_b(j))-x(j) !x_i, y_i�����������I����m
        do k=1,n_loc(j)
          x_loc=x(loc_index(j,k))-x_i  !���p��
          r=abs(x_loc)
          if(r .lt. 1.d-6) then
            nearest_i(j)=loc_index(j,k) !�Yr��G0�A�h���N�Oghost point�ҹ����������I
            goto 3
          end if
        end do
3     end do
!�غc�����x�}
      do j=1,N
  !�p��A�x�}�����e
        do k=1,n_loc(j)
          x_loc=x(loc_index(j,k))-x(j)  !�p��
          call local_poly(x_loc,p,m) !�p��Local Polynomial
          do i=1,m
            a(k,i)=w(j,k)*p(i)
          end do
        end do
  !����A�x�}����m�x�}
        do k=1,n_loc(j)
          do i=1,m
            at(i,k)=a(k,i)
          end do
        end do
  !����Psuudo�x�}
        do i=1,m
          do l=1,m
            ata(i,l)=0
            do k=1,n_loc(j)
              ata(i,l)=ata(i,l)+at(i,k)*a(k,l)    !atxa=at*a
            end do
          end do
        end do
  !����Psuudo�x�}���ϯx�}
        call matrix_inv(m,ata,ata_inv)
  !����lambda�x�}
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


    subroutine og_bc(N,n_near,n_loc,loc_index,kin,x,y,h,uh,vh,ra)  ! not using this function
      implicit none
      integer j,k,N,kin(N),n_loc(1:N),loc_index(1:N,1:n_near*5),n_near,n_count
      real*8 h(1:N),uh(1:N),vh(1:N),x(1:N),y(1:N),ra(1:N),r,x_loc,y_loc
      do j=1,N
        if(kin(j) .eq. 4) then
          h(j)=0
          uh(j)=0
          vh(j)=0           !----vh----
          n_count=0
          do k=1,n_loc(j)
            x_loc=x(loc_index(j,k))-x(j)
            y_loc=y(loc_index(j,k))-y(j)  !----y_loc----
            r=(x_loc**2+y_loc**2)**.5     !----r = x_loc----
            if(r .le. 1.5d0*ra(j)) then !���{
              if(kin(loc_index(j,k)) .ne. 4) then !�Dout going���
                n_count=n_count+1
                h(j)=h(j)+h(loc_index(j,k))
                uh(j)=uh(j)+uh(loc_index(j,k))
                vh(j)=vh(j)+vh(loc_index(j,k))  !----vh----
              end if
            end if
          end do
          h(j)=h(j)/n_count  !out going����I��h�ȵ�����{�Dout going����Ih�Ȥ�����
          uh(j)=uh(j)/n_count !out going����I��uh�ȵ�����{�Dout going����Iuh�Ȥ�����
          vh(j)=vh(j)/n_count !out going����I��vh�ȵ�����{�Dout going����Ivh�Ȥ����� !----vh----
        end if
      end do
    end



    subroutine ghost_BC(N,kin,h,uh,nearest_b,nearest_i,x)
      implicit none
      integer j,N,kin(N),nearest_b(1:N),nearest_i(1:N)
      real*8 h(1:N),uh(1:N),x(1:N),n_x,n_abs,A
      do j=1,N
        if(kin(j) .eq. -1) then

          h(j)=h(nearest_i(j))
          n_x=x(j)-x(nearest_b(j))
          n_abs=(n_x**2)**.5d0
          n_x=n_x/n_abs
          A = -n_x**2
          uh(j) = A * uh(nearest_i(j))

        end if
      end do
    end


    subroutine calculate_q(total_nodes,step,dt,h,w_velocity,q)
      implicit none
      integer j,total_nodes,step
      real*8 h(1:total_nodes,0:2),w_velocity(1:total_nodes,0:2),q(1:total_nodes),dt,W_t
      do j=1,total_nodes
	    W_t = (w_velocity(j,step)-w_velocity(j,step-1)) / dt

        q(j) = 0.5 * (h(j,step)+h(j,step-1)) * W_t

      end do

    end

    subroutine MLS_q(N,x,m,n_near,w,lambda,n_loc,loc_index,q,q_x, h,h_min,zb,kin)
      implicit none
      integer i,j,k,l,m,N,n_near,n_loc(1:N),loc_index(1:N,1:n_near*5),n_smooth,kin(1:N)
      real*8 x(1:N),q(1:N),q_x(1:N),lambda(1:N,1:m,1:n_near*5),w(1:N,1:n_near*5)
      real*8 alpha(1:m),beta(1:n_near*5),smooth_q(1:N),aa
      ! part 1
      real*8 wet(1:N),h_min,h(1:N)
      ! part 2
      real*8 shore(1:N),r_nearest_dry(1:N),x_loc,zb(1:N), r, ra(1:N)
      ! part 3
      real*8 wet_neighbors(1:N),theta(1:N),eta(1:N),ww
      real*8 lambda_e(1:m,1:n_near*5),w_e(1:n_near*5)
      real*8 p(1:m),a(1:n_near*5,1:m),at(1:m,1:n_near*5),ata(1:m,1:m),ata_inv(1:m,1:m)

	  n_smooth=8

	! part 1
      do j=1,N
        if(h(j) .le. h_min) then   !�ʲ��P�_����
          wet(j)=0
        else
          wet(j)=1
        end if
      end do

    ! part 2
      !��X��F���I�����I�A�íp��Z�̪��I���h��
      do j=1,N
        shore(j)=0  !�����]��0
        r_nearest_dry(j)=1.d99 !�����]���ܤj
        do k=1,n_loc(j)
          x_loc=x(loc_index(j,k))-x(j)
          r=(x_loc**2)**.5d0
          if(r .lt. ra(j)*4.d0/3.d0) then                                      !��F���I�����A
            if(wet(j) .eq. 0 .and. wet(loc_index(j,k)) .eq. 1) then            !�����I�A�B�o�@�I�����O���I
              if(zb(j)+h(j) .lt. zb(loc_index(j,k))+h(loc_index(j,k))) then    !�ӥ������I�����{���F���I�������٧C
                shore(j)=1                                                     !�h���O��F���I�����I�A�n�����I�B�z
              end if
            end if
          end if
          if(r .lt. r_nearest_dry(j) .and. wet(loc_index(j,k)) .eq. 0) then    !��F���I�����A�O���I�B�Z�����̪��I�٤p
            r_nearest_dry(j)=r                                                 !�h��s��̪��I���Z��
          end if
        end do
      end do

    !���s�w�q���I
      do j=1,N
        wet(j)=wet(j)+shore(j) !���O���I�A�h����1�F�Y�O���I�A���n�����I�B�z�A�h0+1=1�F���I�Z���I�ܻ��A�h��0
      end do

    ! part 3
      !�p��ƶq�Τ�v
      do j=1,N
        wet_neighbors(j)=0
        do k=1,n_loc(j)
          if(wet(loc_index(j,k)) .eq. 1) then
            wet_neighbors(j)=wet_neighbors(j)+1
          end if
        end do
        theta(j)=1.d0*wet_neighbors(j)/n_loc(j) ! theta �N��P�����I����v�A�Y��1�h��ܩP�򳣬O���I�A�Y��0�h��ܩP�򳣬O���I
      end do
      !�����I�w�qeta�ȡC���I��eta�ä��s�b�A�γ̪����I��eta�ȷ���eta�C�Y�̪����I����ӥH�W�A�h�������C�Y�d�򤺧����L���I�A�h�Ω��ɰ��{��eta�ȡC
      do j=1,N
        if(theta(j) .eq. 0) then
          eta(j)=1.d99  !�Y�d�򤺧����L���I�A�h�]eta=1.99�A�b�p��eta���ɼƮɡA�o�����I�n�ư��A�Y�n���s�غc�����x�}
          goto 11
        end if
        if(wet(j) .eq. 1) then
          eta(j)=zb(j)+h(j)  !�Y�����O���I�A�h��eta=zb+h
          goto 11
        end if
        eta(j)=0   !�ѤU���N�O���I
        ww=0
        do k=1,n_loc(j)
          if(wet(loc_index(j,k)) .eq. 1) then
            ww=ww+w(j,k)**2
            eta(j)=eta(j)+(h(loc_index(j,k))+zb(loc_index(j,k)))*w(j,k)**2
          end if
        end do
        eta(j)=eta(j)/ww
11    end do


    ! Q
    !����Q���Ƥ@�U
    do i=1,n_smooth
      do j=1,N
	    aa=0
		smooth_q(j)=0
        do k=1,n_loc(j)
		  aa=aa+w(j,k)**2
          smooth_q(j)=smooth_q(j)+w(j,k)**2*q(loc_index(j,k))
        end do
        smooth_q(j)=smooth_q(j)/aa
	  end do
      do j=1,N
        q(j)=smooth_q(j)
	  end do
	end do

    ! Q������
      do j=1,N
        if(kin(j) .eq. -1) goto 21            !ghost point�Τ���U���z�q�����ɼơC
        if(wet(j) .eq. 0 .or. theta(j) .lt. .15d0) then !�Y�o�@�I�O���I�A�νd�����I�ƶq�ܤ֡A�h�]�����y�ʡC
          q(j) = 0
          q_x(j) = 0
          goto 21
        end if
      !�p��beta_q
        do k=1,n_loc(j)
          if(eta(loc_index(j,k)) .gt. 1.d98) then
            goto 22 !�Y�d�򤺦�eta=1.d99���I�A�h�n���s�غc�����x�}
          end if
          beta(k)=w(j,k)*q(loc_index(j,k))
        end do
      !�p��alpha_q
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda(j,i,k)*beta(k)
          end do
        end do
      !��X q,q_x
        !q(j)=alpha(1)  ! h=eta-zb
        q_x(j)=alpha(2)
        goto 21
      !���s�غc�����x�}
      !�p��A�x�}�����e
22      do k=1,n_loc(j)
          if(eta(loc_index(j,k)) .gt. 1.d98) then
            w_e(k)=0 !�Yeta=1.d99�A�h�ư��A�]�v����0
          else
            w_e(k)=w(j,k) !���M�N�ӭ쥻���v��
          end if
          x_loc=x(loc_index(j,k))-x(j)
          call local_poly(x_loc,p,m) !�p��Local Polynomial
          do i=1,m
            a(k,i)=w_e(k)*p(i)  ! ��eta�ɡA�u��eta������1.d99���I�Ӻ�C
          end do
        end do
      !����A�x�}����m�x�}
        do k=1,n_loc(j)
          do i=1,m
            at(i,k)=a(k,i)
          end do
        end do
      !����Psuudo�x�}
        do i=1,m
          do l=1,m
            ata(i,l)=0
            do k=1,n_loc(j)
              ata(i,l)=ata(i,l)+at(i,k)*a(k,l)    !atxa=at*a
            end do
          end do
        end do
        call matrix_inv(m,ata,ata_inv)
      !����lambda�x�}
        do i=1,m
          do l=1,n_loc(j)
            lambda_e(i,l)=0
            do k=1,m
              lambda_e(i,l)=lambda_e(i,l)+ata_inv(i,k)*at(k,l)
            end do
          end do
        end do
      !�p��beta_eta
        do k=1,n_loc(j)
          beta(k)=w_e(k)*q(loc_index(j,k))
        end do
      !�p��alpha_eta
        do i=1,m
          alpha(i)=0.d0
          do k=1,n_loc(j)
            alpha(i)=alpha(i)+lambda_e(i,k)*beta(k)
          end do
        end do
      !��X h,h_x,h_y
        !q(j)=alpha(1) ! h=eta-zb
        q_x(j)=alpha(2)
21    end do

    end

